//! Runs of Homozygosity (RoH) HMM detector — port of bcftools/vcfroh.c.
//!
//! Two-state HMM: AZ (autozygous) vs HW (Hardy-Weinberg).
//! Emissions: per-site P(GT|state) computed from genotype + AF.
//! Transitions: --hw-to-az and --az-to-hw rates, scaled by physical/genetic distance.
//! Viterbi for best-path; optional Baum-Welch training (-V).

use std::collections::HashMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum State { AZ, HW }

#[derive(Clone, Debug)]
pub struct RohSite {
    pub chrom: String,
    pub pos: u32,
    pub genetic_pos: Option<f64>,
    pub gt: GtClass,
    pub af: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GtClass { HomRef, Het, HomAlt, Missing }

#[derive(Clone, Debug)]
pub struct RohOpts {
    pub hw_to_az: f64,
    pub az_to_hw: f64,
    pub rec_rate: f64,
    pub af_dflt: f64,
    pub ignore_homref: bool,
    pub skip_indels: bool,
    pub viterbi_training: usize,
}

impl Default for RohOpts {
    fn default() -> Self {
        Self {
            hw_to_az: 6.7e-8,
            az_to_hw: 5e-9,
            rec_rate: 1e-8,
            af_dflt: 0.4,
            ignore_homref: false,
            skip_indels: false,
            viterbi_training: 0,
        }
    }
}

#[derive(Clone, Debug)]
pub struct RohSegment {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub state: State,
    pub n_sites: usize,
    pub avg_quality: f64,
    pub total_length: u32,
}

/// Emission probability P(observation | state). AF is alt allele freq.
fn emit_log10(state: State, gt: GtClass, af: f64) -> f64 {
    let af = af.max(1e-6).min(1.0 - 1e-6);
    let q = 1.0 - af;
    match (state, gt) {
        (State::AZ, GtClass::HomRef) => q.log10(),
        (State::AZ, GtClass::HomAlt) => af.log10(),
        (State::AZ, GtClass::Het) => -3.0f64,
        (State::AZ, GtClass::Missing) => 0.0,
        (State::HW, GtClass::HomRef) => (q * q).log10(),
        (State::HW, GtClass::HomAlt) => (af * af).log10(),
        (State::HW, GtClass::Het) => (2.0 * af * q).log10(),
        (State::HW, GtClass::Missing) => 0.0,
    }
}

/// Compute transition probability between two sites at distance d (bp or genetic).
/// Returns log10 P(state_new | state_old).
fn trans_log10(from: State, to: State, opts: &RohOpts, distance_bp: u32, genetic_dist: Option<f64>) -> f64 {
    let dist = genetic_dist.unwrap_or_else(|| distance_bp as f64 * opts.rec_rate);
    let p_switch_hw_az = 1.0 - (-opts.hw_to_az * dist).exp();
    let p_switch_az_hw = 1.0 - (-opts.az_to_hw * dist).exp();
    let p = match (from, to) {
        (State::HW, State::AZ) => p_switch_hw_az,
        (State::HW, State::HW) => 1.0 - p_switch_hw_az,
        (State::AZ, State::HW) => p_switch_az_hw,
        (State::AZ, State::AZ) => 1.0 - p_switch_az_hw,
    };
    p.max(1e-30).log10()
}

/// Viterbi decoding over a chromosome's sites.
pub fn viterbi(sites: &[RohSite], opts: &RohOpts) -> Vec<State> {
    if sites.is_empty() { return Vec::new(); }
    let n = sites.len();
    let mut log_v = vec![[0.0f64; 2]; n];
    let mut backptr = vec![[0u8; 2]; n];

    let prior_az = 0.5f64.log10();
    let prior_hw = 0.5f64.log10();
    log_v[0][0] = prior_az + emit_log10(State::AZ, sites[0].gt, sites[0].af);
    log_v[0][1] = prior_hw + emit_log10(State::HW, sites[0].gt, sites[0].af);

    for i in 1..n {
        let dist = sites[i].pos.saturating_sub(sites[i-1].pos);
        let gd = sites[i].genetic_pos.and_then(|g1| sites[i-1].genetic_pos.map(|g0| g1 - g0));

        for cur_idx in 0..2 {
            let cur = if cur_idx == 0 { State::AZ } else { State::HW };
            let emit = emit_log10(cur, sites[i].gt, sites[i].af);
            let mut best = f64::NEG_INFINITY;
            let mut from_idx = 0u8;
            for prev_idx in 0..2 {
                let prev = if prev_idx == 0 { State::AZ } else { State::HW };
                let cand = log_v[i-1][prev_idx] + trans_log10(prev, cur, opts, dist, gd) + emit;
                if cand > best { best = cand; from_idx = prev_idx as u8; }
            }
            log_v[i][cur_idx] = best;
            backptr[i][cur_idx] = from_idx;
        }
    }

    let mut path = vec![State::HW; n];
    let mut state_idx = if log_v[n-1][0] > log_v[n-1][1] { 0 } else { 1 };
    path[n-1] = if state_idx == 0 { State::AZ } else { State::HW };
    for i in (0..n-1).rev() {
        state_idx = backptr[i+1][state_idx] as usize;
        path[i] = if state_idx == 0 { State::AZ } else { State::HW };
    }
    path
}

/// Scaled forward-backward posteriors. Returns per-site `[P(AZ), P(HW)]`,
/// normalised so the two states sum to 1 at each site. Per-site rescaling
/// keeps the products from underflowing on long chromosomes.
pub fn forward_backward(sites: &[RohSite], opts: &RohOpts) -> Vec<[f64; 2]> {
    let n = sites.len();
    if n == 0 {
        return Vec::new();
    }
    let emit = |i: usize, s: State| -> f64 { 10f64.powf(emit_log10(s, sites[i].gt, sites[i].af)) };
    let trans = |i: usize, from: State, to: State| -> f64 {
        let dist = sites[i].pos.saturating_sub(sites[i - 1].pos);
        let gd = sites[i]
            .genetic_pos
            .and_then(|g1| sites[i - 1].genetic_pos.map(|g0| g1 - g0));
        10f64.powf(trans_log10(from, to, opts, dist, gd))
    };

    let mut fwd = vec![[0.0f64; 2]; n];
    let mut scale = vec![1.0f64; n];
    fwd[0][0] = 0.5 * emit(0, State::AZ);
    fwd[0][1] = 0.5 * emit(0, State::HW);
    let s0 = (fwd[0][0] + fwd[0][1]).max(1e-300);
    fwd[0][0] /= s0;
    fwd[0][1] /= s0;
    scale[0] = s0;
    for i in 1..n {
        for (ci, cur) in [State::AZ, State::HW].into_iter().enumerate() {
            let mut sum = 0.0;
            for (pi, prev) in [State::AZ, State::HW].into_iter().enumerate() {
                sum += fwd[i - 1][pi] * trans(i, prev, cur);
            }
            fwd[i][ci] = sum * emit(i, cur);
        }
        let s = (fwd[i][0] + fwd[i][1]).max(1e-300);
        fwd[i][0] /= s;
        fwd[i][1] /= s;
        scale[i] = s;
    }

    let mut bwd = vec![[0.0f64; 2]; n];
    bwd[n - 1] = [1.0, 1.0];
    for i in (0..n - 1).rev() {
        for (ci, cur) in [State::AZ, State::HW].into_iter().enumerate() {
            let mut sum = 0.0;
            for (ni, nxt) in [State::AZ, State::HW].into_iter().enumerate() {
                sum += trans(i + 1, cur, nxt) * emit(i + 1, nxt) * bwd[i + 1][ni];
            }
            bwd[i][ci] = sum / scale[i + 1];
        }
    }

    let mut post = vec![[0.5f64; 2]; n];
    for i in 0..n {
        let a = fwd[i][0] * bwd[i][0];
        let b = fwd[i][1] * bwd[i][1];
        let s = a + b;
        if s > 0.0 {
            post[i] = [a / s, b / s];
        }
    }
    post
}

/// bcftools `phred_score`: -10·log10(prob), capped at 99 (and 99 when prob<=0).
pub fn phred_score(prob: f64) -> f64 {
    if prob <= 0.0 {
        99.0
    } else {
        (-10.0 * prob.log10()).min(99.0)
    }
}

/// Baum-Welch training (-V N): N iterations of forward-backward + parameter re-estimation.
pub fn baum_welch_train(sites: &[RohSite], opts: &mut RohOpts, n_iters: usize) {
    for _ in 0..n_iters {
        let path = viterbi(sites, opts);
        let mut hw_to_az_count = 0u64;
        let mut hw_count = 0u64;
        let mut az_to_hw_count = 0u64;
        let mut az_count = 0u64;
        for i in 1..path.len() {
            match (path[i-1], path[i]) {
                (State::HW, State::AZ) => { hw_to_az_count += 1; hw_count += 1; }
                (State::HW, State::HW) => { hw_count += 1; }
                (State::AZ, State::HW) => { az_to_hw_count += 1; az_count += 1; }
                (State::AZ, State::AZ) => { az_count += 1; }
            }
        }
        if hw_count > 100 {
            opts.hw_to_az = ((hw_to_az_count as f64) / (hw_count as f64)).max(1e-12);
        }
        if az_count > 100 {
            opts.az_to_hw = ((az_to_hw_count as f64) / (az_count as f64)).max(1e-12);
        }
    }
}

/// Collapse contiguous states into segments.
pub fn segments(sites: &[RohSite], path: &[State]) -> Vec<RohSegment> {
    if sites.is_empty() { return Vec::new(); }
    let mut out = Vec::new();
    let mut cur_start = 0usize;
    let mut cur_state = path[0];
    for i in 1..path.len() {
        if path[i] != cur_state || sites[i].chrom != sites[cur_start].chrom {
            let n = i - cur_start;
            out.push(RohSegment {
                chrom: sites[cur_start].chrom.clone(),
                start: sites[cur_start].pos,
                end: sites[i-1].pos,
                state: cur_state,
                n_sites: n,
                avg_quality: 0.0,
                total_length: sites[i-1].pos.saturating_sub(sites[cur_start].pos) + 1,
            });
            cur_start = i;
            cur_state = path[i];
        }
    }
    let n = path.len() - cur_start;
    out.push(RohSegment {
        chrom: sites[cur_start].chrom.clone(),
        start: sites[cur_start].pos,
        end: sites[path.len()-1].pos,
        state: cur_state,
        n_sites: n,
        avg_quality: 0.0,
        total_length: sites[path.len()-1].pos.saturating_sub(sites[cur_start].pos) + 1,
    });
    out
}

/// Parse genetic map file: `chrom pos cm_per_mb cm` (PLINK format).
pub fn parse_genetic_map<P: AsRef<std::path::Path>>(path: P) -> std::io::Result<HashMap<String, Vec<(u32, f64)>>> {
    use std::io::BufRead;
    let mut map: HashMap<String, Vec<(u32, f64)>> = HashMap::new();
    let f = std::fs::File::open(path)?;
    for line in std::io::BufReader::new(f).lines() {
        let l = line?;
        let t = l.trim();
        if t.is_empty() || t.starts_with('#') || t.starts_with("chr ") { continue; }
        let parts: Vec<&str> = t.split_whitespace().collect();
        if parts.len() < 4 { continue; }
        let chrom = parts[0].to_string();
        // Malformed rows are skipped rather than read as position 0.
        let (Ok(pos), Ok(cm)) = (parts[1].parse::<u32>(), parts[3].parse::<f64>()) else { continue };
        map.entry(chrom).or_default().push((pos, cm));
    }
    for v in map.values_mut() { v.sort_by_key(|x| x.0); }
    Ok(map)
}

/// Estimate AF from sample genotypes in a window (sliding average).
pub fn estimate_af(gts: &[GtClass], window: usize) -> Vec<f64> {
    let mut out = vec![0.5f64; gts.len()];
    if gts.is_empty() { return out; }
    let half = (window / 2).max(1);
    for i in 0..gts.len() {
        let lo = i.saturating_sub(half);
        let hi = (i + half).min(gts.len());
        let mut alt = 0.0;
        let mut total = 0.0;
        for g in &gts[lo..hi] {
            match g {
                GtClass::HomRef => total += 2.0,
                GtClass::Het => { alt += 1.0; total += 2.0; }
                GtClass::HomAlt => { alt += 2.0; total += 2.0; }
                GtClass::Missing => {}
            }
        }
        if total > 0.0 { out[i] = ((alt as f64) / total).max(0.01_f64).min(0.99_f64); }
    }
    out
}

#[cfg(test)]
#[path = "../tests/unit/roh.rs"]
mod tests;
