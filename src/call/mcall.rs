//! Multi-allelic Bayesian variant caller.
//!
//! Algorithm (matches bcftools `-m`):
//! 1. Convert per-sample FORMAT/PL → P(D|G) via 10^(-PL/10) lookup
//! 2. Enumerate allele subsets (singleton, pair, triple)
//! 3. For each subset: compute QS-weighted AF, sum log-likelihoods under HWE,
//!    add theta prior for each non-ref allele in subset
//! 4. Pick max-likelihood subset
//! 5. QUAL = -10*log10(P(ref_only|D) / P(best|D))  ≈  -4.343 * (ref_lk - max_lk)
//! 6. Per-sample GT call: argmax over genotypes consistent with chosen subset

use super::math::*;

#[derive(Clone, Debug)]
pub struct CallerOpts {
    pub theta: f64,
    pub indel_theta: f64,
    pub keep_alts: bool,
    pub variants_only: bool,
    pub min_ac: u32,
    pub ploidy: u8,
    pub per_sample_ploidy: Option<Vec<u8>>,
    pub prior_af: Option<Vec<f64>>,
    pub gvcf: Option<GvcfOpts>,
    pub groups: Option<Vec<SampleGroup>>,
    pub constrain: ConstrainMode,
    pub family: Option<TrioFamily>,
}

#[derive(Clone, Debug)]
pub struct GvcfOpts {
    pub min_dp: u32,
    pub min_qual: f64,
    pub blocks: Vec<u32>,
}

#[derive(Clone, Debug)]
pub struct SampleGroup {
    pub name: String,
    pub sample_idxs: Vec<usize>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConstrainMode {
    None,
    Alleles,
    Trio,
}

#[derive(Clone, Debug)]
pub struct TrioFamily {
    pub father: Option<usize>,
    pub mother: Option<usize>,
    pub child: Option<usize>,
    pub is_son: bool,
}

impl Default for CallerOpts {
    fn default() -> Self {
        Self {
            theta: 1.1e-3,
            indel_theta: 1.5e-4,
            keep_alts: false,
            variants_only: false,
            min_ac: 0,
            ploidy: 2,
            per_sample_ploidy: None,
            prior_af: None,
            gvcf: None,
            groups: None,
            constrain: ConstrainMode::None,
            family: None,
        }
    }
}

pub struct Caller {
    pub opts: CallerOpts,
    pl2p: [f64; PL2P_SIZE],
    log_theta_snv: f64,
    log_theta_indel: f64,
}

impl Caller {
    pub fn new(opts: CallerOpts, n_samples: usize) -> Self {
        let n_total = total_alleles(n_samples, &opts);
        let a_m = watterson_factor(n_total);
        let scaled_snv = (opts.theta * a_m).min(0.99);
        let scaled_indel = (opts.indel_theta * a_m).min(0.99);
        let log_theta_snv = if scaled_snv > 0.0 { scaled_snv.ln() } else { f64::NEG_INFINITY };
        let log_theta_indel = if scaled_indel > 0.0 { scaled_indel.ln() } else { f64::NEG_INFINITY };
        Self { opts, pl2p: init_pl2p(), log_theta_snv, log_theta_indel }
    }

    fn ploidy_for(&self, sample_idx: usize) -> u8 {
        if let Some(p) = &self.opts.per_sample_ploidy {
            return p.get(sample_idx).copied().unwrap_or(self.opts.ploidy);
        }
        self.opts.ploidy
    }

    fn log_theta_for(&self, is_indel: bool) -> f64 {
        if is_indel { self.log_theta_indel } else { self.log_theta_snv }
    }

    /// Call a single site given per-sample PLs.
    pub fn call_site(&self, site: &mut CallSite) -> CallResult {
        if let Some(groups) = &self.opts.groups {
            return self.call_site_grouped(site, groups);
        }
        self.call_site_single_group(site, None)
    }

    fn call_site_grouped(&self, site: &CallSite, groups: &[SampleGroup]) -> CallResult {
        let mut merged_alleles: Vec<u32> = vec![0];
        let mut per_sample_gts: Vec<(u32, u32)> = vec![(0, 0); site.n_samples];
        let mut per_sample_gqs: Vec<u32> = vec![0; site.n_samples];
        let mut total_qual = 0.0f64;
        let mut any_variant = false;
        for g in groups {
            let sub = self.call_site_single_group(site, Some(&g.sample_idxs));
            if let CallResult::Called { alleles_kept, qual, gts, gqs, .. } = sub {
                for a in &alleles_kept { if !merged_alleles.contains(a) { merged_alleles.push(*a); } }
                for (i, &si) in g.sample_idxs.iter().enumerate() {
                    if si < per_sample_gts.len() && i < gts.len() {
                        per_sample_gts[si] = gts[i];
                        per_sample_gqs[si] = gqs[i];
                    }
                }
                if alleles_kept.len() > 1 { any_variant = true; }
                if qual > total_qual { total_qual = qual; }
            }
        }
        merged_alleles.sort();
        if self.opts.variants_only && !any_variant { return CallResult::Skip; }
        let all_ploidies: Vec<u8> = (0..site.n_samples).map(|si| self.ploidy_for(si)).collect();
        let (ac, an) = compute_ac_an(&per_sample_gts, &merged_alleles, &all_ploidies);
        CallResult::Called {
            alleles_kept: merged_alleles, qual: total_qual,
            gts: per_sample_gts, gqs: per_sample_gqs,
            pls: Vec::new(), ac, an,
        }
    }

    fn call_site_single_group(&self, site: &CallSite, sample_subset: Option<&[usize]>) -> CallResult {
        let n_smpl_orig = site.n_samples;
        let group_idxs: Vec<usize> = sample_subset.map(|s| s.to_vec())
            .unwrap_or_else(|| (0..n_smpl_orig).collect());
        let n_smpl = group_idxs.len();
        let n_als = site.n_alleles;
        let n_gt = n_genotypes(n_als);

        let mut pdg: Vec<f64> = vec![0.0; n_smpl * n_gt];
        for (i, &si) in group_idxs.iter().enumerate() {
            let pls = &site.pls[si * n_gt..(si + 1) * n_gt];
            let row = &mut pdg[i * n_gt..(i + 1) * n_gt];
            init_pdg_row(pls, row, &self.pl2p);
        }

        let qsum = if let Some(af) = &self.opts.prior_af {
            af.iter().take(n_als).map(|f| *f * (n_smpl as f64 * self.opts.ploidy as f64)).collect()
        } else {
            compute_qsum(&pdg, n_smpl, n_als)
        };
        let log_theta = self.log_theta_for(site.is_indel);
        let group_ploidies: Vec<u8> = group_idxs.iter().map(|&si| self.ploidy_for(si)).collect();

        let mut max_lk = f64::NEG_INFINITY;
        let mut max_als: u32 = 0;
        let mut ref_lk = 0.0f64;
        let mut lk_sum = f64::NEG_INFINITY;

        for ia in 0..n_als {
            let lk_tot = single_allele_lk(&pdg, n_smpl, n_gt, ia);
            if ia == 0 { ref_lk = lk_tot; }
            let prior_adj = if ia == 0 { 0.0 } else { log_theta };
            let lk = lk_tot + prior_adj;
            lk_sum = ln_sum_exp(lk_sum, lk);
            if lk > max_lk { max_lk = lk; max_als = 1 << ia; }
        }

        if n_als >= 2 {
            for ia in 0..n_als {
                if qsum[ia] == 0.0 { continue; }
                for ib in 0..ia {
                    if qsum[ib] == 0.0 { continue; }
                    let lk_tot = pair_allele_lk_per_ploidy(&pdg, n_smpl, n_gt, ia, ib, &qsum, &group_ploidies);
                    let mut prior_adj = 0.0;
                    if ia != 0 { prior_adj += log_theta; }
                    if ib != 0 { prior_adj += log_theta; }
                    let lk = lk_tot + prior_adj;
                    lk_sum = ln_sum_exp(lk_sum, lk);
                    if lk > max_lk { max_lk = lk; max_als = (1 << ia) | (1 << ib); }
                }
            }
        }

        if n_als >= 3 {
            for ia in 0..n_als {
                if qsum[ia] == 0.0 { continue; }
                for ib in 0..ia {
                    if qsum[ib] == 0.0 { continue; }
                    for ic in 0..ib {
                        if qsum[ic] == 0.0 { continue; }
                        let lk_tot = triple_allele_lk_per_ploidy(&pdg, n_smpl, n_gt, ia, ib, ic, &qsum, &group_ploidies);
                        let mut prior_adj = 0.0;
                        if ia != 0 { prior_adj += log_theta; }
                        if ib != 0 { prior_adj += log_theta; }
                        if ic != 0 { prior_adj += log_theta; }
                        let lk = lk_tot + prior_adj;
                        lk_sum = ln_sum_exp(lk_sum, lk);
                        if lk > max_lk { max_lk = lk; max_als = (1 << ia) | (1 << ib) | (1 << ic); }
                    }
                }
            }
        }

        let qual = if max_als == 1 && !self.opts.keep_alts {
            if log_theta.is_finite() { M10_OVER_LN10 * log_theta * -1.0 } else { 0.0 }
        } else {
            M10_OVER_LN10 * (ref_lk - max_lk)
        };

        let mut alleles_kept: Vec<u32> = vec![0];
        for i in 1..n_als as u32 {
            if max_als & (1 << i) != 0 { alleles_kept.push(i); }
        }
        if self.opts.keep_alts {
            for i in 1..n_als as u32 {
                if !alleles_kept.contains(&i) { alleles_kept.push(i); }
            }
            alleles_kept.sort();
        }
        let is_pure_ref_site = max_als == 1;
        if is_pure_ref_site && !self.opts.keep_alts {
            alleles_kept = vec![0];
        }
        let _ = lk_sum;

        let mut gts: Vec<(u32, u32)> = Vec::with_capacity(n_smpl);
        let mut gqs: Vec<u32> = Vec::with_capacity(n_smpl);
        let mut new_pls: Vec<i32> = Vec::with_capacity(n_smpl * n_genotypes(alleles_kept.len()));
        let af = compute_af(&qsum, &alleles_kept);

        for i in 0..n_smpl {
            let row = &pdg[i * n_gt..(i + 1) * n_gt];
            let ploidy = group_ploidies[i];
            let (gt, gq, pls_out) = best_gt_for_sample(row, &alleles_kept, &af, ploidy, n_als);
            gts.push(gt);
            gqs.push(gq);
            new_pls.extend(pls_out);
        }

        if self.opts.constrain == ConstrainMode::Trio {
            if let Some(fam) = &self.opts.family {
                apply_trio_constraint(&mut gts, &mut gqs, fam, &group_idxs);
            }
        }

        let (ac, an) = compute_ac_an(&gts, &alleles_kept, &group_ploidies);

        let is_variant = alleles_kept.len() > 1;
        if self.opts.variants_only && !is_variant {
            return CallResult::Skip;
        }

        CallResult::Called {
            alleles_kept,
            qual: qual.max(0.0),
            gts,
            gqs,
            pls: new_pls,
            ac,
            an,
        }
    }
}

pub struct CallSite {
    pub n_samples: usize,
    pub n_alleles: usize,
    pub pls: Vec<i32>,
    pub is_indel: bool,
    pub depths: Option<Vec<u32>>,
}

impl CallSite {
    pub fn new(n_samples: usize, n_alleles: usize, pls: Vec<i32>) -> Self {
        Self { n_samples, n_alleles, pls, is_indel: false, depths: None }
    }
}

fn total_alleles(n_samples: usize, opts: &CallerOpts) -> usize {
    if let Some(p) = &opts.per_sample_ploidy {
        p.iter().take(n_samples).map(|x| *x as usize).sum()
    } else {
        n_samples * opts.ploidy as usize
    }
}

pub enum CallResult {
    Called {
        alleles_kept: Vec<u32>,
        qual: f64,
        gts: Vec<(u32, u32)>,
        gqs: Vec<u32>,
        pls: Vec<i32>,
        ac: Vec<u32>,
        an: u32,
    },
    Skip,
}

fn init_pdg_row(pls: &[i32], pdg: &mut [f64], pl2p: &[f64; PL2P_SIZE]) {
    let n_gt = pdg.len();
    let mut sum = 0.0;
    let mut all_missing = true;
    for j in 0..n_gt {
        let v = pls[j];
        if v == i32::MIN || v == i32::MIN + 1 { pdg[j] = 0.0; continue; }
        pdg[j] = pl_to_prob(v, pl2p);
        sum += pdg[j];
        if v != 0 { all_missing = false; }
    }
    if all_missing && pls[0] == 0 && n_gt > 1 {
        let v = pls.iter().skip(1).all(|&p| p == 0);
        if v { for p in pdg.iter_mut() { *p = 0.0; } return; }
    }
    if sum > 0.0 {
        for p in pdg.iter_mut() { *p /= sum; }
    }
}

fn compute_qsum(pdg: &[f64], n_smpl: usize, n_als: usize) -> Vec<f64> {
    let n_gt = n_genotypes(n_als);
    let mut qsum = vec![0.0f64; n_als];
    for s in 0..n_smpl {
        let row = &pdg[s * n_gt..(s + 1) * n_gt];
        for a in 0..n_als {
            let ii = gt_index(a, a);
            qsum[a] += row[ii];
            for b in 0..n_als {
                if b == a { continue; }
                let ab = gt_index(a, b);
                qsum[a] += 0.5 * row[ab];
            }
        }
    }
    qsum
}

fn single_allele_lk(pdg: &[f64], n_smpl: usize, n_gt: usize, ia: usize) -> f64 {
    let iaa = gt_index(ia, ia);
    let mut lk = 0.0f64;
    for s in 0..n_smpl {
        let v = pdg[s * n_gt + iaa];
        if v > 0.0 { lk += v.ln(); }
    }
    lk
}

fn pair_allele_lk(pdg: &[f64], n_smpl: usize, n_gt: usize, ia: usize, ib: usize, qsum: &[f64], ploidy: u8) -> f64 {
    pair_allele_lk_per_ploidy(pdg, n_smpl, n_gt, ia, ib, qsum, &vec![ploidy; n_smpl])
}

fn pair_allele_lk_per_ploidy(pdg: &[f64], n_smpl: usize, n_gt: usize, ia: usize, ib: usize, qsum: &[f64], ploidies: &[u8]) -> f64 {
    let total = qsum[ia] + qsum[ib];
    if total == 0.0 { return f64::NEG_INFINITY; }
    let fa = qsum[ia] / total;
    let fb = qsum[ib] / total;
    let (fa2, fb2, fab) = (fa * fa, fb * fb, 2.0 * fa * fb);
    let iaa = gt_index(ia, ia);
    let ibb = gt_index(ib, ib);
    let iab = gt_index(ia, ib);
    let mut lk = 0.0f64;
    for s in 0..n_smpl {
        let row = &pdg[s * n_gt..];
        let p = ploidies.get(s).copied().unwrap_or(2);
        let v = if p == 2 {
            fa2 * row[iaa] + fb2 * row[ibb] + fab * row[iab]
        } else if p == 1 {
            fa * row[iaa] + fb * row[ibb]
        } else { 0.0 };
        if v > 0.0 { lk += v.ln(); }
    }
    lk
}

fn triple_allele_lk(pdg: &[f64], n_smpl: usize, n_gt: usize, ia: usize, ib: usize, ic: usize, qsum: &[f64], ploidy: u8) -> f64 {
    triple_allele_lk_per_ploidy(pdg, n_smpl, n_gt, ia, ib, ic, qsum, &vec![ploidy; n_smpl])
}

fn triple_allele_lk_per_ploidy(pdg: &[f64], n_smpl: usize, n_gt: usize, ia: usize, ib: usize, ic: usize, qsum: &[f64], ploidies: &[u8]) -> f64 {
    let total = qsum[ia] + qsum[ib] + qsum[ic];
    if total == 0.0 { return f64::NEG_INFINITY; }
    let fa = qsum[ia] / total;
    let fb = qsum[ib] / total;
    let fc = qsum[ic] / total;
    let (fa2, fb2, fc2) = (fa * fa, fb * fb, fc * fc);
    let (fab, fac, fbc) = (2.0 * fa * fb, 2.0 * fa * fc, 2.0 * fb * fc);
    let iaa = gt_index(ia, ia);
    let ibb = gt_index(ib, ib);
    let icc = gt_index(ic, ic);
    let iab = gt_index(ia, ib);
    let iac = gt_index(ia, ic);
    let ibc = gt_index(ib, ic);
    let mut lk = 0.0f64;
    for s in 0..n_smpl {
        let row = &pdg[s * n_gt..];
        let p = ploidies.get(s).copied().unwrap_or(2);
        let v = if p == 2 {
            fa2 * row[iaa] + fb2 * row[ibb] + fc2 * row[icc]
                + fab * row[iab] + fac * row[iac] + fbc * row[ibc]
        } else if p == 1 {
            fa * row[iaa] + fb * row[ibb] + fc * row[icc]
        } else { 0.0 };
        if v > 0.0 { lk += v.ln(); }
    }
    lk
}

/// Apply Mendelian trio constraint: child's alleles must come one from
/// each parent. If posterior assigns impossible inheritance, swap to
/// best-mendelian-consistent GT.
fn apply_trio_constraint(
    gts: &mut [(u32, u32)], gqs: &mut [u32],
    fam: &TrioFamily, group_idxs: &[usize]
) {
    let pos = |abs_idx: Option<usize>| -> Option<usize> {
        abs_idx.and_then(|s| group_idxs.iter().position(|x| *x == s))
    };
    let f_pos = pos(fam.father);
    let m_pos = pos(fam.mother);
    let c_pos = pos(fam.child);
    let (Some(fi), Some(mi), Some(ci)) = (f_pos, m_pos, c_pos) else { return; };
    let f = gts[fi];
    let m = gts[mi];
    let child = gts[ci];
    if mendelian_ok(f, m, child) { return; }
    let f_alleles = [f.0, f.1];
    let m_alleles = [m.0, m.1];
    let mut best: Option<(u32, u32)> = None;
    for &a in &f_alleles {
        for &b in &m_alleles {
            let cand = (a.min(b), a.max(b));
            best = Some(match best { None => cand, Some(p) => if cand < p { cand } else { p } });
        }
    }
    if let Some(c) = best {
        gts[ci] = c;
        if gqs[ci] > 5 { gqs[ci] -= 5; }
    }
}

fn mendelian_ok(f: (u32, u32), m: (u32, u32), c: (u32, u32)) -> bool {
    let f_set = [f.0, f.1];
    let m_set = [m.0, m.1];
    (f_set.contains(&c.0) && m_set.contains(&c.1)) || (f_set.contains(&c.1) && m_set.contains(&c.0))
}

fn compute_af(qsum: &[f64], kept: &[u32]) -> Vec<f64> {
    let total: f64 = kept.iter().map(|&i| qsum[i as usize]).sum();
    if total == 0.0 { return vec![1.0 / kept.len() as f64; kept.len()]; }
    kept.iter().map(|&i| qsum[i as usize] / total).collect()
}

fn best_gt_for_sample(
    pdg_row: &[f64],
    kept: &[u32],
    af: &[f64],
    ploidy: u8,
    _n_als_orig: usize,
) -> ((u32, u32), u32, Vec<i32>) {
    let nk = kept.len();
    let n_gt_new = n_genotypes(nk);
    let mut posts = vec![f64::NEG_INFINITY; n_gt_new];
    let mut idx = 0;
    for j in 0..nk {
        for i in 0..=j {
            let orig_i = kept[i] as usize;
            let orig_j = kept[j] as usize;
            let prior = if ploidy == 2 {
                let fa = af[i]; let fb = af[j];
                if i == j { (fa * fa).max(1e-30) } else { (2.0 * fa * fb).max(1e-30) }
            } else {
                if i == j { af[i].max(1e-30) } else { 0.0 }
            };
            let pdg_val = pdg_row[gt_index(orig_i, orig_j)];
            posts[idx] = if pdg_val > 0.0 { pdg_val.ln() + prior.ln() } else { f64::NEG_INFINITY };
            idx += 1;
        }
    }
    let max_p = posts.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut best_idx = 0;
    for (i, p) in posts.iter().enumerate() {
        if *p == max_p { best_idx = i; break; }
    }
    let pls: Vec<i32> = posts.iter().map(|p| {
        let phred = if p.is_finite() { (M10_OVER_LN10 * (p - max_p)).round() } else { 255.0 };
        phred.max(0.0).min(255.0) as i32
    }).collect();

    let mut gq: u32 = 0;
    let mut second_best = f64::NEG_INFINITY;
    for (i, p) in posts.iter().enumerate() {
        if i != best_idx && *p > second_best { second_best = *p; }
    }
    if second_best.is_finite() {
        gq = (M10_OVER_LN10 * (second_best - max_p)).max(0.0).min(99.0) as u32;
    }

    let (gi, gj) = idx_to_pair(best_idx);
    let gt_alleles = (kept[gi], kept[gj]);
    (gt_alleles, gq, pls)
}

fn idx_to_pair(idx: usize) -> (usize, usize) {
    let mut k = idx;
    let mut j = 0;
    while k > j { k -= j + 1; j += 1; }
    (k, j)
}

/// AC/AN over the called genotypes; a sample contributes `ploidy` alleles.
fn compute_ac_an(gts: &[(u32, u32)], kept: &[u32], ploidies: &[u8]) -> (Vec<u32>, u32) {
    let mut ac = vec![0u32; kept.len().saturating_sub(1)];
    let mut an = 0u32;
    for (i, &(a, b)) in gts.iter().enumerate() {
        let n = ploidies.get(i).copied().unwrap_or(2).min(2) as usize;
        for &x in [a, b].iter().take(n) {
            an += 1;
            if let Some(pos) = kept.iter().position(|k| *k == x) {
                if pos > 0 { ac[pos - 1] += 1; }
            }
        }
    }
    (ac, an)
}

fn ln_sum_exp(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY { return b; }
    if b == f64::NEG_INFINITY { return a; }
    let m = a.max(b);
    m + ((a - m).exp() + (b - m).exp()).ln()
}

#[cfg(test)]
#[path = "../../tests/unit/call_mcall.rs"]
mod tests;
