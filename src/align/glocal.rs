//! Glocal pair-HMM: port of htslib `probaln_glocal` (Li 2011). Global in the
//! query, local in the reference, with explicit begin/end states and per-row
//! scaling. One kernel serves three consumers:
//! - BAQ: forward-backward gives the MAP state and its posterior per base;
//! - haplotype scoring: the product of the row scaling factors is
//!   P(query | reference), so the forward pass alone is a likelihood;
//! - indel discovery: the MAP state path maps each query base to a reference
//!   column or an insertion.

/// Emission for an inserted base and for a mismatch (`EI`, `EM` in probaln.c).
const EI: f64 = 0.25;
const EM: f64 = 0.333_333_333_33;

#[derive(Clone, Copy, Debug)]
pub struct GlocalParams {
    /// Gap open probability.
    pub d: f64,
    /// Gap extension probability.
    pub e: f64,
    /// Band width.
    pub bw: usize,
}

impl Default for GlocalParams {
    /// `sam_prob_realn` parameters.
    fn default() -> Self {
        Self { d: 0.001, e: 0.1, bw: 10 }
    }
}

pub struct GlocalResult {
    /// log P(query | reference): Σ ln s\[i\] plus ln(l_ref·l_query), the
    /// begin-state normalisation `probaln_glocal` folds into its score.
    pub loglik: f64,
    /// Per query base: `(ref_index << 2) | state`, state 0 = match, 1 = insertion.
    pub state: Vec<i32>,
    /// Per query base: phred of P(MAP state), capped at 99.
    pub q: Vec<u8>,
}

/// Nucleotide code: A/C/G/T = 0..3, anything else 4 (treated as N).
#[inline]
pub fn encode_nt(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

static QUAL2PROB: std::sync::LazyLock<[f64; 256]> = std::sync::LazyLock::new(|| {
    let mut t = [0.0; 256];
    for (i, v) in t.iter_mut().enumerate() {
        *v = 10f64.powf(-(i as f64) / 10.0);
    }
    t
});

/// Band-relative cell offset of reference column `k` in query row `i`
/// (`set_u` in probaln.c). May fall outside `0..bw2*3`; the matrices carry
/// zero padding on both sides of every row so such cells read as 0.
#[inline]
fn set_u(bw: usize, i: usize, k: usize) -> isize {
    let x = i.saturating_sub(bw) as isize;
    (k as isize - x) * 3
}

/// Forward-only likelihood: log P(query | reference).
pub fn glocal_loglik(reference: &[u8], query: &[u8], iqual: Option<&[u8]>, c: &GlocalParams) -> f64 {
    glocal(reference, query, iqual, c, false).map(|r| r.loglik).unwrap_or(f64::NEG_INFINITY)
}

/// `reference` and `query` are nucleotide codes ([`encode_nt`]); `iqual` are
/// phred base qualities (30 when absent). With `posterior`, the backward pass
/// fills `state` and `q`.
pub fn glocal(reference: &[u8], query: &[u8], iqual: Option<&[u8]>, c: &GlocalParams, posterior: bool) -> Option<GlocalResult> {
    let l_ref = reference.len();
    let l_query = query.len();
    if l_ref == 0 || l_query == 0 {
        return None;
    }
    let mut bw = l_ref.max(l_query);
    if bw > c.bw {
        bw = c.bw;
    }
    let diff = l_ref.abs_diff(l_query);
    if bw < diff {
        bw = diff;
    }
    let bw2 = bw * 2 + 1;
    let row = bw2 * 3;
    // Three padding cells on each side of a row absorb out-of-band reads.
    let stride = row + 6;
    let idx = |i: usize, u: isize| -> usize { i * stride + (u + 3) as usize };

    let mut f = vec![0.0f64; (l_query + 1) * stride];
    let mut b = if posterior { vec![0.0f64; (l_query + 1) * stride] } else { Vec::new() };
    let mut s = vec![0.0f64; l_query + 2];
    let q2p = &*QUAL2PROB;
    // 1-based query index -> error probability (30 when no quality is given).
    let qual = |i: usize| -> f64 { q2p[iqual.and_then(|v| v.get(i - 1)).map(|&q| q as usize).unwrap_or(30)] };
    let query1 = |i: usize| -> u8 { query[i - 1] };
    let ref1 = |k: usize| -> u8 { reference[k - 1] };
    let emit = |rk: u8, qy: u8, e: f64| -> f64 {
        if rk > 3 || qy > 3 {
            1.0
        } else if rk == qy {
            1.0 - e
        } else {
            e * EM
        }
    };

    let s_m = 1.0 / (2.0 * l_query as f64 + 2.0);
    let s_i = s_m;
    let m = [
        (1.0 - c.d - c.d) * (1.0 - s_m),
        c.d * (1.0 - s_m),
        c.d * (1.0 - s_m),
        (1.0 - c.e) * (1.0 - s_i),
        c.e * (1.0 - s_i),
        0.0,
        1.0 - c.e,
        0.0,
        c.e,
    ];
    let b_m = (1.0 - c.d) / l_ref as f64;
    let b_i = c.d / l_ref as f64;

    // Forward: row 0.
    f[idx(0, set_u(bw, 0, 0))] = 1.0;
    s[0] = 1.0;
    {
        // Row 1.
        let beg = 1usize;
        let end = l_ref.min(bw + 1);
        let mut sum = 0.0;
        for k in beg..=end {
            let e = emit(ref1(k), query1(1), qual(1));
            let u = set_u(bw, 1, k);
            f[idx(1, u)] = e * b_m;
            f[idx(1, u + 1)] = EI * b_i;
            sum += f[idx(1, u)] + f[idx(1, u + 1)];
        }
        s[1] = sum;
        let (lo, hi) = (set_u(bw, 1, beg), set_u(bw, 1, end) + 2);
        for u in lo..=hi {
            f[idx(1, u)] /= sum;
        }
    }
    for i in 2..=l_query {
        let qli = qual(i);
        let qyi = query1(i);
        let beg = 1usize.max(i.saturating_sub(bw));
        let end = l_ref.min(i + bw);
        let mut sum = 0.0;
        for k in beg..=end {
            let e = emit(ref1(k), qyi, qli);
            let u = set_u(bw, i, k);
            let v11 = set_u(bw, i - 1, k - 1);
            let v10 = set_u(bw, i - 1, k);
            let v01 = set_u(bw, i, k - 1);
            let fm = e * (m[0] * f[idx(i - 1, v11)] + m[3] * f[idx(i - 1, v11 + 1)] + m[6] * f[idx(i - 1, v11 + 2)]);
            let fi = EI * (m[1] * f[idx(i - 1, v10)] + m[4] * f[idx(i - 1, v10 + 1)]);
            let fd = m[2] * f[idx(i, v01)] + m[8] * f[idx(i, v01 + 2)];
            f[idx(i, u)] = fm;
            f[idx(i, u + 1)] = fi;
            f[idx(i, u + 2)] = fd;
            sum += fm + fi + fd;
        }
        if sum <= 0.0 {
            return None;
        }
        s[i] = sum;
        let inv = 1.0 / sum;
        let (lo, hi) = (set_u(bw, i, beg), set_u(bw, i, end) + 2);
        for u in lo..=hi {
            f[idx(i, u)] *= inv;
        }
    }
    {
        // Row l_query + 1: the end state.
        let mut sum = 0.0;
        for k in 1..=l_ref {
            let u = set_u(bw, l_query, k);
            if u < 3 || u >= (row + 3) as isize {
                continue;
            }
            sum += f[idx(l_query, u)] * s_m + f[idx(l_query, u + 1)] * s_i;
        }
        s[l_query + 1] = sum;
    }
    let loglik: f64 = s[1..=l_query + 1].iter().map(|v| if *v > 0.0 { v.ln() } else { f64::NEG_INFINITY }).sum::<f64>()
        + ((l_ref * l_query) as f64).ln();

    if !posterior {
        return Some(GlocalResult { loglik, state: Vec::new(), q: Vec::new() });
    }

    // Backward: row l_query.
    for k in 1..=l_ref {
        let u = set_u(bw, l_query, k);
        if u < 3 || u >= (row + 3) as isize {
            continue;
        }
        b[idx(l_query, u)] = s_m / s[l_query] / s[l_query + 1];
        b[idx(l_query, u + 1)] = s_i / s[l_query] / s[l_query + 1];
    }
    for i in (1..l_query).rev() {
        let beg = 1usize.max(i.saturating_sub(bw));
        let end = l_ref.min(i + bw);
        let y = if i > 1 { 1.0 } else { 0.0 };
        let qli1 = qual(i + 1);
        let qyi1 = query1(i + 1);
        for k in (beg..=end).rev() {
            let u = set_u(bw, i, k);
            let v11 = set_u(bw, i + 1, k + 1);
            let v10 = set_u(bw, i + 1, k);
            let v01 = set_u(bw, i, k + 1);
            let e = if k >= l_ref { 0.0 } else { emit(ref1(k + 1), qyi1, qli1) } * b[idx(i + 1, v11)];
            let bm = e * m[0] + EI * m[1] * b[idx(i + 1, v10 + 1)] + m[2] * b[idx(i, v01 + 2)];
            let bi = e * m[3] + EI * m[4] * b[idx(i + 1, v10 + 1)];
            let bd = (e * m[6] + m[8] * b[idx(i, v01 + 2)]) * y;
            b[idx(i, u)] = bm;
            b[idx(i, u + 1)] = bi;
            b[idx(i, u + 2)] = bd;
        }
        let inv = 1.0 / s[i];
        let (lo, hi) = (set_u(bw, i, beg), set_u(bw, i, end) + 2);
        for u in lo..=hi {
            b[idx(i, u)] *= inv;
        }
    }

    // MAP state per query base.
    let mut state = vec![0i32; l_query];
    let mut q = vec![0u8; l_query];
    for i in 1..=l_query {
        let beg = 1usize.max(i.saturating_sub(bw));
        let end = l_ref.min(i + bw);
        let mut sum = 0.0;
        let mut max = 0.0;
        let mut max_k: i32 = -1;
        for k in beg..=end {
            let u = set_u(bw, i, k);
            let z = f[idx(i, u)] * b[idx(i, u)];
            if z > max {
                max = z;
                max_k = ((k as i32 - 1) << 2) | 0;
            }
            sum += z;
            let z = f[idx(i, u + 1)] * b[idx(i, u + 1)];
            if z > max {
                max = z;
                max_k = ((k as i32 - 1) << 2) | 1;
            }
            sum += z;
        }
        if sum > 0.0 {
            max /= sum;
        }
        state[i - 1] = max_k;
        let kq = if max >= 1.0 { 101 } else { (-4.343 * (1.0 - max).ln() + 0.499) as i32 };
        q[i - 1] = if kq > 100 { 99 } else { kq.max(0) as u8 };
    }
    Some(GlocalResult { loglik, state, q })
}

#[cfg(test)]
#[path = "../../tests/unit/align_glocal.rs"]
mod tests;
