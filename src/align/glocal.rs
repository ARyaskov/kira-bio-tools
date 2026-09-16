//! Glocal pair-HMM: port of htslib `probaln_glocal` (Li 2011). Global in the
//! query, local in the reference, with explicit begin/end states and per-row
//! scaling. One kernel serves three consumers:
//! - BAQ: forward-backward gives the MAP state and its posterior per base;
//! - haplotype scoring: the product of the row scaling factors is
//!   P(query | reference), so the forward pass alone is a likelihood;
//! - indel discovery: the MAP state path maps each query base to a reference
//!   column or an insertion.
//!
//! The three states of a band row live in separate planes (structure of
//! arrays). Within a row the match and insertion cells depend only on the
//! previous row and are computed in one straight loop the compiler
//! vectorises; the deletion cells, which chain along the row, follow in a
//! short serial pass. Every cell is evaluated with the same operations in the
//! same order as probaln.c, so posteriors are bit-identical to htslib's.

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
        Self {
            d: 0.001,
            e: 0.1,
            bw: 10,
        }
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

/// `g_qual2prob` is a `float` table in probaln.c; the single-precision
/// rounding of the error probabilities is kept so posteriors agree.
static QUAL2PROB: std::sync::LazyLock<[f64; 256]> = std::sync::LazyLock::new(|| {
    let mut t = [0.0; 256];
    for (i, v) in t.iter_mut().enumerate() {
        *v = (10f64.powf(-(i as f64) / 10.0)) as f32 as f64;
    }
    t
});

/// Reusable matrices for [`glocal_with`]: one per thread avoids the
/// allocations per aligned read. Each plane holds `(l_query + 1)` rows of
/// `width` cells: the band plus one zero cell on each side, so a neighbour
/// one column outside the band reads as 0 (probaln.c's padding).
#[derive(Default)]
pub struct GlocalScratch {
    fm: Vec<f64>,
    fi: Vec<f64>,
    fd: Vec<f64>,
    bm: Vec<f64>,
    bi: Vec<f64>,
    bd: Vec<f64>,
    s: Vec<f64>,
    /// Per-row emission and the previous row's neighbours, gathered once.
    em: Vec<f64>,
}

/// Forward-only likelihood: log P(query | reference).
pub fn glocal_loglik(
    reference: &[u8],
    query: &[u8],
    iqual: Option<&[u8]>,
    c: &GlocalParams,
) -> f64 {
    glocal(reference, query, iqual, c, false)
        .map(|r| r.loglik)
        .unwrap_or(f64::NEG_INFINITY)
}

/// `reference` and `query` are nucleotide codes ([`encode_nt`]); `iqual` are
/// phred base qualities (30 when absent). With `posterior`, the backward pass
/// fills `state` and `q`.
pub fn glocal(
    reference: &[u8],
    query: &[u8],
    iqual: Option<&[u8]>,
    c: &GlocalParams,
    posterior: bool,
) -> Option<GlocalResult> {
    glocal_with(
        &mut GlocalScratch::default(),
        reference,
        query,
        iqual,
        c,
        posterior,
    )
}

#[inline]
fn emit(rk: u8, qy: u8, e: f64) -> f64 {
    if rk > 3 || qy > 3 {
        1.0
    } else if rk == qy {
        1.0 - e
    } else {
        e * EM
    }
}

/// [`glocal`] with caller-owned matrices.
pub fn glocal_with(
    scratch: &mut GlocalScratch,
    reference: &[u8],
    query: &[u8],
    iqual: Option<&[u8]>,
    c: &GlocalParams,
    posterior: bool,
) -> Option<GlocalResult> {
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
    // One zero cell on each side of the band.
    let width = bw2 + 2;
    let rows = l_query + 1;

    let GlocalScratch {
        fm,
        fi,
        fd,
        bm,
        bi,
        bd,
        s,
        em,
    } = scratch;
    for v in [&mut *fm, &mut *fi, &mut *fd] {
        v.clear();
        v.resize(rows * width, 0.0);
    }
    for v in [&mut *bm, &mut *bi, &mut *bd] {
        v.clear();
        if posterior {
            v.resize(rows * width, 0.0);
        }
    }
    s.clear();
    s.resize(l_query + 2, 0.0);
    em.clear();
    em.resize(width, 0.0);

    let q2p = &*QUAL2PROB;
    // 1-based query index -> error probability (30 when no quality is given).
    let qual = |i: usize| -> f64 {
        q2p[iqual
            .and_then(|v| v.get(i - 1))
            .map(|&q| q as usize)
            .unwrap_or(30)]
    };
    let query1 = |i: usize| -> u8 { query[i - 1] };
    let ref1 = |k: usize| -> u8 { reference[k - 1] };
    // Band origin of row i (`x` in probaln.c) and the band-relative column
    // of reference column k, offset by the padding cell.
    let x_of = |i: usize| -> usize { i.saturating_sub(bw) };
    let cell = |i: usize, j: usize| -> usize { i * width + j + 1 };

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
    fm[cell(0, 0)] = 1.0;
    s[0] = 1.0;
    {
        // Row 1.
        let beg = 1usize;
        let end = l_ref.min(bw + 1);
        let x = x_of(1);
        let (q1, e1) = (query1(1), qual(1));
        let mut sum = 0.0;
        for k in beg..=end {
            let j = k - x;
            let e = emit(ref1(k), q1, e1);
            fm[cell(1, j)] = e * b_m;
            fi[cell(1, j)] = EI * b_i;
            sum += fm[cell(1, j)] + fi[cell(1, j)];
        }
        s[1] = sum;
        for k in beg..=end {
            let j = k - x;
            fm[cell(1, j)] /= sum;
            fi[cell(1, j)] /= sum;
            fd[cell(1, j)] /= sum;
        }
    }
    for i in 2..=l_query {
        let qli = qual(i);
        let qyi = query1(i);
        let beg = 1usize.max(i.saturating_sub(bw));
        let end = l_ref.min(i + bw);
        let x = x_of(i);
        // Column j of this row reads column j + d - 1 (diagonal) and j + d
        // (above) of the previous row.
        let d = x - x_of(i - 1);
        let (jb, je) = (beg - x, end - x);
        for j in jb..=je {
            em[j] = emit(ref1(j + x), qyi, qli);
        }
        // Match and insertion: previous row only, one straight loop.
        let (prev, cur) = (cell(i - 1, 0), cell(i, 0));
        {
            // The previous row ends where this row starts.
            let (fm_head, fm_tail) = fm.split_at_mut(cur - 1);
            let (fi_head, fi_tail) = fi.split_at_mut(cur - 1);
            let pm = &fm_head[prev - 1..prev + width - 1];
            let pi = &fi_head[prev - 1..prev + width - 1];
            let pd = &fd[prev - 1..prev + width - 1];
            let em = &em[..];
            let fm_row = &mut fm_tail[..width];
            let fi_row = &mut fi_tail[..width];
            for j in jb..=je {
                // Padded index: neighbour column c of the previous row sits at c + 1.
                let diag = j + d;
                let up = j + d + 1;
                fm_row[j + 1] = em[j] * (m[0] * pm[diag] + m[3] * pi[diag] + m[6] * pd[diag]);
                fi_row[j + 1] = EI * (m[1] * pm[up] + m[4] * pi[up]);
            }
        }
        // Deletion: chains along the row.
        let mut sum = 0.0;
        for j in jb..=je {
            let left = cell(i, j) - 1;
            let dcell = m[2] * fm[left] + m[8] * fd[left];
            fd[cell(i, j)] = dcell;
            sum += fm[cell(i, j)] + fi[cell(i, j)] + dcell;
        }
        if sum <= 0.0 {
            return None;
        }
        s[i] = sum;
        let inv = 1.0 / sum;
        for j in jb..=je {
            let c0 = cell(i, j);
            fm[c0] *= inv;
            fi[c0] *= inv;
            fd[c0] *= inv;
        }
    }
    {
        // Row l_query + 1: the end state.
        let x = x_of(l_query);
        let mut sum = 0.0;
        for k in 1..=l_ref {
            if k < x || k >= x + bw2 {
                continue;
            }
            let c0 = cell(l_query, k - x);
            sum += fm[c0] * s_m + fi[c0] * s_i;
        }
        s[l_query + 1] = sum;
    }
    let loglik: f64 = s[1..=l_query + 1]
        .iter()
        .map(|v| if *v > 0.0 { v.ln() } else { f64::NEG_INFINITY })
        .sum::<f64>()
        + ((l_ref * l_query) as f64).ln();

    if !posterior {
        return Some(GlocalResult {
            loglik,
            state: Vec::new(),
            q: Vec::new(),
        });
    }

    // Backward: row l_query.
    {
        let x = x_of(l_query);
        for k in 1..=l_ref {
            if k < x || k >= x + bw2 {
                continue;
            }
            let c0 = cell(l_query, k - x);
            bm[c0] = s_m / s[l_query] / s[l_query + 1];
            bi[c0] = s_i / s[l_query] / s[l_query + 1];
        }
    }
    for i in (1..l_query).rev() {
        let beg = 1usize.max(i.saturating_sub(bw));
        let end = l_ref.min(i + bw);
        let y = if i > 1 { 1.0 } else { 0.0 };
        let qli1 = qual(i + 1);
        let qyi1 = query1(i + 1);
        let x = x_of(i);
        // Column j of this row reads column j - d + 1 (diagonal) and j - d
        // (below) of the next row.
        let d = x_of(i + 1) - x;
        let (jb, je) = (beg - x, end - x);
        // e * B_M of the diagonal cell, then the insertion state.
        let (next, cur) = (cell(i + 1, 0), cell(i, 0));
        {
            // This row ends where the next row starts.
            let nm = &bm[next - 1..next + width - 1];
            let (bi_head, bi_tail) = bi.split_at_mut(next - 1);
            let ni = &bi_tail[..width];
            let em = &mut em[..];
            let bi_row = &mut bi_head[cur - 1..cur + width - 1];
            for j in jb..=je {
                let k = j + x;
                let diag = j + 2 - d;
                let below = j + 1 - d;
                let e = if k >= l_ref {
                    0.0
                } else {
                    emit(ref1(k + 1), qyi1, qli1)
                } * nm[diag];
                em[j] = e;
                bi_row[j + 1] = e * m[3] + EI * m[4] * ni[below];
            }
        }
        // Deletion chains from the right.
        for j in (jb..=je).rev() {
            let right = cell(i, j) + 1;
            bd[cell(i, j)] = (em[j] * m[6] + m[8] * bd[right]) * y;
        }
        {
            let ni = &bi[next - 1..next + width - 1];
            let em = &em[..];
            let (bd_row, bm_row) = (
                &bd[cur - 1..cur + width - 1],
                &mut bm[cur - 1..cur + width - 1],
            );
            for j in jb..=je {
                let below = j + 1 - d;
                bm_row[j + 1] = em[j] * m[0] + EI * m[1] * ni[below] + m[2] * bd_row[j + 2];
            }
        }
        let inv = 1.0 / s[i];
        for j in jb..=je {
            let c0 = cell(i, j);
            bm[c0] *= inv;
            bi[c0] *= inv;
            bd[c0] *= inv;
        }
    }

    // MAP state per query base.
    let mut state = vec![0i32; l_query];
    let mut q = vec![0u8; l_query];
    for i in 1..=l_query {
        let beg = 1usize.max(i.saturating_sub(bw));
        let end = l_ref.min(i + bw);
        let x = x_of(i);
        let mut sum = 0.0;
        let mut max = 0.0;
        let mut max_k: i32 = -1;
        for k in beg..=end {
            let c0 = cell(i, k - x);
            let z = fm[c0] * bm[c0];
            if z > max {
                max = z;
                max_k = (k as i32 - 1) << 2;
            }
            sum += z;
            let z = fi[c0] * bi[c0];
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
        let kq = if max >= 1.0 {
            101
        } else {
            (-4.343 * (1.0 - max).ln() + 0.499) as i32
        };
        q[i - 1] = if kq > 100 { 99 } else { kq.max(0) as u8 };
    }
    Some(GlocalResult { loglik, state, q })
}

#[cfg(test)]
#[path = "../../tests/unit/align_glocal.rs"]
mod tests;
