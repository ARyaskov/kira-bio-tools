//! Base Alignment Quality (BAQ). Port of htslib's sam_prob_realn / probaln_glocal.
//!
//! Reference: Li, H. (2011) "Improving SNP discovery by base alignment quality",
//! Bioinformatics 27(8):1157-1158.
//!
//! Linear-scaled forward-backward (no log/exp in the hot loop) using f32.

use noodles_sam::alignment::record::cigar::op::Kind;

const BAND_WIDTH: usize = 7;
const BAQ_EXTEND_LEN: usize = 50;

/// Apply full HMM BAQ. Caps qualities to BAQ_q computed via banded
/// forward-backward over the read vs the local reference window.
pub fn apply_baq_hmm(
    seq: &[u8],
    qual: &mut [u8],
    cigar: &[(Kind, u32)],
    ref_seq: &[u8],
    bandwidth: usize,
) {
    if ref_seq.is_empty() || qual.is_empty() {
        apply_baq_capping(seq, qual, cigar, 30);
        return;
    }
    let (read_extract, ref_extract, q_to_orig) = extract_aligned(seq, qual, cigar, ref_seq);
    if read_extract.is_empty() || ref_extract.is_empty() { return; }

    let posteriors = banded_forward_backward(&read_extract, &ref_extract, &q_to_orig, qual, bandwidth);

    for (i, p) in posteriors.iter().enumerate() {
        let oi = q_to_orig[i];
        if oi >= qual.len() { continue; }
        let baq_q = if *p > 0.999 { 30 }
            else if *p > 0.0 { (-10.0 * (1.0 - *p).log10()).round() as u32 }
            else { 0 };
        let cap = baq_q.min(255) as u8;
        if qual[oi] > cap { qual[oi] = cap; }
    }
}

fn extract_aligned(
    seq: &[u8],
    qual: &[u8],
    cigar: &[(Kind, u32)],
    ref_seq: &[u8],
) -> (Vec<u8>, Vec<u8>, Vec<usize>) {
    let mut r_pos = 0usize;
    let mut q_pos = 0usize;
    let mut read_extract: Vec<u8> = Vec::new();
    let mut ref_extract: Vec<u8> = Vec::new();
    let mut q_to_orig: Vec<usize> = Vec::new();
    for &(kind, len) in cigar {
        let l = len as usize;
        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for i in 0..l {
                    if q_pos + i < seq.len() && r_pos + i < ref_seq.len() {
                        read_extract.push(seq[q_pos + i].to_ascii_uppercase());
                        ref_extract.push(ref_seq[r_pos + i].to_ascii_uppercase());
                        q_to_orig.push(q_pos + i);
                    }
                }
                q_pos += l;
                r_pos += l;
            }
            Kind::Insertion | Kind::SoftClip => q_pos += l,
            Kind::Deletion | Kind::Skip => r_pos += l,
            _ => {}
        }
    }
    let _ = qual;
    (read_extract, ref_extract, q_to_orig)
}

/// Precomputed (p_match, p_mismatch) per quality, f32. p_match = 1 - e_q,
/// p_mismatch = e_q / 3 where e_q = 10^(-q/10).
static EMIT_TABLE: once_cell::sync::Lazy<[(f32, f32); 65]> =
    once_cell::sync::Lazy::new(|| {
        let mut t = [(0.0f32, 0.0f32); 65];
        for q in 0..=64 {
            let e = (10f64).powf(-(q as f64) / 10.0) as f32;
            t[q] = (1.0 - e, e / 3.0);
        }
        t
    });

#[inline(always)]
fn emit_probs(q: u8) -> (f32, f32) {
    EMIT_TABLE[(q as usize).min(64)]
}

/// Per-thread scratch for forward/backward 3-state matrices in linear space.
#[derive(Default)]
struct BaqScratch {
    // Forward: M, I, D per cell (row-major, width = 2*band+1)
    fm: Vec<f32>,
    fi: Vec<f32>,
    fd: Vec<f32>,
    // Backward
    bm: Vec<f32>,
    bi: Vec<f32>,
    bd: Vec<f32>,
    // Per-row scaling factors c_i so total log-likelihood = -sum(log c_i).
    scale: Vec<f32>,
}

impl BaqScratch {
    #[inline]
    fn resize(&mut self, cells: usize, rows: usize) {
        for v in [&mut self.fm, &mut self.fi, &mut self.fd, &mut self.bm, &mut self.bi, &mut self.bd] {
            v.clear();
            v.resize(cells, 0.0);
        }
        self.scale.clear();
        self.scale.resize(rows, 1.0);
    }
}

thread_local! {
    static BAQ_SCRATCH: std::cell::RefCell<BaqScratch> = std::cell::RefCell::new(BaqScratch::default());
}

/// Banded glocal forward-backward HMM (3-state: M, I, D), linear-scaled f32.
/// Returns posterior P(state=M) for each read position.
fn banded_forward_backward(
    read: &[u8],
    refseq: &[u8],
    q_to_orig: &[usize],
    qual: &[u8],
    bandwidth: usize,
) -> Vec<f64> {
    BAQ_SCRATCH.with(|cell| {
        let mut scratch = cell.borrow_mut();
        banded_forward_backward_inner(read, refseq, q_to_orig, qual, bandwidth, &mut scratch)
    })
}

#[inline(always)]
fn cell_idx(i: usize, j: usize, j_center: usize, width: usize, band: usize) -> Option<usize> {
    let off = (j as i64) - (j_center as i64) + (band as i64);
    if off < 0 || off as usize >= width { return None; }
    Some(i * width + off as usize)
}

fn banded_forward_backward_inner(
    read: &[u8],
    refseq: &[u8],
    q_to_orig: &[usize],
    qual: &[u8],
    bandwidth: usize,
    scratch: &mut BaqScratch,
) -> Vec<f64> {
    let n = read.len();
    let m = refseq.len();
    let band = bandwidth.max(1);
    let width = 2 * band + 1;
    let cells = n * width;

    // Linear-space transitions. Tuned to match the old log-space probabilities.
    const T_MM: f32 = 0.999;       // M -> M
    const T_MI: f32 = 0.0005;      // M -> I
    const T_MD: f32 = 0.0005;      // M -> D
    const T_IM: f32 = 0.5;         // I -> M
    const T_II: f32 = 0.5;         // I -> I
    const T_DM: f32 = 0.5;         // D -> M
    const T_DD: f32 = 0.5;         // D -> D

    scratch.resize(cells, n);
    let BaqScratch { fm, fi, fd, bm, bi, bd, scale } = scratch;

    // Convenience: returns Some(k) if cell (i, j) lies inside the band centered at i.
    // Band is centered on the diagonal i ≈ j.
    let in_band = |i: usize, j: usize| -> Option<usize> {
        let center = i; // j is measured around center = i
        cell_idx(i, j, center, width, band)
    };

    let j_range = |i: usize| -> (usize, usize) {
        let lo = (i as i64 - band as i64).max(0) as usize;
        let hi = ((i as i64 + band as i64) as usize).min(m.saturating_sub(1));
        (lo, hi)
    };

    // ---- Forward pass (linear, per-row scaled) ----
    let (lo0, hi0) = j_range(0);
    let q0 = qual.get(q_to_orig.first().copied().unwrap_or(0)).copied().unwrap_or(20);
    let (pm0, pe0) = emit_probs(q0);
    let mut row_sum: f32 = 0.0;
    for j in lo0..=hi0 {
        let e = if read[0] == refseq[j] { pm0 } else { pe0 };
        if let Some(k) = in_band(0, j) {
            fm[k] = e;
            row_sum += e;
        }
    }
    if row_sum == 0.0 { return vec![1.0; n]; }
    let inv = 1.0 / row_sum;
    for k in 0..width { fm[k] *= inv; }
    scale[0] = row_sum;

    for i in 1..n {
        let (lo, hi) = j_range(i);
        let q = qual.get(q_to_orig.get(i).copied().unwrap_or(0)).copied().unwrap_or(20);
        let (pm, pe) = emit_probs(q);
        let mut row_sum: f32 = 0.0;
        for j in lo..=hi {
            let e_m = if read[i] == refseq[j] { pm } else { pe };

            // From (i-1, j-1) — M state new value.
            let mut m_acc: f32 = 0.0;
            if j > 0 {
                if let Some(k) = in_band(i - 1, j - 1) {
                    m_acc = (fm[k] * T_MM + fi[k] * T_IM + fd[k] * T_DM) * e_m;
                }
            }
            // From (i-1, j) — I state.
            let mut i_acc: f32 = 0.0;
            if let Some(k) = in_band(i - 1, j) {
                i_acc = fm[k] * T_MI + fi[k] * T_II;
            }
            // From (i, j-1) — D state.
            let mut d_acc: f32 = 0.0;
            if j > 0 {
                if let Some(k) = in_band(i, j - 1) {
                    d_acc = fm[k] * T_MD + fd[k] * T_DD;
                }
            }
            if let Some(k) = in_band(i, j) {
                fm[k] = m_acc;
                fi[k] = i_acc;
                fd[k] = d_acc;
                row_sum += m_acc + i_acc + d_acc;
            }
        }
        if row_sum == 0.0 {
            // Underflow — fall back to uniform posterior.
            return vec![1.0; n];
        }
        let inv = 1.0 / row_sum;
        let base = i * width;
        for k in 0..width {
            fm[base + k] *= inv;
            fi[base + k] *= inv;
            fd[base + k] *= inv;
        }
        scale[i] = row_sum;
    }

    // ---- Backward pass (linear, scaled by same per-row factors) ----
    let (lo_last, hi_last) = j_range(n - 1);
    for j in lo_last..=hi_last {
        if let Some(k) = in_band(n - 1, j) {
            bm[k] = 1.0;
            bi[k] = 1.0;
            bd[k] = 1.0;
        }
    }
    // Apply final row scale to maintain alpha * beta sum invariant.
    {
        let base = (n - 1) * width;
        let inv = 1.0 / scale[n - 1];
        for k in 0..width {
            bm[base + k] *= inv;
            bi[base + k] *= inv;
            bd[base + k] *= inv;
        }
    }
    for i in (0..n - 1).rev() {
        let (lo, hi) = j_range(i);
        let q_next = qual.get(q_to_orig.get(i + 1).copied().unwrap_or(0)).copied().unwrap_or(20);
        let (pm_n, pe_n) = emit_probs(q_next);
        for j in lo..=hi {
            let next_m_emit = if j + 1 < m {
                if read[i + 1] == refseq[j + 1] { pm_n } else { pe_n }
            } else { 0.0 };

            let bm_next = if j + 1 < m { in_band(i + 1, j + 1).map(|k| bm[k]).unwrap_or(0.0) } else { 0.0 };
            let bi_next = in_band(i + 1, j).map(|k| bi[k]).unwrap_or(0.0);
            let bd_next = in_band(i, j + 1).map(|k| bd[k]).unwrap_or(0.0);

            let m_back = T_MM * next_m_emit * bm_next + T_MI * bi_next + T_MD * bd_next;
            let i_back = T_IM * next_m_emit * bm_next + T_II * bi_next;
            let d_back = T_DM * next_m_emit * bm_next + T_DD * bd_next;
            if let Some(k) = in_band(i, j) {
                bm[k] = m_back;
                bi[k] = i_back;
                bd[k] = d_back;
            }
        }
        let inv = 1.0 / scale[i];
        let base = i * width;
        for k in 0..width {
            bm[base + k] *= inv;
            bi[base + k] *= inv;
            bd[base + k] *= inv;
        }
    }

    // ---- Posterior P(M | read) per read position ----
    // With linear scaling, gamma_M(i) = sum_j fm[i,j] * bm[i,j] * scale[i]
    // (the per-row scale unwinds the normalization for the joint product).
    let mut post = vec![1.0f64; n];
    for i in 0..n {
        let base = i * width;
        let mut acc: f32 = 0.0;
        for k in 0..width {
            acc += fm[base + k] * bm[base + k];
        }
        post[i] = (acc * scale[i]).clamp(0.0, 1.0) as f64;
    }
    post
}

/// Capping fallback: cap quality in ±5 window around indels.
pub fn apply_baq_capping(seq: &[u8], qual: &mut [u8], cigar: &[(Kind, u32)], cap: u8) {
    let _ = seq;
    if qual.is_empty() { return; }
    let mut q_pos: usize = 0;
    let mut indel_positions: Vec<usize> = Vec::new();
    for &(k, l) in cigar {
        match k {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => q_pos += l as usize,
            Kind::Insertion => { indel_positions.push(q_pos); q_pos += l as usize; }
            Kind::Deletion | Kind::Skip => { indel_positions.push(q_pos); }
            Kind::SoftClip => q_pos += l as usize,
            _ => {}
        }
    }
    let window: usize = BAQ_EXTEND_LEN.min(5);
    for ip in &indel_positions {
        let lo = ip.saturating_sub(window).min(qual.len());
        let hi = (ip + window).min(qual.len());
        if hi <= lo { continue; }
        for q in &mut qual[lo..hi] {
            if *q > cap { *q = cap; }
        }
    }
}

/// Per-platform NM down-weighting profile. Rate-based (mismatch fraction, not absolute count) so a
/// single profile auto-scales to any read length.
#[derive(Clone, Copy)]
pub struct NmProfile {
    /// mismatch rate (NM / aligned length) tolerated at full weight
    pub full_rate: f32,
    /// weight falloff per unit of excess mismatch rate
    pub slope: f32,
    /// minimum weight floor
    pub floor: f32,
}

static NM_PROFILE: std::sync::OnceLock<Option<NmProfile>> = std::sync::OnceLock::new();

/// Select the NM-weighting profile ONCE, from the observed median read length (a robust platform
/// proxy) with an optional `KIRA_NM_WEIGHT` override. Must be called before the BAQ pass.
///   - unset / "auto": read-length preset. Short reads (<=320bp, Illumina) use the validated 2%-rate
///     profile; long reads (ONT/PacBio) disable it — their NM is dominated by sequencing error, so it
///     no longer distinguishes paralogs. (HiFi/long-but-accurate users can force a profile via the env.)
///   - "off" | "0": disabled.
///   - "F,S": manual absolute "full,slope" (e.g. "3,0.12"), converted to a rate via `median_len`
///     (back-compat with the hand-tuned value).
pub fn init_nm_profile(median_len: usize) {
    NM_PROFILE.get_or_init(|| {
        match std::env::var("KIRA_NM_WEIGHT").ok().as_deref() {
            Some("off") | Some("0") => None,
            Some(s) if s != "auto" && s.contains(',') => {
                let mut it = s.split(',');
                let full: f32 = it.next().and_then(|x| x.trim().parse().ok()).unwrap_or(3.0);
                let slope: f32 = it.next().and_then(|x| x.trim().parse().ok()).unwrap_or(0.12);
                let rl = median_len.max(1) as f32;
                Some(NmProfile { full_rate: full / rl, slope: slope * rl, floor: 0.12 })
            }
            _ => {
                if median_len == 0 || median_len > 320 {
                    None // long-read platform: NM ≈ error rate, not a paralog signal
                } else {
                    // Illumina short-read — validated on 150bp PE as "3,0.12" (= 2% rate, slope 18).
                    Some(NmProfile { full_rate: 0.02, slope: 18.0, floor: 0.12 })
                }
            }
        }
    });
}

/// NM-aware quality down-weighting (rate-based). Reads whose mismatch rate vs the local reference
/// exceeds the active profile threshold are likely paralog mismaps; scaling their base qualities down
/// suppresses their spurious variant evidence so the caller can use a relaxed QUAL cutoff without
/// admitting paralog FP. No-op until `init_nm_profile` selects a profile (None = disabled).
pub fn nm_weight_qual(seq: &[u8], qual: &mut [u8], cigar: &[(Kind, u32)], ref_seq: &[u8]) {
    let Some(prof) = NM_PROFILE.get().copied().flatten() else { return };
    if qual.is_empty() || ref_seq.is_empty() {
        return;
    }
    let mut q_pos = 0usize;
    let mut r_pos = 0usize;
    let mut mism = 0u32;
    let mut alen = 0u32;
    for &(k, l) in cigar {
        let l = l as usize;
        match k {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for i in 0..l {
                    if q_pos + i < seq.len() && r_pos + i < ref_seq.len() {
                        alen += 1;
                        if seq[q_pos + i].to_ascii_uppercase() != ref_seq[r_pos + i].to_ascii_uppercase() {
                            mism += 1;
                        }
                    }
                }
                q_pos += l;
                r_pos += l;
            }
            Kind::Insertion | Kind::SoftClip => q_pos += l,
            Kind::Deletion | Kind::Skip => r_pos += l,
            _ => {}
        }
    }
    if alen == 0 {
        return;
    }
    let rate = mism as f32 / alen as f32;
    if rate <= prof.full_rate {
        return;
    }
    let w = (1.0 - (rate - prof.full_rate) * prof.slope).max(prof.floor);
    for q in qual.iter_mut() {
        *q = ((*q as f32) * w) as u8;
    }
}

#[cfg(test)]
#[path = "../../tests/unit/bam_baq.rs"]
mod tests;
