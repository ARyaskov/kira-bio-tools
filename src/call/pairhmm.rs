//! Pair-HMM read-vs-haplotype log-likelihood — the HaplotypeCaller scoring core.
//!
//! `read_vs_hap_loglik` returns log P(read | haplotype) under a 3-state (Match/
//! Insert/Delete) pair-HMM with quality-driven match emissions and a fixed affine
//! indel transition model. The read aligns locally within the haplotype (free
//! start/end on the haplotype). Genotyping uses the per-read likelihood ratio
//! between candidate haplotypes; that ratio cancels the start prior, so absolute
//! values are unitless.

const P_GAP_OPEN: f64 = 3.16e-5; // M->I and M->D (Phred ~45)
const P_GAP_CONT: f64 = 0.1; // I->I and D->D (Phred ~10)

#[inline]
fn lse2(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let m = a.max(b);
    m + ((a - m).exp() + (b - m).exp()).ln()
}

#[inline]
fn lse3(a: f64, b: f64, c: f64) -> f64 {
    lse2(lse2(a, b), c)
}

/// log P(read | haplotype). `quals` are Phred base qualities (one per read base).
pub fn read_vs_hap_loglik(read: &[u8], quals: &[u8], hap: &[u8]) -> f64 {
    let n = read.len();
    let m = hap.len();
    if n == 0 || m == 0 {
        return f64::NEG_INFINITY;
    }
    let l_mm = (1.0 - P_GAP_OPEN - P_GAP_OPEN).ln();
    let l_mi = P_GAP_OPEN.ln();
    let l_md = P_GAP_OPEN.ln();
    let l_im = (1.0 - P_GAP_CONT).ln();
    let l_ii = P_GAP_CONT.ln();
    let l_dm = (1.0 - P_GAP_CONT).ln();
    let l_dd = P_GAP_CONT.ln();

    let neg = f64::NEG_INFINITY;
    // Row 0: M[0][j] = ln(1) = 0 (read may begin matching at any haplotype base).
    let mut pm = vec![0.0f64; m + 1];
    let mut pi = vec![neg; m + 1];
    let mut pd = vec![neg; m + 1];
    let mut cm = vec![neg; m + 1];
    let mut ci = vec![neg; m + 1];
    let mut cd = vec![neg; m + 1];

    for i in 1..=n {
        let q = quals.get(i - 1).copied().unwrap_or(30).max(1);
        let eps = 10f64.powf(-(q as f64) / 10.0);
        let l_match = (1.0 - eps).ln();
        let l_mis = (eps / 3.0).ln();
        let rb = read[i - 1];
        cm[0] = neg;
        ci[0] = lse2(pm[0] + l_mi, pi[0] + l_ii);
        cd[0] = neg;
        for j in 1..=m {
            let emit = if rb.eq_ignore_ascii_case(&hap[j - 1]) { l_match } else { l_mis };
            cm[j] = emit + lse3(pm[j - 1] + l_mm, pi[j - 1] + l_im, pd[j - 1] + l_dm);
            ci[j] = lse2(pm[j] + l_mi, pi[j] + l_ii);
            cd[j] = lse2(cm[j - 1] + l_md, cd[j - 1] + l_dd);
        }
        std::mem::swap(&mut pm, &mut cm);
        std::mem::swap(&mut pi, &mut ci);
        std::mem::swap(&mut pd, &mut cd);
    }

    // Read fully consumed; it may end (match or insert) at any haplotype position.
    let mut acc = neg;
    for j in 0..=m {
        acc = lse3(acc, pm[j], pi[j]);
    }
    acc
}

/// Per-read log-likelihood ratio between two haplotypes: positive favours `alt`.
#[inline]
pub fn loglik_ratio(read: &[u8], quals: &[u8], hap_ref: &[u8], hap_alt: &[u8]) -> f64 {
    read_vs_hap_loglik(read, quals, hap_alt) - read_vs_hap_loglik(read, quals, hap_ref)
}

#[cfg(test)]
#[path = "../../tests/unit/call_pairhmm.rs"]
mod tests;
