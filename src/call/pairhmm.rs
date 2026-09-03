//! Read-vs-haplotype log-likelihood on the shared glocal pair-HMM
//! (`crate::align`). The read aligns globally, the haplotype locally, so
//! `read_vs_hap_loglik` is log P(read | haplotype) up to a start prior that
//! cancels in likelihood ratios between haplotypes of similar length.

use crate::align::{GlocalParams, encode_nt, glocal_loglik};

fn params(read_len: usize, hap_len: usize) -> GlocalParams {
    // Band wide enough for the read to sit anywhere in the haplotype.
    GlocalParams { d: 0.001, e: 0.1, bw: hap_len.abs_diff(read_len) + 10 }
}

/// log P(read | haplotype). `quals` are phred base qualities (30 when short).
pub fn read_vs_hap_loglik(read: &[u8], quals: &[u8], hap: &[u8]) -> f64 {
    if read.is_empty() || hap.is_empty() {
        return f64::NEG_INFINITY;
    }
    let r: Vec<u8> = hap.iter().map(|&b| encode_nt(b)).collect();
    let q: Vec<u8> = read.iter().map(|&b| encode_nt(b)).collect();
    let iq: Vec<u8> = (0..read.len()).map(|i| quals.get(i).copied().unwrap_or(30).max(1)).collect();
    glocal_loglik(&r, &q, Some(&iq), &params(read.len(), hap.len()))
}

/// Per-read log-likelihood ratio between two haplotypes: positive favours `alt`.
#[inline]
pub fn loglik_ratio(read: &[u8], quals: &[u8], hap_ref: &[u8], hap_alt: &[u8]) -> f64 {
    read_vs_hap_loglik(read, quals, hap_alt) - read_vs_hap_loglik(read, quals, hap_ref)
}

#[cfg(test)]
#[path = "../../tests/unit/call_pairhmm.rs"]
mod tests;
