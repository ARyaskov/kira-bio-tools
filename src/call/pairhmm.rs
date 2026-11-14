//! Read-vs-haplotype log-likelihood on the shared glocal pair-HMM
//! (`crate::align`). The read aligns globally, the haplotype locally, so
//! `read_vs_hap_loglik` is log P(read | haplotype) up to a start prior that
//! cancels in likelihood ratios between haplotypes of similar length.

use crate::align::glocal::glocal_loglik_with;
use crate::align::{GlocalParams, GlocalScratch, encode_nt};

#[derive(Default)]
struct PairScratch {
    reference: Vec<u8>,
    query: Vec<u8>,
    qualities: Vec<u8>,
    hmm: GlocalScratch,
}

pub(crate) struct HaplotypeScorer {
    haplotypes: Vec<Vec<u8>>,
    scratch: PairScratch,
}

impl HaplotypeScorer {
    pub(crate) fn new<T: AsRef<[u8]>>(haplotypes: &[T]) -> Self {
        Self {
            haplotypes: haplotypes
                .iter()
                .map(|h| h.as_ref().iter().map(|&b| encode_nt(b)).collect())
                .collect(),
            scratch: PairScratch::default(),
        }
    }

    pub(crate) fn score(&mut self, read: &[u8], quals: &[u8]) -> Vec<f64> {
        let PairScratch {
            query,
            qualities,
            hmm,
            ..
        } = &mut self.scratch;
        query.clear();
        query.extend(read.iter().map(|&b| encode_nt(b)));
        qualities.clear();
        qualities.extend((0..read.len()).map(|i| quals.get(i).copied().unwrap_or(30).max(1)));
        self.haplotypes
            .iter()
            .map(|h| {
                glocal_loglik_with(hmm, h, query, Some(qualities), &params(read.len(), h.len()))
            })
            .collect()
    }
}

fn params(read_len: usize, hap_len: usize) -> GlocalParams {
    // Band wide enough for the read to sit anywhere in the haplotype.
    GlocalParams {
        d: 0.001,
        e: 0.1,
        bw: hap_len.abs_diff(read_len) + 10,
    }
}

/// log P(read | haplotype). `quals` are phred base qualities (30 when short).
pub fn read_vs_hap_loglik(read: &[u8], quals: &[u8], hap: &[u8]) -> f64 {
    if read.is_empty() || hap.is_empty() {
        return f64::NEG_INFINITY;
    }
    thread_local! {
        static SCRATCH: std::cell::RefCell<PairScratch> = std::cell::RefCell::new(PairScratch::default());
    }
    SCRATCH.with(|scratch| {
        let mut scratch = scratch.borrow_mut();
        let PairScratch {
            reference,
            query,
            qualities,
            hmm,
        } = &mut *scratch;
        reference.clear();
        reference.extend(hap.iter().map(|&b| encode_nt(b)));
        query.clear();
        query.extend(read.iter().map(|&b| encode_nt(b)));
        qualities.clear();
        qualities.extend((0..read.len()).map(|i| quals.get(i).copied().unwrap_or(30).max(1)));
        glocal_loglik_with(
            hmm,
            reference,
            query,
            Some(qualities),
            &params(read.len(), hap.len()),
        )
    })
}

/// Per-read log-likelihood ratio between two haplotypes: positive favours `alt`.
#[inline]
pub fn loglik_ratio(read: &[u8], quals: &[u8], hap_ref: &[u8], hap_alt: &[u8]) -> f64 {
    read_vs_hap_loglik(read, quals, hap_alt) - read_vs_hap_loglik(read, quals, hap_ref)
}

#[cfg(test)]
#[path = "../../tests/unit/call_pairhmm.rs"]
mod tests;
