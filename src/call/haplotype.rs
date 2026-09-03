//! Local indel recovery in active regions. Each read is realigned to the
//! reference window with the glocal pair-HMM; the MAP state path yields the
//! indel it implies plus the residual mismatch count, and indels supported by
//! enough clean-fitting reads are called. Genotype likelihoods come from
//! read-vs-haplotype likelihoods on the same kernel.

use crate::align::{GlocalParams, encode_nt, glocal};
use crate::bam::pileup::LiveRead;
use crate::call::pairhmm::read_vs_hap_loglik;

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Indel {
    pub is_ins: bool,
    /// Reference-window offset just past the anchor base (the indel sits here).
    pub ref_off: usize,
    /// Indel length in bases.
    pub len: usize,
    /// Inserted bases (empty for deletions).
    pub bases: Vec<u8>,
}

/// Align `read` to `refw` and return the single largest indel on the MAP
/// path plus the number of residual mismatches (the clean-fit signal).
pub fn discover_indel(read: &[u8], refw: &[u8]) -> Option<(Indel, u32)> {
    let n = read.len();
    let m = refw.len();
    if n < 8 || m < 8 {
        return None;
    }
    let r: Vec<u8> = refw.iter().map(|&b| encode_nt(b)).collect();
    let q: Vec<u8> = read.iter().map(|&b| encode_nt(b)).collect();
    let par = GlocalParams { d: 0.001, e: 0.1, bw: n.abs_diff(m) + 10 };
    let res = glocal(&r, &q, None, &par, true)?;
    let mut best: Option<Indel> = None;
    let mut mmc = 0u32;
    let mut prev_col: Option<usize> = None;
    let mut i = 0usize;
    let consider = |cand: Indel, best: &mut Option<Indel>| {
        if best.as_ref().is_none_or(|b| cand.len > b.len) {
            *best = Some(cand);
        }
    };
    while i < n {
        let st = res.state[i];
        if st < 0 {
            return None;
        }
        if st & 3 == 1 {
            let start = i;
            while i < n && res.state[i] & 3 == 1 {
                i += 1;
            }
            if let Some(pc) = prev_col {
                consider(Indel { is_ins: true, ref_off: pc + 1, len: i - start, bases: read[start..i].to_vec() }, &mut best);
            }
            continue;
        }
        let col = (st >> 2) as usize;
        if col >= m {
            return None;
        }
        if let Some(pc) = prev_col {
            if col <= pc {
                // The posterior path is not monotone: no clean alignment.
                return None;
            }
            if col > pc + 1 {
                consider(Indel { is_ins: false, ref_off: pc + 1, len: col - pc - 1, bases: Vec::new() }, &mut best);
            }
        }
        if !read[i].eq_ignore_ascii_case(&refw[col]) {
            mmc += 1;
        }
        prev_col = Some(col);
        i += 1;
    }
    // Indels flush against the window ends are likely artifacts.
    match best {
        Some(b) if b.ref_off > 0 && b.ref_off < m => Some((b, mmc)),
        _ => None,
    }
}

/// Apply `ind` to the reference window, returning the alt haplotype.
pub fn apply_indel(refw: &[u8], ind: &Indel) -> Vec<u8> {
    let mut h = Vec::with_capacity(refw.len() + ind.len);
    h.extend_from_slice(&refw[..ind.ref_off]);
    if ind.is_ins {
        h.extend_from_slice(&ind.bases);
        h.extend_from_slice(&refw[ind.ref_off..]);
    } else {
        let end = (ind.ref_off + ind.len).min(refw.len());
        h.extend_from_slice(&refw[end..]);
    }
    h
}

pub struct AssembledCall {
    pub pos1: u64,
    pub ref_str: String,
    pub alt_str: String,
    pub support: u32,
    pub total: u32,
    /// Window and the two haplotypes the reads are scored against.
    pub win_lo: u32,
    pub win_hi: u32,
    pub hap_ref: Vec<u8>,
    pub hap_alt: Vec<u8>,
}

/// Recover the best clean indel near a window. Each read votes for the indel its
/// realignment implies, but only if it fits cleanly (`<= max_mm` residual
/// mismatches over the window). The indel with the most clean votes wins if it
/// clears `min_support`. `win_lo` is the 0-based genomic start of `refw`.
pub fn assemble_indel(reads: &[LiveRead], win_lo: u32, refw: &[u8], min_support: u32, max_mm: u32) -> Option<AssembledCall> {
    let hi = win_lo + refw.len() as u32;
    let mut votes: std::collections::HashMap<Indel, u32> = std::collections::HashMap::new();
    let mut total = 0u32;
    for lr in reads {
        let Some(sub) = lr.query_window(win_lo, hi) else { continue };
        if sub.len() < 12 {
            continue;
        }
        total += 1;
        if let Some((ind, mm)) = discover_indel(&sub, refw) {
            if mm <= max_mm {
                *votes.entry(ind).or_insert(0) += 1;
            }
        }
    }
    if total == 0 || votes.is_empty() {
        return None;
    }
    let (ind, support) = votes.into_iter().max_by_key(|(_, c)| *c)?;
    if support < min_support || ind.ref_off == 0 {
        return None;
    }
    let anchor_idx = ind.ref_off - 1;
    let anchor_base = refw[anchor_idx] as char;
    let pos1 = (win_lo as u64) + anchor_idx as u64 + 1;
    let (ref_str, alt_str) = if ind.is_ins {
        let mut a = String::with_capacity(1 + ind.len);
        a.push(anchor_base);
        for &b in &ind.bases {
            a.push(b as char);
        }
        (anchor_base.to_string(), a)
    } else {
        let end = (ind.ref_off + ind.len).min(refw.len());
        let mut r = String::with_capacity(1 + ind.len);
        r.push(anchor_base);
        for &b in &refw[ind.ref_off..end] {
            r.push(b as char);
        }
        (r, anchor_base.to_string())
    };
    let hap_alt = apply_indel(refw, &ind);
    Some(AssembledCall { pos1, ref_str, alt_str, support, total, win_lo, win_hi: hi, hap_ref: refw.to_vec(), hap_alt })
}

/// Per-sample evidence for a two-haplotype site.
#[derive(Clone, Copy, Default, Debug)]
pub struct HapSample {
    /// PL for 0/0, 0/1, 1/1 (capped at 255).
    pub pl: [u32; 3],
    pub n_ref: u32,
    pub n_alt: u32,
}

/// Genotype likelihoods from read-vs-haplotype likelihoods: P(read | 0/1)
/// is the mean of the two haplotype likelihoods.
pub fn haplotype_pls(reads: &[LiveRead], n_samples: usize, call: &AssembledCall) -> Vec<HapSample> {
    let mut ll = vec![[0.0f64; 3]; n_samples];
    let mut out = vec![HapSample::default(); n_samples];
    let ln2 = std::f64::consts::LN_2;
    for lr in reads {
        let Some((bases, quals)) = lr.query_window_qual(call.win_lo, call.win_hi) else { continue };
        if bases.len() < 12 || lr.sample_idx >= n_samples {
            continue;
        }
        let lr_ref = read_vs_hap_loglik(&bases, &quals, &call.hap_ref);
        let lr_alt = read_vs_hap_loglik(&bases, &quals, &call.hap_alt);
        if !lr_ref.is_finite() || !lr_alt.is_finite() {
            continue;
        }
        let s = &mut ll[lr.sample_idx];
        s[0] += lr_ref;
        s[2] += lr_alt;
        let m = lr_ref.max(lr_alt);
        s[1] += m + ((lr_ref - m).exp() + (lr_alt - m).exp()).ln() - ln2;
        if lr_alt > lr_ref {
            out[lr.sample_idx].n_alt += 1;
        } else {
            out[lr.sample_idx].n_ref += 1;
        }
    }
    for (o, l) in out.iter_mut().zip(ll.iter()) {
        let max = l.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        for k in 0..3 {
            o.pl[k] = ((-4.343 * (l[k] - max)) + 0.499).clamp(0.0, 255.0) as u32;
        }
    }
    out
}

#[cfg(test)]
#[path = "../../tests/unit/call_haplotype.rs"]
mod tests;
