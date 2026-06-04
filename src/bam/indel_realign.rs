//! Indel realignment. Per indel candidate, score reads against REF
//! and ALT (indel-applied) haplotypes via banded Smith-Waterman;
//! emit refined per-allele likelihoods.
//!
//! Inspired by bcftools/bam2bcf_indel.c; simplified — uses fixed
//! gap/mismatch penalties and a +-50bp window around the indel.

use crate::bam::pileup::LiveRead;

const WINDOW: u32 = 50;

const GAP_OPEN: i32 = -6;
const GAP_EXT: i32 = -1;
const MATCH: i32 = 2;
const MISMATCH: i32 = -4;

#[derive(Debug, Clone)]
pub struct IndelCandidate {
    pub chr: String,
    pub pos: u32,
    pub ref_base: u8,
    pub kind: IndelKind,
    pub seq: Vec<u8>,
    pub length: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndelKind { Insertion, Deletion }

#[derive(Debug, Clone)]
pub struct RealignScore {
    pub allele_idx: usize,
    pub log_prob: f64,
}

pub fn realign_reads_at_indel(
    candidate: &IndelCandidate,
    reads: &[&LiveRead],
    ref_window: &[u8],
    window_start: u32,
) -> Vec<RealignScore> {
    if reads.is_empty() || ref_window.is_empty() { return Vec::new(); }
    let local_indel_pos = (candidate.pos - window_start) as usize;
    let alt_ref = match candidate.kind {
        IndelKind::Insertion => apply_insertion(ref_window, local_indel_pos + 1, &candidate.seq),
        IndelKind::Deletion => apply_deletion(ref_window, local_indel_pos + 1, candidate.length as usize),
    };

    let mut scores = Vec::with_capacity(reads.len() * 2);
    for (i, lr) in reads.iter().enumerate() {
        let read_seq = read_window_seq(lr, candidate.pos, WINDOW);
        if read_seq.is_empty() { continue; }
        let s_ref = align_score(&read_seq, ref_window);
        let s_alt = align_score(&read_seq, &alt_ref);
        let log_p_ref = score_to_logp(s_ref, read_seq.len());
        let log_p_alt = score_to_logp(s_alt, read_seq.len());
        scores.push(RealignScore { allele_idx: i * 2, log_prob: log_p_ref });
        scores.push(RealignScore { allele_idx: i * 2 + 1, log_prob: log_p_alt });
    }
    scores
}

fn apply_insertion(refseq: &[u8], pos: usize, ins: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(refseq.len() + ins.len());
    out.extend_from_slice(&refseq[..pos.min(refseq.len())]);
    out.extend_from_slice(ins);
    if pos < refseq.len() { out.extend_from_slice(&refseq[pos..]); }
    out
}

fn apply_deletion(refseq: &[u8], pos: usize, len: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(refseq.len().saturating_sub(len));
    out.extend_from_slice(&refseq[..pos.min(refseq.len())]);
    let after = (pos + len).min(refseq.len());
    if after < refseq.len() { out.extend_from_slice(&refseq[after..]); }
    out
}

fn read_window_seq(lr: &LiveRead, anchor_pos: u32, window: u32) -> Vec<u8> {
    use noodles_sam::alignment::record::cigar::op::Kind;
    let win_start = anchor_pos.saturating_sub(window);
    let win_end = anchor_pos + window;
    let mut r_pos = lr.ref_start;
    let mut q_pos: u32 = 0;
    let mut out: Vec<u8> = Vec::with_capacity(2 * window as usize + 1);
    for &(kind, len) in &lr.cigar_pairs {
        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for i in 0..len {
                    let r = r_pos + i;
                    if r >= win_start && r < win_end {
                        let qi = (q_pos + i) as usize;
                        if qi < lr.seq.len() { out.push(lr.seq[qi].to_ascii_uppercase()); }
                    }
                }
                r_pos += len; q_pos += len;
            }
            Kind::Insertion => {
                if r_pos >= win_start && r_pos < win_end {
                    for i in 0..len {
                        let qi = (q_pos + i) as usize;
                        if qi < lr.seq.len() { out.push(lr.seq[qi].to_ascii_uppercase()); }
                    }
                }
                q_pos += len;
            }
            Kind::Deletion | Kind::Skip => { r_pos += len; }
            Kind::SoftClip => q_pos += len,
            _ => {}
        }
        if r_pos > win_end { break; }
    }
    out
}

/// Banded Smith-Waterman score. Returns max alignment score.
fn align_score(read: &[u8], refseq: &[u8]) -> i32 {
    let n = read.len();
    let m = refseq.len();
    if n == 0 || m == 0 { return 0; }
    let band: i32 = ((n as i32 - m as i32).abs() + 16).min(64);
    let mut prev_h: Vec<i32> = vec![0; m + 1];
    let mut cur_h: Vec<i32> = vec![0; m + 1];
    let mut prev_e: Vec<i32> = vec![i32::MIN / 4; m + 1];
    let mut cur_e: Vec<i32> = vec![i32::MIN / 4; m + 1];
    let mut cur_f: Vec<i32> = vec![i32::MIN / 4; m + 1];
    let mut best: i32 = 0;

    for i in 1..=n {
        let lo = (i as i32 - band).max(1) as usize;
        let hi = ((i as i32 + band) as usize).min(m);
        for j in lo..=hi {
            let m_ij = if read[i - 1].eq_ignore_ascii_case(&refseq[j - 1]) { MATCH } else { MISMATCH };
            cur_e[j] = (cur_h[j - 1] + GAP_OPEN).max(prev_e[j - 1] + GAP_EXT);
            cur_f[j] = (prev_h[j] + GAP_OPEN).max(cur_f[j].max(prev_h[j] + GAP_EXT));
            let s_diag = prev_h[j - 1] + m_ij;
            cur_h[j] = s_diag.max(cur_e[j]).max(cur_f[j]).max(0);
            if cur_h[j] > best { best = cur_h[j]; }
        }
        std::mem::swap(&mut prev_h, &mut cur_h);
        std::mem::swap(&mut prev_e, &mut cur_e);
        for v in cur_h.iter_mut() { *v = 0; }
        for v in cur_e.iter_mut() { *v = i32::MIN / 4; }
        for v in cur_f.iter_mut() { *v = i32::MIN / 4; }
    }
    best
}

fn score_to_logp(score: i32, read_len: usize) -> f64 {
    let max_possible = (read_len as i32) * MATCH;
    if max_possible <= 0 { return -100.0; }
    let frac = score.max(0) as f64 / max_possible as f64;
    let p = frac.max(0.001).min(0.999);
    p.ln()
}

/// Aggregate realign-scores: returns per-allele likelihood ratio that
/// can be added to errmod-based PL.
pub fn aggregate_scores(scores: &[RealignScore], n_alleles: usize) -> Vec<f64> {
    let mut log_lik = vec![0.0f64; n_alleles];
    for s in scores {
        if s.allele_idx < n_alleles { log_lik[s.allele_idx] += s.log_prob; }
    }
    log_lik
}

#[cfg(test)]
#[path = "../../tests/unit/bam_indel_realign.rs"]
mod tests;
