//! Base Alignment Quality: port of htslib `sam_prob_realn` with
//! `BAQ_APPLY | BAQ_EXTEND` (the `bcftools mpileup` default) on the shared
//! glocal pair-HMM. Li H. (2011) Bioinformatics 27(8):1157.

use noodles_sam::alignment::record::cigar::op::Kind;

use crate::align::{GlocalParams, encode_nt, glocal};

/// Apply BAQ in place: each aligned base keeps min(qual, extended BAQ);
/// inserted and soft-clipped bases get 0, as in htslib. `ref_window` holds
/// the reference for `[win_start, win_start + len)` (0-based) and should
/// reach a read length past both ends of the alignment. Returns false when
/// nothing was applied (reference skip, no aligned bases, no reference).
pub fn apply_baq_hmm(seq: &[u8], qual: &mut [u8], cigar: &[(Kind, u32)], ref_window: &[u8], win_start: u32, read_pos: u32) -> bool {
    let l_qseq = seq.len();
    if l_qseq == 0 || qual.len() != l_qseq || ref_window.is_empty() {
        return false;
    }
    // Alignment extent: x over the reference, y over the read.
    let (mut x, mut y) = (read_pos as i64, 0i64);
    let (mut yb, mut ye, mut xb, mut xe) = (-1i64, -1i64, -1i64, -1i64);
    for &(op, l) in cigar {
        let l = l as i64;
        match op {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                if yb < 0 {
                    yb = y;
                }
                if xb < 0 {
                    xb = x;
                }
                ye = y + l;
                xe = x + l;
                x += l;
                y += l;
            }
            Kind::SoftClip | Kind::Insertion => y += l,
            Kind::Deletion => x += l,
            Kind::Skip => return false,
            _ => {}
        }
    }
    if xb < 0 {
        return false;
    }
    let mut bw: i64 = 7;
    let skew = ((xe - xb) - (ye - yb)).abs();
    if skew > bw {
        bw = skew + 3;
    }
    xb -= yb + bw / 2;
    if xb < 0 {
        xb = 0;
    }
    xe += l_qseq as i64 - ye + bw / 2;
    if xe - xb - l_qseq as i64 > bw {
        let adj = (xe - xb - l_qseq as i64 - bw) / 2;
        xb += adj;
        xe -= adj;
    }
    // Reference codes for [xb, xe), truncated where the window ends.
    let ws = win_start as i64;
    let mut r: Vec<u8> = Vec::with_capacity((xe - xb).max(0) as usize);
    let mut i = xb;
    while i < xe {
        let off = i - ws;
        if off < 0 || off as usize >= ref_window.len() || ref_window[off as usize] == 0 {
            break;
        }
        r.push(encode_nt(ref_window[off as usize]));
        i += 1;
    }
    if r.is_empty() {
        return false;
    }
    let s: Vec<u8> = seq.iter().map(|&b| encode_nt(b)).collect();
    let conf = GlocalParams { d: 0.001, e: 0.1, bw: bw as usize };
    let Some(res) = glocal(&r, &s, Some(qual), &conf, true) else { return false };

    // Extended BAQ: within a match run, a base's BAQ is the smaller of the
    // running maxima from the left and from the right.
    let mut bq = vec![0u8; l_qseq];
    let mut left = vec![0u8; l_qseq];
    let mut rght = vec![0u8; l_qseq];
    let (mut x, mut y) = (read_pos as i64, 0usize);
    for &(op, l) in cigar {
        let l = l as usize;
        match op {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                if l == 0 {
                    continue;
                }
                for i in y..y + l {
                    let st = res.state[i];
                    let col = x - xb + (i - y) as i64;
                    bq[i] = if (st & 3) != 0 || (st >> 2) as i64 != col { 0 } else { res.q[i] };
                }
                left[y] = bq[y];
                for i in y + 1..y + l {
                    left[i] = bq[i].max(left[i - 1]);
                }
                rght[y + l - 1] = bq[y + l - 1];
                for i in (y..y + l - 1).rev() {
                    rght[i] = bq[i].max(rght[i + 1]);
                }
                for i in y..y + l {
                    bq[i] = left[i].min(rght[i]);
                }
                x += l as i64;
                y += l;
            }
            Kind::SoftClip | Kind::Insertion => y += l,
            Kind::Deletion => x += l as i64,
            _ => {}
        }
    }
    for (q, &b) in qual.iter_mut().zip(bq.iter()) {
        if *q > b {
            *q = b;
        }
    }
    true
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

/// Select the NM-weighting profile once (`mpileup --nm-weight`), from the
/// observed median read length:
///   - "off" (default): disabled, bcftools behaviour;
///   - "auto": short reads (<=320 bp) use the validated 2%-rate profile,
///     long reads disable it (their NM is sequencing error, not a paralog signal);
///   - "F,S": absolute full/slope (e.g. "3,0.12"), converted to a rate via `median_len`.
pub fn init_nm_profile(spec: &str, median_len: usize) {
    NM_PROFILE.get_or_init(|| match spec.trim() {
        "" | "off" | "0" => None,
        s if s != "auto" && s.contains(',') => {
            let mut it = s.split(',');
            let full: f32 = it.next().and_then(|x| x.trim().parse().ok()).unwrap_or(3.0);
            let slope: f32 = it.next().and_then(|x| x.trim().parse().ok()).unwrap_or(0.12);
            let rl = median_len.max(1) as f32;
            Some(NmProfile { full_rate: full / rl, slope: slope * rl, floor: 0.12 })
        }
        _ => {
            if median_len == 0 || median_len > 320 {
                None
            } else {
                Some(NmProfile { full_rate: 0.02, slope: 18.0, floor: 0.12 })
            }
        }
    });
}

/// NM-aware quality down-weighting (rate-based). Reads whose mismatch rate vs the local reference
/// exceeds the active profile threshold are likely paralog mismaps; scaling their base qualities down
/// suppresses their spurious variant evidence. No-op until `init_nm_profile` selects a profile.
/// `ref_seq` must start at the read's first reference base.
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
