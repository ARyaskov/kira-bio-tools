//! Local indel recovery: realign each read to the reference window (affine), read
//! off the indel it implies plus the residual mismatch count, and accept indels
//! supported by enough reads that fit *cleanly* (few residual mismatches) over a
//! wide window. A spurious indel betrays itself with downstream mismatches the
//! true indel does not have; the wide window + mismatch cap is what gives this
//! ~2% false-call rate where a bare edit-distance scan floods. Recovers indels
//! the aligner modelled as mismatches and so never placed in a CIGAR.

use crate::bam::pileup::LiveRead;

const NEG: f64 = f64::NEG_INFINITY;
// affine alignment scoring for discovery
const MA: f64 = 1.0; // match
const MI: f64 = -4.0; // mismatch
const GO: f64 = -6.0; // gap open
const GE: f64 = -1.0; // gap extend

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Indel {
    pub is_ins: bool,
    pub ref_off: usize,   // ref-window offset just past the anchor base (indel sits here)
    pub len: usize,       // indel length in bases
    pub bases: Vec<u8>,   // inserted bases (empty for deletions)
}

/// Affine align `read` to `refw`; return the single largest indel plus the number
/// of residual mismatches in the alignment (the clean-fit signal).
pub fn discover_indel(read: &[u8], refw: &[u8]) -> Option<(Indel, u32)> {
    let n = read.len();
    let m = refw.len();
    if n < 8 || m < 8 { return None; }
    let mut mm = vec![vec![NEG; m + 1]; n + 1];
    let mut ix = vec![vec![NEG; m + 1]; n + 1]; // gap in ref (insertion in read)
    let mut iy = vec![vec![NEG; m + 1]; n + 1]; // gap in read (deletion)
    let mut bm = vec![vec![0u8; m + 1]; n + 1];
    let mut bx = vec![vec![0u8; m + 1]; n + 1];
    let mut by = vec![vec![0u8; m + 1]; n + 1];
    mm[0][0] = 0.0;
    for i in 1..=n { ix[i][0] = GO + (i as f64 - 1.0) * GE; }
    for j in 1..=m { iy[0][j] = GO + (j as f64 - 1.0) * GE; }
    for i in 1..=n {
        for j in 1..=m {
            let sc = if read[i - 1].eq_ignore_ascii_case(&refw[j - 1]) { MA } else { MI };
            let mut best = mm[i - 1][j - 1]; let mut s = 0u8;
            if ix[i - 1][j - 1] > best { best = ix[i - 1][j - 1]; s = 1; }
            if iy[i - 1][j - 1] > best { best = iy[i - 1][j - 1]; s = 2; }
            mm[i][j] = best + sc; bm[i][j] = s;
            let (o, e) = (mm[i - 1][j] + GO, ix[i - 1][j] + GE);
            if o >= e { ix[i][j] = o; bx[i][j] = 0; } else { ix[i][j] = e; bx[i][j] = 1; }
            let (o, e) = (mm[i][j - 1] + GO, iy[i][j - 1] + GE);
            if o >= e { iy[i][j] = o; by[i][j] = 0; } else { iy[i][j] = e; by[i][j] = 2; }
        }
    }
    let (mut i, mut j) = (n, m);
    let mut cur = {
        let mut c = 0u8; let mut v = mm[n][m];
        if ix[n][m] > v { v = ix[n][m]; c = 1; }
        if iy[n][m] > v { c = 2; }
        c
    };
    let mut ops: Vec<(u8, usize)> = Vec::new(); // (0=M,1=I,2=D)
    let push = |op: u8, ops: &mut Vec<(u8, usize)>| {
        if let Some(last) = ops.last_mut() {
            if last.0 == op { last.1 += 1; return; }
        }
        ops.push((op, 1));
    };
    while i > 0 || j > 0 {
        if i == 0 { push(2, &mut ops); j -= 1; continue; }
        if j == 0 { push(1, &mut ops); i -= 1; continue; }
        match cur {
            0 => { push(0, &mut ops); let s = bm[i][j]; i -= 1; j -= 1; cur = s; }
            1 => { push(1, &mut ops); let s = bx[i][j]; i -= 1; cur = s; }
            _ => { push(2, &mut ops); let s = by[i][j]; j -= 1; cur = s; }
        }
    }
    ops.reverse();
    // walk forward: biggest indel + its offset/bases, and count residual mismatches
    let (mut ref_pos, mut read_pos) = (0usize, 0usize);
    let mut best: Option<Indel> = None;
    let mut mmc = 0u32;
    for &(op, len) in &ops {
        match op {
            0 => {
                for t in 0..len {
                    if !read[read_pos + t].eq_ignore_ascii_case(&refw[ref_pos + t]) { mmc += 1; }
                }
                ref_pos += len; read_pos += len;
            }
            1 => {
                if best.as_ref().map_or(true, |b| len > b.len) {
                    best = Some(Indel { is_ins: true, ref_off: ref_pos, len,
                                        bases: read[read_pos..read_pos + len].to_vec() });
                }
                read_pos += len;
            }
            _ => {
                if best.as_ref().map_or(true, |b| len > b.len) {
                    best = Some(Indel { is_ins: false, ref_off: ref_pos, len, bases: Vec::new() });
                }
                ref_pos += len;
            }
        }
    }
    // ignore indels flush against the window ends (likely artifacts)
    match best {
        Some(b) if b.ref_off > 0 && b.ref_off < m => Some((b, mmc)),
        _ => None,
    }
}

/// Apply `ind` to the reference window, returning the alt haplotype.
fn apply_indel(refw: &[u8], ind: &Indel) -> Vec<u8> {
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
}

/// Recover the best clean indel near a window. Each read votes for the indel its
/// realignment implies, but only if it fits cleanly (`<= max_mm` residual
/// mismatches over the window). The indel with the most clean votes wins if it
/// clears `min_support`. `win_lo` is the 0-based genomic start of `refw`.
/// Validated operating point: ~50bp window, max_mm=3, min_support=3.
pub fn assemble_indel(
    reads: &[&LiveRead],
    win_lo: u32,
    refw: &[u8],
    min_support: u32,
    max_mm: u32,
) -> Option<AssembledCall> {
    let hi = win_lo + refw.len() as u32;
    let mut votes: std::collections::HashMap<Indel, u32> = std::collections::HashMap::new();
    let mut total = 0u32;
    for lr in reads {
        let Some(sub) = lr.query_window(win_lo, hi) else { continue };
        if sub.len() < 12 { continue; }
        total += 1;
        if let Some((ind, mm)) = discover_indel(&sub, refw) {
            if mm <= max_mm {
                *votes.entry(ind).or_insert(0) += 1;
            }
        }
    }
    if total == 0 || votes.is_empty() { return None; }
    let (ind, support) = votes.into_iter().max_by_key(|(_, c)| *c)?;
    if support < min_support { return None; }
    // build the VCF record (anchor = base just before the indel)
    if ind.ref_off == 0 { return None; }
    let anchor_idx = ind.ref_off - 1;
    let anchor_base = refw[anchor_idx] as char;
    let pos1 = (win_lo as u64) + anchor_idx as u64 + 1;
    let (ref_str, alt_str) = if ind.is_ins {
        let mut a = String::with_capacity(1 + ind.len);
        a.push(anchor_base);
        for &b in &ind.bases { a.push(b as char); }
        (anchor_base.to_string(), a)
    } else {
        let end = (ind.ref_off + ind.len).min(refw.len());
        let mut r = String::with_capacity(1 + ind.len);
        r.push(anchor_base);
        for &b in &refw[ind.ref_off..end] { r.push(b as char); }
        (r, anchor_base.to_string())
    };
    Some(AssembledCall { pos1, ref_str, alt_str, support, total })
}

#[cfg(test)]
#[path = "../../tests/unit/call_haplotype.rs"]
mod tests;
