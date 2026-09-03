//! Binning arithmetic shared by CSI and TBI (`CSIv1.pdf`, htslib `hts.h`).

/// First bin id of `level`.
#[inline]
pub fn bin_first(level: u8) -> u32 {
    (((1u64 << (3 * level as u64)) - 1) / 7) as u32
}

/// Number of regular bins for `depth` levels (metadata pseudo-bin is `n_bins + 1`).
#[inline]
pub fn n_bins(depth: u8) -> u32 {
    bin_first(depth + 1)
}

#[inline]
pub fn metadata_bin(depth: u8) -> u32 {
    n_bins(depth) + 1
}

#[inline]
pub fn bin_parent(bin: u32) -> u32 {
    (bin.saturating_sub(1)) >> 3
}

/// Level of `bin` within a tree of `depth` levels.
pub fn bin_level(bin: u32, depth: u8) -> u8 {
    let mut l = depth;
    loop {
        if bin >= bin_first(l) {
            return l;
        }
        if l == 0 {
            return 0;
        }
        l -= 1;
    }
}

/// Leftmost leaf-window index covered by `bin`.
pub fn bin_bot(bin: u32, depth: u8) -> u64 {
    let l = bin_level(bin, depth);
    ((bin - bin_first(l)) as u64) << (3 * (depth - l) as u64)
}

/// Largest representable position (exclusive).
#[inline]
pub fn max_pos(min_shift: u8, depth: u8) -> u64 {
    1u64 << (min_shift as u64 + 3 * depth as u64)
}

/// Bin containing `[beg, end)` (0-based, half-open). htslib `hts_reg2bin`.
pub fn reg2bin(beg: u64, end: u64, min_shift: u8, depth: u8) -> u32 {
    let end = end.saturating_sub(1).max(beg);
    let mut s = min_shift as u64;
    let mut t: u64 = ((1u64 << (3 * depth as u64)) - 1) / 7;
    let mut l = depth;
    while l > 0 {
        if beg >> s == end >> s {
            return (t + (beg >> s)) as u32;
        }
        l -= 1;
        s += 3;
        t -= 1u64 << (3 * l as u64);
    }
    0
}

/// All bins overlapping `[beg, end)` (0-based, half-open). htslib `hts_reg2bins`.
pub fn reg2bins(beg: u64, end: u64, min_shift: u8, depth: u8, out: &mut Vec<u32>) {
    out.clear();
    if end <= beg {
        return;
    }
    let end = end - 1;
    let mut t: u64 = 0;
    let mut s = min_shift as u64 + 3 * depth as u64;
    for l in 0..=depth {
        let b = t + (beg >> s);
        let e = t + (end >> s);
        for i in b..=e {
            out.push(i as u32);
        }
        s = s.saturating_sub(3);
        t += 1u64 << (3 * l as u64);
    }
}

/// Number of levels for a CSI of `min_shift` covering contigs up to `max_len`
/// bases (htslib `bcf_index`).
pub fn depth_for(min_shift: u8, max_len: u64) -> u8 {
    let max_len = max_len.saturating_add(256);
    let mut n = 0u8;
    let mut s = 1u64 << min_shift;
    while max_len > s && n < 10 {
        n += 1;
        s <<= 3;
    }
    n
}

/// Levels used by tabix for a CSI of `min_shift` (htslib `tbx_index`).
pub fn tabix_depth_for(min_shift: u8) -> u8 {
    const TBX_MAX_SHIFT: u8 = 31;
    (TBX_MAX_SHIFT.saturating_sub(min_shift) + 2) / 3
}

#[cfg(test)]
#[path = "../../tests/unit/csi_utils.rs"]
mod tests;
