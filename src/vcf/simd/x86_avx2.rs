#![allow(unsafe_op_in_unsafe_fn)]

use crate::util::chr_name_to_id;
use crate::vcf::structs::{VcfFields, VcfFieldsFull};

#[cfg(target_arch = "x86_64")]
use core::arch::x86_64::*;

const MAX_TABS_FULL: usize = 138;
const MAX_TABS_MIN: usize = 8;

#[target_feature(enable = "avx2")]
pub unsafe fn parse_vcf_line_simd(line: &[u8]) -> Option<VcfFieldsFull<'_>> {
    if line.len() < 10 {
        return None;
    }

    let mut tabs = [0usize; MAX_TABS_FULL];
    let found = find_tabs_avx2(line, &mut tabs);

    // VCF requires at least 8 columns (INFO included)
    if found < 7 {
        return None;
    }

    let line_str = std::str::from_utf8_unchecked(line);

    // Helper: safe slice
    #[inline(always)]
    fn slice<'a>(s: &'a str, start: usize, end: usize) -> &'a str {
        if start <= end && end <= s.len() {
            &s[start..end]
        } else {
            ""
        }
    }

    let chrom = slice(line_str, 0, tabs[0]);
    let pos = slice(line_str, tabs[0] + 1, tabs[1]);
    let id = slice(line_str, tabs[1] + 1, tabs[2]);
    let ref_allele = slice(line_str, tabs[2] + 1, tabs[3]);
    let alt = slice(line_str, tabs[3] + 1, tabs[4]);
    let qual = slice(line_str, tabs[4] + 1, tabs[5]);
    let filter = slice(line_str, tabs[5] + 1, tabs[6]);

    // INFO: may be empty
    let info_start = tabs[6] + 1;
    let info_end = if found >= 8 { tabs[7] } else { line_str.len() };
    let info = slice(line_str, info_start, info_end);

    // FORMAT + samples
    let (format, samples) = if found >= 8 {
        let fmt_start = tabs[7] + 1;
        let fmt_end = tabs.get(8).copied().unwrap_or(line_str.len());
        let format = slice(line_str, fmt_start, fmt_end);

        let mut sample_vec = Vec::new();
        for i in 8..found {
            let s_start = tabs[i] + 1;
            let s_end = tabs.get(i + 1).copied().unwrap_or(line_str.len());
            let s = slice(line_str, s_start, s_end);
            sample_vec.push(s);
        }

        (Some(format), sample_vec)
    } else {
        (None, Vec::new())
    };

    Some(VcfFieldsFull {
        chrom,
        pos,
        id,
        ref_allele,
        alt,
        qual,
        filter,
        info,
        format,
        samples,
    })
}

#[target_feature(enable = "avx2")]
pub unsafe fn parse_vcf_fields_avx2(line: &[u8]) -> Option<VcfFields<'_>> {
    if line.len() < 10 {
        return None;
    }

    let mut tabs = [0usize; MAX_TABS_MIN];
    let found = find_tabs_avx2(line, &mut tabs);

    if found < 7 {
        return None;
    }

    let line_str = std::str::from_utf8_unchecked(line);

    #[inline(always)]
    fn slice<'a>(s: &'a str, start: usize, end: usize) -> &'a str {
        if start <= end && end <= s.len() {
            &s[start..end]
        } else {
            ""
        }
    }

    let info_start = tabs[6] + 1;
    let info_end = if found >= 8 { tabs[7] } else { line_str.len() };

    Some(VcfFields {
        chrom: slice(line_str, 0, tabs[0]),
        pos: slice(line_str, tabs[0] + 1, tabs[1]),
        id: slice(line_str, tabs[1] + 1, tabs[2]),
        ref_allele: slice(line_str, tabs[2] + 1, tabs[3]),
        alt: slice(line_str, tabs[3] + 1, tabs[4]),
        qual: slice(line_str, tabs[4] + 1, tabs[5]),
        filter: slice(line_str, tabs[5] + 1, tabs[6]),
        info: slice(line_str, info_start, info_end),
    })
}

#[target_feature(enable = "avx2")]
pub unsafe fn parse_chr_pos_avx2(line: &[u8]) -> Option<(u8, u32)> {
    let mut tabs = [0usize; 2];
    let found = find_tabs_avx2(line, &mut tabs);

    if found < 1 {
        return None;
    }

    let chrom_bytes = &line[..tabs[0]];
    let chrom_str = std::str::from_utf8_unchecked(chrom_bytes);
    let chr_id = chr_name_to_id(chrom_str)?;

    let pos_start = tabs[0] + 1;
    let pos_end = tabs.get(1).copied().unwrap_or(line.len());

    if pos_start >= pos_end {
        return None;
    }

    let pos = parse_u32_fast(line.as_ptr().add(pos_start), pos_end - pos_start);
    Some((chr_id, pos))
}

#[target_feature(enable = "avx2")]
unsafe fn find_tabs_avx2(buf: &[u8], out: &mut [usize]) -> usize {
    let max = out.len();
    let mut found = 0;
    let mut pos = 0;
    let tab = _mm256_set1_epi8(b'\t' as i8);

    while pos + 32 <= buf.len() && found < max {
        let chunk = _mm256_loadu_si256(buf.as_ptr().add(pos) as *const __m256i);
        let cmp = _mm256_cmpeq_epi8(chunk, tab);
        let mask = _mm256_movemask_epi8(cmp) as u32;

        let mut m = mask;
        while m != 0 && found < max {
            let tz = m.trailing_zeros() as usize;
            out[found] = pos + tz;
            found += 1;
            m &= m - 1;
        }

        pos += 32;
    }

    while pos < buf.len() && found < max {
        if buf[pos] == b'\t' {
            out[found] = pos;
            found += 1;
        }
        pos += 1;
    }

    found
}

#[target_feature(enable = "avx2")]
pub unsafe fn parse_u32_fast(ptr: *const u8, len: usize) -> u32 {
    let mut x = 0u32;
    let mut i = 0;

    while i < len && i < 10 {
        let c = *ptr.add(i);
        if c < b'0' || c > b'9' {
            break;
        }
        x = x * 10 + (c - b'0') as u32;
        i += 1;
    }

    x
}

#[target_feature(enable = "avx2")]
pub unsafe fn find_char_avx2(buf: &[u8], ch: u8) -> Option<usize> {
    let target = _mm256_set1_epi8(ch as i8);
    let mut pos = 0;

    while pos + 32 <= buf.len() {
        let chunk = _mm256_loadu_si256(buf.as_ptr().add(pos) as *const __m256i);
        let cmp = _mm256_cmpeq_epi8(chunk, target);
        let mask = _mm256_movemask_epi8(cmp);

        if mask != 0 {
            return Some(pos + mask.trailing_zeros() as usize);
        }
        pos += 32;
    }

    while pos < buf.len() {
        if buf[pos] == ch {
            return Some(pos);
        }
        pos += 1;
    }

    None
}
