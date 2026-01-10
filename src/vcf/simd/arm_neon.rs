use crate::util::chr_name_to_id;
use crate::vcf::structs::{VcfFields, VcfFieldsFull};

#[cfg(target_arch = "aarch64")]
use core::arch::aarch64::*;

const MAX_TABS: usize = 138;

#[target_feature(enable = "neon")]
pub unsafe fn parse_vcf_line_simd(line: &[u8]) -> Option<VcfFieldsFull<'_>> {
    if line.len() < 20 {
        return None;
    }

    let mut tabs = [0usize; MAX_TABS];
    let found = find_tabs_neon(line, &mut tabs);

    if found < 7 {
        return None;
    }

    let line_str = std::str::from_utf8_unchecked(line);

    let chrom = &line_str[0..tabs[0]];
    let pos = &line_str[tabs[0] + 1..tabs[1]];
    let id = &line_str[tabs[1] + 1..tabs[2]];
    let ref_allele = &line_str[tabs[2] + 1..tabs[3]];
    let alt = &line_str[tabs[3] + 1..tabs[4]];
    let qual = &line_str[tabs[4] + 1..tabs[5]];
    let filter = &line_str[tabs[5] + 1..tabs[6]];

    let info_end = if found >= 8 { tabs[7] } else { line_len };
    let info_start = tabs[6] + 1;
    let info = if info_start < info_end {
        &line_str[info_start..info_end]
    } else {
        ""
    };

    let (format, samples) = if found >= 8 {
        let format_start = tabs[7] + 1;
        let format_end = tabs.get(8).copied().unwrap_or(line_len);
        let format_field = if format_start < format_end {
            &line_str[format_start..format_end]
        } else {
            ""
        };

        let mut sample_vec = Vec::with_capacity(found.saturating_sub(8));
        for i in 8..found {
            let sample_start = tabs[i] + 1;
            let sample_end = tabs.get(i + 1).copied().unwrap_or(line_len);
            if sample_start < sample_end {
                sample_vec.push(&line_str[sample_start..sample_end]);
            }
        }

        (Some(format_field), sample_vec)
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

#[target_feature(enable = "neon")]
pub unsafe fn parse_vcf_fields_neon(line: &[u8]) -> Option<VcfFields<'_>> {
    if line.len() < 20 {
        return None;
    }

    let mut tabs = [0usize; 8];
    let n = find_tabs_neon(line, &mut tabs);

    if n < 7 {
        return None;
    }

    let line_str = std::str::from_utf8_unchecked(line);
    let line_len = line.len();

    let info_start = tabs[6] + 1;
    let info_end = if n >= 8 { tabs[7] } else { line_str.len() };

    if info_start > info_end || info_end > line_str.len() {
        return None;
    }

    Some(VcfFields {
        chrom: &line_str[0..tabs[0]],
        pos: &line_str[tabs[0] + 1..tabs[1]],
        id: &line_str[tabs[1] + 1..tabs[2]],
        ref_allele: &line_str[tabs[2] + 1..tabs[3]],
        alt: &line_str[tabs[3] + 1..tabs[4]],
        qual: &line_str[tabs[4] + 1..tabs[5]],
        filter: &line_str[tabs[5] + 1..tabs[6]],
        info: &line_str[info_start..info_end],
    })
}

#[target_feature(enable = "neon")]
pub unsafe fn parse_chr_pos_neon(line: &[u8]) -> Option<(u8, u32)> {
    let mut tabs = [0usize; 2];
    let n = find_tabs_neon(line, &mut tabs);

    if n < 1 {
        return None;
    }

    let chrom_end = tabs[0];
    let chrom_bytes = &line[..chrom_end];
    let chrom_str = std::str::from_utf8_unchecked(chrom_bytes);
    let chr_id = chr_name_to_id(chrom_str)?;

    let pos_start = tabs[0] + 1;
    let pos_end = if n >= 2 { tabs[1] } else { line.len() };

    if pos_start >= pos_end {
        return None;
    }

    let pos = parse_u32_fast(line.as_ptr().add(pos_start), pos_end - pos_start);

    Some((chr_id, pos))
}

#[target_feature(enable = "neon")]
unsafe fn find_tabs_neon(buf: &[u8], out: &mut [usize]) -> usize {
    let max_tabs = out.len();
    let mut found = 0;
    let mut pos = 0;
    let tab = vdupq_n_u8(b'\t');

    while pos + 16 <= buf.len() && found < max_tabs {
        let chunk = vld1q_u8(buf.as_ptr().add(pos));
        let cmp = vceqq_u8(chunk, tab);

        let mask_parts: [u8; 16] = std::mem::transmute(cmp);
        for (i, &byte) in mask_parts.iter().enumerate() {
            if byte == 0xFF && found < max_tabs {
                out[found] = pos + i;
                found += 1;
            }
        }

        pos += 16;
    }

    while pos < buf.len() && found < max_tabs {
        if buf[pos] == b'\t' {
            out[found] = pos;
            found += 1;
        }
        pos += 1;
    }

    found
}

#[target_feature(enable = "neon")]
pub unsafe fn parse_u32_fast(ptr: *const u8, len: usize) -> u32 {
    let mut x: u32 = 0;
    for i in 0..len.min(10) {
        let c = *ptr.add(i) as u32;
        if c < b'0' as u32 || c > b'9' as u32 {
            break;
        }
        x = x * 10 + (c - b'0' as u32);
    }
    x
}

#[target_feature(enable = "neon")]
pub unsafe fn find_char_neon(buf: &[u8], ch: u8) -> Option<usize> {
    let target = vdupq_n_u8(ch);
    let mut pos = 0;

    while pos + 16 <= buf.len() {
        let chunk = vld1q_u8(buf.as_ptr().add(pos));
        let cmp = vceqq_u8(chunk, target);

        let mask_parts: [u8; 16] = std::mem::transmute(cmp);
        for (i, &byte) in mask_parts.iter().enumerate() {
            if byte == 0xFF {
                return Some(pos + i);
            }
        }

        pos += 16;
    }

    while pos < buf.len() {
        if buf[pos] == ch {
            return Some(pos);
        }
        pos += 1;
    }

    None
}
