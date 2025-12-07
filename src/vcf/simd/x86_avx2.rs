use crate::util::chr_name_to_id;
use crate::vcf::structs::{VcfFields, VcfFieldsFull};

#[cfg(target_arch = "x86_64")]
use core::arch::x86_64::*;

const MAX_TABS: usize = 138;

#[target_feature(enable = "avx2")]
pub unsafe fn parse_vcf_line_simd(line: &[u8]) -> Option<VcfFieldsFull<'_>> {
    if line.len() < 20 {
        return None;
    }

    let mut tabs = [0usize; MAX_TABS];
    let found = find_tabs_avx2(line, &mut tabs);

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
    let info = &line_str[tabs[6] + 1..tabs.get(7).copied().unwrap_or(line.len())];

    let (format, samples) = if found >= 8 {
        let format_field = &line_str[tabs[7] + 1..tabs.get(8).copied().unwrap_or(line.len())];

        let mut sample_vec = Vec::with_capacity(found.saturating_sub(8));
        for i in 8..found {
            let sample_start = tabs[i] + 1;
            let sample_end = tabs.get(i + 1).copied().unwrap_or(line.len());
            sample_vec.push(&line_str[sample_start..sample_end]);
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

#[target_feature(enable = "avx2")]
pub unsafe fn parse_vcf_fields_avx2(line: &[u8]) -> Option<VcfFields<'_>> {
    if line.len() < 20 {
        return None;
    }

    let mut tabs = [0usize; 8];
    let n = find_tabs_avx2(line, &mut tabs);

    if n < 7 {
        return None;
    }

    let line_str = std::str::from_utf8_unchecked(line);

    Some(VcfFields {
        chrom: &line_str[0..tabs[0]],
        pos: &line_str[tabs[0] + 1..tabs[1]],
        id: &line_str[tabs[1] + 1..tabs[2]],
        ref_allele: &line_str[tabs[2] + 1..tabs[3]],
        alt: &line_str[tabs[3] + 1..tabs[4]],
        qual: &line_str[tabs[4] + 1..tabs[5]],
        filter: &line_str[tabs[5] + 1..tabs[6]],
        info: &line_str[tabs[6] + 1..tabs.get(7).copied().unwrap_or(line.len())],
    })
}

#[target_feature(enable = "avx2")]
pub unsafe fn parse_chr_pos_avx2(line: &[u8]) -> Option<(u8, u32)> {
    let mut tabs = [0usize; 2];
    let n = find_tabs_avx2(line, &mut tabs);

    if n < 1 {
        return None;
    }

    let chrom_end = tabs[0];
    let chrom_bytes = &line[..chrom_end];
    let chrom_str = std::str::from_utf8_unchecked(chrom_bytes);
    let chr_id = chr_name_to_id(chrom_str)?;

    let pos_start = tabs[0] + 1;
    let pos_end = if n >= 2 { tabs[1] } else { line.len() };
    let pos = parse_u32_fast(line.as_ptr().add(pos_start), pos_end - pos_start);

    Some((chr_id, pos))
}

#[target_feature(enable = "avx2")]
unsafe fn find_tabs_avx2(buf: &[u8], out: &mut [usize]) -> usize {
    let max_tabs = out.len();
    let mut found = 0;
    let mut pos = 0;
    let tab = _mm256_set1_epi8(b'\t' as i8);

    while pos + 32 <= buf.len() && found < max_tabs {
        let chunk = _mm256_loadu_si256(buf.as_ptr().add(pos) as *const __m256i);
        let cmp = _mm256_cmpeq_epi8(chunk, tab);
        let mask = _mm256_movemask_epi8(cmp);

        if mask != 0 {
            let mut m = mask as u32;
            while m != 0 && found < max_tabs {
                let tz = m.trailing_zeros() as usize;
                out[found] = pos + tz;
                found += 1;
                m &= m - 1;
            }
        }

        pos += 32;
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

#[target_feature(enable = "avx2")]
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
