use crate::vcf::structs::{VcfFields, VcfFieldsFull};

const MAX_TABS: usize = 138;

pub fn parse_vcf_line_scalar(line: &[u8]) -> Option<VcfFieldsFull<'_>> {
    if line.len() < 20 {
        return None;
    }

    let mut tabs = [0usize; MAX_TABS];
    let found = find_tabs_scalar(line, &mut tabs);

    if found < 7 {
        return None;
    }

    let line_str = std::str::from_utf8(line).ok()?;

    let chrom = &line_str[0..tabs[0]];
    let pos = &line_str[tabs[0] + 1..tabs[1]];
    let id = &line_str[tabs[1] + 1..tabs[2]];
    let ref_allele = &line_str[tabs[2] + 1..tabs[3]];
    let alt = &line_str[tabs[3] + 1..tabs[4]];
    let qual = &line_str[tabs[4] + 1..tabs[5]];
    let filter = &line_str[tabs[5] + 1..tabs[6]];
    let info_end = if found >= 8 { tabs[7] } else { line_str.len() };
    let info = &line_str[tabs[6] + 1..info_end];

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

pub fn parse_vcf_fields_scalar(line: &[u8]) -> Option<VcfFields<'_>> {
    if line.len() < 20 {
        return None;
    }

    let mut tabs = [0usize; 8];
    let n = find_tabs_scalar(line, &mut tabs);

    if n < 7 {
        return None;
    }

    let line_str = std::str::from_utf8(line).ok()?;

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

pub fn parse_chr_pos_scalar(line: &[u8]) -> Option<(u8, u32)> {
    let chrom_end = memchr::memchr(b'\t', line)?;
    let chr_id = parse_chromosome_scalar(&line[..chrom_end])?;

    let pos_start = chrom_end + 1;
    let remaining = &line[pos_start..];
    let pos_end = memchr::memchr(b'\t', remaining).unwrap_or(remaining.len());
    let pos = parse_u32_scalar(&remaining[..pos_end])?;

    Some((chr_id, pos))
}

fn find_tabs_scalar(buf: &[u8], out: &mut [usize]) -> usize {
    let max_tabs = out.len();
    let mut found = 0;

    for (i, &byte) in buf.iter().enumerate() {
        if byte == b'\t' {
            out[found] = i;
            found += 1;
            if found >= max_tabs {
                break;
            }
        }
    }

    found
}

pub fn parse_u32_scalar(bytes: &[u8]) -> Option<u32> {
    if bytes.is_empty() || bytes.len() > 10 {
        return None;
    }

    let mut result = 0u32;
    for &byte in bytes {
        if !byte.is_ascii_digit() {
            return None;
        }
        result = result.wrapping_mul(10).wrapping_add((byte - b'0') as u32);
    }

    Some(result)
}

fn parse_chromosome_scalar(bytes: &[u8]) -> Option<u8> {
    let bytes = if bytes.len() > 3 && &bytes[..3] == b"chr" {
        &bytes[3..]
    } else {
        bytes
    };

    match bytes {
        b"X" | b"x" => Some(23),
        b"Y" | b"y" => Some(24),
        b"M" | b"MT" | b"m" | b"mt" => Some(25),
        _ if bytes.len() <= 2 && bytes.iter().all(|b| b.is_ascii_digit()) => {
            let mut num = 0u8;
            for &byte in bytes {
                num = num.wrapping_mul(10).wrapping_add(byte - b'0');
            }
            if (1..=22).contains(&num) {
                Some(num)
            } else {
                None
            }
        }
        _ => None,
    }
}
