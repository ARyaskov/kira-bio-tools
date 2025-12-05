use anyhow::Result;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use super::structs::*;

pub fn annotate_vcf_ani(db: &Path, input: &Path, output: &Path) -> Result<()> {
    let ani = AniIndex::open(db)?;

    let fin = File::open(input)?;
    let rdr = BufReader::new(fin);

    let fout = File::create(output)?;
    let mut bw = BufWriter::new(fout);

    for line in rdr.lines() {
        let line = line?;

        if line.starts_with('#') {
            bw.write_all(line.as_bytes())?;
            bw.write_all(b"\n")?;
            continue;
        }

        if let Some(row) = parse_vcf_record(&line) {
            let (chr, pos, id, r, alt_raw, qual, filter, info, rest) = row;
            let alt_list: Vec<&str> = alt_raw.split(',').collect();

            if let Some((ann, _ann_alt_list)) = ani.lookup_full(chr, pos, r, alt_raw) {
                let merged =
                    merge_record(chr, pos, id, r, &alt_list, qual, filter, info, rest, ann);

                bw.write_all(merged.as_bytes())?;
                bw.write_all(b"\n")?;
                continue;
            }
        }

        bw.write_all(line.as_bytes())?;
        bw.write_all(b"\n")?;
    }

    Ok(())
}

fn parse_vcf_record<'a>(
    line: &'a str,
) -> Option<(
    &'a str,
    u32,
    &'a str,
    &'a str,
    &'a str,
    &'a str,
    &'a str,
    &'a str,
    Vec<&'a str>,
)> {
    let mut c = line.split('\t');
    let chr = c.next()?;
    let pos = c.next()?.parse().ok()?;
    let id = c.next()?;
    let r = c.next()?;
    let alt = c.next()?;
    let qual = c.next()?;
    let filter = c.next()?;
    let info = c.next()?;
    let rest: Vec<&str> = c.collect();
    Some((chr, pos, id, r, alt, qual, filter, info, rest))
}

fn merge_record(
    chr: &str,
    pos: u32,
    id: &str,
    r: &str,
    alt_list: &[&str],
    qual: &str,
    filter: &str,
    info: &str,
    rest: Vec<&str>,
    ann: AnnotationBundle,
) -> String {
    // Merge ID
    let id2 = match ann.id {
        Some(ann_id) if ann_id != "." => ann_id,
        _ => id,
    };

    // Merge QUAL
    let qual2 = match ann.qual {
        Some(ann_qual) if ann_qual != "." => ann_qual,
        _ => qual,
    };

    // Merge FILTER
    let filter2 = match ann.filter {
        Some(ann_filter) if ann_filter != "." => ann_filter,
        _ => filter,
    };

    // Merge INFO with bcftools-compatible sorting
    let merged_info = merge_info_bcftools_order(info, alt_list, &ann.info, r);

    let mut out = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        chr,
        pos,
        id2,
        r,
        alt_list.join(","),
        qual2,
        filter2,
        merged_info,
    );

    for f in rest {
        out.push('\t');
        out.push_str(f);
    }

    out
}

/// Custom sorting key to match bcftools INFO field order
/// bcftools uses case-insensitive alphabetical sort, but with stable ordering
fn bcftools_sort_key(key: &str) -> String {
    // Convert to uppercase for case-insensitive comparison
    // This ensures: AC, AN, DP4, INDEL, STR, TEST, T_FLOAT, T_INT
    key.to_uppercase()
}

/// Merge INFO fields with bcftools-compatible alphabetical sorting
fn merge_info_bcftools_order(
    base: &str,
    alt_list: &[&str],
    ann_fields: &[StructuredInfoField],
    r: &str,
) -> String {
    use std::collections::HashMap;

    let mut map = HashMap::new();

    // Parse base INFO
    if !base.is_empty() && base != "." {
        for s in base.split(';') {
            if s.is_empty() {
                continue;
            }
            if let Some((k, v)) = s.split_once('=') {
                map.insert(k.to_string(), v.to_string());
            } else {
                map.insert(s.to_string(), String::new());
            }
        }
    }

    // Add annotation fields
    for f in ann_fields {
        let value = match f.number {
            FieldNumber::Zero => {
                map.insert(f.key.to_string(), String::new());
                continue;
            }

            FieldNumber::One => {
                if let Some(v) = f.values.first() {
                    v.to_string()
                } else {
                    continue;
                }
            }

            FieldNumber::Many => f.values.join(","),

            FieldNumber::A => {
                // Expand to match number of ALT alleles
                let mut vals = Vec::with_capacity(alt_list.len());
                for i in 0..alt_list.len() {
                    vals.push(if i < f.values.len() {
                        f.values[i]
                    } else {
                        f.values[0]
                    });
                }
                vals.join(",")
            }

            _ => f.values.join(","),
        };

        map.insert(f.key.to_string(), value);
    }

    // Add INDEL flag if applicable
    let is_indel = r.len() != 1 || alt_list.iter().any(|a| a.len() != 1);
    if is_indel {
        map.insert("INDEL".to_string(), String::new());
    }

    // Sort keys using bcftools ordering (case-insensitive alphabetical)
    let mut keys: Vec<String> = map.keys().cloned().collect();
    keys.sort_by(|a, b| bcftools_sort_key(a).cmp(&bcftools_sort_key(b)));

    // Build output in sorted order
    let mut parts = Vec::new();
    for key in keys {
        if let Some(value) = map.get(&key) {
            if value.is_empty() {
                parts.push(key);
            } else {
                parts.push(format!("{}={}", key, value));
            }
        }
    }

    if parts.is_empty() {
        ".".to_string()
    } else {
        parts.join(";")
    }
}
