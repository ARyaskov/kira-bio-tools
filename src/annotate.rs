// Full annotate.rs adapted for ANI v2 structured index
// CPU implementation

use crate::annotate_index::{AniIndex, AnnotationBundle, FieldNumber, StructuredInfoField};
use anyhow::Result;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

/// ------------------------------------------------------------
/// Main annotate function (CPU)
/// ------------------------------------------------------------
pub fn annotate_vcf_ani(
    db: &std::path::Path,
    input: &std::path::Path,
    output: &std::path::Path,
) -> Result<()> {
    let ani = AniIndex::open(db)?;

    let fin = File::open(input)?;
    let rdr = BufReader::new(fin);

    let fout = File::create(output)?;
    let mut bw = BufWriter::new(fout);

    for line in rdr.lines() {
        let line = line?;

        // Header lines are passed through unchanged
        if line.starts_with('#') {
            bw.write_all(line.as_bytes())?;
            bw.write_all(b"\n")?;
            continue;
        }

        // Parse data line
        if let Some(row) = parse_vcf_record(&line) {
            let (chr, pos, id, r, alt_raw, qual, filter, info, rest) = row;
            let alt_list: Vec<&str> = alt_raw.split(',').collect();

            // Now lookup returns: (bundle, ann_alt_list)
            if let Some((ann, ann_alt_list)) = ani.lookup_full(chr, pos, r, alt_raw) {
                let merged = merge_record(
                    chr,
                    pos,
                    id,
                    r,
                    &alt_list,
                    qual,
                    filter,
                    info,
                    rest,
                    ann,
                    &ann_alt_list,
                );

                bw.write_all(merged.as_bytes())?;
                bw.write_all(b"\n")?;
                continue;
            }
        }

        // If no annotation found, write original
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

/// Merge full VCF record with ANI annotation.
pub fn merge_record(
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
    ann_alt_list: &[&str], // required for Number=A mapping
) -> String {
    // -------------------------
    // Merge ID
    // -------------------------
    let id2 = match ann.id {
        Some(".") => id,
        Some(v) => v,
        None => id,
    };

    // -------------------------
    // Merge QUAL
    // -------------------------
    let qual2 = match ann.qual {
        Some(".") => qual,
        Some(v) => v,
        None => qual,
    };

    // -------------------------
    // Merge FILTER
    // -------------------------
    let filter2 = match ann.filter {
        Some(".") => filter,
        Some(v) => v,
        None => filter,
    };

    // -------------------------
    // Merge INFO (delegated)
    // -------------------------
    let merged_info = merge_info(info, alt_list, &ann.info, r);

    // -------------------------
    // Build final VCF line
    // -------------------------
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

fn merge_info(
    base: &str,
    alt_list: &[&str],
    ann_fields: &[StructuredInfoField],
    r: &str,
) -> String {
    let mut out: Vec<String> = Vec::new();

    // 1) Base INFO
    if !base.is_empty() && base != "." {
        for s in base.split(';') {
            if !s.is_empty() {
                out.push(s.to_string());
            }
        }
    }

    // 2) ANI INFO
    for f in ann_fields {
        match f.number {
            FieldNumber::Zero => {
                out.push(f.key.to_string());
            }

            FieldNumber::One => {
                out.push(format!("{}={}", f.key, f.values[0]));
            }

            FieldNumber::Many => {
                out.push(format!("{}={}", f.key, f.values.join(",")));
            }

            FieldNumber::A => {
                // map allele values
                let mut vals = Vec::with_capacity(alt_list.len());
                for i in 0..alt_list.len() {
                    vals.push(if i < f.values.len() {
                        f.values[i]
                    } else {
                        f.values[0]
                    });
                }
                out.push(format!("{}={}", f.key, vals.join(",")));
            }

            _ => {
                out.push(format!("{}={}", f.key, f.values.join(",")));
            }
        }
    }

    // 3) INDEL tag
    let is_indel = r.len() != 1 || alt_list.iter().any(|a| a.len() != 1);

    if is_indel {
        out.push("INDEL".to_string());
    }

    // 4) Deduplicate keys
    let mut seen = std::collections::HashSet::new();
    let mut dedup = Vec::new();

    for item in out {
        let key = item.split('=').next().unwrap().to_string();
        if seen.insert(key) {
            dedup.push(item);
        }
    }

    // 5) Sort keys
    dedup.sort();

    // 6) Join
    if dedup.is_empty() {
        ".".to_string()
    } else {
        dedup.join(";")
    }
}
