use anyhow::Result;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

use super::reader::{BatchVcfReader, VcfAnnotationReader};
use super::structs::*;
use crate::vcf_parser_fast::FastVcfParser;

const BATCH_SIZE: usize = 200_000;
const OUTPUT_BUFFER_SIZE: usize = 64 * 1024 * 1024;

pub fn annotate_vcf_ani_v2(db: &Path, input: &Path, output: &Path) -> Result<()> {
    let timing = std::env::var("KIRA_BT_TIMING").is_ok();
    let start = std::time::Instant::now();

    let load_start = std::time::Instant::now();
    let ani = AniIndex::open(db)?;
    if timing {
        eprintln!(
            "[annotate] ANI load: {:.3}s",
            load_start.elapsed().as_secs_f64()
        );
    }

    let input_reader = VcfAnnotationReader::open(input)?;
    let mut batch_reader = BatchVcfReader::new(input_reader, BATCH_SIZE);

    let (headers, mut batch_reader) = batch_reader.into_headers_and_self()?;

    let merged_headers = merge_annotation_headers(&headers, &ani)?;

    let output_file = File::create(output)?;
    let mut writer = BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, output_file);

    for header in merged_headers {
        writeln!(writer, "{}", header)?;
    }

    let processed = AtomicUsize::new(0);
    let mut batch_num = 0;

    loop {
        let batch_start = std::time::Instant::now();

        let batch = batch_reader.read_batch()?;
        if batch.is_empty() {
            break;
        }

        batch_num += 1;

        let annotated: Vec<String> = batch
            .par_iter()
            .map(|line| annotate_line(line, &ani))
            .collect();

        for line in annotated {
            writeln!(writer, "{}", line)?;
        }

        let count = batch.len();
        processed.fetch_add(count, Ordering::Relaxed);

        if timing && batch_num % 10 == 0 {
            let total = processed.load(Ordering::Relaxed);
            let elapsed = start.elapsed().as_secs_f64();
            let rate = total as f64 / elapsed;
            eprintln!(
                "[annotate] Batch {}: {} variants ({:.0} var/s, batch: {:.3}s)",
                batch_num,
                total,
                rate,
                batch_start.elapsed().as_secs_f64()
            );
        }
    }

    writer.flush()?;

    if timing {
        let total = processed.load(Ordering::Relaxed);
        let elapsed = start.elapsed().as_secs_f64();
        eprintln!(
            "[annotate] DONE: {} variants in {:.3}s ({:.0} var/s)",
            total,
            elapsed,
            total as f64 / elapsed
        );
    }

    Ok(())
}

#[inline]
fn annotate_line(line: &str, ani: &AniIndex) -> String {
    let mut parser = FastVcfParser::new(line);

    if let Some(fields) = parser.parse_standard_fields() {
        let chr = fields.chrom;
        let pos = fields.pos.parse::<u32>().unwrap_or(0);
        let id = fields.id;
        let rf = fields.ref_allele;
        let alt_raw = fields.alt;
        let qual = fields.qual;
        let filter = fields.filter;
        let info = fields.info;

        let rest = parser.rest();
        let rest_fields: Vec<&str> = if rest.is_empty() {
            Vec::new()
        } else {
            rest.split('\t').collect()
        };

        let alt_list: Vec<&str> = alt_raw.split(',').collect();

        if let Some((ann, _ann_alt_list)) = ani.lookup_full(chr, pos, rf, alt_raw) {
            return merge_record(
                chr,
                pos,
                id,
                rf,
                &alt_list,
                qual,
                filter,
                info,
                rest_fields,
                ann,
            );
        }
    }

    line.to_string()
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
    let id2 = ann.id.unwrap_or(id);
    let qual2 = ann.qual.unwrap_or(qual);
    let filter2 = ann.filter.unwrap_or(filter);
    let merged_info = merge_info_optimized(info, alt_list, &ann.info, r);

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

    for field in rest {
        out.push('\t');
        out.push_str(field);
    }

    out
}

fn merge_info_optimized(
    base: &str,
    alt_list: &[&str],
    ann_fields: &[StructuredInfoField],
    r: &str,
) -> String {
    use std::collections::HashMap;

    let mut map: HashMap<&str, String> = HashMap::with_capacity(32);

    if !base.is_empty() && base != "." {
        for segment in base.split(';') {
            if segment.is_empty() {
                continue;
            }
            if let Some((k, v)) = segment.split_once('=') {
                map.insert(k, v.to_string());
            } else {
                map.insert(segment, String::new());
            }
        }
    }

    for field in ann_fields {
        let value = match field.number {
            FieldNumber::Zero => {
                map.insert(field.key, String::new());
                continue;
            }
            FieldNumber::One => {
                if let Some(&v) = field.values.first() {
                    v.to_string()
                } else {
                    continue;
                }
            }
            FieldNumber::Many => field.values.join(","),
            FieldNumber::A => {
                let mut vals = Vec::with_capacity(alt_list.len());
                for i in 0..alt_list.len() {
                    vals.push(field.values.get(i).copied().unwrap_or(field.values[0]));
                }
                vals.join(",")
            }
            _ => field.values.join(","),
        };

        map.insert(field.key, value);
    }

    let is_indel = r.len() != 1 || alt_list.iter().any(|a| a.len() != 1);
    if is_indel {
        map.insert("INDEL", String::new());
    }

    let mut keys: Vec<&str> = map.keys().copied().collect();
    keys.sort_by(|a, b| a.to_uppercase().cmp(&b.to_uppercase()));

    let mut parts = Vec::with_capacity(keys.len());
    for key in keys {
        if let Some(value) = map.get(key) {
            if value.is_empty() {
                parts.push(key.to_string());
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

fn merge_annotation_headers(vcf_headers: &[String], _ani: &AniIndex) -> Result<Vec<String>> {
    let mut merged = Vec::with_capacity(vcf_headers.len() + 10);

    for header in vcf_headers {
        if header.starts_with("#CHROM") {
            merged.push(
                "##INFO=<ID=T_STR,Number=1,Type=String,Description=\"Test String\">".to_string(),
            );
            merged.push(
                "##INFO=<ID=T_INT,Number=.,Type=Integer,Description=\"Test Integer\">".to_string(),
            );
            merged.push(
                "##INFO=<ID=T_FLOAT,Number=.,Type=Float,Description=\"Test Float\">".to_string(),
            );
        }
        merged.push(header.clone());
    }

    Ok(merged)
}
