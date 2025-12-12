use anyhow::Result;
use crossbeam_channel::{bounded, Receiver, Sender};
use indexmap::IndexMap;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::thread;
use std::time::Instant;

use super::constants::*;
use super::reader::{StreamingVcfReader, VcfAnnotationReader};
use super::structs::bundle::{parse_info_field, AnnotationBundle, FieldNumber};
use super::structs::*;
use crate::bgzf::BgzfWriter;
use crate::util::{
    choose_best_number, detect_format, extract_info_key, extract_info_number, read_cstring,
    url_decode_info_value, VcfFormat,
};
use crate::vcf::VcfParser;

pub fn annotate_vcf_ani_v2(db: &Path, input: &Path, output: &Path) -> Result<()> {
    let timing = std::env::var("KIRA_BT_TIMING").is_ok();
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok() || timing;
    let start = Instant::now();

    let num_threads = rayon::current_num_threads() / 4;
    if debug {
        eprintln!("[annotate] Using {} CPU threads", num_threads);
        eprintln!("[annotate] Batch size: {} lines", BATCH_SIZE);
    }

    let load_start = Instant::now();
    let ani = AniIndex::open(db)?;
    if debug {
        eprintln!(
            "[annotate] ANI load: {:.3}s",
            load_start.elapsed().as_secs_f64()
        );
    }

    let mut field_meta = load_field_metadata(&ani)?;

    if field_meta.is_empty() {
        if debug {
            eprintln!("[annotate] No metadata in ANI headers, inferring from data...");
        }

        let field_names: Vec<String> = ani
            .entries
            .iter()
            .take(10)
            .flat_map(|e| {
                let info_str = read_cstring(&ani.strings, e.info_ofs as usize);
                let decoded = url_decode_info_value(info_str);
                decoded
                    .split(';')
                    .filter_map(|kv| kv.split('=').next().map(|k| k.to_string()))
                    .collect::<Vec<_>>()
            })
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();

        field_meta = infer_field_metadata_from_data(&ani, &field_names);
    }

    let input_format = detect_format(input)?;
    let use_bgzf = matches!(input_format, VcfFormat::Bgzf);
    if debug {
        eprintln!(
            "[annotate] Input: {:?}, Output BGZF: {}",
            input_format, use_bgzf
        );
    }

    let reader_start = Instant::now();
    let input_reader = VcfAnnotationReader::open(input)?;
    let streaming_reader = StreamingVcfReader::new(input_reader);

    let (headers, mut reader) = streaming_reader.into_headers_and_self()?;
    if debug {
        eprintln!(
            "[annotate] Reader init: {:.3}s, headers: {}",
            reader_start.elapsed().as_secs_f64(),
            headers.len()
        );
    }

    let merged_headers = merge_annotation_headers(&headers, &ani)?;

    let (read_tx, read_rx) = bounded::<Vec<String>>(CHANNEL_DEPTH);
    let (work_tx, work_rx) = bounded::<Vec<String>>(CHANNEL_DEPTH);

    let ani_clone = std::sync::Arc::new(ani);
    let field_meta_clone = std::sync::Arc::new(field_meta);
    let ani_worker = ani_clone.clone();
    let field_meta_worker = field_meta_clone.clone();

    let worker = thread::spawn(move || {
        worker_thread(read_rx, work_tx, ani_worker, field_meta_worker, num_threads)
    });

    let output_clone = output.to_path_buf();
    let writer = thread::spawn(move || {
        writer_thread(work_rx, merged_headers, &output_clone, use_bgzf, timing)
    });

    let mut batch = Vec::with_capacity(BATCH_SIZE);
    let mut total_lines = 0usize;
    let mut last_report = Instant::now();

    while let Some(line) = reader.read_line()? {
        batch.push(line);

        if batch.len() >= BATCH_SIZE {
            total_lines += batch.len();
            if let Err(_) = read_tx.send(std::mem::replace(
                &mut batch,
                Vec::with_capacity(BATCH_SIZE),
            )) {
                break;
            }

            if timing && last_report.elapsed().as_secs() >= 2 {
                eprintln!("[annotate] Progress: {} lines read", total_lines);
                last_report = Instant::now();
            }
        }
    }

    if !batch.is_empty() {
        total_lines += batch.len();
        let _ = read_tx.send(batch);
    }

    drop(read_tx);

    worker.join().unwrap()?;
    writer.join().unwrap()?;

    if debug {
        eprintln!(
            "[annotate] Total time: {:.3}s for {} lines",
            start.elapsed().as_secs_f64(),
            total_lines
        );
    }

    Ok(())
}

fn load_field_metadata(ani: &AniIndex) -> Result<HashMap<String, FieldNumber>> {
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok();
    let mut metadata = HashMap::new();

    let strings_str = std::str::from_utf8(&ani.strings).unwrap_or("");

    if debug {
        eprintln!("[DEBUG] ANI strings total length: {}", ani.strings.len());
        let preview_len = strings_str.len().min(500);
        eprintln!(
            "[DEBUG] First {} chars of strings: {:?}",
            preview_len,
            &strings_str[..preview_len]
        );

        let header_count = strings_str
            .split('\0')
            .filter(|s| s.starts_with("##INFO="))
            .count();
        eprintln!(
            "[DEBUG] Found {} ##INFO headers in ANI strings",
            header_count
        );
    }

    for line in strings_str.split('\0') {
        if !line.starts_with("##INFO=") {
            continue;
        }

        if let Some(key) = extract_info_key(line) {
            if let Some(number) = extract_info_number(line) {
                metadata.insert(key.clone(), number);

                if debug {
                    eprintln!(
                        "[DEBUG] Loaded metadata from header: {} -> {:?}",
                        key, number
                    );
                }
            }
        }
    }

    if debug {
        eprintln!(
            "[DEBUG] Total metadata entries loaded from headers: {}",
            metadata.len()
        );
    }

    Ok(metadata)
}

fn infer_field_metadata_from_data(
    ani: &AniIndex,
    field_names: &[String],
) -> HashMap<String, FieldNumber> {
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok();
    let mut candidates: HashMap<String, Vec<FieldNumber>> = HashMap::new();

    if debug {
        eprintln!(
            "[DEBUG] Inferring metadata from data for {} fields",
            field_names.len()
        );
    }

    for entry_idx in 0..ani.entries.len().min(100) {
        let entry = &ani.entries[entry_idx];
        let info_str = read_cstring(&ani.strings, entry.info_ofs as usize);
        let decoded_info = url_decode_info_value(info_str);

        for field_name in field_names {
            let pattern = format!("{}=", field_name);
            if let Some(start) = decoded_info.find(&pattern) {
                let value_start = start + pattern.len();
                let rest = &decoded_info[value_start..];
                let value_end = rest.find(';').unwrap_or(rest.len());
                let value = &rest[..value_end];

                let alt_str = read_cstring(&ani.strings, entry.alt_ofs as usize);
                let num_alts = alt_str.split(',').count();
                let num_values = value.split(',').count();

                let number = if num_values == num_alts {
                    FieldNumber::A
                } else if num_values == num_alts + 1 {
                    FieldNumber::R
                } else if num_values == 1 {
                    FieldNumber::One
                } else {
                    FieldNumber::Many
                };

                candidates
                    .entry(field_name.clone())
                    .or_insert_with(Vec::new)
                    .push(number);

                if debug && entry_idx < 3 {
                    eprintln!(
                        "[DEBUG] Entry {}: {} -> {:?} (alts={}, values={}, raw={:?})",
                        entry_idx,
                        field_name,
                        number,
                        num_alts,
                        num_values,
                        &value[..value.len().min(30)]
                    );
                }
            }
        }
    }

    let mut metadata = HashMap::new();
    for (field_name, numbers) in candidates.into_iter() {
        let best = choose_best_number(&numbers);

        if debug {
            eprintln!(
                "[DEBUG] Final inference: {} -> {:?} (from {} samples)",
                field_name,
                best,
                numbers.len()
            );
        }

        metadata.insert(field_name, best);
    }

    if debug {
        eprintln!(
            "[DEBUG] Total metadata inferred from data: {}",
            metadata.len()
        );
    }

    metadata
}

fn merge_allele_values(
    existing: &[&str],
    database: &[String],
    number: FieldNumber,
    field_type: &str,
) -> Vec<String> {
    let mut result = Vec::new();

    let all_missing = existing.iter().all(|v| *v == "." || v.is_empty());

    if field_type == "Integer" && !all_missing {
        return database.to_vec();
    }

    match number {
        FieldNumber::A => {
            for i in 0..database.len() {
                if let Some(&ex) = existing.get(i) {
                    if ex == "." || ex.is_empty() {
                        result.push(database[i].clone());
                    } else {
                        result.push(ex.to_string());
                    }
                } else {
                    result.push(database[i].clone());
                }
            }
        }
        FieldNumber::R => {
            for i in 0..database.len() {
                if let Some(&ex) = existing.get(i) {
                    if ex == "." || ex.is_empty() {
                        result.push(database[i].clone());
                    } else {
                        result.push(ex.to_string());
                    }
                } else {
                    result.push(database[i].clone());
                }
            }
        }
        _ => result = database.to_vec(),
    }

    result
}

fn infer_field_type(key: &str) -> &'static str {
    if key.starts_with('I') || key.ends_with("INT") {
        return "Integer";
    }
    if key.starts_with('F') || key.ends_with("FLT") || key.ends_with("FLOAT") {
        return "Float";
    }
    if key.starts_with('S') || key.ends_with("STR") || key.ends_with("STRING") {
        return "String";
    }
    "String"
}

fn worker_thread(
    rx: Receiver<Vec<String>>,
    tx: Sender<Vec<String>>,
    ani: std::sync::Arc<AniIndex>,
    field_meta: std::sync::Arc<HashMap<String, FieldNumber>>,
    num_threads: usize,
) -> Result<()> {
    use rayon::prelude::*;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    while let Ok(batch) = rx.recv() {
        let annotated: Vec<String> = pool.install(|| {
            batch
                .par_iter()
                .map(|line| annotate_line(line, &ani, &field_meta))
                .collect()
        });

        if tx.send(annotated).is_err() {
            break;
        }
    }

    Ok(())
}

fn writer_thread(
    rx: Receiver<Vec<String>>,
    headers: Vec<String>,
    output: &Path,
    use_bgzf: bool,
    timing: bool,
) -> Result<()> {
    let start = Instant::now();
    let mut lines_written = 0usize;
    let mut bytes_written = 0usize;

    if use_bgzf {
        let mut writer = BgzfWriter::create(output)?;

        for h in headers {
            writeln!(writer, "{}", h)?;
            bytes_written += h.len() + 1;
        }

        while let Ok(batch) = rx.recv() {
            for line in batch {
                writeln!(writer, "{}", line)?;
                bytes_written += line.len() + 1;
                lines_written += 1;
            }
        }

        writer.finish()?;
    } else {
        let file = File::create(output)?;
        let mut writer = BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, file);

        for h in headers {
            writeln!(writer, "{}", h)?;
            bytes_written += h.len() + 1;
        }

        while let Ok(batch) = rx.recv() {
            for line in batch {
                writeln!(writer, "{}", line)?;
                bytes_written += line.len() + 1;
                lines_written += 1;
            }
        }

        writer.flush()?;
    }

    if timing {
        let elapsed = start.elapsed().as_secs_f64();
        let mb_sec = (bytes_written as f64 / 1_048_576.0) / elapsed;
        eprintln!(
            "[annotate] Write complete: {} lines, {:.1} MB/s",
            lines_written, mb_sec
        );
    }

    Ok(())
}

fn annotate_line(line: &str, ani: &AniIndex, field_meta: &HashMap<String, FieldNumber>) -> String {
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok();

    if debug && field_meta.is_empty() {
        eprintln!("[DEBUG] WARNING: field_meta is EMPTY! No metadata loaded from headers!");
    }

    let mut parser = VcfParser::new(line);
    let Some(rec) = parser.parse_standard_fields() else {
        return line.to_string();
    };

    let pos = rec.pos.parse::<u32>().unwrap_or(0);

    let vcf_alt_alleles: Vec<&str> = rec.alt.split(',').collect();

    let mut updated_id = rec.id.to_string();
    let mut updated_filter = rec.filter.to_string();
    let mut info_map: IndexMap<String, String> = IndexMap::new();

    for kv in rec.info.split(';') {
        if kv.is_empty() || kv == "." {
            continue;
        }
        let mut parts = kv.splitn(2, '=');
        let k = parts.next().unwrap();
        let v = parts.next().unwrap_or("");

        info_map.insert(k.to_string(), v.to_string());
    }

    for vcf_alt in &vcf_alt_alleles {
        let vcf_alt_trimmed = vcf_alt.trim();
        let bundle = match ani.lookup(rec.chrom, pos, rec.ref_allele.trim(), vcf_alt_trimmed) {
            Some(b) => b,
            None => continue,
        };

        if debug {
            eprintln!(
                "[DEBUG] VCF ALTs: {:?}, DB ALTs: {:?}",
                vcf_alt_alleles,
                bundle.alt.split(',').collect::<Vec<_>>()
            );
        }

        let allele_map = bundle.match_alleles(&vcf_alt_alleles);

        if debug {
            eprintln!("[DEBUG] Allele map: {:?}", allele_map);
        }

        if let Some(id) = &bundle.id {
            if updated_id == "." || updated_id.is_empty() {
                updated_id = id.clone();
            } else if !rec.id.contains(id.as_str()) {
                updated_id = format!("{};{}", updated_id, id);
            }
        }

        if let Some(filt) = &bundle.filter {
            if updated_filter == "." || updated_filter.is_empty() {
                updated_filter = filt.clone();
            }
        }

        for field in &bundle.info {
            let key = &field.key;

            let Some(field_number) = field_meta.get(key.as_str()).copied() else {
                if debug {
                    eprintln!(
                        "[DEBUG] No metadata for field {}, using fallback (values.len={})",
                        key,
                        field.values.len()
                    );
                }

                if field.values.is_empty() {
                    info_map.insert(key.clone(), String::new());
                } else {
                    info_map.insert(key.clone(), field.values.join(","));
                }
                continue;
            };

            if debug {
                eprintln!(
                    "[DEBUG] Processing field {} with Number={:?}, values.len={}, allele_map.len={}",
                    key,
                    field_number,
                    field.values.len(),
                    allele_map.len()
                );
            }

            if field.values.is_empty() {
                match field_number {
                    FieldNumber::Zero => {
                        info_map.insert(key.clone(), String::new());
                    }
                    FieldNumber::A => {
                        let missing: Vec<String> =
                            vcf_alt_alleles.iter().map(|_| ".".to_string()).collect();
                        info_map.insert(key.clone(), missing.join(","));
                    }
                    FieldNumber::R => {
                        let mut missing = vec![".".to_string()];
                        for _ in &vcf_alt_alleles {
                            missing.push(".".to_string());
                        }
                        info_map.insert(key.clone(), missing.join(","));
                    }
                    _ => {}
                }
                continue;
            }

            match field_number {
                FieldNumber::Zero => {
                    info_map.insert(key.clone(), String::new());
                }

                FieldNumber::One => {
                    if !field.values.is_empty() {
                        info_map.insert(key.clone(), field.values[0].clone());
                    }
                }

                FieldNumber::Many => {
                    info_map.insert(key.clone(), field.values.join(","));
                }

                FieldNumber::A | FieldNumber::R => {
                    let remapped = bundle.remap_field_values(&field, &allele_map);

                    if info_map.contains_key(key.as_str()) {
                        let existing_val = info_map.get(key.as_str()).unwrap();
                        let existing_parts: Vec<&str> = existing_val.split(',').collect();
                        let field_type = infer_field_type(key);
                        let merged = merge_allele_values(
                            &existing_parts,
                            &remapped,
                            field_number,
                            field_type,
                        );
                        info_map.insert(key.clone(), merged.join(","));
                    } else {
                        info_map.insert(key.clone(), remapped.join(","));
                    }
                }

                FieldNumber::G => {
                    info_map.insert(key.clone(), field.values.join(","));
                }
            }
        }
    }

    let info_str = if info_map.is_empty() {
        ".".to_string()
    } else {
        let parts: Vec<String> = info_map
            .into_iter()
            .filter_map(|(k, v)| {
                if v.is_empty() {
                    Some(k)
                } else {
                    Some(format!("{}={}", k, v))
                }
            })
            .collect();
        if parts.is_empty() {
            ".".to_string()
        } else {
            parts.join(";")
        }
    };

    let rest = parser.rest();
    let mut fields = vec![
        rec.chrom.to_string(),
        rec.pos.to_string(),
        updated_id,
        rec.ref_allele.to_string(),
        rec.alt.to_string(),
        rec.qual.to_string(),
        updated_filter,
        info_str,
    ];

    if !rest.is_empty() {
        fields.push(rest.to_string());
    }

    fields.join("\t")
}

fn merge_annotation_headers(vcf_headers: &[String], _ani: &AniIndex) -> Result<Vec<String>> {
    Ok(vcf_headers.to_vec())
}
