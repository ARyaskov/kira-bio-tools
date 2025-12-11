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
use super::structs::bundle::{parse_info_field, AnnotationBundle};
use super::structs::*;
use crate::bgzf::BgzfWriter;
use crate::util::{detect_format, url_decode_info_value, VcfFormat};
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
    let mut field_order: Vec<String> = Vec::new();

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
                info_str
                    .split(';')
                    .filter_map(|kv| kv.split('=').next().map(|k| k.to_string()))
                    .collect::<Vec<_>>()
            })
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();

        field_meta = infer_field_metadata_from_data(&ani, &field_names);

        field_order = field_names;
    } else {
        field_order = field_meta.keys().cloned().collect();
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

        for field_name in field_names {
            let pattern = format!("{}=", field_name);
            if let Some(start) = info_str.find(&pattern) {
                let value_start = start + pattern.len();
                let rest = &info_str[value_start..];
                let value_end = rest.find(';').unwrap_or(rest.len());
                let value = &rest[..value_end];

                let decoded_value = url_decode_info_value(value);

                let alt_str = read_cstring(&ani.strings, entry.alt_ofs as usize);
                let num_alts = alt_str.split(',').count();
                let num_values = decoded_value.split(',').count();

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

fn choose_best_number(numbers: &[FieldNumber]) -> FieldNumber {
    let mut has_r = false;
    let mut has_a = false;
    let mut has_many = false;
    let mut has_one = false;

    for n in numbers {
        match n {
            FieldNumber::R => has_r = true,
            FieldNumber::A => has_a = true,
            FieldNumber::Many => has_many = true,
            FieldNumber::One => has_one = true,
            _ => {}
        }
    }

    if has_r {
        return FieldNumber::R;
    }
    if has_a {
        return FieldNumber::A;
    }
    if has_many {
        return FieldNumber::Many;
    }
    if has_one {
        return FieldNumber::One;
    }

    FieldNumber::One
}

fn read_cstring<'a>(data: &'a [u8], pos: usize) -> &'a str {
    if pos >= data.len() {
        return "";
    }
    let mut end = pos;
    while end < data.len() && data[end] != 0 {
        end += 1;
    }
    std::str::from_utf8(&data[pos..end]).unwrap_or("")
}

fn linear_search_ani<'a>(
    ani: &'a AniIndex,
    chr: &str,
    pos: u32,
    rf: &str,
) -> Option<AnnotationBundle<'a>> {
    use crate::chr_name_to_id;

    let chr_id = chr_name_to_id(chr)? as u8;

    for entry in &ani.entries {
        if entry.chr_id == chr_id && entry.pos == pos {
            let rf_str = read_cstring(&ani.strings, entry.ref_ofs as usize);

            if rf_str == rf {
                let alt_str = read_cstring(&ani.strings, entry.alt_ofs as usize);
                let id_str = read_cstring(&ani.strings, entry.id_ofs as usize);
                let qual_str = read_cstring(&ani.strings, entry.qual_ofs as usize);
                let filter_str = read_cstring(&ani.strings, entry.filter_ofs as usize);
                let info_str = read_cstring(&ani.strings, entry.info_ofs as usize);

                let info_fields = parse_info_field(info_str);

                return Some(AnnotationBundle {
                    alt: alt_str,
                    id: if id_str == "." || id_str.is_empty() {
                        None
                    } else {
                        Some(id_str)
                    },
                    qual: if qual_str == "." || qual_str.is_empty() {
                        None
                    } else {
                        Some(qual_str)
                    },
                    filter: if filter_str == "." || filter_str.is_empty() {
                        None
                    } else {
                        Some(filter_str)
                    },
                    info: info_fields,
                });
            }
        }
    }

    None
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

fn extract_info_key(line: &str) -> Option<String> {
    if let Some(start) = line.find("ID=") {
        let rest = &line[start + 3..];
        if let Some(end) = rest.find(',') {
            return Some(rest[..end].to_string());
        }
    }
    None
}

fn extract_info_number(line: &str) -> Option<FieldNumber> {
    if let Some(start) = line.find("Number=") {
        let rest = &line[start + 7..];
        if let Some(end) = rest.find(',') {
            let number_str = &rest[..end];
            return match number_str {
                "0" => Some(FieldNumber::Zero),
                "1" => Some(FieldNumber::One),
                "A" => Some(FieldNumber::A),
                "R" => Some(FieldNumber::R),
                "G" => Some(FieldNumber::G),
                "." => Some(FieldNumber::Many),
                _ => {
                    if number_str.parse::<i32>().is_ok() {
                        Some(FieldNumber::Many)
                    } else {
                        None
                    }
                }
            };
        }
    }
    None
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

    let bundle = match ani.lookup(rec.chrom, pos, rec.ref_allele) {
        Some(b) => b,
        None => {
            if debug {
                eprintln!(
                    "[DEBUG] Primary lookup failed for {}:{} {}, trying linear search...",
                    rec.chrom, pos, rec.ref_allele
                );
            }

            match linear_search_ani(ani, rec.chrom, pos, rec.ref_allele) {
                Some(b) => {
                    if debug {
                        eprintln!(
                            "[DEBUG] Linear search found match for {}:{} {} - DB ALT: {}",
                            rec.chrom, pos, rec.ref_allele, b.alt
                        );
                    }
                    b
                }
                None => {
                    if debug {
                        eprintln!(
                            "[DEBUG] No lookup match for {}:{} {}",
                            rec.chrom, pos, rec.ref_allele
                        );
                    }
                    return line.to_string();
                }
            }
        }
    };

    if debug && ani.lookup(rec.chrom, pos, rec.ref_allele).is_some() {
        eprintln!(
            "[DEBUG] Found match for {}:{} {} - DB ALT: {}",
            rec.chrom, pos, rec.ref_allele, bundle.alt
        );
        eprintln!("[DEBUG] Metadata contains {} fields", field_meta.len());
    }

    let vcf_alt_alleles: Vec<&str> = rec.alt.split(',').collect();
    let db_alt_alleles: Vec<&str> = bundle.alt.split(',').collect();

    let allele_map = match_alleles(&db_alt_alleles, &vcf_alt_alleles);

    if debug {
        eprintln!(
            "[DEBUG] VCF ALTs: {:?}, DB ALTs: {:?}, Map: {:?}",
            vcf_alt_alleles, db_alt_alleles, allele_map
        );
    }

    let mut updated_id = rec.id.to_string();
    let mut updated_qual = rec.qual.to_string();
    let mut updated_filter = rec.filter.to_string();
    let mut info_map: IndexMap<String, String> = IndexMap::new();

    if let Some(id) = bundle.id {
        if id != "." && !id.is_empty() {
            if rec.id == "." || rec.id.is_empty() {
                updated_id = id.to_string();
            } else if !rec.id.contains(id) {
                updated_id = format!("{};{}", rec.id, id);
            }
        }
    }

    if let Some(qual) = bundle.qual {
        if qual != "." && !qual.is_empty() && (rec.qual == "." || rec.qual.is_empty()) {
            updated_qual = qual.to_string();
        }
    }

    if let Some(filt) = bundle.filter {
        if filt != "." && !filt.is_empty() && (rec.filter == "." || rec.filter.is_empty()) {
            updated_filter = filt.to_string();
        }
    }

    for kv in rec.info.split(';') {
        if kv.is_empty() || kv == "." {
            continue;
        }
        let mut parts = kv.splitn(2, '=');
        let k = parts.next().unwrap();
        let v = parts.next().unwrap_or("");

        info_map.insert(k.to_string(), v.to_string());
    }

    for field in &bundle.info {
        let key = field.key;

        let Some(field_number) = field_meta.get(key).copied() else {
            if debug {
                eprintln!(
                    "[DEBUG] No metadata for field {}, using fallback (values.len={})",
                    key,
                    field.values.len()
                );
            }

            if field.values.is_empty() {
                info_map.insert(key.to_string(), String::new());
            } else {
                let decoded: Vec<String> = field
                    .values
                    .iter()
                    .map(|v| url_decode_info_value(v))
                    .collect();
                info_map.insert(key.to_string(), decoded.join(","));
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
                    info_map.insert(key.to_string(), String::new());
                }
                FieldNumber::A => {
                    let missing: Vec<String> =
                        vcf_alt_alleles.iter().map(|_| ".".to_string()).collect();
                    info_map.insert(key.to_string(), missing.join(","));
                }
                FieldNumber::R => {
                    let mut missing = vec![".".to_string()];
                    for _ in &vcf_alt_alleles {
                        missing.push(".".to_string());
                    }
                    info_map.insert(key.to_string(), missing.join(","));
                }
                _ => {}
            }
            continue;
        }

        match field_number {
            FieldNumber::Zero => {
                info_map.insert(key.to_string(), String::new());
            }

            FieldNumber::One => {
                if !field.values.is_empty() {
                    let decoded = url_decode_info_value(field.values[0]);
                    info_map.insert(key.to_string(), decoded);
                }
            }

            FieldNumber::Many => {
                let decoded: Vec<String> = field
                    .values
                    .iter()
                    .map(|v| url_decode_info_value(v))
                    .collect();
                info_map.insert(key.to_string(), decoded.join(","));
            }

            FieldNumber::A | FieldNumber::R => {
                let remapped = remap_field_values(&field.values, field_number, &allele_map);

                let final_values = if let Some(existing) = info_map.shift_remove(key) {
                    let existing_vals: Vec<&str> = existing.split(',').collect();
                    let field_type = infer_field_type(key);
                    merge_allele_values(&existing_vals, &remapped, field_number, field_type)
                } else {
                    remapped
                };

                let joined = final_values.join(",");

                if debug {
                    eprintln!(
                        "[DEBUG] Field {} remapped: {:?} -> {:?}",
                        key, field.values, final_values
                    );
                }

                info_map.insert(key.to_string(), joined);
            }

            FieldNumber::G => {
                let decoded: Vec<String> = field
                    .values
                    .iter()
                    .map(|v| url_decode_info_value(v))
                    .collect();
                info_map.insert(key.to_string(), decoded.join(","));
            }
        }
    }

    let info_str = if info_map.is_empty() {
        ".".to_string()
    } else {
        let mut parts: Vec<String> = Vec::new();
        for (k, v) in info_map {
            if v.is_empty() {
                parts.push(k);
            } else {
                parts.push(format!("{}={}", k, v));
            }
        }
        parts.join(";")
    };

    let rest = parser.rest();
    let mut fields = vec![
        rec.chrom.to_string(),
        rec.pos.to_string(),
        updated_id,
        rec.ref_allele.to_string(),
        rec.alt.to_string(),
        updated_qual,
        updated_filter,
        info_str,
    ];

    if !rest.is_empty() && rest != "." {
        fields.push(rest.to_string());
    }

    fields.join("\t")
}

fn match_alleles(db_alts: &[&str], vcf_alts: &[&str]) -> Vec<Option<usize>> {
    let mut db_map: HashMap<&str, Vec<usize>> = HashMap::new();
    for (idx, alt) in db_alts.iter().enumerate() {
        db_map.entry(*alt).or_insert_with(Vec::new).push(idx);
    }

    vcf_alts
        .iter()
        .map(|vcf_alt| {
            db_map
                .get(vcf_alt)
                .and_then(|indices| indices.first().copied())
        })
        .collect()
}

fn remap_field_values(
    values: &[&str],
    number: FieldNumber,
    allele_map: &[Option<usize>],
) -> Vec<String> {
    let decoded_values: Vec<String> = values
        .iter()
        .flat_map(|v| {
            let decoded = url_decode_info_value(v);
            decoded
                .split(',')
                .map(|s| s.to_string())
                .collect::<Vec<_>>()
        })
        .collect();

    match number {
        FieldNumber::A => allele_map
            .iter()
            .map(|opt_idx| {
                opt_idx
                    .and_then(|idx| decoded_values.get(idx))
                    .cloned()
                    .unwrap_or_else(|| ".".to_string())
            })
            .collect(),
        FieldNumber::R => {
            let mut result = Vec::new();

            if let Some(ref_val) = decoded_values.first() {
                result.push(ref_val.clone());
            } else {
                result.push(".".to_string());
            }

            for opt_idx in allele_map {
                if let Some(idx) = opt_idx {
                    if let Some(val) = decoded_values.get(idx + 1) {
                        result.push(val.clone());
                    } else {
                        result.push(".".to_string());
                    }
                } else {
                    result.push(".".to_string());
                }
            }

            result
        }
        _ => decoded_values,
    }
}

fn merge_annotation_headers(vcf_headers: &[String], _ani: &AniIndex) -> Result<Vec<String>> {
    Ok(vcf_headers.to_vec())
}
