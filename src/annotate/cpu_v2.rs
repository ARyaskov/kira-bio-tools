use anyhow::Result;
use crossbeam_channel::{bounded, Receiver, Sender};
use indexmap::IndexMap;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::Arc;
use std::thread;
use std::time::Instant;

use super::constants::*;
use super::reader::{StreamingVcfReader, VcfAnnotationReader};
use super::structs::annotate_mode::AnnotateMode;
use super::structs::bundle::{parse_info_field, AnnotationBundle, FieldNumber};
use super::structs::*;
use crate::bgzf::BgzfWriter;
use crate::util::{
    choose_best_number, detect_format, extract_info_key, extract_info_number, read_cstring,
    url_decode_info_value, VcfFormat,
};
use crate::vcf::VcfParser;

#[derive(Debug, Clone)]
pub struct ColumnSpec {
    pub key: String,
    pub dst_key: String,
    pub mode: AnnotateMode,
}

impl ColumnSpec {
    pub fn parse(spec: &str) -> Self {
        let (mode, rest) = AnnotateMode::parse(spec);

        let (src_key, dst_key) = if rest.contains(":=") {
            let parts: Vec<&str> = rest.splitn(2, ":=").collect();
            if parts.len() == 2 {
                let dst = parts[0].strip_prefix("INFO/").unwrap_or(parts[0]);
                let src = parts[1].strip_prefix("INFO/").unwrap_or(parts[1]);
                (src.to_string(), dst.to_string())
            } else {
                let key = rest.strip_prefix("INFO/").unwrap_or(rest).to_string();
                (key.clone(), key)
            }
        } else {
            let key = rest.strip_prefix("INFO/").unwrap_or(rest).to_string();
            (key.clone(), key)
        };

        Self {
            key: src_key,
            dst_key,
            mode,
        }
    }

    pub fn parse_all(columns: &[String]) -> Vec<Self> {
        columns
            .iter()
            .filter(|c| {
                let upper = c.to_uppercase();
                !upper.starts_with("CHROM")
                    && !upper.starts_with("POS")
                    && !upper.starts_with("REF")
                    && !upper.starts_with("ALT")
                    && !upper.starts_with("FROM")
                    && !upper.starts_with("TO")
                    && !upper.starts_with("BEG")
                    && !upper.starts_with("END")
                    && *c != "-"
            })
            .map(|c| Self::parse(c))
            .collect()
    }
}

pub fn annotate_vcf_ani_v2(
    db: &Path,
    input: &Path,
    output: &Path,
    columns: &[String],
) -> Result<()> {
    let timing = std::env::var("KIRA_BT_TIMING").is_ok();
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok() || timing;
    let start = Instant::now();

    let column_specs = ColumnSpec::parse_all(columns);

    let field_order: Vec<String> = column_specs.iter().map(|c| c.dst_key.clone()).collect();

    let num_threads = rayon::current_num_threads() / 4;
    if debug {
        eprintln!("[annotate] Using {} CPU threads", num_threads);
        eprintln!("[annotate] Batch size: {} lines", BATCH_SIZE);
        eprintln!(
            "[annotate] Column specs: {:?}",
            column_specs
                .iter()
                .map(|c| format!("{}{}", c.mode, c.key))
                .collect::<Vec<_>>()
        );
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
            .collect::<HashSet<_>>()
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

    let ani_clone = Arc::new(ani);
    let field_meta_clone = Arc::new(field_meta);
    let field_order_arc = Arc::new(field_order);
    let column_specs_arc = Arc::new(column_specs);

    let ani_worker = ani_clone.clone();
    let field_meta_worker = field_meta_clone.clone();
    let field_order_worker = field_order_arc.clone();
    let column_specs_worker = column_specs_arc.clone();

    let worker = thread::spawn(move || {
        worker_thread(
            read_rx,
            work_tx,
            ani_worker,
            field_meta_worker,
            field_order_worker,
            column_specs_worker,
            num_threads,
        )
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
            if read_tx
                .send(std::mem::replace(
                    &mut batch,
                    Vec::with_capacity(BATCH_SIZE),
                ))
                .is_err()
            {
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
                    .or_default()
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
    for (field_name, numbers) in candidates {
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

fn is_missing_value(val: &str) -> bool {
    val.is_empty() || val == "."
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
    ani: Arc<AniIndex>,
    field_meta: Arc<HashMap<String, FieldNumber>>,
    field_order: Arc<Vec<String>>,
    column_specs: Arc<Vec<ColumnSpec>>,
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
                .map(|line| annotate_line(line, &ani, &field_meta, &field_order, &column_specs))
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

fn annotate_line(
    line: &str,
    ani: &AniIndex,
    field_meta: &HashMap<String, FieldNumber>,
    field_order: &[String],
    column_specs: &[ColumnSpec],
) -> String {
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok();

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

    let mut found_bundles: Vec<(usize, AnnotationBundle)> = Vec::new();

    for (vcf_idx, vcf_alt) in vcf_alt_alleles.iter().enumerate() {
        let vcf_alt_trimmed = vcf_alt.trim();
        if let Some(bundle) = ani.lookup(rec.chrom, pos, rec.ref_allele.trim(), vcf_alt_trimmed) {
            found_bundles.push((vcf_idx, bundle));
        }
    }

    let multiallelic_bundle: Option<AnnotationBundle> = if found_bundles.is_empty() {
        ani.lookup_any_alt(rec.chrom, pos, rec.ref_allele.trim())
    } else {
        None
    };

    if found_bundles.is_empty() && multiallelic_bundle.is_none() {
        return line.to_string();
    }

    let column_spec_map: HashMap<&str, &ColumnSpec> =
        column_specs.iter().map(|c| (c.key.as_str(), c)).collect();

    for (_vcf_idx, bundle) in &found_bundles {
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
    }

    let mut per_field_values: HashMap<String, Vec<Option<String>>> = HashMap::new();
    let mut per_field_ref: HashMap<String, Option<String>> = HashMap::new();

    for (vcf_idx, bundle) in &found_bundles {
        for field in &bundle.info {
            let key = &field.key;

            let Some(field_number) = field_meta.get(key.as_str()).copied() else {
                if debug {
                    eprintln!("[DEBUG] No metadata for field {}, using fallback", key);
                }
                continue;
            };

            let entry = per_field_values
                .entry(key.clone())
                .or_insert_with(|| vec![None; vcf_alt_alleles.len()]);

            match field_number {
                FieldNumber::A => {
                    if let Some(val) = field.values.first() {
                        entry[*vcf_idx] = Some(val.clone());
                    }
                }
                FieldNumber::R => {
                    if let Some(ref_val) = field.values.first() {
                        per_field_ref
                            .entry(key.clone())
                            .or_insert_with(|| Some(ref_val.clone()));
                    }
                    entry[*vcf_idx] = field.values.get(1).cloned();
                }
                _ => {}
            }
        }
    }

    if let Some(ref bundle) = multiallelic_bundle {
        let db_alts: Vec<&str> = bundle.alt.split(',').collect();
        let mut alt_to_db_idx: HashMap<&str, usize> = HashMap::new();
        for (i, alt) in db_alts.iter().enumerate() {
            alt_to_db_idx.insert(*alt, i);
        }

        if let Some(id) = &bundle.id {
            if updated_id == "." || updated_id.is_empty() {
                updated_id = id.clone();
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
                continue;
            };

            let entry = per_field_values
                .entry(key.clone())
                .or_insert_with(|| vec![None; vcf_alt_alleles.len()]);

            match field_number {
                FieldNumber::A => {
                    for (vcf_idx, vcf_alt) in vcf_alt_alleles.iter().enumerate() {
                        if let Some(&db_idx) = alt_to_db_idx.get(*vcf_alt) {
                            if let Some(val) = field.values.get(db_idx) {
                                entry[vcf_idx] = Some(val.clone());
                            }
                        }
                    }
                }
                FieldNumber::R => {
                    if let Some(ref_val) = field.values.first() {
                        per_field_ref
                            .entry(key.clone())
                            .or_insert_with(|| Some(ref_val.clone()));
                    }
                    for (vcf_idx, vcf_alt) in vcf_alt_alleles.iter().enumerate() {
                        if let Some(&db_idx) = alt_to_db_idx.get(*vcf_alt) {
                            if let Some(val) = field.values.get(db_idx + 1) {
                                entry[vcf_idx] = Some(val.clone());
                            }
                        }
                    }
                }
                _ => {}
            }
        }
    }

    for (key, vcf_values) in per_field_values {
        let field_number = match field_meta.get(key.as_str()).copied() {
            Some(n) => n,
            None => continue,
        };

        let spec = column_spec_map.get(key.as_str()).copied();
        let mode = spec
            .map(|s| s.mode)
            .unwrap_or_else(AnnotateMode::default_mode);
        let dst_key = spec.map(|s| s.dst_key.as_str()).unwrap_or(&key);

        let vcf_has_field = info_map.contains_key(dst_key);
        let existing_val = info_map.get(dst_key).map(|s| s.as_str()).unwrap_or("");
        let existing_parts: Vec<&str> = if !existing_val.is_empty() && existing_val != "." {
            existing_val.split(',').collect()
        } else {
            vec![]
        };

        let field_type = infer_field_type(&key);
        let is_integer = field_type == "Integer";

        let final_values: Vec<String> = match field_number {
            FieldNumber::A => {
                if mode.replace_missing && is_integer {
                    vcf_values
                        .iter()
                        .map(|opt| opt.clone().unwrap_or_else(|| ".".to_string()))
                        .collect()
                } else if mode.replace_missing {
                    vcf_values
                        .iter()
                        .enumerate()
                        .map(|(i, opt_val)| {
                            let vcf_val = existing_parts.get(i).copied().unwrap_or(".");
                            if !is_missing_value(vcf_val) {
                                vcf_val.to_string()
                            } else if let Some(val) = opt_val {
                                if !is_missing_value(val) || mode.carry_over_missing {
                                    val.clone()
                                } else {
                                    ".".to_string()
                                }
                            } else {
                                ".".to_string()
                            }
                        })
                        .collect()
                } else if mode.should_append() && vcf_has_field {
                    let mut result = existing_parts
                        .iter()
                        .map(|s| s.to_string())
                        .collect::<Vec<_>>();
                    for (i, opt_val) in vcf_values.iter().enumerate() {
                        if let Some(val) = opt_val {
                            if !is_missing_value(val) {
                                if i < result.len() {
                                    if is_missing_value(&result[i]) {
                                        result[i] = val.clone();
                                    } else {
                                        result[i] = format!("{},{}", result[i], val);
                                    }
                                } else {
                                    result.push(val.clone());
                                }
                            }
                        }
                    }
                    result
                } else if mode.replace_non_missing {
                    if !vcf_has_field {
                        continue;
                    }
                    vcf_values
                        .iter()
                        .enumerate()
                        .map(|(i, opt_val)| {
                            let vcf_val = existing_parts.get(i).copied().unwrap_or(".");
                            if is_missing_value(vcf_val) {
                                vcf_val.to_string()
                            } else if let Some(val) = opt_val {
                                if !is_missing_value(val) || mode.carry_over_missing {
                                    val.clone()
                                } else {
                                    vcf_val.to_string()
                                }
                            } else {
                                vcf_val.to_string()
                            }
                        })
                        .collect()
                } else {
                    vcf_values
                        .iter()
                        .enumerate()
                        .map(|(i, opt_val)| {
                            if let Some(val) = opt_val {
                                if !is_missing_value(val) || mode.carry_over_missing {
                                    val.clone()
                                } else {
                                    existing_parts
                                        .get(i)
                                        .map(|s| s.to_string())
                                        .unwrap_or_else(|| ".".to_string())
                                }
                            } else {
                                existing_parts
                                    .get(i)
                                    .map(|s| s.to_string())
                                    .unwrap_or_else(|| ".".to_string())
                            }
                        })
                        .collect()
                }
            }
            FieldNumber::R => {
                let mut result = vec![];

                let db_ref_val = per_field_ref
                    .get(&key)
                    .and_then(|opt| opt.as_ref())
                    .cloned();
                let vcf_ref_val = existing_parts.first().copied().unwrap_or(".");

                let ref_val = if mode.replace_missing && is_integer {
                    db_ref_val.unwrap_or_else(|| ".".to_string())
                } else if mode.replace_missing {
                    if !is_missing_value(vcf_ref_val) {
                        vcf_ref_val.to_string()
                    } else {
                        db_ref_val.unwrap_or_else(|| ".".to_string())
                    }
                } else if mode.replace_non_missing {
                    if is_missing_value(vcf_ref_val) {
                        vcf_ref_val.to_string()
                    } else {
                        db_ref_val.unwrap_or_else(|| vcf_ref_val.to_string())
                    }
                } else {
                    db_ref_val.unwrap_or_else(|| vcf_ref_val.to_string())
                };
                result.push(ref_val);

                if mode.replace_missing && is_integer {
                    for opt_val in vcf_values.iter() {
                        result.push(opt_val.clone().unwrap_or_else(|| ".".to_string()));
                    }
                } else if mode.replace_missing {
                    for (i, opt_val) in vcf_values.iter().enumerate() {
                        let vcf_val = existing_parts.get(i + 1).copied().unwrap_or(".");
                        if !is_missing_value(vcf_val) {
                            result.push(vcf_val.to_string());
                        } else if let Some(val) = opt_val {
                            if !is_missing_value(val) || mode.carry_over_missing {
                                result.push(val.clone());
                            } else {
                                result.push(".".to_string());
                            }
                        } else {
                            result.push(".".to_string());
                        }
                    }
                } else if mode.should_append() && vcf_has_field {
                    for (i, opt_val) in vcf_values.iter().enumerate() {
                        let existing = existing_parts.get(i + 1).copied().unwrap_or(".");
                        if let Some(val) = opt_val {
                            if !is_missing_value(val) {
                                if is_missing_value(existing) {
                                    result.push(val.clone());
                                } else {
                                    result.push(format!("{},{}", existing, val));
                                }
                            } else {
                                result.push(existing.to_string());
                            }
                        } else {
                            result.push(existing.to_string());
                        }
                    }
                } else if mode.replace_non_missing {
                    if !vcf_has_field {
                        continue;
                    }
                    for (i, opt_val) in vcf_values.iter().enumerate() {
                        let vcf_val = existing_parts.get(i + 1).copied().unwrap_or(".");
                        if is_missing_value(vcf_val) {
                            result.push(vcf_val.to_string());
                        } else if let Some(val) = opt_val {
                            if !is_missing_value(val) || mode.carry_over_missing {
                                result.push(val.clone());
                            } else {
                                result.push(vcf_val.to_string());
                            }
                        } else {
                            result.push(vcf_val.to_string());
                        }
                    }
                } else {
                    for (i, opt_val) in vcf_values.iter().enumerate() {
                        if let Some(val) = opt_val {
                            if !is_missing_value(val) || mode.carry_over_missing {
                                result.push(val.clone());
                            } else {
                                result.push(
                                    existing_parts
                                        .get(i + 1)
                                        .map(|s| s.to_string())
                                        .unwrap_or_else(|| ".".to_string()),
                                );
                            }
                        } else {
                            result.push(
                                existing_parts
                                    .get(i + 1)
                                    .map(|s| s.to_string())
                                    .unwrap_or_else(|| ".".to_string()),
                            );
                        }
                    }
                }

                result
            }
            _ => vcf_values.iter().filter_map(|v| v.clone()).collect(),
        };

        info_map.insert(dst_key.to_string(), final_values.join(","));
    }

    let info_str = if info_map.is_empty() {
        ".".to_string()
    } else {
        let mut ordered_keys: Vec<String> = Vec::new();
        let mut unordered_keys: Vec<String> = Vec::new();

        for key in info_map.keys() {
            if field_order.contains(key) {
                ordered_keys.push(key.clone());
            } else {
                unordered_keys.push(key.clone());
            }
        }

        ordered_keys.sort_by_key(|k| {
            field_order
                .iter()
                .position(|f| f == k)
                .unwrap_or(usize::MAX)
        });
        unordered_keys.sort();

        let mut all_keys = ordered_keys;
        all_keys.extend(unordered_keys);

        let parts: Vec<String> = all_keys
            .into_iter()
            .filter_map(|k| {
                let v = info_map.get(&k)?;
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
