pub mod annotation;
pub mod column_spec;
pub mod field_metadata;
pub mod merge_info;
pub mod merge_info_helpers;
pub mod threads;
pub mod vcf_output;
pub mod vcf_parsing;

pub use annotation::*;
pub use column_spec::*;
pub use field_metadata::*;
pub use merge_info::*;
pub use threads::*;
pub use vcf_output::*;
pub use vcf_parsing::*;

use anyhow::Result;
use crossbeam_channel::bounded;
use std::path::Path;
use std::sync::Arc;
use std::thread;
use std::time::Instant;

use super::constants::*;
use super::reader::{StreamingVcfReader, VcfAnnotationReader};
use crate::annotate::structs::ani::AniIndex;
use crate::annotate::structs::bundle::FieldNumber;
use crate::util::detect_format;

fn split_mapped_ref(raw: &str) -> (&str, &str) {
    if let Some((src, dst)) = raw.split_once("=>") {
        (src, dst)
    } else {
        (raw, raw)
    }
}

fn is_format_ref(raw: &str) -> bool {
    let upper = raw.to_ascii_uppercase();
    upper == "FMT" || upper == "FORMAT" || upper.starts_with("FMT/") || upper.starts_with("FORMAT/")
}

fn info_id_from_ref(raw: &str) -> String {
    raw.strip_prefix("INFO/").unwrap_or(raw).to_string()
}

fn format_id_from_ref(raw: &str) -> String {
    raw.strip_prefix("FMT/")
        .or_else(|| raw.strip_prefix("FORMAT/"))
        .unwrap_or(raw)
        .to_string()
}

fn rewrite_header_id(line: &str, dst_id: &str) -> String {
    if let Some(start) = line.find("ID=") {
        let id_start = start + 3;
        if let Some(rel_end) = line[id_start..].find([',', '>']) {
            let id_end = id_start + rel_end;
            let mut out = String::with_capacity(line.len() + dst_id.len());
            out.push_str(&line[..id_start]);
            out.push_str(dst_id);
            out.push_str(&line[id_end..]);
            return out;
        }
    }
    line.to_string()
}

fn synthetic_info_header_for_fixed_source(src_ref: &str, dst_id: &str) -> Option<String> {
    let src_upper = src_ref.to_ascii_uppercase();
    if src_upper == "ID" {
        return Some(format!(
            "##INFO=<ID={dst_id},Number=1,Type=String,Description=\"Transferred ID column\">"
        ));
    }
    if src_upper == "FILTER" {
        return Some(format!(
            "##INFO=<ID={dst_id},Number=1,Type=String,Description=\"Transferred FILTER column\">"
        ));
    }
    if src_upper == "QUAL" {
        return Some(format!(
            "##INFO=<ID={dst_id},Number=1,Type=Float,Description=\"Transferred QUAL column\">"
        ));
    }
    None
}

pub fn annotate_vcf_ani_v2(
    db: &Path,
    input: &Path,
    output: &Path,
    columns: &[String],
    bgzf_level: Option<u32>,
    mmap_output: bool,
    mmap_no_flush: bool,
    ram_output: bool,
    ram_max_mb: u32,
) -> Result<()> {
    let timing = std::env::var("KIRA_BT_TIMING").is_ok();
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok() || timing;
    let start = Instant::now();

    let mut column_specs = ColumnSpec::parse_all(columns);
    let info_overwrite_all = column_specs
        .iter()
        .any(|c| c.key.eq_ignore_ascii_case("INFO") && c.mode.replace_all);
    let format_overwrite_all = column_specs.iter().any(|c| {
        (c.key.eq_ignore_ascii_case("FMT") || c.key.eq_ignore_ascii_case("FORMAT"))
            && c.mode.replace_all
    });

    let num_threads = rayon::current_num_threads() / 2;
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

    let field_meta = load_and_infer_metadata(&ani, debug)?;
    let ani_headers = iter_ani_header_lines(&ani);

    column_specs = expand_column_specs(&column_specs, &ani_headers, &field_meta);

    if debug {
        eprintln!(
            "[annotate] Column specs: {:?}",
            column_specs
                .iter()
                .map(|c| format!("{}{}", c.mode, c.key))
                .collect::<Vec<_>>()
        );
    }

    let input_format = detect_format(input)?;
    let output_ext = output.extension().and_then(|s| s.to_str()).unwrap_or("");
    let output_wants_bgzf = matches!(output_ext, "gz" | "bgz" | "bgzf");
    let use_bgzf = output_wants_bgzf;
    if debug {
        eprintln!(
            "[annotate] Input: {:?}, Output path: {:?}, ext: {:?}, wants_bgzf: {}, use_bgzf: {}",
            input_format, output, output_ext, output_wants_bgzf, use_bgzf
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

    let merged_headers = merge_annotation_headers(&headers, &ani_headers, &column_specs)?;

    let input_samples = extract_samples_from_headers(&headers);
    let db_samples = extract_samples_from_headers(&ani_headers);
    let sample_map = build_sample_map(&input_samples, &db_samples);

    let (read_tx, read_rx) = bounded::<Vec<String>>(CHANNEL_DEPTH);
    let (work_tx, work_rx) = bounded::<Vec<String>>(CHANNEL_DEPTH);

    let ani_clone = Arc::new(ani);
    let field_meta_clone = Arc::new(field_meta);
    let column_specs_arc = Arc::new(column_specs);
    let sample_map_arc = Arc::new(sample_map);
    let bundle_acc = Arc::new(BundleTimingAccum::new());
    let bundle_acc_worker = bundle_acc.clone();

    let worker = thread::spawn(move || {
        worker_thread(
            read_rx,
            work_tx,
            ani_clone.clone(),
            field_meta_clone.clone(),
            column_specs_arc.clone(),
            sample_map_arc.clone(),
            info_overwrite_all,
            format_overwrite_all,
            num_threads,
            timing,
            bundle_acc_worker,
        )
    });

    let output_clone = output.to_path_buf();
    let writer = thread::spawn(move || {
        writer_thread(
            work_rx,
            merged_headers,
            &output_clone,
            use_bgzf,
            bgzf_level,
            mmap_output,
            mmap_no_flush,
            ram_output,
            ram_max_mb,
            timing,
            "cpu",
        )
    });

    read_batches(&mut reader, read_tx, timing)?;

    worker.join().unwrap()?;
    writer.join().unwrap()?;

    if debug {
        eprintln!(
            "[annotate] Total time: {:.3}s",
            start.elapsed().as_secs_f64()
        );
    }
    if timing {
        let (r, i, o, s) = bundle_acc.snapshot_seconds();
        eprintln!(
            "[annotate] bundle_read: {:.3}s, bundle_info: {:.3}s, bundle_optional: {:.3}s, bundle_samples: {:.3}s",
            r, i, o, s
        );
    }

    Ok(())
}

fn read_batches(
    reader: &mut StreamingVcfReader,
    read_tx: crossbeam_channel::Sender<Vec<String>>,
    timing: bool,
) -> Result<()> {
    let mut batch = Vec::with_capacity(BATCH_SIZE);
    let mut total_lines = 0usize;
    let mut last_report = Instant::now();

    while let Some(line) = reader.read_line()? {
        batch.push(line);

        if batch.len() >= BATCH_SIZE {
            if timing {
                total_lines += batch.len();
            }
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
        let _ = read_tx.send(batch);
    }

    Ok(())
}

pub(crate) fn merge_annotation_headers(
    vcf_headers: &[String],
    ani_headers: &[String],
    column_specs: &[ColumnSpec],
) -> Result<Vec<String>> {
    let mut info_needed = Vec::new();
    let mut info_needed_seen = std::collections::HashSet::new();
    let mut format_needed = Vec::new();
    let mut format_needed_seen = std::collections::HashSet::new();
    let mut info_src_to_dst: Vec<(String, String, String)> = Vec::new();
    let mut format_src_to_dst: Vec<(String, String)> = Vec::new();
    let mut need_format = false;
    let mut need_filter = false;

    for col in column_specs {
        let (src_ref, dst_ref) = split_mapped_ref(&col.key);
        match dst_ref.to_uppercase().as_str() {
            "ID" | "QUAL" => {}
            "FILTER" => need_filter = true,
            "FMT" | "FORMAT" => need_format = true,
            _ => {
                if is_format_ref(dst_ref) {
                    need_format = true;
                    let src_id = format_id_from_ref(src_ref);
                    let dst_id = format_id_from_ref(dst_ref);
                    if format_needed_seen.insert(dst_id.clone()) {
                        format_needed.push(dst_id.clone());
                    }
                    format_src_to_dst.push((src_id, dst_id));
                } else {
                    let src_id = info_id_from_ref(src_ref);
                    let dst_id = info_id_from_ref(dst_ref);
                    if info_needed_seen.insert(dst_id.clone()) {
                        info_needed.push(dst_id.clone());
                    }
                    info_src_to_dst.push((src_ref.to_string(), src_id, dst_id));
                }
            }
        }
    }

    let mut input_info_ids = std::collections::HashSet::new();
    let mut input_format_ids = std::collections::HashSet::new();
    let mut input_filter_ids = std::collections::HashSet::new();

    for h in vcf_headers {
        if h.starts_with("##INFO=") {
            if let Some(id) = extract_header_id(h) {
                input_info_ids.insert(id);
            }
        } else if h.starts_with("##FORMAT=") {
            if let Some(id) = extract_header_id(h) {
                input_format_ids.insert(id);
            }
        } else if h.starts_with("##FILTER=") {
            if let Some(id) = extract_header_id(h) {
                input_filter_ids.insert(id);
            }
        }
    }

    let mut info_headers = Vec::new();
    let mut format_headers = Vec::new();
    let mut filter_headers = Vec::new();

    for h in ani_headers {
        if h.starts_with("##INFO=") {
            info_headers.push(h.clone());
        } else if h.starts_with("##FORMAT=") {
            format_headers.push(h.clone());
        } else if h.starts_with("##FILTER=") {
            filter_headers.push(h.clone());
        }
    }

    let mut extra = Vec::new();

    if !info_needed.is_empty() {
        let mut info_header_by_id = std::collections::HashMap::new();
        for h in &info_headers {
            if let Some(id) = extract_header_id(h) {
                info_header_by_id.insert(id, h.clone());
            }
        }

        for dst_id in &info_needed {
            if input_info_ids.contains(dst_id) {
                continue;
            }
            let mapped_src = info_src_to_dst
                .iter()
                .find(|(_, _, dst)| dst == dst_id)
                .map(|(src_ref, src_id, _)| (src_ref.clone(), src_id.clone()));
            if let Some((src_ref, src_id)) = mapped_src {
                if let Some(synth) = synthetic_info_header_for_fixed_source(&src_ref, dst_id) {
                    extra.push(synth);
                    continue;
                }
                if let Some(src_line) = info_header_by_id.get(&src_id) {
                    extra.push(rewrite_header_id(src_line, dst_id));
                    continue;
                }
            }
            if let Some(line) = info_header_by_id.get(dst_id) {
                extra.push(line.clone());
            }
        }
    }

    if need_format {
        let mut format_header_by_id = std::collections::HashMap::new();
        for h in &format_headers {
            if let Some(id) = extract_header_id(h) {
                format_header_by_id.insert(id, h.clone());
            }
        }

        for dst_id in &format_needed {
            if input_format_ids.contains(dst_id) {
                continue;
            }
            let mapped_src = format_src_to_dst
                .iter()
                .find(|(_, dst)| dst == dst_id)
                .map(|(src, _)| src.clone());
            if let Some(src_id) = mapped_src {
                if let Some(src_line) = format_header_by_id.get(&src_id) {
                    extra.push(rewrite_header_id(src_line, dst_id));
                    continue;
                }
            }
            if let Some(line) = format_header_by_id.get(dst_id) {
                extra.push(line.clone());
            }
        }

        for h in format_headers {
            if let Some(id) = extract_header_id(&h) {
                if format_needed_seen.contains(&id) {
                    continue;
                }
                if !input_format_ids.contains(&id)
                    && column_specs.iter().any(|c| {
                        split_mapped_ref(&c.key).1.eq_ignore_ascii_case("FMT")
                            || split_mapped_ref(&c.key).1.eq_ignore_ascii_case("FORMAT")
                    })
                {
                    extra.push(h);
                }
            }
        }
    }

    if need_filter {
        for h in filter_headers {
            if let Some(id) = extract_header_id(&h) {
                if !input_filter_ids.contains(&id) {
                    extra.push(h);
                }
            }
        }

        if !input_filter_ids.contains("PASS") {
            extra.push("##FILTER=<ID=PASS,Description=\"All filters passed\">".to_string());
        }
    }

    let mut merged = Vec::new();
    let mut chrom_line = None;
    for h in vcf_headers {
        if h.starts_with("#CHROM") {
            chrom_line = Some(h.clone());
        } else {
            merged.push(h.clone());
        }
    }

    merged.extend(extra);

    if let Some(chrom) = chrom_line {
        merged.push(chrom);
    }

    Ok(merged)
}

pub(crate) fn expand_column_specs(
    column_specs: &[ColumnSpec],
    ani_headers: &[String],
    field_meta: &std::collections::HashMap<String, FieldNumber>,
) -> Vec<ColumnSpec> {
    let mut expanded = Vec::new();
    let mut seen_info = std::collections::HashSet::new();
    let mut info_all_mode = None;

    for col in column_specs {
        if col.key.eq_ignore_ascii_case("INFO") {
            info_all_mode = Some(col.mode);
            continue;
        }

        expanded.push(col.clone());

        let (_, dst_ref) = split_mapped_ref(&col.key);
        match dst_ref.to_uppercase().as_str() {
            "ID" | "QUAL" | "FILTER" | "FMT" | "FORMAT" => {}
            _ => {
                if !is_format_ref(dst_ref) {
                    seen_info.insert(info_id_from_ref(dst_ref));
                }
            }
        }
    }

    if let Some(mode) = info_all_mode {
        let mut info_ids = Vec::new();
        for h in ani_headers {
            if h.starts_with("##INFO=") {
                if let Some(id) = extract_header_id(h) {
                    info_ids.push(id);
                }
            }
        }

        if info_ids.is_empty() {
            info_ids = field_meta.keys().cloned().collect();
            info_ids.sort();
        }

        for id in info_ids {
            if seen_info.insert(id.clone()) {
                expanded.push(ColumnSpec {
                    key: id.clone(),
                    dst_key: id,
                    mode,
                });
            }
        }
    }

    expanded
}

fn extract_header_id(line: &str) -> Option<String> {
    crate::util::extract_info_key(line)
}

pub(crate) fn extract_samples_from_headers(headers: &[String]) -> Vec<String> {
    for h in headers {
        if h.starts_with("#CHROM") {
            let parts: Vec<&str> = h.trim().split('\t').collect();
            if parts.len() > 9 {
                return parts[9..].iter().map(|s| s.to_string()).collect();
            }
            break;
        }
    }
    Vec::new()
}

pub(crate) fn build_sample_map(
    input_samples: &[String],
    db_samples: &[String],
) -> Vec<Option<usize>> {
    let mut map = Vec::with_capacity(input_samples.len());
    let mut db_index = std::collections::HashMap::new();
    for (i, name) in db_samples.iter().enumerate() {
        db_index.insert(name, i);
    }

    for name in input_samples {
        map.push(db_index.get(name).copied());
    }

    map
}

#[cfg(test)]
mod tests {
    use super::build_sample_map;

    #[test]
    fn test_build_sample_map_by_name() {
        let input = vec!["A".to_string(), "B".to_string(), "C".to_string()];
        let db = vec!["B".to_string(), "A".to_string()];
        let map = build_sample_map(&input, &db);
        assert_eq!(map, vec![Some(1), Some(0), None]);
    }
}
