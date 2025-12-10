use anyhow::Result;
use crossbeam_channel::{bounded, Receiver, Sender};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;
use std::time::Instant;

use super::constants::*;
use super::reader::{StreamingVcfReader, VcfAnnotationReader};
use super::structs::*;
use crate::bgzf::BgzfWriter;
use crate::util::{clean_info_values, detect_format, url_decode_info_value, VcfFormat};
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
    let ani_worker = ani_clone.clone();
    let worker = thread::spawn(move || worker_thread(read_rx, work_tx, ani_worker, num_threads));

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

fn worker_thread(
    rx: Receiver<Vec<String>>,
    tx: Sender<Vec<String>>,
    ani: std::sync::Arc<AniIndex>,
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
                .map(|line| annotate_line(line, &ani))
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

fn annotate_line(line: &str, ani: &AniIndex) -> String {
    let mut parser = VcfParser::new(line);
    let Some(rec) = parser.parse_standard_fields() else {
        return line.to_string();
    };

    let pos = rec.pos.parse::<u32>().unwrap_or(0);

    let Some(bundle) = ani.lookup(rec.chrom, pos, rec.ref_allele) else {
        return line.to_string();
    };

    let src_alt_alleles: Vec<&str> = rec.alt.split(',').collect();

    let mut updated_id = rec.id.to_string();
    let mut updated_qual = rec.qual.to_string();
    let mut updated_filter = rec.filter.to_string();
    let mut info_map: HashMap<String, String> = HashMap::new();

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
        if kv.is_empty() {
            continue;
        }
        let mut parts = kv.splitn(2, '=');
        let k = parts.next().unwrap();
        let v = parts.next().unwrap_or("");
        info_map.insert(k.to_string(), v.to_string());
    }

    for field in &bundle.info {
        let key = field.key;

        if field.values.is_empty() {
            info_map.insert(key.to_string(), String::new());
            continue;
        }

        let vals = &field.values;

        match field.number {
            FieldNumber::Zero => {
                info_map.insert(key.to_string(), String::new());
            }

            FieldNumber::One => {
                if !vals.is_empty() {
                    let decoded = url_decode_info_value(vals[0]);
                    info_map.insert(key.to_string(), decoded);
                }
            }

            FieldNumber::Many => {
                let decoded: Vec<String> = vals.iter().map(|v| url_decode_info_value(v)).collect();
                info_map.insert(key.to_string(), decoded.join(","));
            }

            FieldNumber::A => {
                if vals.len() == src_alt_alleles.len() {
                    let decoded: Vec<String> =
                        vals.iter().map(|v| url_decode_info_value(v)).collect();
                    let joined = decoded.join(",");
                    let cleaned = clean_info_values(&joined);
                    if !cleaned.is_empty() {
                        info_map.insert(key.to_string(), cleaned);
                    }
                } else if vals.len() > src_alt_alleles.len() {
                    let taken: Vec<String> = vals
                        .iter()
                        .take(src_alt_alleles.len())
                        .map(|v| url_decode_info_value(v))
                        .collect();
                    let joined = taken.join(",");
                    let cleaned = clean_info_values(&joined);
                    if !cleaned.is_empty() {
                        info_map.insert(key.to_string(), cleaned);
                    }
                } else {
                    let mut extended = vals.to_vec();
                    while extended.len() < src_alt_alleles.len() {
                        extended.push(".");
                    }
                    let decoded: Vec<String> =
                        extended.iter().map(|v| url_decode_info_value(v)).collect();
                    let joined = decoded.join(",");
                    let cleaned = clean_info_values(&joined);
                    if !cleaned.is_empty() {
                        info_map.insert(key.to_string(), cleaned);
                    }
                }
            }

            FieldNumber::R => {
                let expected = src_alt_alleles.len() + 1;

                if vals.len() == expected {
                    let decoded: Vec<String> = vals
                        .iter()
                        .skip(1)
                        .map(|v| url_decode_info_value(v))
                        .collect();
                    let joined = decoded.join(",");
                    let cleaned = clean_info_values(&joined);
                    if !cleaned.is_empty() {
                        info_map.insert(key.to_string(), cleaned);
                    }
                } else if vals.len() > expected {
                    let decoded: Vec<String> = vals
                        .iter()
                        .skip(1)
                        .take(src_alt_alleles.len())
                        .map(|v| url_decode_info_value(v))
                        .collect();
                    let joined = decoded.join(",");
                    let cleaned = clean_info_values(&joined);
                    if !cleaned.is_empty() {
                        info_map.insert(key.to_string(), cleaned);
                    }
                } else {
                    let mut extended = vals.to_vec();
                    while extended.len() < expected {
                        extended.push(".");
                    }
                    let decoded: Vec<String> = extended
                        .iter()
                        .skip(1)
                        .map(|v| url_decode_info_value(v))
                        .collect();
                    let joined = decoded.join(",");
                    let cleaned = clean_info_values(&joined);
                    if !cleaned.is_empty() {
                        info_map.insert(key.to_string(), cleaned);
                    }
                }
            }

            FieldNumber::G => {
                let decoded: Vec<String> = vals.iter().map(|v| url_decode_info_value(v)).collect();
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

    if !rest.is_empty() {
        fields.push(rest.to_string());
    }

    fields.join("\t")
}

fn merge_annotation_headers(vcf_headers: &[String], _ani: &AniIndex) -> Result<Vec<String>> {
    Ok(vcf_headers.to_vec())
}
