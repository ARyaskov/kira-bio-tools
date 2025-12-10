use anyhow::Result;
use crossbeam_channel::{bounded, Receiver, Sender};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;
use std::time::Instant;

use super::reader::{StreamingVcfReader, VcfAnnotationReader};
use super::structs::*;
use crate::bgzf::BgzfWriter;
use crate::util::{detect_format, VcfFormat};
use crate::vcf::VcfParser;

const BATCH_SIZE: usize = 100_000;
const CHANNEL_DEPTH: usize = 32;
const OUTPUT_BUFFER_SIZE: usize = 256 * 1024 * 1024;
const ESTIMATE_BYTES_PER_LINE: usize = 200;

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

    if use_bgzf {
        run_optimized_bgzf_pipeline(merged_headers, reader, ani, output, timing, debug)
    } else {
        run_plain_pipeline_parallel(merged_headers, reader, ani, output, timing, debug)
    }
}

fn run_optimized_bgzf_pipeline(
    headers: Vec<String>,
    mut reader: StreamingVcfReader,
    ani: AniIndex,
    output: &Path,
    timing: bool,
    debug: bool,
) -> Result<()> {
    let start = Instant::now();

    let writer_start = Instant::now();
    let mut writer = BgzfWriter::create(output)?;
    if debug {
        eprintln!(
            "[annotate] BGZF writer ready: {:.3}s",
            writer_start.elapsed().as_secs_f64()
        );
    }

    for header in &headers {
        writer.write_all(header.as_bytes())?;
        writer.write_all(b"\n")?;
    }
    if debug {
        eprintln!("[annotate] Headers written");
    }

    let (line_tx, line_rx): (Sender<Vec<String>>, Receiver<Vec<String>>) = bounded(CHANNEL_DEPTH);
    let (result_tx, result_rx): (Sender<Vec<u8>>, Receiver<Vec<u8>>) = bounded(CHANNEL_DEPTH);

    let reader_thread = thread::spawn(move || -> Result<usize> {
        let mut lines_read = 0;
        let mut buffer = Vec::with_capacity(BATCH_SIZE);

        loop {
            for _ in 0..BATCH_SIZE {
                match reader.read_line()? {
                    Some(line) => buffer.push(line),
                    None => break,
                }
            }

            if buffer.is_empty() {
                break;
            }

            lines_read += buffer.len();
            let owned_batch = std::mem::replace(&mut buffer, Vec::with_capacity(BATCH_SIZE));

            if line_tx.send(owned_batch).is_err() {
                break;
            }
        }

        drop(line_tx);
        Ok(lines_read)
    });

    let num_workers = rayon::current_num_threads() / 4;
    let ani_arc = std::sync::Arc::new(ani);

    let worker_threads: Vec<_> = (0..num_workers)
        .map(|_worker_id| {
            let rx = line_rx.clone();
            let tx = result_tx.clone();
            let ani = ani_arc.clone();

            thread::spawn(move || {
                while let Ok(batch) = rx.recv() {
                    let mut mega_buffer =
                        Vec::with_capacity(batch.len() * ESTIMATE_BYTES_PER_LINE * 2);

                    for line in batch {
                        if annotate_line_to_buffer(&line, &ani, &mut mega_buffer).is_ok() {
                            mega_buffer.push(b'\n');
                        } else {
                            mega_buffer.extend_from_slice(line.as_bytes());
                            mega_buffer.push(b'\n');
                        }
                    }

                    if tx.send(mega_buffer).is_err() {
                        break;
                    }
                }
            })
        })
        .collect();

    drop(line_rx);
    drop(result_tx);

    let processed = AtomicUsize::new(0);
    let mut batch_num = 0;

    for chunk in result_rx.iter() {
        writer.write_all(&chunk)?;

        batch_num += 1;
        let chunk_lines = chunk.iter().filter(|&&b| b == b'\n').count();
        processed.fetch_add(chunk_lines, Ordering::Relaxed);

        if timing && batch_num % 10 == 0 {
            let total = processed.load(Ordering::Relaxed);
            let elapsed = start.elapsed().as_secs_f64();
            eprintln!(
                "[annotate] Batch {}: {} total variants ({:.0} var/s)",
                batch_num,
                total,
                total as f64 / elapsed
            );
        }
    }

    writer.finish()?;

    for handle in worker_threads {
        handle.join().ok();
    }

    let lines_read = reader_thread.join().unwrap()?;

    if timing || debug {
        let elapsed = start.elapsed().as_secs_f64();
        eprintln!(
            "[annotate] DONE: {} variants in {:.3}s ({:.0} var/s)",
            lines_read,
            elapsed,
            lines_read as f64 / elapsed
        );
    }

    Ok(())
}

fn run_plain_pipeline_parallel(
    headers: Vec<String>,
    mut reader: StreamingVcfReader,
    ani: AniIndex,
    output: &Path,
    timing: bool,
    debug: bool,
) -> Result<()> {
    let start = Instant::now();

    let output_file = File::create(output)?;
    let mut writer = BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, output_file);

    for header in &headers {
        writer.write_all(header.as_bytes())?;
        writer.write_all(b"\n")?;
    }

    if debug {
        eprintln!("[annotate] Headers written");
    }

    let (line_tx, line_rx): (Sender<Vec<String>>, Receiver<Vec<String>>) = bounded(CHANNEL_DEPTH);
    let (result_tx, result_rx): (Sender<Vec<u8>>, Receiver<Vec<u8>>) = bounded(CHANNEL_DEPTH);

    let reader_thread = thread::spawn(move || -> Result<usize> {
        let mut lines_read = 0;
        let mut buffer = Vec::with_capacity(BATCH_SIZE);

        loop {
            for _ in 0..BATCH_SIZE {
                match reader.read_line()? {
                    Some(line) => buffer.push(line),
                    None => break,
                }
            }

            if buffer.is_empty() {
                break;
            }

            lines_read += buffer.len();
            let owned_batch = std::mem::replace(&mut buffer, Vec::with_capacity(BATCH_SIZE));

            if line_tx.send(owned_batch).is_err() {
                break;
            }
        }

        drop(line_tx);
        Ok(lines_read)
    });

    let num_workers = rayon::current_num_threads() / 4;
    let ani_arc = std::sync::Arc::new(ani);

    let worker_threads: Vec<_> = (0..num_workers)
        .map(|_worker_id| {
            let rx = line_rx.clone();
            let tx = result_tx.clone();
            let ani = ani_arc.clone();

            thread::spawn(move || {
                while let Ok(batch) = rx.recv() {
                    let mut mega_buffer =
                        Vec::with_capacity(batch.len() * ESTIMATE_BYTES_PER_LINE * 2);

                    for line in batch {
                        if annotate_line_to_buffer(&line, &ani, &mut mega_buffer).is_ok() {
                            mega_buffer.push(b'\n');
                        } else {
                            mega_buffer.extend_from_slice(line.as_bytes());
                            mega_buffer.push(b'\n');
                        }
                    }

                    if tx.send(mega_buffer).is_err() {
                        break;
                    }
                }
            })
        })
        .collect();

    drop(line_rx);
    drop(result_tx);

    let processed = AtomicUsize::new(0);

    for chunk in result_rx.iter() {
        writer.write_all(&chunk)?;

        let chunk_lines = chunk.iter().filter(|&&b| b == b'\n').count();
        let total = processed.fetch_add(chunk_lines, Ordering::Relaxed) + chunk_lines;

        if timing && total % 100_000 == 0
            || (timing && chunk_lines > 0 && total / 100_000 != (total - chunk_lines) / 100_000)
        {
            let elapsed = start.elapsed().as_secs_f64();
            eprintln!(
                "[annotate] {} variants ({:.0} var/s)",
                total,
                total as f64 / elapsed
            );
        }
    }

    writer.flush()?;

    for handle in worker_threads {
        handle.join().ok();
    }

    let lines_read = reader_thread.join().unwrap()?;

    if timing || debug {
        let elapsed = start.elapsed().as_secs_f64();
        eprintln!(
            "[annotate] DONE: {} variants in {:.3}s ({:.0} var/s)",
            lines_read,
            elapsed,
            lines_read as f64 / elapsed
        );
    }

    Ok(())
}

#[allow(dead_code)]
fn run_plain_pipeline_single_threaded(
    headers: Vec<String>,
    mut reader: StreamingVcfReader,
    ani: AniIndex,
    output: &Path,
    timing: bool,
    _debug: bool,
) -> Result<()> {
    let start = Instant::now();

    let output_file = File::create(output)?;
    let mut writer = BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, output_file);

    for header in &headers {
        writer.write_all(header.as_bytes())?;
        writer.write_all(b"\n")?;
    }

    let mut buffer = Vec::with_capacity(ESTIMATE_BYTES_PER_LINE);
    let mut lines_processed = 0;

    while let Some(line) = reader.read_line()? {
        buffer.clear();

        if annotate_line_to_buffer(&line, &ani, &mut buffer).is_err() {
            buffer.extend_from_slice(line.as_bytes());
        }
        buffer.push(b'\n');
        writer.write_all(&buffer)?;

        lines_processed += 1;

        if timing && lines_processed % 100_000 == 0 {
            let elapsed = start.elapsed().as_secs_f64();
            eprintln!(
                "[annotate] {} variants ({:.0} var/s)",
                lines_processed,
                lines_processed as f64 / elapsed
            );
        }
    }

    writer.flush()?;

    if timing {
        let elapsed = start.elapsed().as_secs_f64();
        eprintln!(
            "[annotate] DONE: {} variants in {:.3}s ({:.0} var/s)",
            lines_processed,
            elapsed,
            lines_processed as f64 / elapsed
        );
    }

    Ok(())
}

#[inline]
fn annotate_line_to_buffer(line: &str, ani: &AniIndex, out: &mut Vec<u8>) -> Result<()> {
    use crate::vcf::simd::SimdVcfParser;

    let line_bytes = line.as_bytes();

    let fields = match SimdVcfParser::parse_fields(line_bytes) {
        Some(f) => f,
        None => {
            out.extend_from_slice(line_bytes);
            return Ok(());
        }
    };

    let chr = fields.chrom;
    let pos = fields.position().unwrap_or(0);
    let id = fields.id;
    let rf = fields.ref_allele;
    let alt_raw = fields.alt;
    let qual = fields.qual;
    let filter = fields.filter;
    let info = fields.info;

    let rest_offset = info.as_ptr() as usize - line.as_ptr() as usize + info.len();
    let rest = if rest_offset < line.len() {
        &line[rest_offset..]
    } else {
        ""
    };

    if let Some(ann) = ani.lookup(chr, pos, rf) {
        write_merged_record(
            out, chr, fields.pos, &ann, id, rf, alt_raw, qual, filter, info, rest,
        )?;
        return Ok(());
    }

    out.extend_from_slice(line_bytes);
    Ok(())
}

#[inline]
fn parse_u32_fast(bytes: &[u8]) -> Option<u32> {
    let mut result = 0u32;
    for &byte in bytes {
        if !byte.is_ascii_digit() {
            return None;
        }
        result = result.wrapping_mul(10).wrapping_add((byte - b'0') as u32);
    }
    Some(result)
}

fn write_merged_record(
    out: &mut Vec<u8>,
    chr: &str,
    pos: &str,
    ann: &AnnotationBundle,
    id: &str,
    rf: &str,
    alt_raw: &str,
    qual: &str,
    filter: &str,
    info: &str,
    rest: &str,
) -> Result<()> {
    out.extend_from_slice(chr.as_bytes());
    out.push(b'\t');
    out.extend_from_slice(pos.as_bytes());
    out.push(b'\t');

    let id2 = ann.id.unwrap_or(id);
    out.extend_from_slice(id2.as_bytes());
    out.push(b'\t');

    out.extend_from_slice(rf.as_bytes());
    out.push(b'\t');
    out.extend_from_slice(alt_raw.as_bytes());
    out.push(b'\t');

    let qual2 = ann.qual.unwrap_or(qual);
    out.extend_from_slice(qual2.as_bytes());
    out.push(b'\t');

    let filter2 = ann.filter.unwrap_or(filter);
    out.extend_from_slice(filter2.as_bytes());
    out.push(b'\t');

    write_merged_info_bcftools_compat(out, info, alt_raw, ann)?;

    if !rest.is_empty() {
        out.push(b'\t');
        out.extend_from_slice(rest.as_bytes());
    }

    Ok(())
}

fn write_merged_info_bcftools_compat(
    out: &mut Vec<u8>,
    base_info: &str,
    alt_raw: &str,
    ann: &AnnotationBundle,
) -> Result<()> {
    let alt_count = alt_raw.split(',').count();

    let mut merged_fields: HashMap<String, InfoFieldMerged> = HashMap::new();

    if !base_info.is_empty() && base_info != "." {
        for segment in base_info.split(';') {
            if segment.is_empty() {
                continue;
            }

            if !segment.contains('=') {
                merged_fields.insert(
                    segment.to_string(),
                    InfoFieldMerged {
                        values: vec![],
                        is_flag: true,
                        number: FieldNumber::Zero,
                        ty: FieldType::Flag,
                    },
                );
                continue;
            }

            let mut kv = segment.splitn(2, '=');
            let key = kv.next().unwrap().to_string();
            let vals = kv.next().unwrap_or("");
            let values: Vec<String> = vals.split(',').map(|s| s.to_string()).collect();

            merged_fields.insert(
                key,
                InfoFieldMerged {
                    values,
                    is_flag: false,
                    number: FieldNumber::Many,
                    ty: FieldType::Str,
                },
            );
        }
    }

    for field in &ann.info {
        let key = field.key.to_string();
        let is_append = key.starts_with('+');
        let clean_key = if is_append { &key[1..] } else { &key };

        let entry = merged_fields
            .entry(clean_key.to_string())
            .or_insert_with(|| InfoFieldMerged {
                values: vec![],
                is_flag: false,
                number: field.number,
                ty: field.ty,
            });

        entry.number = field.number;
        entry.ty = field.ty;

        if matches!(field.ty, FieldType::Flag) {
            entry.is_flag = true;
            continue;
        }

        match field.number {
            FieldNumber::A => {
                let mut new_values = vec![".".to_string(); alt_count];

                for i in 0..alt_count.min(field.values.len()) {
                    let decoded = url_decode_info_value(field.values[i]);
                    if !decoded.is_empty() && decoded != "." {
                        new_values[i] = decoded;
                    }
                }

                entry.values = new_values;
            }

            FieldNumber::R => {
                let mut new_values = vec![".".to_string(); alt_count + 1];

                for i in 0..(alt_count + 1).min(field.values.len()) {
                    let decoded = url_decode_info_value(field.values[i]);
                    if !decoded.is_empty() && decoded != "." {
                        new_values[i] = decoded;
                    }
                }

                entry.values = new_values;
            }

            FieldNumber::One | FieldNumber::Many => {
                let decoded_values: Vec<String> = field
                    .values
                    .iter()
                    .map(|v| url_decode_info_value(v))
                    .collect();

                entry.values = decoded_values;
            }

            _ => {}
        }
    }

    let mut keys: Vec<String> = merged_fields.keys().cloned().collect();
    keys.sort();

    let mut first = true;

    for key in keys {
        let f = &merged_fields[&key];

        if f.is_flag {
            if !first {
                out.push(b';');
            }
            out.extend_from_slice(key.as_bytes());
            first = false;
            continue;
        }

        let cleaned = clean_info_values(&f.values, f.number);

        if cleaned.is_empty() {
            continue;
        }

        if !first {
            out.push(b';');
        }
        out.extend_from_slice(key.as_bytes());
        out.push(b'=');
        out.extend_from_slice(cleaned.as_bytes());
        first = false;
    }

    if first {
        out.push(b'.');
    }

    Ok(())
}

struct InfoFieldMerged {
    values: Vec<String>,
    is_flag: bool,
    number: FieldNumber,
    ty: FieldType,
}

fn clean_info_values(values: &[String], number: FieldNumber) -> String {
    let mut v = values.to_vec();

    match number {
        FieldNumber::A => {
            while v.last().map_or(false, |s| s == ".") {
                v.pop();
            }
        }

        FieldNumber::R => {
            while v.len() > 1 && v.last().map_or(false, |s| s == ".") {
                v.pop();
            }
        }

        _ => {}
    }

    if v.is_empty() || v.iter().all(|s| s == ".") {
        return String::new();
    }

    v.join(",")
}

fn url_decode_info_value(val: &str) -> String {
    let mut result = String::with_capacity(val.len());
    let mut chars = val.chars().peekable();

    while let Some(c) = chars.next() {
        if c == '%' {
            let hex1 = chars.next();
            let hex2 = chars.next();

            if let (Some(h1), Some(h2)) = (hex1, hex2) {
                let hex_str = format!("{}{}", h1, h2);
                if let Ok(byte) = u8::from_str_radix(&hex_str, 16) {
                    result.push(byte as char);
                    continue;
                }
            }

            result.push(c);
        } else {
            result.push(c);
        }
    }

    result
}

fn merge_annotation_headers(vcf_headers: &[String], _ani: &AniIndex) -> Result<Vec<String>> {
    Ok(vcf_headers.to_vec())
}
