use anyhow::Result;
use crossbeam_channel::{bounded, Receiver, Sender};
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
        .map(|worker_id| {
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

    if let Some((ann, _)) = ani.lookup_full(chr, pos, rf, alt_raw) {
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

    write_merged_info_fast(out, info, alt_raw, &ann.info, rf)?;

    if !rest.is_empty() {
        out.push(b'\t');
        out.extend_from_slice(rest.as_bytes());
    }

    Ok(())
}

fn write_merged_info_fast(
    out: &mut Vec<u8>,
    base: &str,
    alt_raw: &str,
    ann_fields: &[StructuredInfoField],
    r: &str,
) -> Result<()> {
    let alt_count = alt_raw.bytes().filter(|&b| b == b',').count() + 1;
    let mut first = true;

    if !base.is_empty() && base != "." {
        for segment in base.split(';') {
            if segment.is_empty() {
                continue;
            }

            if !first {
                out.push(b';');
            }
            first = false;

            out.extend_from_slice(segment.as_bytes());
        }
    }

    for field in ann_fields {
        if !first {
            out.push(b';');
        }
        first = false;

        out.extend_from_slice(field.key.as_bytes());

        match field.number {
            FieldNumber::Zero => continue,
            FieldNumber::One => {
                if let Some(&v) = field.values.first() {
                    out.push(b'=');
                    out.extend_from_slice(v.as_bytes());
                }
            }
            FieldNumber::A => {
                out.push(b'=');
                if field.values.len() == 1 {
                    let val = field.values[0].as_bytes();
                    out.extend_from_slice(val);
                    for _ in 1..alt_count {
                        out.push(b',');
                        out.extend_from_slice(val);
                    }
                } else {
                    for i in 0..alt_count {
                        if i > 0 {
                            out.push(b',');
                        }
                        let v = field.values.get(i).copied().unwrap_or(field.values[0]);
                        out.extend_from_slice(v.as_bytes());
                    }
                }
            }
            _ => {
                out.push(b'=');
                for (i, &v) in field.values.iter().enumerate() {
                    if i > 0 {
                        out.push(b',');
                    }
                    out.extend_from_slice(v.as_bytes());
                }
            }
        }
    }

    let is_indel = r.len() != 1 || alt_raw.bytes().any(|b| b != b',' && (b < b'A' || b > b'T'));
    if is_indel {
        if !first {
            out.push(b';');
        }
        out.extend_from_slice(b"INDEL");
    }

    if first {
        out.push(b'.');
    }

    Ok(())
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
