use anyhow::Result;
use crossbeam_channel::{Receiver, Sender};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

use crate::annotate::constants::OUTPUT_BUFFER_SIZE;
use crate::annotate::cpu_v2::annotation::{
    BundleTimingAccum, LineKeys, annotate_line, annotate_line_with_mph_hints,
    annotate_line_with_timing, column_targets_format, extract_line_keys,
};
use crate::annotate::cpu_v2::column_spec::ColumnSpec;
use crate::annotate::cpu_v2::read_batch::ReadBatch;
use crate::annotate::structs::ani::AniIndex;
use crate::annotate::structs::annotate_mode::AnnotateMode;
use crate::annotate::structs::bundle::FieldNumber;
use crate::bgzf::BgzfWriter;

/// Batch of input lines + their annotated outputs. `output[i] == Some(s)`
/// means line `i` changed and the writer emits `s`; `output[i] == None`
/// means unchanged and the writer takes the bytes from `input.line(i)` —
/// no allocation, no clone.
pub struct AnnotatedBatch {
    pub input: ReadBatch,
    pub output: Vec<Option<String>>,
}

pub fn worker_thread(
    rx: Receiver<ReadBatch>,
    tx: Sender<AnnotatedBatch>,
    ani: Arc<AniIndex>,
    field_meta: Arc<HashMap<String, FieldNumber>>,
    column_specs: Arc<Vec<ColumnSpec>>,
    sample_map: Arc<Vec<Option<usize>>>,
    info_overwrite_all: bool,
    format_overwrite_all: bool,
    num_threads: usize,
    timing: bool,
    bundle_acc: Arc<BundleTimingAccum>,
) -> Result<()> {
    let column_modes: Arc<Vec<(String, AnnotateMode)>> = Arc::new(
        column_specs
            .iter()
            .map(|c| (c.key.clone(), c.mode))
            .collect(),
    );

    // Whether the column spec actually touches FORMAT / sample fields. Same
    // condition `annotate_line` derives internally — surfaced here so the
    // batched key-extractor knows whether to parse the (expensive) sample
    // block. For pure INFO annotations (the dominant case) we skip it.
    let want_format = format_overwrite_all
        || column_modes
            .iter()
            .any(|(k, _)| column_targets_format(k));

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    while let Ok(batch) = rx.recv() {
        // Timing mode keeps the legacy per-line path so the bundle-timing
        // counters stay meaningful (they require the entry-by-entry call
        // graph). Non-timing mode uses the batched-SIMD pipeline below.
        let annotated: Vec<Option<String>> = if timing {
            pool.install(|| {
                (0..batch.len())
                    .into_par_iter()
                    .map(|i| {
                        annotate_line_with_timing(
                            batch.line(i),
                            &ani,
                            &field_meta,
                            &column_modes,
                            &sample_map,
                            info_overwrite_all,
                            format_overwrite_all,
                            &bundle_acc,
                        )
                    })
                    .collect()
            })
        } else {
            pool.install(|| {
                annotate_batch_simd(
                    &batch,
                    &ani,
                    &field_meta,
                    &column_modes,
                    &sample_map,
                    info_overwrite_all,
                    format_overwrite_all,
                    want_format,
                )
            })
        };

        if tx
            .send(AnnotatedBatch {
                input: batch,
                output: annotated,
            })
            .is_err()
        {
            break;
        }
    }

    Ok(())
}

/// Batched annotation pipeline.
///
/// 1. **Parallel parse + key extraction.** Each line is parsed once into a
///    borrowed `ParsedVcfRecord<'_>`, its contig is resolved, and the
///    MPH keys for all its ALT alleles are computed. No DB access happens
///    in this phase — purely CPU on the input line bytes.
/// 2. **Sequential key flatten.** All per-line key SmallVecs are dumped into
///    a single contiguous `Vec<u64>` plus a parallel `line_ranges` table.
/// 3. **Single SIMD batch MPH probe.** `Index::lookup_batch_u64_simd_into`
///    runs an AVX2 4-wide hash + 16-deep Bloom prefetch chain over the
///    entire flat-key array. Caller-owned scratch (`canon`, `out`) means
///    zero allocation in the lookup itself.
/// 4. **Parallel merge.** Each line consumes its precomputed MPH hits as
///    fast-path hints; misses or verification failures fall back to the
///    pos-index + vcmp pathway.
///
/// The win comes from collapsing N per-line scalar `lookup_u64` calls into
/// one wide SIMD probe — documented as 2-3× faster than scalar at our
/// batch sizes (`BATCH_SIZE` = 400 K).
#[allow(clippy::too_many_arguments)]
fn annotate_batch_simd(
    batch: &ReadBatch,
    ani: &Arc<AniIndex>,
    field_meta: &Arc<HashMap<String, FieldNumber>>,
    column_modes: &Arc<Vec<(String, AnnotateMode)>>,
    sample_map: &Arc<Vec<Option<usize>>>,
    info_overwrite_all: bool,
    format_overwrite_all: bool,
    want_format: bool,
) -> Vec<Option<String>> {
    // Phase 1: parallel parse + per-line key extraction.
    let line_keys: Vec<Option<LineKeys<'_>>> = (0..batch.len())
        .into_par_iter()
        .map(|i| extract_line_keys(batch.line(i), ani, want_format))
        .collect();

    // Phase 2: flatten keys into a single contiguous buffer.
    let total_keys: usize = line_keys
        .iter()
        .map(|lk| lk.as_ref().map_or(0, |x| x.keys.len()))
        .sum();
    let mut flat_keys: Vec<u64> = Vec::with_capacity(total_keys);
    let mut line_ranges: Vec<(usize, usize)> = Vec::with_capacity(batch.len());
    for lk in &line_keys {
        let start = flat_keys.len();
        if let Some(lk) = lk {
            flat_keys.extend_from_slice(&lk.keys);
        }
        line_ranges.push((start, flat_keys.len()));
    }

    // Phase 3: SIMD batch MPH probe (zero-alloc with caller-owned scratch).
    let mut canon = vec![0u64; flat_keys.len()];
    let mut mph_out = vec![None; flat_keys.len()];
    if !flat_keys.is_empty() {
        ani.index
            .lookup_batch_u64_simd_into(&flat_keys, &mut canon, &mut mph_out);
    }

    // Phase 4: parallel merge with MPH hints. Each line owns its slice of
    // `mph_out` so the per-line closures don't need synchronisation.
    (0..batch.len())
        .into_par_iter()
        .map(|i| {
            let Some(line_keys) = line_keys[i].as_ref() else {
                return None;
            };
            let (s, e) = line_ranges[i];
            annotate_line_with_mph_hints(
                batch.line(i),
                line_keys,
                &mph_out[s..e],
                ani,
                field_meta,
                column_modes,
                sample_map,
                info_overwrite_all,
                format_overwrite_all,
            )
        })
        .collect()
}

// Suppress the never-read warning for the path-import we keep as the
// reference scalar variant — useful for benchmarks and as a fallback if
// the SIMD path ever needs to be disabled at compile time.
#[allow(unused_imports)]
use annotate_line as _scalar_annotate_line;

fn annotated_line<'a>(batch: &'a AnnotatedBatch, idx: usize) -> &'a str {
    if let Some(Some(line)) = batch.output.get(idx) {
        return line.as_str();
    }
    batch.input.line(idx)
}

pub fn writer_thread_annotated(
    rx: Receiver<AnnotatedBatch>,
    headers: Vec<String>,
    output: &Path,
    use_bgzf: bool,
    bgzf_level: Option<u32>,
    mmap_output: bool,
    mmap_no_flush: bool,
    ram_output: bool,
    _ram_max_mb: u32,
    timing: bool,
    log_tag: &'static str,
) -> Result<()> {
    let start = Instant::now();
    let mut lines_written = 0usize;
    let mut bytes_written = 0usize;

    if timing {
        eprintln!("[{}] writer start", log_tag);
    }
    if ram_output {
        for h in headers {
            bytes_written += h.len() + 1;
        }

        while let Ok(batch) = rx.recv() {
            for idx in 0..batch.input.len() {
                let line = annotated_line(&batch, idx);
                bytes_written += line.len() + 1;
                lines_written += 1;
            }
        }
    } else if use_bgzf {
        let level = bgzf_level.unwrap_or(1).min(9);
        let mut writer = BgzfWriter::with_compression(output, flate2::Compression::new(level))?;

        for h in headers {
            writeln!(writer, "{}", h)?;
            bytes_written += h.len() + 1;
        }

        while let Ok(batch) = rx.recv() {
            for idx in 0..batch.input.len() {
                let line = annotated_line(&batch, idx);
                writeln!(writer, "{}", line)?;
                bytes_written += line.len() + 1;
                lines_written += 1;
            }
        }

        writer.finish()?;
    } else if mmap_output {
        let mut writer = crate::util::MmapWriter::create(output, OUTPUT_BUFFER_SIZE)?;

        for h in headers {
            writer.write_all(h.as_bytes())?;
            writer.write_all(b"\n")?;
            bytes_written += h.len() + 1;
        }

        while let Ok(batch) = rx.recv() {
            for idx in 0..batch.input.len() {
                let line = annotated_line(&batch, idx);
                writer.write_all(line.as_bytes())?;
                writer.write_all(b"\n")?;
                bytes_written += line.len() + 1;
                lines_written += 1;
            }
        }

        writer.finish(!mmap_no_flush)?;
    } else {
        let file = File::create(output)?;
        let mut writer = BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, file);

        for h in headers {
            writeln!(writer, "{}", h)?;
            bytes_written += h.len() + 1;
        }

        while let Ok(batch) = rx.recv() {
            for idx in 0..batch.input.len() {
                let line = annotated_line(&batch, idx);
                writeln!(writer, "{}", line)?;
                bytes_written += line.len() + 1;
                lines_written += 1;
            }
        }

        writer.flush()?;
    }

    if timing {
        let elapsed = start.elapsed().as_secs_f64();
        if ram_output {
            eprintln!(
                "[annotate] Write complete: {} lines (RAM sink), {:.3}s",
                lines_written, elapsed
            );
        } else {
            let mb_sec = (bytes_written as f64 / 1_048_576.0) / elapsed;
            eprintln!(
                "[annotate] Write complete: {} lines, {:.1} MB/s",
                lines_written, mb_sec
            );
        }
        eprintln!("[{}] writer end", log_tag);
    }

    Ok(())
}

pub fn writer_thread(
    rx: Receiver<Vec<String>>,
    headers: Vec<String>,
    output: &Path,
    use_bgzf: bool,
    bgzf_level: Option<u32>,
    mmap_output: bool,
    mmap_no_flush: bool,
    ram_output: bool,
    _ram_max_mb: u32,
    timing: bool,
    log_tag: &'static str,
) -> Result<()> {
    let start = Instant::now();
    let mut lines_written = 0usize;
    let mut bytes_written = 0usize;

    if timing {
        eprintln!("[{}] writer start", log_tag);
    }
    if ram_output {
        for h in headers {
            bytes_written += h.len() + 1;
        }

        while let Ok(batch) = rx.recv() {
            for line in batch {
                bytes_written += line.len() + 1;
                lines_written += 1;
            }
        }
    } else if use_bgzf {
        let level = bgzf_level.unwrap_or(1).min(9);
        let mut writer = BgzfWriter::with_compression(output, flate2::Compression::new(level))?;

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
    } else if mmap_output {
        let mut writer = crate::util::MmapWriter::create(output, OUTPUT_BUFFER_SIZE)?;

        for h in headers {
            writer.write_all(h.as_bytes())?;
            writer.write_all(b"\n")?;
            bytes_written += h.len() + 1;
        }

        while let Ok(batch) = rx.recv() {
            for line in batch {
                writer.write_all(line.as_bytes())?;
                writer.write_all(b"\n")?;
                bytes_written += line.len() + 1;
                lines_written += 1;
            }
        }

        writer.finish(!mmap_no_flush)?;
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
        if ram_output {
            eprintln!(
                "[annotate] Write complete: {} lines (RAM sink), {:.3}s",
                lines_written, elapsed
            );
        } else {
            let mb_sec = (bytes_written as f64 / 1_048_576.0) / elapsed;
            eprintln!(
                "[annotate] Write complete: {} lines, {:.1} MB/s",
                lines_written, mb_sec
            );
        }
        eprintln!("[{}] writer end", log_tag);
    }

    Ok(())
}
