use anyhow::Result;
use crossbeam_channel::{bounded, Receiver, Sender};
use kira_kv_engine::HybridIndex;
use rayon::prelude::*;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::thread;
use std::time::Instant;

use crate::annotate::constants::{BATCH_SIZE, CHANNEL_DEPTH};
use crate::annotate::cpu_v2::annotation::annotate_line;
use crate::annotate::cpu_v2::field_metadata::{iter_ani_header_lines, load_and_infer_metadata};
use crate::annotate::cpu_v2::{
    annotate_record_with_bundles, build_sample_map, expand_column_specs,
    extract_samples_from_headers, merge_annotation_headers, ColumnSpec,
};
use crate::annotate::cpu_v2::{parse_vcf_record_simd, patch_samples_from_line, writer_thread};
use crate::annotate::reader::{StreamingVcfReader, VcfAnnotationReader};
use crate::annotate::structs::ani::AniIndex;
use crate::annotate::structs::annotate_mode::AnnotateMode;
use crate::annotate::structs::bundle::{AnnotationBundle, FieldNumber};
use crate::util::{chr_name_to_id, detect_format, fast_hash64, VcfFormat};
use std::collections::HashMap;

pub struct GpuAni {
    index: HybridIndex,
    entry_keys: Vec<u64>,
}

struct GpuLookupBuffers {
    out: Vec<u32>,
    capacity: usize,
}

pub struct GpuAnnotator {
    gpu: GpuAni,
    buffers: GpuLookupBuffers,
}

impl GpuAnnotator {
    pub fn new(ani: &AniIndex) -> Result<Self> {
        let gpu = GpuAni::load(ani)?;
        let buffers = GpuLookupBuffers::new(BATCH_SIZE)?;
        Ok(Self { gpu, buffers })
    }
}

impl GpuAni {
    pub fn load(ani: &AniIndex) -> Result<Self> {
        let index_bytes = ani.index.to_bytes()?;
        let index = HybridIndex::from_bytes(&index_bytes)?;
        let entry_keys = build_entry_keys(ani);
        Ok(Self { index, entry_keys })
    }

    pub fn lookup_batch(&self, keys: &[u64]) -> Result<Vec<u32>> {
        Ok(self.lookup_batch_cpu(keys))
    }

    fn lookup_batch_with_buffers(
        &self,
        keys: &[u64],
        buffers: &mut GpuLookupBuffers,
    ) -> Result<Vec<u32>> {
        Ok(self.lookup_batch_cpu_with_buffers(keys, buffers))
    }

    pub fn lookup_batch_from_strings(
        &self,
        ref_pool: &[u8],
        ref_offsets: &[u32],
        ref_lens: &[u32],
        alt_pool: &[u8],
        alt_offsets: &[u32],
        alt_lens: &[u32],
        key_ref_idx: &[u32],
        key_chr: &[u8],
        key_pos: &[u32],
    ) -> Result<Vec<u32>> {
        Ok(self.lookup_batch_from_strings_cpu(
            ref_pool,
            ref_offsets,
            ref_lens,
            alt_pool,
            alt_offsets,
            alt_lens,
            key_ref_idx,
            key_chr,
            key_pos,
        ))
    }

    fn lookup_batch_cpu(&self, keys: &[u64]) -> Vec<u32> {
        let mut out = vec![u32::MAX; keys.len()];
        for (i, &key) in keys.iter().enumerate() {
            if let Ok(idx) = self.index.lookup_u64(key) {
                if idx < self.entry_keys.len() && self.entry_keys[idx] == key {
                    out[i] = idx as u32;
                }
            }
        }
        out
    }

    fn lookup_batch_cpu_with_buffers(
        &self,
        keys: &[u64],
        buffers: &mut GpuLookupBuffers,
    ) -> Vec<u32> {
        if keys.is_empty() {
            return Vec::new();
        }
        buffers.ensure_capacity(keys.len());
        for v in buffers.out.iter_mut().take(keys.len()) {
            *v = u32::MAX;
        }
        for (i, &key) in keys.iter().enumerate() {
            if let Ok(idx) = self.index.lookup_u64(key) {
                if idx < self.entry_keys.len() && self.entry_keys[idx] == key {
                    buffers.out[i] = idx as u32;
                }
            }
        }
        buffers.out[..keys.len()].to_vec()
    }

    fn lookup_batch_from_strings_cpu(
        &self,
        ref_pool: &[u8],
        ref_offsets: &[u32],
        ref_lens: &[u32],
        alt_pool: &[u8],
        alt_offsets: &[u32],
        alt_lens: &[u32],
        key_ref_idx: &[u32],
        key_chr: &[u8],
        key_pos: &[u32],
    ) -> Vec<u32> {
        let n = alt_offsets.len();
        let mut out = vec![u32::MAX; n];
        for i in 0..n {
            let ref_idx = key_ref_idx[i] as usize;
            if ref_idx >= ref_offsets.len() || ref_idx >= ref_lens.len() {
                continue;
            }
            let ref_start = ref_offsets[ref_idx] as usize;
            let ref_len = ref_lens[ref_idx] as usize;
            if ref_start + ref_len > ref_pool.len() {
                continue;
            }
            let alt_start = alt_offsets[i] as usize;
            let alt_len = alt_lens[i] as usize;
            if alt_start + alt_len > alt_pool.len() {
                continue;
            }
            let chr_id = key_chr[i];
            let pos = key_pos[i];
            let mut key = (chr_id as u64) << 32 | pos as u64;
            key ^= fast_hash64(&ref_pool[ref_start..ref_start + ref_len]);
            key ^= fast_hash64(&alt_pool[alt_start..alt_start + alt_len]);
            if let Ok(idx) = self.index.lookup_u64(key) {
                if idx < self.entry_keys.len() && self.entry_keys[idx] == key {
                    out[i] = idx as u32;
                }
            }
        }
        out
    }
}

impl GpuLookupBuffers {
    fn new(capacity: usize) -> Result<Self> {
        Ok(Self {
            out: vec![0u32; capacity],
            capacity,
        })
    }

    fn ensure_capacity(&mut self, capacity: usize) {
        if capacity <= self.capacity {
            return;
        }
        self.out = vec![0u32; capacity];
        self.capacity = capacity;
    }
}

pub fn annotate_vcf_ani_gpu(
    ani: &AniIndex,
    input: &Path,
    output: &Path,
    columns: &[String],
    bgzf_level: Option<u32>,
    mmap_output: bool,
    mmap_no_flush: bool,
    ram_output: bool,
    ram_max_mb: u32,
) -> Result<()> {
    let total_start = Instant::now();
    let timing = std::env::var("KIRA_BT_TIMING").is_ok();
    if timing {
        eprintln!("[gpu] timing: enabled");
    }
    let mut column_specs = ColumnSpec::parse_all(columns);
    let info_overwrite_all = column_specs
        .iter()
        .any(|c| c.key.eq_ignore_ascii_case("INFO"));
    let format_overwrite_all = column_specs
        .iter()
        .any(|c| c.key.eq_ignore_ascii_case("FMT") || c.key.eq_ignore_ascii_case("FORMAT"));

    let meta_start = Instant::now();
    let field_meta = load_and_infer_metadata(ani, false)?;
    let ani_headers = iter_ani_header_lines(ani);
    column_specs = expand_column_specs(&column_specs, &ani_headers, &field_meta);
    let column_modes: Vec<(String, AnnotateMode)> = column_specs
        .iter()
        .map(|c| (c.key.clone(), c.mode))
        .collect();
    eprintln!("[gpu] metadata: {:.3}s", meta_start.elapsed().as_secs_f64());

    let input_format = detect_format(input)?;
    let output_ext = output.extension().and_then(|s| s.to_str()).unwrap_or("");
    let output_wants_bgzf = matches!(output_ext, "gz" | "bgz" | "bgzf");
    let use_bgzf = matches!(input_format, VcfFormat::Bgzf) || output_wants_bgzf;

    let input_reader = VcfAnnotationReader::open(input)?;
    let streaming_reader = StreamingVcfReader::new(input_reader);
    let (headers, mut reader) = streaming_reader.into_headers_and_self()?;

    let header_start = Instant::now();
    let merged_headers = merge_annotation_headers(&headers, &ani_headers, &column_specs)?;
    let input_samples = extract_samples_from_headers(&headers);
    let db_samples = extract_samples_from_headers(&ani_headers);
    let sample_map = build_sample_map(&input_samples, &db_samples);
    eprintln!(
        "[gpu] headers: {:.3}s",
        header_start.elapsed().as_secs_f64()
    );

    let (read_tx, read_rx) = bounded::<Vec<String>>(CHANNEL_DEPTH);
    let (work_tx, work_rx) = bounded::<Vec<String>>(CHANNEL_DEPTH);
    let field_meta = Arc::new(field_meta);
    let column_modes = Arc::new(column_modes);
    let sample_map = Arc::new(sample_map);
    let use_bgzf = use_bgzf;

    thread::scope(|s| -> Result<()> {
        let worker = s.spawn(move || {
            worker_thread_gpu(
                read_rx,
                work_tx,
                ani,
                field_meta,
                column_modes,
                sample_map,
                info_overwrite_all,
                format_overwrite_all,
                timing,
            )
        });

        let output_path = output.to_path_buf();
        let headers = merged_headers;
        let writer = s.spawn(move || {
            writer_thread(
                work_rx,
                headers,
                &output_path,
                use_bgzf,
                bgzf_level,
                mmap_output,
                mmap_no_flush,
                ram_output,
                ram_max_mb,
                timing,
                "gpu",
            )
        });

        let join_start = Instant::now();
        read_batches_gpu(&mut reader, read_tx, timing)?;
        if timing {
            eprintln!("[gpu] waiting worker...");
        }
        worker.join().unwrap()?;
        if timing {
            eprintln!(
                "[gpu] worker done: {:.3}s",
                join_start.elapsed().as_secs_f64()
            );
        }
        if timing {
            eprintln!("[gpu] waiting writer...");
        }
        writer.join().unwrap()?;
        if timing {
            eprintln!(
                "[gpu] writer done: {:.3}s",
                join_start.elapsed().as_secs_f64()
            );
        }
        Ok(())
    })?;

    eprintln!("[gpu] total: {:.3}s", total_start.elapsed().as_secs_f64());
    Ok(())
}

struct MinParsed {
    chr_id: u8,
    pos: u32,
    ref_range: (usize, usize),
    alt_range: (usize, usize),
}

fn make_key(chr_id: u8, pos: u32, ref_allele: &str, alt: &str) -> u64 {
    let mut h = (chr_id as u64) << 32 | pos as u64;
    h ^= fast_hash64(ref_allele.as_bytes());
    h ^= fast_hash64(alt.as_bytes());
    h
}

fn fast_parse_min(line: &str) -> Option<MinParsed> {
    let bytes = line.as_bytes();
    let mut tabs = [0usize; 5];
    let mut count = 0usize;
    for (i, b) in bytes.iter().enumerate() {
        if *b == b'\t' {
            tabs[count] = i;
            count += 1;
            if count == 5 {
                break;
            }
        }
    }
    if count < 4 {
        return None;
    }

    let chrom = &line[..tabs[0]];
    let chr_id = chr_name_to_id(chrom)?;

    let pos_bytes = &bytes[(tabs[0] + 1)..tabs[1]];
    let pos = parse_u32_bytes(pos_bytes)?;

    let ref_start = tabs[2] + 1;
    let ref_end = tabs[3];
    if ref_end <= ref_start || ref_end > bytes.len() {
        return None;
    }

    let alt_start = tabs[3] + 1;
    let alt_end = if count >= 5 { tabs[4] } else { bytes.len() };
    if alt_end <= alt_start || alt_end > bytes.len() {
        return None;
    }

    Some(MinParsed {
        chr_id,
        pos,
        ref_range: (ref_start, ref_end),
        alt_range: (alt_start, alt_end),
    })
}

fn alt_slice_for_idx(bytes: &[u8], alt_range: (usize, usize), alt_idx: usize) -> Option<&[u8]> {
    let mut start = alt_range.0;
    let mut idx = 0usize;
    for i in alt_range.0..alt_range.1 {
        if bytes[i] == b',' {
            if idx == alt_idx {
                return Some(&bytes[start..i]);
            }
            idx += 1;
            start = i + 1;
        }
    }
    if idx == alt_idx && start < alt_range.1 {
        return Some(&bytes[start..alt_range.1]);
    }
    None
}

fn parse_u32_bytes(bytes: &[u8]) -> Option<u32> {
    if bytes.is_empty() || bytes.len() > 10 {
        return None;
    }
    let mut v: u32 = 0;
    for b in bytes {
        if !b.is_ascii_digit() {
            return None;
        }
        v = v.wrapping_mul(10).wrapping_add((b - b'0') as u32);
    }
    Some(v)
}

fn want_format_for_line(line: &str, need_format: bool) -> bool {
    if need_format {
        return true;
    }

    let mut tabs = 0usize;
    for &b in line.as_bytes() {
        if b == b'\t' {
            tabs += 1;
            if tabs >= 8 {
                return true;
            }
        }
    }
    false
}

fn worker_thread_gpu(
    rx: Receiver<Vec<String>>,
    tx: Sender<Vec<String>>,
    ani: &AniIndex,
    field_meta: Arc<HashMap<String, FieldNumber>>,
    column_modes: Arc<Vec<(String, AnnotateMode)>>,
    sample_map: Arc<Vec<Option<usize>>>,
    info_overwrite_all: bool,
    format_overwrite_all: bool,
    timing: bool,
) -> Result<()> {
    let mut gpu = GpuAni::load(ani)?;
    let mut buffers = GpuLookupBuffers::new(BATCH_SIZE)?;
    worker_loop_gpu(
        rx,
        tx,
        ani,
        field_meta,
        column_modes,
        sample_map,
        info_overwrite_all,
        format_overwrite_all,
        timing,
        &mut gpu,
        &mut buffers,
    )
}

fn worker_loop_gpu(
    rx: Receiver<Vec<String>>,
    tx: Sender<Vec<String>>,
    ani: &AniIndex,
    field_meta: Arc<HashMap<String, FieldNumber>>,
    column_modes: Arc<Vec<(String, AnnotateMode)>>,
    sample_map: Arc<Vec<Option<usize>>>,
    info_overwrite_all: bool,
    format_overwrite_all: bool,
    timing: bool,
    gpu: &mut GpuAni,
    buffers: &mut GpuLookupBuffers,
) -> Result<()> {
    let num_threads = rayon::current_num_threads().max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    let mut total_lines = 0usize;
    let mut parse_total = 0f64;
    let mut lookup_total = 0f64;
    let mut annotate_total = 0f64;
    let mut keys_total = 0f64;
    let mut bundle_total = 0f64;
    let mut bundle_read_total = 0f64;
    let mut bundle_info_total = 0f64;
    let mut bundle_optional_total = 0f64;
    let mut bundle_samples_total = 0f64;
    let mut send_total = 0f64;
    let mut last_report = Instant::now();

    let need_format = format_overwrite_all
        || column_modes
            .iter()
            .any(|(k, _)| k.eq_ignore_ascii_case("FMT") || k.eq_ignore_ascii_case("FORMAT"));
    let need_info = info_overwrite_all
        || column_modes.iter().any(|(k, _)| {
            !(k.eq_ignore_ascii_case("ID")
                || k.eq_ignore_ascii_case("QUAL")
                || k.eq_ignore_ascii_case("FILTER")
                || k.eq_ignore_ascii_case("FMT")
                || k.eq_ignore_ascii_case("FORMAT"))
        });
    let use_cpu_annotate = sample_map.iter().any(|v| v.is_none());
    if timing && use_cpu_annotate {
        eprintln!("[gpu] cpu-annotate: enabled");
    }
    while let Ok(mut batch) = rx.recv() {
        if batch.is_empty() {
            continue;
        }

        if use_cpu_annotate {
            let annotated: Vec<String> = pool.install(|| {
                batch
                    .par_iter()
                    .map(|line| {
                        annotate_line(
                            line,
                            ani,
                            &field_meta,
                            &column_modes,
                            &sample_map,
                            info_overwrite_all,
                            format_overwrite_all,
                        )
                    })
                    .collect()
            });
            let send_start = Instant::now();
            if tx.send(annotated).is_err() {
                break;
            }
            if timing {
                send_total += send_start.elapsed().as_secs_f64();
            }
            continue;
        }

        let batch_start = Instant::now();
        let parse_start = Instant::now();
        let mins: Vec<Option<MinParsed>> =
            pool.install(|| batch.par_iter().map(|line| fast_parse_min(line)).collect());
        parse_total += parse_start.elapsed().as_secs_f64();

        let keys_start = Instant::now();
        let mut keys: Vec<u64> = Vec::new();
        let mut key_line_idx: Vec<usize> = Vec::new();
        let mut key_alt_idx: Vec<usize> = Vec::new();

        for (line_idx, line) in batch.iter().enumerate() {
            let Some(ref min) = mins[line_idx] else {
                continue;
            };
            let bytes = line.as_bytes();
            let ref_bytes = &bytes[min.ref_range.0..min.ref_range.1];
            let ref_hash = fast_hash64(ref_bytes);
            let base = (min.chr_id as u64) << 32 | min.pos as u64;

            let (alt_start, alt_end) = min.alt_range;
            let mut alt_idx = 0usize;
            let mut start = alt_start;
            for i in alt_start..alt_end {
                if bytes[i] == b',' {
                    let alt_hash = fast_hash64(&bytes[start..i]);
                    let key = base ^ ref_hash ^ alt_hash;
                    keys.push(key);
                    key_line_idx.push(line_idx);
                    key_alt_idx.push(alt_idx);
                    alt_idx += 1;
                    start = i + 1;
                }
            }
            if start < alt_end {
                let alt_hash = fast_hash64(&bytes[start..alt_end]);
                let key = base ^ ref_hash ^ alt_hash;
                keys.push(key);
                key_line_idx.push(line_idx);
                key_alt_idx.push(alt_idx);
            }
        }
        keys_total += keys_start.elapsed().as_secs_f64();

        if keys.is_empty() {
            if tx.send(batch).is_err() {
                break;
            }
            continue;
        }

        let lookup_start = Instant::now();
        let idxs = gpu.lookup_batch_with_buffers(&keys, buffers)?;
        lookup_total += lookup_start.elapsed().as_secs_f64();
        let bundle_start = Instant::now();
        let mut bundles_per_line: Vec<Vec<(usize, AnnotationBundle)>> =
            vec![Vec::new(); batch.len()];
        for (i, idx) in idxs.iter().enumerate() {
            let line_idx = key_line_idx[i];
            let alt_idx = key_alt_idx[i];
            if *idx == u32::MAX {
                if let Some(min) = mins[line_idx].as_ref() {
                    let bytes = batch[line_idx].as_bytes();
                    let ref_bytes = &bytes[min.ref_range.0..min.ref_range.1];
                    if let Some(alt_bytes) = alt_slice_for_idx(bytes, min.alt_range, alt_idx) {
                        if let (Ok(ref_str), Ok(alt_str)) = (
                            std::str::from_utf8(ref_bytes),
                            std::str::from_utf8(alt_bytes),
                        ) {
                            let ref_hash = fast_hash64(ref_str.as_bytes());
                            if let Some(bundle) = ani.lookup_exact_by_chr_id_opts(
                                min.chr_id,
                                min.pos,
                                ref_str,
                                ref_hash,
                                alt_str,
                                need_info,
                                need_format,
                            ) {
                                bundles_per_line[line_idx].push((alt_idx, bundle));
                            }
                        }
                    }
                }
                continue;
            }
            let bundle = if timing {
                let (bundle, t) = ani.build_bundle_from_entry_timed_opts(
                    &ani.entries[*idx as usize],
                    need_info,
                    need_format,
                );
                bundle_read_total += t.read_s;
                bundle_info_total += t.info_s;
                bundle_optional_total += t.optional_s;
                bundle_samples_total += t.samples_s;
                bundle
            } else {
                ani.build_bundle_from_entry_opts(
                    &ani.entries[*idx as usize],
                    need_info,
                    need_format,
                )
            };
            bundles_per_line[line_idx].push((alt_idx, bundle));
        }
        bundle_total += bundle_start.elapsed().as_secs_f64();

        let annotate_start = Instant::now();
        let annotated: Vec<Option<String>> = pool.install(|| {
            batch
                .par_iter()
                .enumerate()
                .map(|(i, line)| {
                    if bundles_per_line[i].is_empty() {
                        return None;
                    }
                    let want_format = want_format_for_line(line, need_format);
                    parse_vcf_record_simd(line, want_format).map(|mut parsed| {
                        patch_samples_from_line(&mut parsed, line);
                        annotate_record_with_bundles(
                            &parsed,
                            &bundles_per_line[i],
                            &field_meta,
                            &column_modes,
                            &sample_map,
                            info_overwrite_all,
                            format_overwrite_all,
                            false,
                        )
                    })
                })
                .collect()
        });
        annotate_total += annotate_start.elapsed().as_secs_f64();

        for (i, val) in annotated.into_iter().enumerate() {
            if let Some(s) = val {
                batch[i] = s;
            }
        }

        let batch_len = batch.len();
        let send_start = Instant::now();
        if tx.send(batch).is_err() {
            break;
        }
        send_total += send_start.elapsed().as_secs_f64();

        if timing {
            total_lines += key_line_idx.len();
            if last_report.elapsed().as_secs_f64() >= 2.0 {
                eprintln!(
                    "[gpu] lines: {}, parse: {:.3}s, keys: {:.3}s, lookup: {:.3}s, bundle: {:.3}s, annotate: {:.3}s, send: {:.3}s",
                    total_lines,
                    parse_total,
                    keys_total,
                    lookup_total,
                    bundle_total,
                    annotate_total,
                    send_total
                );
                eprintln!(
                    "[gpu] bundle_read: {:.3}s, bundle_info: {:.3}s, bundle_optional: {:.3}s, bundle_samples: {:.3}s",
                    bundle_read_total,
                    bundle_info_total,
                    bundle_optional_total,
                    bundle_samples_total
                );
                last_report = Instant::now();
            }
            if batch_start.elapsed().as_secs_f64() >= 2.0 {
                eprintln!(
                    "[gpu] batch: {:.3}s, keys: {}, lines: {}",
                    batch_start.elapsed().as_secs_f64(),
                    keys.len(),
                    batch_len
                );
            }
            if send_start.elapsed().as_secs_f64() >= 2.0 {
                eprintln!(
                    "[gpu] send blocked: {:.3}s",
                    send_start.elapsed().as_secs_f64()
                );
            }
        }
    }

    if timing {
        eprintln!(
            "[gpu] done: lines: {}, parse: {:.3}s, keys: {:.3}s, lookup: {:.3}s, bundle: {:.3}s, annotate: {:.3}s, send: {:.3}s",
            total_lines,
            parse_total,
            keys_total,
            lookup_total,
            bundle_total,
            annotate_total,
            send_total
        );
        eprintln!("[gpu] bundle: {:.3}s", bundle_total);
        eprintln!(
            "[gpu] bundle_read: {:.3}s, bundle_info: {:.3}s, bundle_optional: {:.3}s, bundle_samples: {:.3}s",
            bundle_read_total,
            bundle_info_total,
            bundle_optional_total,
            bundle_samples_total
        );
    }

    Ok(())
}

pub fn annotate_vcf_ani_gpu_with_state(
    state: &mut GpuAnnotator,
    ani: &AniIndex,
    input: &Path,
    output: &Path,
    columns: &[String],
    bgzf_level: Option<u32>,
    mmap_output: bool,
    mmap_no_flush: bool,
    ram_output: bool,
    ram_max_mb: u32,
) -> Result<()> {
    let total_start = Instant::now();
    let timing = std::env::var("KIRA_BT_TIMING").is_ok();
    if timing {
        eprintln!("[gpu] timing: enabled");
    }
    let mut column_specs = ColumnSpec::parse_all(columns);
    let info_overwrite_all = column_specs
        .iter()
        .any(|c| c.key.eq_ignore_ascii_case("INFO"));
    let format_overwrite_all = column_specs
        .iter()
        .any(|c| c.key.eq_ignore_ascii_case("FMT") || c.key.eq_ignore_ascii_case("FORMAT"));

    let meta_start = Instant::now();
    let field_meta = load_and_infer_metadata(ani, false)?;
    let ani_headers = iter_ani_header_lines(ani);
    column_specs = expand_column_specs(&column_specs, &ani_headers, &field_meta);
    let column_modes: Vec<(String, AnnotateMode)> = column_specs
        .iter()
        .map(|c| (c.key.clone(), c.mode))
        .collect();
    eprintln!("[gpu] metadata: {:.3}s", meta_start.elapsed().as_secs_f64());

    let input_format = detect_format(input)?;
    let output_ext = output.extension().and_then(|s| s.to_str()).unwrap_or("");
    let output_wants_bgzf = matches!(output_ext, "gz" | "bgz" | "bgzf");
    let use_bgzf = matches!(input_format, VcfFormat::Bgzf) || output_wants_bgzf;

    let input_reader = VcfAnnotationReader::open(input)?;
    let streaming_reader = StreamingVcfReader::new(input_reader);
    let (headers, mut reader) = streaming_reader.into_headers_and_self()?;

    let header_start = Instant::now();
    let merged_headers = merge_annotation_headers(&headers, &ani_headers, &column_specs)?;
    let input_samples = extract_samples_from_headers(&headers);
    let db_samples = extract_samples_from_headers(&ani_headers);
    let sample_map = build_sample_map(&input_samples, &db_samples);
    eprintln!(
        "[gpu] headers: {:.3}s",
        header_start.elapsed().as_secs_f64()
    );

    let (read_tx, read_rx) = bounded::<Vec<String>>(CHANNEL_DEPTH);
    let (work_tx, work_rx) = bounded::<Vec<String>>(CHANNEL_DEPTH);
    let field_meta = Arc::new(field_meta);
    let column_modes = Arc::new(column_modes);
    let sample_map = Arc::new(sample_map);
    let use_bgzf = use_bgzf;

    let output_path = output.to_path_buf();
    let headers = merged_headers;
    let writer = thread::spawn(move || {
        writer_thread(
            work_rx,
            headers,
            &output_path,
            use_bgzf,
            bgzf_level,
            mmap_output,
            mmap_no_flush,
            ram_output,
            ram_max_mb,
            timing,
            "gpu",
        )
    });

    let reader_thread = thread::spawn(move || read_batches_gpu(&mut reader, read_tx, timing));

    let join_start = Instant::now();
    worker_loop_gpu(
        read_rx,
        work_tx,
        ani,
        field_meta,
        column_modes,
        sample_map,
        info_overwrite_all,
        format_overwrite_all,
        timing,
        &mut state.gpu,
        &mut state.buffers,
    )?;
    if timing {
        eprintln!(
            "[gpu] worker done: {:.3}s",
            join_start.elapsed().as_secs_f64()
        );
    }
    reader_thread.join().unwrap()?;
    if timing {
        eprintln!("[gpu] waiting writer...");
    }
    writer.join().unwrap()?;
    if timing {
        eprintln!(
            "[gpu] writer done: {:.3}s",
            join_start.elapsed().as_secs_f64()
        );
    }

    eprintln!("[gpu] total: {:.3}s", total_start.elapsed().as_secs_f64());
    Ok(())
}

pub struct GpuBatch {
    pub vcf_id: usize,
    pub lines: Vec<String>,
}

pub struct GpuJob {
    pub input: PathBuf,
    pub output: PathBuf,
    pub use_bgzf: bool,
    pub headers: Vec<String>,
    pub sample_map: Arc<Vec<Option<usize>>>,
}

const MULTI_MAX_LINES: usize = BATCH_SIZE * 4;

pub fn annotate_vcf_ani_gpu_multi_with_state(
    state: &mut GpuAnnotator,
    ani: &AniIndex,
    jobs: Vec<GpuJob>,
    columns: &[String],
    bgzf_level: Option<u32>,
    mmap_output: bool,
    mmap_no_flush: bool,
    ram_output: bool,
    ram_max_mb: u32,
    max_writers: usize,
) -> Result<()> {
    let total_start = Instant::now();
    let timing = std::env::var("KIRA_BT_TIMING").is_ok();
    if timing {
        eprintln!("[gpu] timing: enabled");
    }

    let mut column_specs = ColumnSpec::parse_all(columns);
    let info_overwrite_all = column_specs
        .iter()
        .any(|c| c.key.eq_ignore_ascii_case("INFO"));
    let format_overwrite_all = column_specs
        .iter()
        .any(|c| c.key.eq_ignore_ascii_case("FMT") || c.key.eq_ignore_ascii_case("FORMAT"));

    let meta_start = Instant::now();
    let field_meta = load_and_infer_metadata(ani, false)?;
    let ani_headers = iter_ani_header_lines(ani);
    column_specs = expand_column_specs(&column_specs, &ani_headers, &field_meta);
    let column_modes: Vec<(String, AnnotateMode)> = column_specs
        .iter()
        .map(|c| (c.key.clone(), c.mode))
        .collect();
    eprintln!("[gpu] metadata: {:.3}s", meta_start.elapsed().as_secs_f64());

    let field_meta = Arc::new(field_meta);
    let column_modes = Arc::new(column_modes);

    let (read_tx, read_rx) = bounded::<GpuBatch>(CHANNEL_DEPTH);

    let mut writer_txs: Vec<Sender<Vec<String>>> = Vec::with_capacity(jobs.len());
    let mut writer_threads = Vec::with_capacity(jobs.len());
    let mut reader_threads = Vec::with_capacity(jobs.len());
    let writer_limit =
        std::sync::Arc::new((std::sync::Mutex::new(0usize), std::sync::Condvar::new()));

    for (vcf_id, job) in jobs.iter().enumerate() {
        let (work_tx, work_rx) = bounded::<Vec<String>>(CHANNEL_DEPTH);
        writer_txs.push(work_tx);

        let headers = job.headers.clone();
        let output = job.output.clone();
        let use_bgzf = job.use_bgzf;
        let limit = writer_limit.clone();
        let writer = thread::spawn(move || {
            let (lock, cv) = &*limit;
            let mut active = lock.lock().unwrap();
            while *active >= max_writers {
                active = cv.wait(active).unwrap();
            }
            *active += 1;
            drop(active);

            let result = writer_thread(
                work_rx,
                headers,
                &output,
                use_bgzf,
                bgzf_level,
                mmap_output,
                mmap_no_flush,
                ram_output,
                ram_max_mb,
                timing,
                "gpu",
            );

            let (lock, cv) = &*limit;
            let mut active = lock.lock().unwrap();
            *active = active.saturating_sub(1);
            cv.notify_one();
            result
        });
        writer_threads.push(writer);

        let input = job.input.clone();
        let read_tx = read_tx.clone();
        let reader = thread::spawn(move || read_batches_gpu_multi(&input, vcf_id, read_tx, timing));
        reader_threads.push(reader);
    }
    drop(read_tx);

    let sample_maps: Vec<Arc<Vec<Option<usize>>>> =
        jobs.iter().map(|j| j.sample_map.clone()).collect();

    worker_loop_gpu_multi(
        read_rx,
        writer_txs,
        ani,
        field_meta,
        column_modes,
        sample_maps,
        info_overwrite_all,
        format_overwrite_all,
        timing,
        &mut state.gpu,
        &mut state.buffers,
    )?;

    for t in reader_threads {
        t.join().unwrap()?;
    }
    for t in writer_threads {
        t.join().unwrap()?;
    }

    eprintln!("[gpu] total: {:.3}s", total_start.elapsed().as_secs_f64());
    Ok(())
}

fn read_batches_gpu_multi(
    input: &Path,
    vcf_id: usize,
    read_tx: Sender<GpuBatch>,
    timing: bool,
) -> Result<()> {
    let input_reader = VcfAnnotationReader::open(input)?;
    let streaming_reader = StreamingVcfReader::new(input_reader);
    let (_headers, mut reader) = streaming_reader.into_headers_and_self()?;
    let mut batch = Vec::with_capacity(BATCH_SIZE);
    let mut total_lines = 0usize;
    let mut last_report = Instant::now();

    while let Some(line) = reader.read_line()? {
        if line.starts_with('#') {
            continue;
        }
        batch.push(line);

        if batch.len() >= BATCH_SIZE {
            if timing {
                total_lines += batch.len();
            }
            if read_tx
                .send(GpuBatch {
                    vcf_id,
                    lines: std::mem::replace(&mut batch, Vec::with_capacity(BATCH_SIZE)),
                })
                .is_err()
            {
                break;
            }
            if timing && last_report.elapsed().as_secs() >= 2 {
                eprintln!("[gpu] read: {} lines", total_lines);
                last_report = Instant::now();
            }
        }
    }

    if !batch.is_empty() {
        let final_len = batch.len();
        let _ = read_tx.send(GpuBatch {
            vcf_id,
            lines: batch,
        });
        if timing {
            total_lines += final_len;
        }
    }
    if timing {
        eprintln!("[gpu] read done: {} lines", total_lines);
    }
    Ok(())
}

fn worker_loop_gpu_multi(
    rx: Receiver<GpuBatch>,
    txs: Vec<Sender<Vec<String>>>,
    ani: &AniIndex,
    field_meta: Arc<HashMap<String, FieldNumber>>,
    column_modes: Arc<Vec<(String, AnnotateMode)>>,
    sample_maps: Vec<Arc<Vec<Option<usize>>>>,
    info_overwrite_all: bool,
    format_overwrite_all: bool,
    timing: bool,
    gpu: &mut GpuAni,
    buffers: &mut GpuLookupBuffers,
) -> Result<()> {
    let num_threads = rayon::current_num_threads().max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    let mut total_lines = 0usize;
    let mut parse_total = 0f64;
    let mut lookup_total = 0f64;
    let mut annotate_total = 0f64;
    let mut keys_total = 0f64;
    let mut bundle_total = 0f64;
    let mut bundle_read_total = 0f64;
    let mut bundle_info_total = 0f64;
    let mut bundle_optional_total = 0f64;
    let mut bundle_samples_total = 0f64;
    let mut send_total = 0f64;
    let mut last_report = Instant::now();

    let need_format = format_overwrite_all
        || column_modes
            .iter()
            .any(|(k, _)| k.eq_ignore_ascii_case("FMT") || k.eq_ignore_ascii_case("FORMAT"));
    let need_info = info_overwrite_all
        || column_modes.iter().any(|(k, _)| {
            !(k.eq_ignore_ascii_case("ID")
                || k.eq_ignore_ascii_case("QUAL")
                || k.eq_ignore_ascii_case("FILTER")
                || k.eq_ignore_ascii_case("FMT")
                || k.eq_ignore_ascii_case("FORMAT"))
        });
    let use_cpu_annotate = sample_maps.iter().any(|m| m.iter().any(|v| v.is_none()));
    if timing && use_cpu_annotate {
        eprintln!("[gpu] cpu-annotate: enabled");
    }

    loop {
        let first = match rx.recv() {
            Ok(v) => v,
            Err(_) => break,
        };
        if first.lines.is_empty() {
            continue;
        }
        let mut batches = Vec::new();
        let mut total_batch_lines = 0usize;
        total_batch_lines += first.lines.len();
        batches.push(first);

        while total_batch_lines < MULTI_MAX_LINES {
            match rx.try_recv() {
                Ok(b) => {
                    total_batch_lines += b.lines.len();
                    batches.push(b);
                }
                Err(crossbeam_channel::TryRecvError::Empty) => break,
                Err(crossbeam_channel::TryRecvError::Disconnected) => break,
            }
        }

        if use_cpu_annotate {
            let annotated_batches: Vec<(usize, Vec<String>)> = batches
                .iter()
                .map(|batch| {
                    let sample_map = &sample_maps[batch.vcf_id];
                    let annotated: Vec<String> = pool.install(|| {
                        batch
                            .lines
                            .par_iter()
                            .map(|line| {
                                annotate_line(
                                    line,
                                    ani,
                                    &field_meta,
                                    &column_modes,
                                    sample_map,
                                    info_overwrite_all,
                                    format_overwrite_all,
                                )
                            })
                            .collect()
                    });
                    (batch.vcf_id, annotated)
                })
                .collect();

            let send_start = Instant::now();
            for (vcf_id, lines) in annotated_batches {
                if txs[vcf_id].send(lines).is_err() {
                    return Ok(());
                }
            }
            if timing {
                send_total += send_start.elapsed().as_secs_f64();
            }
            continue;
        }

        let batch_start = Instant::now();
        let parse_start = Instant::now();
        let mins_per_batch: Vec<Vec<Option<MinParsed>>> = batches
            .iter()
            .map(|b| {
                pool.install(|| {
                    b.lines
                        .par_iter()
                        .map(|line| fast_parse_min(line))
                        .collect()
                })
            })
            .collect();
        parse_total += parse_start.elapsed().as_secs_f64();

        let keys_start = Instant::now();
        let mut keys: Vec<u64> = Vec::new();
        let mut key_batch_idx: Vec<usize> = Vec::new();
        let mut key_line_idx: Vec<usize> = Vec::new();
        let mut key_alt_idx: Vec<usize> = Vec::new();

        for (batch_idx, batch) in batches.iter().enumerate() {
            let mins = &mins_per_batch[batch_idx];
            for (line_idx, line) in batch.lines.iter().enumerate() {
                let Some(ref min) = mins[line_idx] else {
                    continue;
                };
                let bytes = line.as_bytes();
                let ref_bytes = &bytes[min.ref_range.0..min.ref_range.1];
                let ref_hash = fast_hash64(ref_bytes);
                let base = (min.chr_id as u64) << 32 | min.pos as u64;

                let (alt_start, alt_end) = min.alt_range;
                let mut alt_idx = 0usize;
                let mut start = alt_start;
                for i in alt_start..alt_end {
                    if bytes[i] == b',' {
                        let alt_hash = fast_hash64(&bytes[start..i]);
                        let key = base ^ ref_hash ^ alt_hash;
                        keys.push(key);
                        key_batch_idx.push(batch_idx);
                        key_line_idx.push(line_idx);
                        key_alt_idx.push(alt_idx);
                        alt_idx += 1;
                        start = i + 1;
                    }
                }
                if start < alt_end {
                    let alt_hash = fast_hash64(&bytes[start..alt_end]);
                    let key = base ^ ref_hash ^ alt_hash;
                    keys.push(key);
                    key_batch_idx.push(batch_idx);
                    key_line_idx.push(line_idx);
                    key_alt_idx.push(alt_idx);
                }
            }
        }
        keys_total += keys_start.elapsed().as_secs_f64();

        if keys.is_empty() {
            for batch in batches {
                if txs[batch.vcf_id].send(batch.lines).is_err() {
                    return Ok(());
                }
            }
            continue;
        }

        let lookup_start = Instant::now();
        let idxs = gpu.lookup_batch_with_buffers(&keys, buffers)?;
        lookup_total += lookup_start.elapsed().as_secs_f64();

        let bundle_start = Instant::now();
        let mut bundles_per_batch: Vec<Vec<Vec<(usize, AnnotationBundle)>>> = batches
            .iter()
            .map(|b| vec![Vec::new(); b.lines.len()])
            .collect();
        let mut bundle_cache: HashMap<u32, AnnotationBundle> = HashMap::new();

        for (i, idx) in idxs.iter().enumerate() {
            let batch_idx = key_batch_idx[i];
            let line_idx = key_line_idx[i];
            let alt_idx = key_alt_idx[i];
            if *idx == u32::MAX {
                if let Some(min) = mins_per_batch[batch_idx][line_idx].as_ref() {
                    let bytes = batches[batch_idx].lines[line_idx].as_bytes();
                    let ref_bytes = &bytes[min.ref_range.0..min.ref_range.1];
                    if let Some(alt_bytes) = alt_slice_for_idx(bytes, min.alt_range, alt_idx) {
                        if let (Ok(ref_str), Ok(alt_str)) = (
                            std::str::from_utf8(ref_bytes),
                            std::str::from_utf8(alt_bytes),
                        ) {
                            let ref_hash = fast_hash64(ref_str.as_bytes());
                            if let Some(bundle) = ani.lookup_exact_by_chr_id_opts(
                                min.chr_id,
                                min.pos,
                                ref_str,
                                ref_hash,
                                alt_str,
                                need_info,
                                need_format,
                            ) {
                                bundles_per_batch[batch_idx][line_idx].push((alt_idx, bundle));
                            }
                        }
                    }
                }
                continue;
            }
            let bundle = if let Some(cached) = bundle_cache.get(idx) {
                cached.clone()
            } else if timing {
                let (bundle, t) = ani.build_bundle_from_entry_timed_opts(
                    &ani.entries[*idx as usize],
                    need_info,
                    need_format,
                );
                bundle_read_total += t.read_s;
                bundle_info_total += t.info_s;
                bundle_optional_total += t.optional_s;
                bundle_samples_total += t.samples_s;
                bundle_cache.insert(*idx, bundle.clone());
                bundle
            } else {
                let bundle = ani.build_bundle_from_entry_opts(
                    &ani.entries[*idx as usize],
                    need_info,
                    need_format,
                );
                bundle_cache.insert(*idx, bundle.clone());
                bundle
            };
            bundles_per_batch[batch_idx][line_idx].push((alt_idx, bundle));
        }
        bundle_total += bundle_start.elapsed().as_secs_f64();

        let annotate_start = Instant::now();
        let annotated_batches: Vec<(usize, Vec<String>)> = batches
            .iter()
            .enumerate()
            .map(|(batch_idx, batch)| {
                let sample_map = &sample_maps[batch.vcf_id];
                let annotated: Vec<Option<String>> = pool.install(|| {
                    batch
                        .lines
                        .par_iter()
                        .enumerate()
                        .map(|(i, line)| {
                            if bundles_per_batch[batch_idx][i].is_empty() {
                                return None;
                            }
                            let want_format = want_format_for_line(line, need_format);
                            parse_vcf_record_simd(line, want_format).map(|mut parsed| {
                                patch_samples_from_line(&mut parsed, line);
                                annotate_record_with_bundles(
                                    &parsed,
                                    &bundles_per_batch[batch_idx][i],
                                    &field_meta,
                                    &column_modes,
                                    sample_map,
                                    info_overwrite_all,
                                    format_overwrite_all,
                                    false,
                                )
                            })
                        })
                        .collect()
                });
                let mut out = batch.lines.clone();
                for (i, val) in annotated.into_iter().enumerate() {
                    if let Some(s) = val {
                        out[i] = s;
                    }
                }
                (batch.vcf_id, out)
            })
            .collect();
        annotate_total += annotate_start.elapsed().as_secs_f64();

        let send_start = Instant::now();
        for (vcf_id, lines) in annotated_batches {
            if txs[vcf_id].send(lines).is_err() {
                return Ok(());
            }
        }
        send_total += send_start.elapsed().as_secs_f64();

        if timing {
            total_lines += key_line_idx.len();
            if last_report.elapsed().as_secs_f64() >= 2.0 {
                eprintln!(
                    "[gpu] lines: {}, parse: {:.3}s, keys: {:.3}s, lookup: {:.3}s, bundle: {:.3}s, annotate: {:.3}s, send: {:.3}s",
                    total_lines,
                    parse_total,
                    keys_total,
                    lookup_total,
                    bundle_total,
                    annotate_total,
                    send_total
                );
                eprintln!(
                    "[gpu] bundle_read: {:.3}s, bundle_info: {:.3}s, bundle_optional: {:.3}s, bundle_samples: {:.3}s",
                    bundle_read_total,
                    bundle_info_total,
                    bundle_optional_total,
                    bundle_samples_total
                );
                last_report = Instant::now();
            }
            if batch_start.elapsed().as_secs_f64() >= 2.0 {
                eprintln!(
                    "[gpu] batch: {:.3}s, keys: {}, lines: {}",
                    batch_start.elapsed().as_secs_f64(),
                    keys.len(),
                    total_batch_lines
                );
            }
            if send_start.elapsed().as_secs_f64() >= 2.0 {
                eprintln!(
                    "[gpu] send blocked: {:.3}s",
                    send_start.elapsed().as_secs_f64()
                );
            }
        }
    }

    if timing {
        eprintln!(
            "[gpu] done: lines: {}, parse: {:.3}s, keys: {:.3}s, lookup: {:.3}s, bundle: {:.3}s, annotate: {:.3}s, send: {:.3}s",
            total_lines,
            parse_total,
            keys_total,
            lookup_total,
            bundle_total,
            annotate_total,
            send_total
        );
        eprintln!("[gpu] bundle: {:.3}s", bundle_total);
        eprintln!(
            "[gpu] bundle_read: {:.3}s, bundle_info: {:.3}s, bundle_optional: {:.3}s, bundle_samples: {:.3}s",
            bundle_read_total,
            bundle_info_total,
            bundle_optional_total,
            bundle_samples_total
        );
    }

    Ok(())
}

fn read_batches_gpu(
    reader: &mut StreamingVcfReader,
    read_tx: Sender<Vec<String>>,
    timing: bool,
) -> Result<()> {
    let mut batch = Vec::with_capacity(BATCH_SIZE);
    let mut total_lines = 0usize;
    let mut last_report = Instant::now();
    let mut blocked = false;
    while let Some(line) = reader.read_line()? {
        if line.starts_with('#') {
            continue;
        }
        batch.push(line);

        if batch.len() >= BATCH_SIZE {
            if timing {
                total_lines += batch.len();
                eprintln!("[gpu] read: {} lines", total_lines);
                last_report = Instant::now();
            }
            if timing && !blocked && read_tx.is_full() {
                eprintln!("[gpu] read: channel full, waiting");
                blocked = true;
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
            blocked = false;
        }
    }

    if !batch.is_empty() {
        if timing {
            total_lines += batch.len();
            eprintln!("[gpu] read: {} lines", total_lines);
        }
        if timing && !blocked && read_tx.is_full() {
            eprintln!("[gpu] read: channel full, waiting");
        }
        let _ = read_tx.send(batch);
    }
    if timing && last_report.elapsed().as_secs_f64() >= 0.0 {
        eprintln!("[gpu] read done: {} lines", total_lines);
    }
    Ok(())
}

fn build_entry_keys(ani: &AniIndex) -> Vec<u64> {
    let mut keys = Vec::with_capacity(ani.entries.len());
    for entry in &ani.entries {
        let ref_str = ani.read_cstring(entry.ref_ofs as usize);
        let alt_str = ani.read_cstring(entry.alt_ofs as usize);
        let key = make_key(entry.chr_id, entry.pos, ref_str.as_ref(), alt_str.as_ref());
        keys.push(key);
    }
    keys
}
