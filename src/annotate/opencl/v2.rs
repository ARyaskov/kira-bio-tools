#![cfg(feature = "opencl")]

use crate::util::fast_hash64;
use anyhow::Result;
use crossbeam_channel::{Receiver, Sender, bounded};
use ocl::core::{ProgramInfo, ProgramInfoResult};
use ocl::{Context as OclContext, Device, DeviceType, Platform, ProQue, Program, Queue};
use rayon::prelude::*;
use std::collections::HashMap;
use std::ffi::CString;
use std::fs;
use std::path::Path;
use std::path::PathBuf;
use std::sync::Arc;
use std::thread;
use std::time::Instant;

use crate::annotate::constants::CHANNEL_DEPTH;
use crate::annotate::cpu_v2::field_metadata::{iter_ani_header_lines, load_and_infer_metadata};
use crate::annotate::cpu_v2::vcf_parsing::{parse_vcf_record, patch_samples_from_line};
use crate::annotate::cpu_v2::{
    ColumnSpec, annotate_record_with_bundles, build_sample_map, expand_column_specs,
    extract_samples_from_headers, merge_annotation_headers, writer_thread,
};
use crate::annotate::reader::{StreamingVcfReader, VcfAnnotationReader};
use crate::annotate::structs::ani::AniIndex;
use crate::annotate::structs::annotate_mode::AnnotateMode;
use crate::annotate::structs::bundle::{AnnotationBundle, FieldNumber};
use crate::util::{VcfFormat, chr_name_to_id, detect_format};
use kira_kv_engine::Index;

const LINE_BATCH: usize = 200_000;

pub struct OpenCLv2 {
    index: Index,
    entry_keys: Vec<u64>,
    batch_cap: usize,
}

impl OpenCLv2 {
    pub fn new(ani: &AniIndex, batch_cap: usize) -> Result<Self> {
        let index_bytes = ani.index.serialize()?;
        let index = Index::deserialize(&index_bytes)?;
        // Prefer the cached entry_keys section in .ani (same fast path as
        // the CUDA loader). Legacy ANIs without the section fall through
        // to the eager scan with a warning.
        let entry_keys = if let Some(cached) = ani.cached_entry_keys() {
            eprintln!(
                "[opencl] entry_keys loaded from .ani cache ({} keys, zero-copy mmap)",
                cached.len()
            );
            cached.to_vec()
        } else {
            eprintln!(
                "[opencl] WARNING — .ani lacks cached entry_keys section; \
                 scanning entries (slow). Rebuild .ani to make this instant."
            );
            build_entry_keys(ani)
        };
        Ok(Self {
            index,
            entry_keys,
            batch_cap,
        })
    }

    #[inline]
    pub fn run_batch(
        &mut self,
        keys: &[u64],
        timing: bool,
        write_total: &mut f64,
        kernel_total: &mut f64,
        read_total: &mut f64,
    ) -> Result<Vec<u32>> {
        let batch = keys.len();
        if batch > self.batch_cap {
            anyhow::bail!("OpenCL batch size {} exceeds cap {}", batch, self.batch_cap);
        }

        if timing {
            let start = Instant::now();
            *write_total += start.elapsed().as_secs_f64();
        }
        if timing {
            let start = Instant::now();
            *kernel_total += start.elapsed().as_secs_f64();
        }
        let out = self.lookup_batch_cpu(keys);
        if timing {
            let start = Instant::now();
            *read_total += start.elapsed().as_secs_f64();
        }
        Ok(out)
    }

    pub fn run_batch_from_strings(
        &self,
        ref_pool: &[u8],
        ref_offsets: &[u32],
        ref_lens: &[u32],
        alt_pool: &[u8],
        alt_offsets: &[u32],
        alt_lens: &[u32],
        key_ref_idx: &[u32],
        key_chr: &[u32],
        key_pos: &[u32],
        timing: bool,
        write_total: &mut f64,
        kernel_total: &mut f64,
        read_total: &mut f64,
    ) -> Result<Vec<u32>> {
        let nkeys = alt_offsets.len();
        if nkeys == 0 {
            return Ok(Vec::new());
        }

        if timing {
            let start = Instant::now();
            *write_total += start.elapsed().as_secs_f64();
        }
        if timing {
            let start = Instant::now();
            *kernel_total += start.elapsed().as_secs_f64();
        }
        let out = self.lookup_batch_from_strings_cpu(
            ref_pool,
            ref_offsets,
            ref_lens,
            alt_pool,
            alt_offsets,
            alt_lens,
            key_ref_idx,
            key_chr,
            key_pos,
        );
        if timing {
            let start = Instant::now();
            *read_total += start.elapsed().as_secs_f64();
        }
        Ok(out)
    }

    pub fn batch_cap(&self) -> usize {
        self.batch_cap
    }

    fn lookup_batch_cpu(&self, keys: &[u64]) -> Vec<u32> {
        // kira_kv_engine 0.6.0: SIMD batch lookup. AVX2 4-wide hash + 16-deep
        // prefetch chain — 2-3× faster than scalar lookup_u64 for batches of
        // 100K+ keys. Always verify the hit against entry_keys to catch lean_mph
        // false positives (we don't use lean_mph today, but the check is cheap).
        let opts = self.index.lookup_batch_u64_simd(keys);
        let mut out = vec![u32::MAX; keys.len()];
        for (i, idx_opt) in opts.into_iter().enumerate() {
            if let Some(idx) = idx_opt {
                if idx < self.entry_keys.len() && self.entry_keys[idx] == keys[i] {
                    out[i] = idx as u32;
                }
            }
        }
        out
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
        key_chr: &[u32],
        key_pos: &[u32],
    ) -> Vec<u32> {
        use crate::annotate::structs::ani::make_variant_key;
        let n = alt_offsets.len();
        // Stage 1: build the u64 key array so we can hit the SIMD batch path
        // in a single pass instead of N scalar lookup_u64 calls.
        let mut keys: Vec<u64> = Vec::with_capacity(n);
        let mut valid: Vec<bool> = Vec::with_capacity(n);
        for i in 0..n {
            let ref_idx = key_ref_idx[i] as usize;
            if ref_idx >= ref_offsets.len() || ref_idx >= ref_lens.len() {
                keys.push(0);
                valid.push(false);
                continue;
            }
            let ref_start = ref_offsets[ref_idx] as usize;
            let ref_len = ref_lens[ref_idx] as usize;
            if ref_start + ref_len > ref_pool.len() {
                keys.push(0);
                valid.push(false);
                continue;
            }
            let alt_start = alt_offsets[i] as usize;
            let alt_len = alt_lens[i] as usize;
            if alt_start + alt_len > alt_pool.len() {
                keys.push(0);
                valid.push(false);
                continue;
            }
            let chr_id = key_chr[i];
            let pos = key_pos[i];
            // Non-commutative key — identical recipe at build + CPU lookup.
            let key = make_variant_key(
                chr_id,
                pos,
                &ref_pool[ref_start..ref_start + ref_len],
                &alt_pool[alt_start..alt_start + alt_len],
            );
            keys.push(key);
            valid.push(true);
        }
        // Stage 2: single batched MPH probe.
        let opts = self.index.lookup_batch_u64_simd(&keys);
        let mut out = vec![u32::MAX; n];
        for i in 0..n {
            if !valid[i] {
                continue;
            }
            if let Some(idx) = opts[i] {
                if idx < self.entry_keys.len() && self.entry_keys[idx] == keys[i] {
                    out[i] = idx as u32;
                }
            }
        }
        out
    }
}

fn pick_opencl_device() -> Result<(Platform, Device)> {
    let platforms = Platform::list();
    let mut fallback: Option<(Platform, Device)> = None;

    for platform in platforms.iter().copied() {
        let devices = Device::list(platform, Some(DeviceType::GPU))?;
        for device in devices {
            let vendor = device.vendor().unwrap_or_default().to_ascii_lowercase();
            let name = device.name().unwrap_or_default().to_ascii_lowercase();
            let is_nvidia = vendor.contains("nvidia") || name.contains("nvidia");
            if is_nvidia {
                return Ok((platform, device));
            }
            if fallback.is_none() {
                fallback = Some((platform, device));
            }
        }
    }

    if let Some(pick) = fallback {
        return Ok(pick);
    }

    for platform in platforms.iter().copied() {
        let devices = Device::list(platform, None)?;
        if let Some(device) = devices.into_iter().next() {
            return Ok((platform, device));
        }
    }

    anyhow::bail!("No OpenCL devices found");
}

fn build_proque_with_cache(src: &str, platform: Platform, device: Device) -> Result<ProQue> {
    let vendor = device.vendor().unwrap_or_default();
    let name = device.name().unwrap_or_default();
    let dtype = device
        .info(ocl::enums::DeviceInfo::Type)
        .ok()
        .and_then(|v| match v {
            ocl::enums::DeviceInfoResult::Type(t) => Some(t),
            _ => None,
        })
        .unwrap_or(DeviceType::empty());
    eprintln!(
        "[opencl] device: {} (vendor: {}, type: {:?})",
        name, vendor, dtype
    );

    let context = OclContext::builder()
        .platform(platform)
        .devices(device.clone())
        .build()?;
    let queue = Queue::new(&context, device.clone(), None)?;
    let program = build_program_cached(&context, &device, src)?;

    Ok(ProQue::new(context, queue, program, Some(1)))
}

fn build_program_cached(context: &OclContext, device: &Device, src: &str) -> Result<Program> {
    let cmplr_opts = CString::new("").unwrap();
    let cache_path = opencl_cache_path(device, src)?;
    if let Ok(bin) = fs::read(&cache_path) {
        eprintln!("[opencl] cache: hit ({})", cache_path.to_string_lossy());
        if let Ok(program) = Program::with_binary(context, &[device.clone()], &[&bin], &cmplr_opts)
        {
            return Ok(program);
        }
        eprintln!("[opencl] cache: hit invalid, rebuilding");
    }

    eprintln!("[opencl] cache: miss ({})", cache_path.to_string_lossy());
    let src_c = CString::new(src)?;
    let program = Program::with_source(context, &[src_c], Some(&[device.clone()]), &cmplr_opts)?;

    if let Ok(ProgramInfoResult::Binaries(bins)) = program.info(ProgramInfo::Binaries) {
        if let Some(bin) = bins.first() {
            if let Some(parent) = cache_path.parent() {
                let _ = fs::create_dir_all(parent);
            }
            if fs::write(&cache_path, bin).is_ok() {
                eprintln!("[opencl] cache: wrote {}", cache_path.to_string_lossy());
            }
        }
    }

    Ok(program)
}

fn opencl_cache_path(device: &Device, src: &str) -> Result<PathBuf> {
    let vendor = device.vendor().unwrap_or_default();
    let name = device.name().unwrap_or_default();
    let driver = device.version().map(|v| v.to_string()).unwrap_or_default();

    let mut h = fast_hash64(src.as_bytes());
    h = h.wrapping_mul(0x9e3779b97f4a7c15) ^ fast_hash64(vendor.as_bytes());
    h = h.wrapping_mul(0x9e3779b97f4a7c15) ^ fast_hash64(name.as_bytes());
    h = h.wrapping_mul(0x9e3779b97f4a7c15) ^ fast_hash64(driver.as_bytes());

    let mut path = PathBuf::from("target");
    path.push("opencl_cache");
    path.push(format!("ocl_{:016x}.bin", h));
    Ok(path)
}

pub fn annotate_vcf_opencl_v2(
    ani: &AniIndex,
    input: &Path,
    output: &Path,
    columns: &[String],
    bgzf_level: Option<u32>,
    mmap_output: bool,
    mmap_no_flush: bool,
    ram_output: bool,
    ram_max_mb: u32,
    batch_cap: usize,
) -> Result<()> {
    let timing = std::env::var("KIRA_BT_TIMING").is_ok();
    let init_start = Instant::now();
    let mut gpu = OpenCLv2::new(ani, batch_cap)?;
    if timing {
        eprintln!("[opencl] init: {:.3}s", init_start.elapsed().as_secs_f64());
    }
    annotate_vcf_opencl_v2_with_gpu(
        &mut gpu,
        ani,
        input,
        output,
        columns,
        bgzf_level,
        mmap_output,
        mmap_no_flush,
        ram_output,
        ram_max_mb,
    )
}

pub fn annotate_vcf_opencl_v2_with_gpu(
    gpu: &mut OpenCLv2,
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
    let timing = std::env::var("KIRA_BT_TIMING").is_ok();
    let mut column_specs = ColumnSpec::parse_all(columns);
    let info_overwrite_all = column_specs
        .iter()
        .any(|c| c.key.eq_ignore_ascii_case("INFO") && c.mode.replace_all);
    let format_overwrite_all = column_specs.iter().any(|c| {
        (c.key.eq_ignore_ascii_case("FMT") || c.key.eq_ignore_ascii_case("FORMAT"))
            && c.mode.replace_all
    });

    let field_meta = load_and_infer_metadata(ani, false)?;
    let ani_headers = iter_ani_header_lines(ani);
    column_specs = expand_column_specs(&column_specs, &ani_headers, &field_meta);
    let column_modes: Vec<(String, AnnotateMode)> = column_specs
        .iter()
        .map(|c| (c.key.clone(), c.mode))
        .collect();

    let input_format = detect_format(input)?;
    let output_ext = output.extension().and_then(|s| s.to_str()).unwrap_or("");
    let output_wants_bgzf = matches!(output_ext, "gz" | "bgz" | "bgzf");
    let use_bgzf = output_wants_bgzf;

    let input_reader = VcfAnnotationReader::open(input)?;
    let streaming_reader = StreamingVcfReader::new(input_reader);
    let (headers, mut reader) = streaming_reader.into_headers_and_self()?;

    let merged_headers = merge_annotation_headers(&headers, &ani_headers, &column_specs)?;
    let input_samples = extract_samples_from_headers(&headers);
    let db_samples = extract_samples_from_headers(&ani_headers);
    let sample_map = build_sample_map(&input_samples, &db_samples);

    let (read_tx, read_rx) = bounded::<Vec<String>>(CHANNEL_DEPTH);
    let (work_tx, work_rx) = bounded::<Vec<String>>(CHANNEL_DEPTH);
    let ani_ref = ani;
    let field_meta = Arc::new(field_meta);
    let column_modes = Arc::new(column_modes);
    let sample_map = Arc::new(sample_map);
    let use_bgzf = use_bgzf;
    thread::scope(|s| -> Result<()> {
        let worker = s.spawn(move || {
            worker_thread_opencl(
                read_rx,
                work_tx,
                gpu,
                ani_ref,
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
                "opencl",
            )
        });

        let wait_start = Instant::now();
        read_batches_opencl(&mut reader, read_tx, timing)?;
        if timing {
            eprintln!("[opencl] waiting worker...");
        }
        worker.join().unwrap()?;
        if timing {
            eprintln!(
                "[opencl] worker done: {:.3}s",
                wait_start.elapsed().as_secs_f64()
            );
        }
        if timing {
            eprintln!("[opencl] waiting writer...");
        }
        writer.join().unwrap()?;
        if timing {
            eprintln!(
                "[opencl] writer done: {:.3}s",
                wait_start.elapsed().as_secs_f64()
            );
        }
        Ok(())
    })?;

    Ok(())
}

/// Thin wrapper over the canonical `make_variant_key` for parity with the
/// CUDA fallback. Identical wire output, used at both build and lookup time.
#[inline]
fn make_key(chr_id: u32, pos: u32, ref_allele: &str, alt: &str) -> u64 {
    crate::annotate::structs::ani::make_variant_key(
        chr_id,
        pos,
        ref_allele.as_bytes(),
        alt.as_bytes(),
    )
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

struct MinParsed {
    chr_id: u32,
    pos: u32,
    ref_range: (usize, usize),
    alt_range: (usize, usize),
}

fn fast_parse_min(line: &str, ani: &AniIndex) -> Option<MinParsed> {
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
    let chr_id = ani
        .contig_id(chrom)
        .or_else(|| chr_name_to_id(chrom).map(u32::from))?;

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

fn build_bundles_pos_index_batch(
    ani: &AniIndex,
    batch: &[String],
    mins: &[Option<MinParsed>],
    field_meta: &HashMap<String, FieldNumber>,
    need_info: bool,
    need_format: bool,
) -> Vec<Vec<(usize, AnnotationBundle)>> {
    let mut bundles_per_line: Vec<Vec<(usize, AnnotationBundle)>> = vec![Vec::new(); batch.len()];
    for (line_idx, line) in batch.iter().enumerate() {
        let Some(ref min) = mins[line_idx] else {
            continue;
        };
        let entry_indices = match ani.lookup_pos_index(min.chr_id, min.pos) {
            Some(v) => v,
            None => continue,
        };
        let bytes = line.as_bytes();
        let ref_bytes = &bytes[min.ref_range.0..min.ref_range.1];
        let ref_str = match std::str::from_utf8(ref_bytes) {
            Ok(v) => v,
            Err(_) => continue,
        };
        let alt_range = min.alt_range;
        let mut alt_idx = 0usize;
        let mut start = alt_range.0;
        for i in alt_range.0..alt_range.1 {
            if bytes[i] == b',' {
                if let Ok(alt_str) = std::str::from_utf8(&bytes[start..i]) {
                    if let Some(bundle) = find_bundle_pos_index(
                        ani,
                        entry_indices,
                        min.chr_id,
                        min.pos,
                        ref_str,
                        alt_str,
                        field_meta,
                        need_info,
                        need_format,
                    ) {
                        bundles_per_line[line_idx].push((alt_idx, bundle));
                    }
                }
                alt_idx += 1;
                start = i + 1;
            }
        }
        if start < alt_range.1 {
            if let Ok(alt_str) = std::str::from_utf8(&bytes[start..alt_range.1]) {
                if let Some(bundle) = find_bundle_pos_index(
                    ani,
                    entry_indices,
                    min.chr_id,
                    min.pos,
                    ref_str,
                    alt_str,
                    field_meta,
                    need_info,
                    need_format,
                ) {
                    bundles_per_line[line_idx].push((alt_idx, bundle));
                }
            }
        }
    }
    bundles_per_line
}

fn find_bundle_pos_index(
    ani: &AniIndex,
    entry_indices: &[u32],
    chr_id: u32,
    pos: u32,
    rf: &str,
    alt: &str,
    field_meta: &HashMap<String, FieldNumber>,
    need_info: bool,
    need_format: bool,
) -> Option<AnnotationBundle> {
    for &idx in entry_indices {
        let e = &ani.entries[idx as usize];
        if e.chr_id != chr_id || e.pos != pos {
            continue;
        }
        let rf_str = ani.read_cstring(e.ref_ofs as usize);
        if rf_str.as_ref() != rf {
            continue;
        }
        let alt_str = ani.read_cstring(e.alt_ofs as usize);
        if alt_str.as_ref() != alt {
            continue;
        }
        return Some(ani.build_bundle_from_entry_idx_opts_with_meta(
            idx as usize,
            field_meta,
            need_info,
            need_format,
        ));
    }
    None
}

fn worker_thread_opencl(
    rx: Receiver<Vec<String>>,
    tx: Sender<Vec<String>>,
    gpu: &mut OpenCLv2,
    ani: &AniIndex,
    field_meta: Arc<HashMap<String, FieldNumber>>,
    column_modes: Arc<Vec<(String, AnnotateMode)>>,
    sample_map: Arc<Vec<Option<usize>>>,
    info_overwrite_all: bool,
    format_overwrite_all: bool,
    timing: bool,
) -> Result<()> {
    let num_threads = rayon::current_num_threads().max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    let mut total_lines = 0usize;
    let mut parse_total = 0f64;
    let mut key_total = 0f64;
    let mut pack_total = 0f64;
    let mut lookup_total = 0f64;
    let mut write_total = 0f64;
    let mut kernel_total = 0f64;
    let mut read_total = 0f64;
    let mut bundle_total = 0f64;
    let mut bundle_read_total = 0f64;
    let mut bundle_info_total = 0f64;
    let mut bundle_optional_total = 0f64;
    let mut bundle_samples_total = 0f64;
    let mut annotate_total = 0f64;
    let mut send_total = 0f64;
    let mut send_max = 0f64;
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

    while let Ok(mut batch) = rx.recv() {
        if batch.is_empty() {
            continue;
        }

        let parse_start = Instant::now();
        let mins: Vec<Option<MinParsed>> =
            pool.install(|| batch.par_iter().map(|line| fast_parse_min(line, ani)).collect());
        parse_total += parse_start.elapsed().as_secs_f64();

        if ani.has_pos_index() {
            let bundle_start = Instant::now();
            let bundles_per_line = build_bundles_pos_index_batch(
                ani,
                &batch,
                &mins,
                &field_meta,
                need_info,
                need_format,
            );
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
                        parse_vcf_record(line).map(|mut parsed| {
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

            let send_start = Instant::now();
            if tx.send(batch).is_err() {
                break;
            }
            let send_elapsed = send_start.elapsed().as_secs_f64();
            send_total += send_elapsed;
            if send_elapsed > send_max {
                send_max = send_elapsed;
            }

            if timing {
                total_lines += bundles_per_line.iter().map(|v| v.len()).sum::<usize>();
                if last_report.elapsed().as_secs_f64() >= 2.0 {
                    eprintln!(
                        "[opencl] lines: {}, parse: {:.3}s, key: {:.3}s, pack: {:.3}s, lookup: {:.3}s, annotate: {:.3}s",
                        total_lines,
                        parse_total,
                        key_total,
                        pack_total,
                        lookup_total,
                        annotate_total
                    );
                    eprintln!(
                        "[opencl] kernel: {:.3}s, upload: {:.3}s, download: {:.3}s, bundle: {:.3}s, send: {:.3}s (max {:.3}s)",
                        kernel_total, write_total, read_total, bundle_total, send_total, send_max
                    );
                    last_report = Instant::now();
                }
            }
            continue;
        }

        let key_start = Instant::now();
        let mut key_count = 0usize;
        for (line_idx, line) in batch.iter().enumerate() {
            let Some(ref min) = mins[line_idx] else {
                continue;
            };
            let bytes = line.as_bytes();
            let (alt_start, alt_end) = min.alt_range;
            if alt_end <= alt_start || alt_end > bytes.len() {
                continue;
            }
            let mut count = 1usize;
            for i in alt_start..alt_end {
                if bytes[i] == b',' {
                    count += 1;
                }
            }
            key_count += count;
        }

        let use_gpu_hash = key_count > 0 && key_count <= gpu.batch_cap();

        let mut keys: Vec<u64> = Vec::new();
        let mut key_line_idx: Vec<usize> = Vec::new();
        let mut key_alt_idx: Vec<usize> = Vec::new();

        let mut ref_pool: Vec<u8> = Vec::new();
        let mut ref_offsets: Vec<u32> = vec![0; batch.len()];
        let mut ref_lens: Vec<u32> = vec![0; batch.len()];
        let mut alt_pool: Vec<u8> = Vec::new();
        let mut alt_offsets: Vec<u32> = Vec::new();
        let mut alt_lens: Vec<u32> = Vec::new();
        let mut key_ref_idx: Vec<u32> = Vec::new();
        let mut key_chr: Vec<u32> = Vec::new();
        let mut key_pos: Vec<u32> = Vec::new();

        if use_gpu_hash {
            let pack_start = Instant::now();
            for (line_idx, line) in batch.iter().enumerate() {
                let Some(ref min) = mins[line_idx] else {
                    continue;
                };
                let bytes = line.as_bytes();
                let ref_start = min.ref_range.0;
                let ref_end = min.ref_range.1;
                if ref_end <= ref_start || ref_end > bytes.len() {
                    continue;
                }
                let ref_off = ref_pool.len();
                ref_pool.extend_from_slice(&bytes[ref_start..ref_end]);
                ref_offsets[line_idx] = ref_off as u32;
                ref_lens[line_idx] = (ref_end - ref_start) as u32;

                let (alt_start, alt_end) = min.alt_range;
                let mut alt_idx = 0usize;
                let mut start = alt_start;
                for i in alt_start..alt_end {
                    if bytes[i] == b',' {
                        let alt_off = alt_pool.len();
                        alt_pool.extend_from_slice(&bytes[start..i]);
                        alt_offsets.push(alt_off as u32);
                        alt_lens.push((i - start) as u32);
                        key_ref_idx.push(line_idx as u32);
                        key_chr.push(min.chr_id);
                        key_pos.push(min.pos);
                        key_line_idx.push(line_idx);
                        key_alt_idx.push(alt_idx);
                        alt_idx += 1;
                        start = i + 1;
                    }
                }
                if start < alt_end {
                    let alt_off = alt_pool.len();
                    alt_pool.extend_from_slice(&bytes[start..alt_end]);
                    alt_offsets.push(alt_off as u32);
                    alt_lens.push((alt_end - start) as u32);
                    key_ref_idx.push(line_idx as u32);
                    key_chr.push(min.chr_id);
                    key_pos.push(min.pos);
                    key_line_idx.push(line_idx);
                    key_alt_idx.push(alt_idx);
                }
            }
            pack_total += pack_start.elapsed().as_secs_f64();
        } else {
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
        }

        key_total += key_start.elapsed().as_secs_f64();

        if (!use_gpu_hash && keys.is_empty()) || (use_gpu_hash && alt_offsets.is_empty()) {
            if tx.send(batch).is_err() {
                break;
            }
            continue;
        }

        let lookup_start = Instant::now();
        let mut bundles_per_line: Vec<Vec<(usize, AnnotationBundle)>> =
            vec![Vec::new(); batch.len()];
        if use_gpu_hash {
            let idxs = gpu.run_batch_from_strings(
                &ref_pool,
                &ref_offsets,
                &ref_lens,
                &alt_pool,
                &alt_offsets,
                &alt_lens,
                &key_ref_idx,
                &key_chr,
                &key_pos,
                timing,
                &mut write_total,
                &mut kernel_total,
                &mut read_total,
            )?;
            for (i, idx) in idxs.iter().enumerate() {
                if *idx == u32::MAX {
                    continue;
                }
                let line_idx = key_line_idx[i];
                let alt_idx = key_alt_idx[i];
                let bundle_start = Instant::now();
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
                    ani.build_bundle_from_entry_idx_opts_with_meta(
                        *idx as usize,
                        &field_meta,
                        need_info,
                        need_format,
                    )
                };
                bundle_total += bundle_start.elapsed().as_secs_f64();
                bundles_per_line[line_idx].push((alt_idx, bundle));
            }
        } else {
            let mut offset = 0usize;
            while offset < keys.len() {
                let end = (offset + gpu.batch_cap()).min(keys.len());
                let idxs = gpu.run_batch(
                    &keys[offset..end],
                    timing,
                    &mut write_total,
                    &mut kernel_total,
                    &mut read_total,
                )?;
                for (i, idx) in idxs.iter().enumerate() {
                    if *idx == u32::MAX {
                        continue;
                    }
                    let global = offset + i;
                    let line_idx = key_line_idx[global];
                    let alt_idx = key_alt_idx[global];
                    let bundle_start = Instant::now();
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
                        ani.build_bundle_from_entry_idx_opts_with_meta(
                            *idx as usize,
                            &field_meta,
                            need_info,
                            need_format,
                        )
                    };
                    bundle_total += bundle_start.elapsed().as_secs_f64();
                    bundles_per_line[line_idx].push((alt_idx, bundle));
                }
                offset = end;
            }
        }
        lookup_total += lookup_start.elapsed().as_secs_f64();

        let annotate_start = Instant::now();
        let annotated: Vec<Option<String>> = pool.install(|| {
            batch
                .par_iter()
                .enumerate()
                .map(|(i, line)| {
                    if bundles_per_line[i].is_empty() {
                        return None;
                    }
                    parse_vcf_record(line).map(|mut parsed| {
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

        let send_start = Instant::now();
        if tx.send(batch).is_err() {
            break;
        }
        let send_elapsed = send_start.elapsed().as_secs_f64();
        send_total += send_elapsed;
        if send_elapsed > send_max {
            send_max = send_elapsed;
        }

        if timing {
            total_lines += key_line_idx.len();
            if last_report.elapsed().as_secs_f64() >= 2.0 {
                eprintln!(
                    "[opencl] lines: {}, parse: {:.3}s, key: {:.3}s, pack: {:.3}s, lookup: {:.3}s, annotate: {:.3}s",
                    total_lines, parse_total, key_total, pack_total, lookup_total, annotate_total
                );
                eprintln!(
                    "[opencl] kernel: {:.3}s, upload: {:.3}s, download: {:.3}s, bundle: {:.3}s, send: {:.3}s (max {:.3}s)",
                    kernel_total, write_total, read_total, bundle_total, send_total, send_max
                );
                eprintln!(
                    "[opencl] bundle_read: {:.3}s, bundle_info: {:.3}s, bundle_optional: {:.3}s, bundle_samples: {:.3}s",
                    bundle_read_total,
                    bundle_info_total,
                    bundle_optional_total,
                    bundle_samples_total
                );
                last_report = Instant::now();
            }
        }
    }

    if timing {
        eprintln!(
            "[opencl] done: lines: {}, parse: {:.3}s, key: {:.3}s, pack: {:.3}s, lookup: {:.3}s, annotate: {:.3}s",
            total_lines, parse_total, key_total, pack_total, lookup_total, annotate_total
        );
        eprintln!(
            "[opencl] kernel: {:.3}s, upload: {:.3}s, download: {:.3}s, bundle: {:.3}s, send: {:.3}s (max {:.3}s)",
            kernel_total, write_total, read_total, bundle_total, send_total, send_max
        );
        eprintln!(
            "[opencl] bundle_read: {:.3}s, bundle_info: {:.3}s, bundle_optional: {:.3}s, bundle_samples: {:.3}s",
            bundle_read_total, bundle_info_total, bundle_optional_total, bundle_samples_total
        );
    }

    Ok(())
}

fn read_batches_opencl(
    reader: &mut StreamingVcfReader,
    read_tx: Sender<Vec<String>>,
    timing: bool,
) -> Result<()> {
    use crate::annotate::constants::{BATCH_MAX_LINES, BATCH_MIN_LINES, batch_target_bytes};

    let byte_target = batch_target_bytes();
    let mut batch: Vec<String> = Vec::with_capacity(BATCH_MIN_LINES.max(1024));
    let mut batch_bytes: usize = 0;
    let mut total_lines = 0usize;
    let mut send_total = 0f64;
    let mut send_max = 0f64;
    if timing {
        eprintln!(
            "[opencl] reader: byte-bounded batches, target {} MB",
            byte_target / (1024 * 1024)
        );
    }
    while let Some(line) = reader.read_line()? {
        if line.starts_with('#') {
            continue;
        }
        batch_bytes += line.len() + 1;
        batch.push(line);

        let oversized = batch_bytes >= byte_target && batch.len() >= BATCH_MIN_LINES;
        let too_many = batch.len() >= BATCH_MAX_LINES;
        if oversized || too_many {
            if timing {
                total_lines += batch.len();
                eprintln!("[opencl] read: {total_lines} lines");
            }
            let next_capacity = batch.len().clamp(BATCH_MIN_LINES, BATCH_MAX_LINES);
            let send_start = Instant::now();
            if read_tx
                .send(std::mem::replace(&mut batch, Vec::with_capacity(next_capacity)))
                .is_err()
            {
                break;
            }
            batch_bytes = 0;
            let send_elapsed = send_start.elapsed().as_secs_f64();
            send_total += send_elapsed;
            if send_elapsed > send_max {
                send_max = send_elapsed;
            }
        }
    }

    if !batch.is_empty() {
        if timing {
            total_lines += batch.len();
            eprintln!("[opencl] read: {} lines", total_lines);
        }
        let send_start = Instant::now();
        let _ = read_tx.send(batch);
        let send_elapsed = send_start.elapsed().as_secs_f64();
        send_total += send_elapsed;
        if send_elapsed > send_max {
            send_max = send_elapsed;
        }
    }
    if timing {
        eprintln!("[opencl] read done: {} lines", total_lines);
        eprintln!(
            "[opencl] read send: {:.3}s (max {:.3}s)",
            send_total, send_max
        );
    }
    Ok(())
}
