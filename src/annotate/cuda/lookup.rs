use anyhow::Result;
use crossbeam_channel::{Receiver, Sender, bounded};
use cust::prelude::*;
use kira_kv_engine::Index;
use rayon::prelude::*;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::thread;
use std::time::Instant;

use crate::annotate::constants::{BATCH_SIZE, CHANNEL_DEPTH};
use crate::annotate::cpu_v2::annotation::{annotate_line, annotate_record_with_bundles_and_info};
use crate::annotate::cpu_v2::field_metadata::{iter_ani_header_lines, load_and_infer_metadata};
use crate::annotate::cpu_v2::{
    ColumnSpec, annotate_record_with_bundles, build_sample_map, expand_column_specs,
    extract_samples_from_headers, merge_annotation_headers,
};
use crate::annotate::cpu_v2::{
    AnnotatedBatch, ParsedVcfRecord, ReadBatch, parse_vcf_record_simd,
    patch_samples_from_line, writer_thread, writer_thread_annotated,
};
use crate::annotate::reader::{StreamingVcfReader, VcfAnnotationReader};
use crate::annotate::structs::ani::{AniIndex, AniPosBlock, AniPosContig};
use crate::annotate::structs::annotate_mode::AnnotateMode;
use crate::annotate::structs::bundle::{AnnotationBundle, FieldNumber, FieldType};
use crate::util::{VcfFormat, chr_name_to_id, detect_format, fast_hash64};
use std::collections::HashMap;

pub struct GpuAni {
    _ctx: cust::context::Context,
    module: Module,
    stream: Stream,
    /// CPU-side `Index` clone used only when `mph_state` is `None` (the
    /// GPU export isn't available). Skipped on the happy path.
    index: Option<Index>,
    entry_keys: Vec<u64>,
    pos_index: Option<GpuPosIndex>,
    info_cache: Option<GpuInfoCache>,
    /// PtrHash25 MPH state as POD constants + device buffers. `None` when
    /// the engine isn't a PtrHash25 (e.g. AES-NI builds).
    mph_state: Option<GpuMphState>,
}

/// Device-resident PtrHash25 lookup state.
struct GpuMphState {
    prehash_seed: u64,
    mph_salt: u64,
    num_buckets: u32,
    num_slots: u64,
    prerotate: u8,
    /// One u8 per bucket (flat layout).
    pilots: DeviceBuffer<u8>,
    /// Optional block-Bloom filter (8 × u64 per block).
    bloom_words: Option<DeviceBuffer<u64>>,
    bloom_blocks: u32,
    /// Optional u16 per slot fingerprints. `None` when built with `lean_mph`.
    fingerprints: Option<DeviceBuffer<u16>>,
}

struct GpuPosIndex {
    contigs: DeviceBuffer<AniPosContig>,
    blocks: DeviceBuffer<AniPosBlock>,
    pos_offsets: DeviceBuffer<u32>,
    pos_counts: DeviceBuffer<u16>,
    contig_count: u32,
}

struct GpuInfoCache {
    tag_types: DeviceBuffer<u8>,
    entry_offsets: DeviceBuffer<u32>,
    entry_counts: DeviceBuffer<u16>,
    pair_tag_ids: DeviceBuffer<u32>,
    pair_value_off: DeviceBuffer<u32>,
    pair_value_len: DeviceBuffer<u32>,
    raw_values: DeviceBuffer<u8>,
    pair_offsets: DeviceBuffer<u32>,
    pair_counts: DeviceBuffer<u16>,
    int_values: DeviceBuffer<i32>,
    float_values: DeviceBuffer<f32>,
    str_offsets: DeviceBuffer<u32>,
    str_lens: DeviceBuffer<u32>,
    str_data: DeviceBuffer<u8>,
}

struct GpuLookupBuffers {
    out: Vec<u32>,
    capacity: usize,
    d_chr: DeviceBuffer<u32>,
    d_pos: DeviceBuffer<u32>,
    d_out_offsets: DeviceBuffer<u32>,
    d_out_counts: DeviceBuffer<u16>,
    host_offsets: Vec<u32>,
    host_counts: Vec<u16>,
    lookup_canon: Vec<u64>,
    lookup_out: Vec<Option<usize>>,
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
        let init_start = Instant::now();
        eprintln!("[gpu] init: starting (CUDA context + PTX load)...");
        let _ctx = cust::quick_init()?;
        let ptx_path = std::path::Path::new(env!("OUT_DIR")).join("ani_kernel.ptx");
        let ptx = std::fs::read_to_string(&ptx_path)
            .or_else(|_| std::fs::read_to_string("ani_kernel.ptx"))?;
        let module = Module::from_ptx(&ptx, &[])?;
        let stream = Stream::new(StreamFlags::NON_BLOCKING, None)?;
        eprintln!(
            "[gpu] init: CUDA context + PTX ready ({:.2}s)",
            init_start.elapsed().as_secs_f64()
        );

        let t = Instant::now();
        let entry_keys = if let Some(cached) = ani.cached_entry_keys() {
            eprintln!(
                "[gpu] init: entry_keys loaded from .ani cache ({:.3}s, {} keys)",
                t.elapsed().as_secs_f64(),
                cached.len()
            );
            cached.to_vec()
        } else {
            eprintln!(
                "[gpu] init: WARNING — .ani lacks cached entry_keys section, \
                 rebuilding from scratch. Rebuild the .ani to skip this on future runs."
            );
            eprintln!(
                "[gpu] init: building entry_keys ({} entries)...",
                ani.entries.len()
            );
            let keys = build_entry_keys(ani);
            eprintln!(
                "[gpu] init: entry_keys built ({:.2}s, {} keys)",
                t.elapsed().as_secs_f64(),
                keys.len()
            );
            keys
        };

        let t = Instant::now();
        eprintln!("[gpu] init: uploading PtrHash25 (pilots+bloom+fingerprints) to device...");
        let mph_state = match ani.index.gpu_export() {
            Some(exp) => Some(GpuMphState {
                prehash_seed: exp.prehash_seed,
                mph_salt: exp.mph_salt,
                num_buckets: exp.num_buckets,
                num_slots: exp.num_slots,
                prerotate: exp.prerotate,
                pilots: DeviceBuffer::from_slice(&exp.pilots)?,
                bloom_words: match &exp.bloom {
                    Some(b) => Some(DeviceBuffer::from_slice(&b.words)?),
                    None => None,
                },
                bloom_blocks: exp.bloom.as_ref().map(|b| b.blocks as u32).unwrap_or(0),
                fingerprints: match &exp.fingerprints {
                    Some(fp) => Some(DeviceBuffer::from_slice(fp)?),
                    None => None,
                },
            }),
            None => None,
        };
        eprintln!(
            "[gpu] init: PtrHash25 uploaded ({:.2}s, mph_state: {})",
            t.elapsed().as_secs_f64(),
            if mph_state.is_some() { "yes" } else { "no" }
        );

        let t = Instant::now();
        let pos_index = if let Some(pos) = ani.pos_index() {
            eprintln!(
                "[gpu] init: uploading pos-index ({} contigs)...",
                pos.contigs.len()
            );
            Some(GpuPosIndex {
                contigs: DeviceBuffer::from_slice(&pos.contigs)?,
                blocks: DeviceBuffer::from_slice(&pos.blocks)?,
                pos_offsets: DeviceBuffer::from_slice(&pos.pos_offsets)?,
                pos_counts: DeviceBuffer::from_slice(&pos.pos_counts)?,
                contig_count: pos.contigs.len() as u32,
            })
        } else {
            None
        };
        eprintln!(
            "[gpu] init: pos-index ready ({:.2}s)",
            t.elapsed().as_secs_f64()
        );

        let t = Instant::now();
        eprintln!("[gpu] init: building info_cache (one-time)...");
        let info_cache = if let Some(cache) = ani.info_cache() {
            eprintln!(
                "[gpu] init: info_cache built in {:.2}s; uploading to device...",
                t.elapsed().as_secs_f64()
            );
            let upload_t = Instant::now();
            let blob = ani.info_blob().unwrap();
            let pair_value_off: Vec<u32> = blob.pairs.iter().map(|p| p.value_off).collect();
            let pair_value_len: Vec<u32> = blob.pairs.iter().map(|p| p.value_len).collect();
            let tag_types: Vec<u8> = cache.tag_types.iter().map(field_type_to_u8).collect();
            let approx_mb = (blob.values.len()
                + cache.str_data.len()
                + cache.int_values.len() * 4
                + cache.float_values.len() * 4)
                / (1024 * 1024);
            eprintln!(
                "[gpu] init: uploading info_cache to device (~{} MB raw values + strings)...",
                approx_mb
            );
            // Pin host pages for DMA; failure is non-fatal.
            let pinned_values = pin_for_dma(&blob.values);
            let pinned_str_data = pin_for_dma(&cache.str_data);
            if pinned_values.is_some() || pinned_str_data.is_some() {
                eprintln!(
                    "[gpu] init: pinned host pages for DMA (values: {}, str_data: {})",
                    if pinned_values.is_some() { "yes" } else { "no" },
                    if pinned_str_data.is_some() { "yes" } else { "no" }
                );
            }

            let result = Some(GpuInfoCache {
                tag_types: DeviceBuffer::from_slice(&tag_types)?,
                entry_offsets: DeviceBuffer::from_slice(&cache.entry_offsets)?,
                entry_counts: DeviceBuffer::from_slice(&cache.entry_counts)?,
                pair_tag_ids: DeviceBuffer::from_slice(&cache.pair_tag_ids)?,
                pair_value_off: DeviceBuffer::from_slice(&pair_value_off)?,
                pair_value_len: DeviceBuffer::from_slice(&pair_value_len)?,
                raw_values: DeviceBuffer::from_slice(&blob.values)?,
                pair_offsets: DeviceBuffer::from_slice(&cache.pair_offsets)?,
                pair_counts: DeviceBuffer::from_slice(&cache.pair_counts)?,
                int_values: DeviceBuffer::from_slice(&cache.int_values)?,
                float_values: DeviceBuffer::from_slice(&cache.float_values)?,
                str_offsets: DeviceBuffer::from_slice(&cache.str_offsets)?,
                str_lens: DeviceBuffer::from_slice(&cache.str_lens)?,
                str_data: DeviceBuffer::from_slice(&cache.str_data)?,
            });

            if let Some(guard) = pinned_values {
                drop(guard);
            }
            if let Some(guard) = pinned_str_data {
                drop(guard);
            }

            eprintln!(
                "[gpu] init: info_cache device upload done ({:.2}s)",
                upload_t.elapsed().as_secs_f64()
            );
            result
        } else {
            eprintln!("[gpu] init: no info_blob in .ani, skipping info_cache");
            None
        };

        eprintln!(
            "[gpu] init: TOTAL {:.2}s (worker can now drain reader channel)",
            init_start.elapsed().as_secs_f64()
        );

        Ok(Self {
            _ctx,
            module,
            stream,
            index: if mph_state.is_some() {
                None
            } else {
                eprintln!(
                    "[gpu] init: no GPU MPH state (likely AES-NI build); \
                     building CPU index for fallback..."
                );
                let cpu_t = Instant::now();
                let bytes = ani.index.serialize()?;
                let idx = Index::deserialize(&bytes)?;
                drop(bytes);
                eprintln!(
                    "[gpu] init: CPU fallback index ready ({:.2}s)",
                    cpu_t.elapsed().as_secs_f64()
                );
                Some(idx)
            },
            entry_keys,
            pos_index,
            info_cache,
            mph_state,
        })
    }

    /// Whether the device-side PtrHash25 state was populated.
    pub fn has_gpu_mph(&self) -> bool {
        self.mph_state.is_some()
    }

    pub fn lookup_batch(&self, keys: &[u64]) -> Result<Vec<u32>> {
        if self.mph_state.is_some() {
            match self.lookup_batch_gpu(keys) {
                Ok(v) => return Ok(v),
                Err(e) => {
                    eprintln!(
                        "[gpu] kernel launch failed ({e}); falling back to CPU SIMD batch"
                    );
                }
            }
        }
        Ok(self.lookup_batch_cpu(keys))
    }

    fn lookup_batch_with_buffers(
        &self,
        keys: &[u64],
        buffers: &mut GpuLookupBuffers,
    ) -> Result<Vec<u32>> {
        if self.mph_state.is_some() {
            if let Ok(v) = self.lookup_batch_gpu(keys) {
                return Ok(v);
            }
        }
        Ok(self.lookup_batch_cpu_with_buffers(keys, buffers))
    }

    /// Device-side PtrHash25 lookup using `ani_ptrhash25_lookup_kernel`.
    /// Caller must validate hits against `entry_keys` (handled here).
    pub fn lookup_batch_gpu(&self, keys: &[u64]) -> Result<Vec<u32>> {
        let Some(state) = self.mph_state.as_ref() else {
            anyhow::bail!("GPU MPH state not available (AES-NI-hashed index?)");
        };
        if keys.is_empty() {
            return Ok(Vec::new());
        }
        let n = keys.len();
        let d_keys = DeviceBuffer::from_slice(keys)?;
        let mut d_out: DeviceBuffer<u32> = unsafe { DeviceBuffer::uninitialized(n)? };
        let func = self.module.get_function("ani_ptrhash25_lookup_kernel")?;
        let threads_per_block: u32 = 256;
        let blocks: u32 =
            ((n as u32) + threads_per_block - 1) / threads_per_block;
        let stream = &self.stream;
        let pilots_ptr = state.pilots.as_device_ptr();
        let bloom_ptr = state
            .bloom_words
            .as_ref()
            .map(|b| b.as_device_ptr().as_raw())
            .unwrap_or(0);
        let fingerprints_ptr = state
            .fingerprints
            .as_ref()
            .map(|fp| fp.as_device_ptr().as_raw())
            .unwrap_or(0);
        unsafe {
            cust::launch!(func<<<blocks, threads_per_block, 0, stream>>>(
                d_keys.as_device_ptr(),
                state.prehash_seed,
                state.mph_salt,
                state.num_buckets,
                state.num_slots,
                state.prerotate as u32,
                pilots_ptr,
                bloom_ptr,
                state.bloom_blocks,
                fingerprints_ptr,
                d_out.as_device_ptr(),
                n as i32,
            ))?;
        }
        stream.synchronize()?;
        let mut out_host: Vec<u32> = vec![0u32; n];
        d_out.copy_to(&mut out_host)?;
        let mut out = vec![u32::MAX; n];
        for (i, &slot) in out_host.iter().enumerate() {
            if slot != u32::MAX
                && (slot as usize) < self.entry_keys.len()
                && self.entry_keys[slot as usize] == keys[i]
            {
                out[i] = slot;
            }
        }
        Ok(out)
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
        key_chr: &[u32],
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
        let Some(index) = self.index.as_ref() else {
            eprintln!(
                "[gpu] WARN: lookup_batch_cpu invoked but CPU index not loaded; \
                 returning all-miss"
            );
            return vec![u32::MAX; keys.len()];
        };
        let opts = index.lookup_batch_u64_simd(keys);
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

    fn lookup_batch_cpu_with_buffers(
        &self,
        keys: &[u64],
        buffers: &mut GpuLookupBuffers,
    ) -> Vec<u32> {
        if keys.is_empty() {
            return Vec::new();
        }
        let _ = buffers.ensure_capacity(keys.len());
        for v in buffers.out.iter_mut().take(keys.len()) {
            *v = u32::MAX;
        }
        let n = keys.len();
        let Some(index) = self.index.as_ref() else {
            eprintln!(
                "[gpu] WARN: lookup_batch_cpu_with_buffers invoked but CPU index not loaded; \
                 returning all-miss"
            );
            return vec![u32::MAX; n];
        };
        buffers.ensure_lookup_scratch(n);
        let canon = &mut buffers.lookup_canon[..n];
        let opt_slice = &mut buffers.lookup_out[..n];
        index.lookup_batch_u64_simd_into(keys, canon, opt_slice);
        for (i, idx_opt) in opt_slice.iter().enumerate() {
            if let Some(idx) = idx_opt {
                if *idx < self.entry_keys.len() && self.entry_keys[*idx] == keys[i] {
                    buffers.out[i] = *idx as u32;
                }
            }
        }
        buffers.out[..n].to_vec()
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
            let key = make_variant_key(
                chr_id,
                pos,
                &ref_pool[ref_start..ref_start + ref_len],
                &alt_pool[alt_start..alt_start + alt_len],
            );
            keys.push(key);
            valid.push(true);
        }
        let Some(index) = self.index.as_ref() else {
            eprintln!(
                "[gpu] WARN: lookup_batch_from_strings_cpu invoked but CPU index not loaded; \
                 returning all-miss"
            );
            return vec![u32::MAX; n];
        };
        let opts = index.lookup_batch_u64_simd(&keys);
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

    pub fn has_pos_index(&self) -> bool {
        self.pos_index.is_some()
    }

    pub fn lookup_pos_batch(
        &self,
        chr_ids: &[u32],
        positions: &[u32],
        buffers: &mut GpuLookupBuffers,
    ) -> Result<(Vec<u32>, Vec<u16>)> {
        let Some(pos_index) = &self.pos_index else {
            return Ok((Vec::new(), Vec::new()));
        };
        let n = chr_ids.len();
        if n == 0 {
            return Ok((Vec::new(), Vec::new()));
        }
        let _ = buffers.ensure_capacity(n)?;

        buffers.d_chr.copy_from(chr_ids)?;
        buffers.d_pos.copy_from(positions)?;

        let threads = 256;
        let blocks = ((n as u32) + threads - 1) / threads;
        let func = self.module.get_function("ani_pos_lookup_kernel")?;

        let stream = &self.stream;
        unsafe {
            launch!(
                func<<<blocks, threads, 0, stream>>>(
                    buffers.d_chr.as_device_ptr(),
                    buffers.d_pos.as_device_ptr(),
                    pos_index.contigs.as_device_ptr(),
                    pos_index.contig_count,
                    pos_index.blocks.as_device_ptr(),
                    pos_index.pos_offsets.as_device_ptr(),
                    pos_index.pos_counts.as_device_ptr(),
                    buffers.d_out_offsets.as_device_ptr(),
                    buffers.d_out_counts.as_device_ptr(),
                    n as i32
                )
            )?;
        }
        self.stream.synchronize()?;

        buffers.host_offsets.resize(n, 0);
        buffers.host_counts.resize(n, 0);
        buffers.d_out_offsets.copy_to(&mut buffers.host_offsets)?;
        buffers.d_out_counts.copy_to(&mut buffers.host_counts)?;

        Ok((buffers.host_offsets.clone(), buffers.host_counts.clone()))
    }
}

fn field_type_to_u8(ty: &FieldType) -> u8 {
    match ty {
        FieldType::Integer => 0,
        FieldType::Float => 1,
        FieldType::String => 2,
        FieldType::Flag => 3,
    }
}

impl GpuLookupBuffers {
    fn new(capacity: usize) -> Result<Self> {
        Ok(Self {
            out: vec![0u32; capacity],
            capacity,
            d_chr: unsafe { DeviceBuffer::uninitialized(capacity)? },
            d_pos: unsafe { DeviceBuffer::uninitialized(capacity)? },
            d_out_offsets: unsafe { DeviceBuffer::uninitialized(capacity)? },
            d_out_counts: unsafe { DeviceBuffer::uninitialized(capacity)? },
            host_offsets: vec![0u32; capacity],
            host_counts: vec![0u16; capacity],
            lookup_canon: vec![0u64; capacity],
            lookup_out: vec![None; capacity],
        })
    }

    fn ensure_capacity(&mut self, capacity: usize) -> Result<()> {
        if capacity <= self.capacity {
            return Ok(());
        }
        self.out = vec![0u32; capacity];
        self.d_chr = unsafe { DeviceBuffer::uninitialized(capacity)? };
        self.d_pos = unsafe { DeviceBuffer::uninitialized(capacity)? };
        self.d_out_offsets = unsafe { DeviceBuffer::uninitialized(capacity)? };
        self.d_out_counts = unsafe { DeviceBuffer::uninitialized(capacity)? };
        self.host_offsets.resize(capacity, 0);
        self.host_counts.resize(capacity, 0);
        self.lookup_canon.resize(capacity, 0);
        self.lookup_out.resize(capacity, None);
        self.capacity = capacity;
        Ok(())
    }

    /// Grow lookup-scratch buffers (canon + out) without touching device buffers.
    fn ensure_lookup_scratch(&mut self, n: usize) {
        if self.lookup_canon.len() < n {
            self.lookup_canon.resize(n, 0);
        }
        if self.lookup_out.len() < n {
            self.lookup_out.resize(n, None);
        }
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
        .any(|c| c.key.eq_ignore_ascii_case("INFO") && c.mode.replace_all);
    let format_overwrite_all = column_specs.iter().any(|c| {
        (c.key.eq_ignore_ascii_case("FMT") || c.key.eq_ignore_ascii_case("FORMAT"))
            && c.mode.replace_all
    });

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
    let use_bgzf = output_wants_bgzf;

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

    let (read_tx, read_rx) = bounded::<ReadBatch>(CHANNEL_DEPTH);
    let (work_tx, work_rx) = bounded::<AnnotatedBatch>(CHANNEL_DEPTH);
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
            writer_thread_annotated(
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
    chr_id: u32,
    pos: u32,
    ref_range: (usize, usize),
    alt_range: (usize, usize),
}

struct GpuInfoTagSpec {
    key: String,
    tag_id: u32,
    number: FieldNumber,
    ty: FieldType,
}

struct GpuInfoMergeState {
    tags: Vec<GpuInfoTagSpec>,
    d_tag_ids: DeviceBuffer<u32>,
    d_tag_types: DeviceBuffer<u8>,
}

/// Thin wrapper over `make_variant_key`.
#[inline]
fn make_key(chr_id: u32, pos: u32, ref_allele: &str, alt: &str) -> u64 {
    crate::annotate::structs::ani::make_variant_key(
        chr_id,
        pos,
        ref_allele.as_bytes(),
        alt.as_bytes(),
    )
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

fn build_gpu_info_merge_state(
    ani: &AniIndex,
    field_meta: &HashMap<String, FieldNumber>,
    column_modes: &[(String, AnnotateMode)],
) -> Option<GpuInfoMergeState> {
    if ani.info_cache().is_none() {
        return None;
    }
    let blob = ani.info_blob()?;
    let mut id_map: HashMap<&str, u32> = HashMap::new();
    for (idx, key) in blob.dict_strings.iter().enumerate() {
        id_map.insert(key.as_str(), idx as u32);
    }
    let mut tags = Vec::new();
    for (key_raw, mode) in column_modes {
        let key = key_raw.strip_prefix("INFO/").unwrap_or(key_raw);
        if key.eq_ignore_ascii_case("ID")
            || key.eq_ignore_ascii_case("QUAL")
            || key.eq_ignore_ascii_case("FILTER")
            || key.eq_ignore_ascii_case("FMT")
            || key.eq_ignore_ascii_case("FORMAT")
        {
            continue;
        }
        if mode.replace_missing
            || mode.replace_non_missing
            || mode.set_or_append
            || mode.match_value
        {
            return None;
        }
        let number = field_meta.get(key).copied().unwrap_or(FieldNumber::One);
        if matches!(number, FieldNumber::A | FieldNumber::R | FieldNumber::G) {
            return None;
        }
        let tag_id = id_map.get(key).copied().unwrap_or(u32::MAX);
        let ty = ani
            .info_cache()
            .and_then(|cache| cache.tag_types.get(tag_id as usize).copied())
            .unwrap_or(FieldType::String);
        tags.push(GpuInfoTagSpec {
            key: key.to_string(),
            tag_id,
            number,
            ty,
        });
    }
    if tags.is_empty() {
        return None;
    }
    let tag_ids: Vec<u32> = tags.iter().map(|t| t.tag_id).collect();
    let tag_types: Vec<u8> = tags.iter().map(|t| field_type_to_u8(&t.ty)).collect();
    let d_tag_ids = DeviceBuffer::from_slice(&tag_ids).ok()?;
    let d_tag_types = DeviceBuffer::from_slice(&tag_types).ok()?;
    Some(GpuInfoMergeState {
        tags,
        d_tag_ids,
        d_tag_types,
    })
}

fn build_bundles_pos_index_batch(
    ani: &AniIndex,
    batch: &[&str],
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

fn build_bundles_from_pos_results(
    ani: &AniIndex,
    batch: &[&str],
    mins: &[Option<MinParsed>],
    offsets: &[u32],
    counts: &[u16],
    entry_indices: &[u32],
    field_meta: &HashMap<String, FieldNumber>,
    need_info: bool,
    need_format: bool,
) -> Vec<Vec<(usize, AnnotationBundle)>> {
    let mut bundles_per_line: Vec<Vec<(usize, AnnotationBundle)>> = vec![Vec::new(); batch.len()];
    for (line_idx, line) in batch.iter().enumerate() {
        let Some(ref min) = mins[line_idx] else {
            continue;
        };
        let count = *counts.get(line_idx).unwrap_or(&0) as usize;
        if count == 0 {
            continue;
        }
        let offset = *offsets.get(line_idx).unwrap_or(&0) as usize;
        let end = offset + count;
        if end > entry_indices.len() {
            continue;
        }
        let entries = &entry_indices[offset..end];
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
                        entries,
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
                    entries,
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

fn build_entry_matches_from_pos_results(
    ani: &AniIndex,
    batch: &[&str],
    mins: &[Option<MinParsed>],
    offsets: &[u32],
    counts: &[u16],
    entry_indices: &[u32],
) -> Vec<Vec<Option<u32>>> {
    let mut matches: Vec<Vec<Option<u32>>> = vec![Vec::new(); batch.len()];
    for (line_idx, line) in batch.iter().enumerate() {
        let Some(ref min) = mins[line_idx] else {
            continue;
        };
        let count = *counts.get(line_idx).unwrap_or(&0) as usize;
        if count == 0 {
            continue;
        }
        let offset = *offsets.get(line_idx).unwrap_or(&0) as usize;
        let end = offset + count;
        if end > entry_indices.len() {
            continue;
        }
        let entries = &entry_indices[offset..end];
        let bytes = line.as_bytes();
        let ref_bytes = &bytes[min.ref_range.0..min.ref_range.1];
        let ref_str = match std::str::from_utf8(ref_bytes) {
            Ok(v) => v,
            Err(_) => continue,
        };
        let alt_range = min.alt_range;
        let mut start = alt_range.0;
        for i in alt_range.0..alt_range.1 {
            if bytes[i] == b',' {
                let entry = if let Ok(alt_str) = std::str::from_utf8(&bytes[start..i]) {
                    find_entry_pos_index(ani, entries, min.chr_id, min.pos, ref_str, alt_str)
                } else {
                    None
                };
                matches[line_idx].push(entry);
                start = i + 1;
            }
        }
        if start < alt_range.1 {
            let entry = if let Ok(alt_str) = std::str::from_utf8(&bytes[start..alt_range.1]) {
                find_entry_pos_index(ani, entries, min.chr_id, min.pos, ref_str, alt_str)
            } else {
                None
            };
            matches[line_idx].push(entry);
        }
    }
    matches
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

fn find_entry_pos_index(
    ani: &AniIndex,
    entry_indices: &[u32],
    chr_id: u32,
    pos: u32,
    rf: &str,
    alt: &str,
) -> Option<u32> {
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
        return Some(idx);
    }
    None
}

fn gpu_merge_info_batch(
    gpu: &GpuAni,
    state: &GpuInfoMergeState,
    alt_entry_idx: &[u32],
    alt_offsets: &[u32],
    alt_counts: &[u16],
) -> Result<Vec<u32>> {
    let Some(cache) = &gpu.info_cache else {
        return Ok(Vec::new());
    };
    let n_records = alt_offsets.len();
    let n_tags = state.tags.len();
    if n_records == 0 || n_tags == 0 {
        return Ok(Vec::new());
    }
    let total = n_records * n_tags;
    let d_alt_entry_idx = DeviceBuffer::from_slice(alt_entry_idx)?;
    let d_alt_offsets = DeviceBuffer::from_slice(alt_offsets)?;
    let d_alt_counts = DeviceBuffer::from_slice(alt_counts)?;
    let d_out = unsafe { DeviceBuffer::uninitialized(total)? };
    let func = gpu.module.get_function("ani_info_merge_kernel")?;
    let threads = 256;
    let blocks = ((total as u32) + threads - 1) / threads;
    let stream = &gpu.stream;
    unsafe {
        launch!(
            func<<<blocks, threads, 0, stream>>>(
                d_alt_entry_idx.as_device_ptr(),
                d_alt_offsets.as_device_ptr(),
                d_alt_counts.as_device_ptr(),
                cache.entry_offsets.as_device_ptr(),
                cache.entry_counts.as_device_ptr(),
                cache.pair_tag_ids.as_device_ptr(),
                cache.pair_value_off.as_device_ptr(),
                cache.pair_value_len.as_device_ptr(),
                cache.raw_values.as_device_ptr(),
                cache.tag_types.as_device_ptr(),
                state.d_tag_ids.as_device_ptr(),
                d_out.as_device_ptr(),
                n_records as i32,
                n_tags as i32
            )
        )?;
    }
    gpu.stream.synchronize()?;
    let mut out = vec![0u32; total];
    d_out.copy_to(&mut out)?;
    Ok(out)
}

fn format_info_from_pairs(
    ani: &AniIndex,
    state: &GpuInfoMergeState,
    pair_idx: &[u32],
    n_records: usize,
) -> Vec<Option<String>> {
    let n_tags = state.tags.len();
    let mut out = vec![None; n_records];
    for rec in 0..n_records {
        let mut parts: Vec<String> = Vec::new();
        for (tag_i, tag) in state.tags.iter().enumerate() {
            let idx = pair_idx[rec * n_tags + tag_i];
            if idx == u32::MAX {
                continue;
            }
            let val = if tag.number == FieldNumber::Zero || tag.ty == FieldType::Flag {
                Some(String::new())
            } else {
                format_pair_value(ani, tag.ty, idx as usize)
            };
            if let Some(v) = val {
                if v.is_empty() {
                    parts.push(tag.key.clone());
                } else {
                    parts.push(format!("{}={}", tag.key, v));
                }
            }
        }
        if !parts.is_empty() {
            out[rec] = Some(parts.join(";"));
        }
    }
    out
}

fn format_pair_value(ani: &AniIndex, _ty: FieldType, idx: usize) -> Option<String> {
    let blob = ani.info_blob()?;
    if idx >= blob.pairs.len() {
        return None;
    }
    let pair = &blob.pairs[idx];
    if pair.value_len == 0 {
        return None;
    }
    let start = pair.value_off as usize;
    let end = start + pair.value_len as usize;
    if end > blob.values.len() {
        return None;
    }
    let raw = std::str::from_utf8(&blob.values[start..end]).unwrap_or("");
    if raw.is_empty() || raw == "." {
        return None;
    }
    Some(raw.to_string())
}

fn worker_thread_gpu(
    rx: Receiver<ReadBatch>,
    tx: Sender<AnnotatedBatch>,
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
    rx: Receiver<ReadBatch>,
    tx: Sender<AnnotatedBatch>,
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
    let info_merge_state = if need_info {
        build_gpu_info_merge_state(ani, &field_meta, column_modes.as_ref())
    } else {
        None
    };
    // CPU-annotate fallback needed only when FORMAT must merge AND
    // sample_map has gaps. For INFO-only annotation we stay on GPU.
    let use_cpu_annotate =
        need_format && sample_map.iter().any(|v| v.is_none());
    if timing {
        if use_cpu_annotate {
            eprintln!(
                "[gpu] cpu-annotate: enabled (need_format + sample_map gaps)"
            );
        } else if need_format && !sample_map.iter().any(|v| v.is_none()) {
            eprintln!("[gpu] cpu-annotate: not needed (format + samples align)");
        } else {
            eprintln!("[gpu] cpu-annotate: not needed (INFO-only annotation)");
        }
    }
    let mut first_batch_logged = false;
    while let Ok(mut batch) = rx.recv() {
        if batch.is_empty() {
            continue;
        }
        if !first_batch_logged {
            eprintln!(
                "[gpu] first batch arrived ({} lines), entering pipeline...",
                batch.len()
            );
            first_batch_logged = true;
        }

        let lines: Vec<&str> = (0..batch.len()).map(|i| batch.line(i)).collect();

        if use_cpu_annotate {
            let output: Vec<Option<String>> = pool.install(|| {
                lines
                    .par_iter()
                    .map(|line| {
                        let out = annotate_line(
                            line,
                            ani,
                            &field_meta,
                            &column_modes,
                            &sample_map,
                            info_overwrite_all,
                            format_overwrite_all,
                        );
                        if out == *line { None } else { Some(out) }
                    })
                    .collect()
            });
            let send_start = Instant::now();
            if tx
                .send(AnnotatedBatch {
                    input: batch,
                    output,
                })
                .is_err()
            {
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
            pool.install(|| lines.par_iter().map(|line| fast_parse_min(line, ani)).collect());
        parse_total += parse_start.elapsed().as_secs_f64();

        // `KIRA_BT_GPU_FORCE_MPH=1` skips pos-index, forcing the MPH kernel.
        let force_mph = std::env::var("KIRA_BT_GPU_FORCE_MPH").is_ok();
        if ani.has_pos_index() && !force_mph {
            let mut chr_ids = Vec::with_capacity(lines.len());
            let mut positions = Vec::with_capacity(lines.len());
            for min in &mins {
                if let Some(m) = min {
                    chr_ids.push(m.chr_id);
                    positions.push(m.pos);
                } else {
                    chr_ids.push(u32::MAX);
                    positions.push(0);
                }
            }
            let lookup_start = Instant::now();
            let (offsets, counts) = if gpu.has_pos_index() {
                gpu.lookup_pos_batch(&chr_ids, &positions, buffers)?
            } else {
                (Vec::new(), Vec::new())
            };
            lookup_total += lookup_start.elapsed().as_secs_f64();

            let parse_records_start = Instant::now();
            let all_info_empty = if info_overwrite_all {
                true
            } else {
                lines.par_iter().all(|line| {
                    let bytes = line.as_bytes();
                    let mut tabs = 0usize;
                    let mut info_start = 0usize;
                    let mut info_end = bytes.len();
                    for (i, &b) in bytes.iter().enumerate() {
                        if b == b'\t' {
                            tabs += 1;
                            if tabs == 7 {
                                info_start = i + 1;
                            } else if tabs == 8 {
                                info_end = i;
                                break;
                            }
                        }
                    }
                    let info = &bytes[info_start..info_end];
                    info.is_empty() || info == b"."
                })
            };
            parse_total += parse_records_start.elapsed().as_secs_f64();

            let use_gpu_info_merge = info_merge_state.is_some()
                && gpu.info_cache.is_some()
                && gpu.has_pos_index()
                && (info_overwrite_all || all_info_empty);

            let entry_indices = ani.pos_index().unwrap().entry_indices.as_slice();
            let bundle_start = Instant::now();
            let entry_matches = if use_gpu_info_merge {
                build_entry_matches_from_pos_results(
                    ani,
                    &lines,
                    &mins,
                    &offsets,
                    &counts,
                    entry_indices,
                )
            } else {
                Vec::new()
            };
            let bundles_per_line = if use_gpu_info_merge {
                build_bundles_from_pos_results(
                    ani,
                    &lines,
                    &mins,
                    &offsets,
                    &counts,
                    entry_indices,
                    &field_meta,
                    false,
                    need_format,
                )
            } else if gpu.has_pos_index() {
                build_bundles_from_pos_results(
                    ani,
                    &lines,
                    &mins,
                    &offsets,
                    &counts,
                    entry_indices,
                    &field_meta,
                    need_info,
                    need_format,
                )
            } else {
                build_bundles_pos_index_batch(
                    ani,
                    &lines,
                    &mins,
                    &field_meta,
                    need_info,
                    need_format,
                )
            };
            bundle_total += bundle_start.elapsed().as_secs_f64();

            let merged_info = if use_gpu_info_merge {
                let mut alt_offsets = Vec::with_capacity(lines.len());
                let mut alt_counts = Vec::with_capacity(lines.len());
                let mut flat_entries = Vec::new();
                for (i, line) in lines.iter().enumerate() {
                    let alt_count = mins[i]
                        .as_ref()
                        .map(|m| {
                            let bytes = line.as_bytes();
                            let (s, e) = m.alt_range;
                            if s >= e {
                                0
                            } else {
                                bytes[s..e].iter().filter(|&&b| b == b',').count() + 1
                            }
                        })
                        .unwrap_or(0);
                    alt_offsets.push(flat_entries.len() as u32);
                    alt_counts.push(alt_count as u16);
                    let matches = entry_matches.get(i);
                    for a in 0..alt_count {
                        let v = matches
                            .and_then(|m| m.get(a).and_then(|v| *v))
                            .unwrap_or(u32::MAX);
                        flat_entries.push(v);
                    }
                }
                if let Some(state) = &info_merge_state {
                    let pair_idx =
                        gpu_merge_info_batch(gpu, state, &flat_entries, &alt_offsets, &alt_counts)?;
                    format_info_from_pairs(ani, state, &pair_idx, lines.len())
                } else {
                    vec![None; lines.len()]
                }
            } else {
                vec![None; lines.len()]
            };

            let annotate_start = Instant::now();
            // Output: Some(s) = annotated/changed; None = unchanged
            // (writer reuses input bytes). Matches the CPU pipeline
            // shape.
            let output: Vec<Option<String>> = pool.install(|| {
                lines
                    .par_iter()
                    .enumerate()
                    .map(|(i, line)| {
                        if bundles_per_line[i].is_empty() {
                            return None;
                        }
                        let want_format = want_format_for_line(line, need_format);
                        let mut parsed = parse_vcf_record_simd(line, want_format)?;
                        patch_samples_from_line(&mut parsed, line);
                        if use_gpu_info_merge {
                            let info = merged_info[i].as_deref().unwrap_or("");
                            Some(annotate_record_with_bundles_and_info(
                                &parsed,
                                &bundles_per_line[i],
                                &field_meta,
                                &column_modes,
                                &sample_map,
                                info_overwrite_all,
                                format_overwrite_all,
                                false,
                                Some(info),
                            ))
                        } else {
                            Some(annotate_record_with_bundles(
                                &parsed,
                                &bundles_per_line[i],
                                &field_meta,
                                &column_modes,
                                &sample_map,
                                info_overwrite_all,
                                format_overwrite_all,
                                false,
                            ))
                        }
                    })
                    .collect()
            });
            annotate_total += annotate_start.elapsed().as_secs_f64();

            let batch_len = batch.len();

            // Drop the `lines` view *before* moving `batch` into the
            // channel — `lines` borrows from `batch`, so we need the
            // borrow to end first.
            drop(lines);

            let send_start = Instant::now();
            if tx
                .send(AnnotatedBatch {
                    input: batch,
                    output,
                })
                .is_err()
            {
                break;
            }
            let send_elapsed = send_start.elapsed().as_secs_f64();
            send_total += send_elapsed;

            if timing {
                total_lines += bundles_per_line.iter().map(|v| v.len()).sum::<usize>();
                // Force a first-batch summary so the user sees throughput
                // immediately rather than after 2s of silence. Also print
                // batch wall-clock so multi-second per-batch processing
                // is visible.
                let force_first = total_lines > 0 && parse_total + lookup_total + bundle_total + annotate_total > 0.0
                    && last_report.elapsed().as_secs_f64() >= 2.0;
                if force_first {
                    eprintln!(
                        "[gpu] lines: {}, parse: {:.3}s, lookup: {:.3}s, bundle: {:.3}s, annotate: {:.3}s, send: {:.3}s",
                        total_lines,
                        parse_total,
                        lookup_total,
                        bundle_total,
                        annotate_total,
                        send_total
                    );
                    last_report = Instant::now();
                }
                // First batch unconditionally: 1 line of stats so the user
                // can immediately see what fraction of the wall-clock
                // budget goes where.
                if total_lines == batch_len {
                    eprintln!(
                        "[gpu] first batch processed: {} lines, parse: {:.3}s, lookup: {:.3}s, bundle: {:.3}s, annotate: {:.3}s",
                        batch_len,
                        parse_total,
                        lookup_total,
                        bundle_total,
                        annotate_total,
                    );
                }
            }
            continue;
        }

        let keys_start = Instant::now();
        let mut keys: Vec<u64> = Vec::new();
        let mut key_line_idx: Vec<usize> = Vec::new();
        let mut key_alt_idx: Vec<usize> = Vec::new();

        // **Bugfix.** Previously this site used the legacy XOR-commutative
        // recipe `(chr<<32 | pos) ^ fxhash(ref) ^ fxhash(alt)`. The MPH was
        // built with the new non-commutative `make_variant_key`, so XOR keys
        // would miss 100 % of the time on GPU → every lookup fell back to
        // pos_index, silently turning the GPU kernel into dead code. Use
        // the same key constructor the builder + CPU lookup use.
        use crate::annotate::structs::ani::make_variant_key;
        for (line_idx, line) in batch.iter().enumerate() {
            let Some(ref min) = mins[line_idx] else {
                continue;
            };
            let bytes = line.as_bytes();
            let ref_bytes = &bytes[min.ref_range.0..min.ref_range.1];

            let (alt_start, alt_end) = min.alt_range;
            let mut alt_idx = 0usize;
            let mut start = alt_start;
            for i in alt_start..alt_end {
                if bytes[i] == b',' {
                    let key = make_variant_key(min.chr_id, min.pos, ref_bytes, &bytes[start..i]);
                    keys.push(key);
                    key_line_idx.push(line_idx);
                    key_alt_idx.push(alt_idx);
                    alt_idx += 1;
                    start = i + 1;
                }
            }
            if start < alt_end {
                let key =
                    make_variant_key(min.chr_id, min.pos, ref_bytes, &bytes[start..alt_end]);
                keys.push(key);
                key_line_idx.push(line_idx);
                key_alt_idx.push(alt_idx);
            }
        }
        keys_total += keys_start.elapsed().as_secs_f64();

        if keys.is_empty() {
            // No keys generated → no annotations → pass batch through
            // unchanged. AnnotatedBatch with all-None outputs means the
            // writer just emits the original input bytes verbatim.
            let n = batch.len();
            drop(lines);
            if tx
                .send(AnnotatedBatch {
                    input: batch,
                    output: vec![None; n],
                })
                .is_err()
            {
                break;
            }
            continue;
        }

        let lookup_start = Instant::now();
        let idxs = gpu.lookup_batch_with_buffers(&keys, buffers)?;
        lookup_total += lookup_start.elapsed().as_secs_f64();
        let bundle_start = Instant::now();
        let mut bundles_per_line: Vec<Vec<(usize, AnnotationBundle)>> =
            vec![Vec::new(); lines.len()];
        for (i, idx) in idxs.iter().enumerate() {
            let line_idx = key_line_idx[i];
            let alt_idx = key_alt_idx[i];
            if *idx == u32::MAX {
                if let Some(min) = mins[line_idx].as_ref() {
                    let bytes = lines[line_idx].as_bytes();
                    let ref_bytes = &bytes[min.ref_range.0..min.ref_range.1];
                    if let Some(alt_bytes) = alt_slice_for_idx(bytes, min.alt_range, alt_idx) {
                        if let (Ok(ref_str), Ok(alt_str)) = (
                            std::str::from_utf8(ref_bytes),
                            std::str::from_utf8(alt_bytes),
                        ) {
                            let ref_hash = fast_hash64(ref_str.as_bytes());
                            if let Some(bundle) = ani.lookup_exact_by_chr_id_pos_index_opts(
                                min.chr_id,
                                min.pos,
                                ref_str,
                                ref_hash,
                                alt_str,
                                &field_meta,
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
                ani.build_bundle_from_entry_idx_opts_with_meta(
                    *idx as usize,
                    &field_meta,
                    need_info,
                    need_format,
                )
            };
            bundles_per_line[line_idx].push((alt_idx, bundle));
        }
        bundle_total += bundle_start.elapsed().as_secs_f64();

        let annotate_start = Instant::now();
        let output: Vec<Option<String>> = pool.install(|| {
            lines
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

        let batch_len = batch.len();
        drop(lines);
        let send_start = Instant::now();
        if tx
            .send(AnnotatedBatch {
                input: batch,
                output,
            })
            .is_err()
        {
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
            bundle_read_total, bundle_info_total, bundle_optional_total, bundle_samples_total
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
        .any(|c| c.key.eq_ignore_ascii_case("INFO") && c.mode.replace_all);
    let format_overwrite_all = column_specs.iter().any(|c| {
        (c.key.eq_ignore_ascii_case("FMT") || c.key.eq_ignore_ascii_case("FORMAT"))
            && c.mode.replace_all
    });

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
    let use_bgzf = output_wants_bgzf;

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

    let (read_tx, read_rx) = bounded::<ReadBatch>(CHANNEL_DEPTH);
    let (work_tx, work_rx) = bounded::<AnnotatedBatch>(CHANNEL_DEPTH);
    let field_meta = Arc::new(field_meta);
    let column_modes = Arc::new(column_modes);
    let sample_map = Arc::new(sample_map);
    let use_bgzf = use_bgzf;

    let output_path = output.to_path_buf();
    let headers = merged_headers;
    let writer = thread::spawn(move || {
        writer_thread_annotated(
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
        .any(|c| c.key.eq_ignore_ascii_case("INFO") && c.mode.replace_all);
    let format_overwrite_all = column_specs.iter().any(|c| {
        (c.key.eq_ignore_ascii_case("FMT") || c.key.eq_ignore_ascii_case("FORMAT"))
            && c.mode.replace_all
    });

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
    use crate::annotate::constants::{BATCH_MAX_LINES, BATCH_MIN_LINES, batch_target_bytes};

    let input_reader = VcfAnnotationReader::open(input)?;
    let streaming_reader = StreamingVcfReader::new(input_reader);
    let (_headers, mut reader) = streaming_reader.into_headers_and_self()?;
    let byte_target = batch_target_bytes();
    let mut batch: Vec<String> = Vec::with_capacity(BATCH_MIN_LINES.max(1024));
    let mut batch_bytes: usize = 0;
    let mut total_lines = 0usize;
    let mut last_report = Instant::now();

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
            }
            let next_capacity = batch.len().clamp(BATCH_MIN_LINES, BATCH_MAX_LINES);
            if read_tx
                .send(GpuBatch {
                    vcf_id,
                    lines: std::mem::replace(&mut batch, Vec::with_capacity(next_capacity)),
                })
                .is_err()
            {
                break;
            }
            batch_bytes = 0;
            if timing && last_report.elapsed().as_secs() >= 2 {
                eprintln!("[gpu] read: {total_lines} lines (vcf_id={vcf_id})");
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
        eprintln!("[gpu] read done: {total_lines} lines (vcf_id={vcf_id})");
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
    let use_cpu_annotate =
        need_format && sample_maps.iter().any(|m| m.iter().any(|v| v.is_none()));
    if timing && use_cpu_annotate {
        eprintln!("[gpu] cpu-annotate: enabled (multi-job, need_format + gaps)");
    } else if timing {
        eprintln!("[gpu] cpu-annotate: not needed (multi-job)");
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
                        .map(|line| fast_parse_min(line, ani))
                        .collect()
                })
            })
            .collect();
        parse_total += parse_start.elapsed().as_secs_f64();

        if ani.has_pos_index() {
            let bundle_start = Instant::now();
            let mut bundles_per_batch: Vec<Vec<Vec<(usize, AnnotationBundle)>>> = Vec::new();
            for (batch_idx, batch) in batches.iter().enumerate() {
                let mins = &mins_per_batch[batch_idx];
                let lines_view: Vec<&str> = batch.lines.iter().map(|s| s.as_str()).collect();
                let bundles = if gpu.has_pos_index() {
                    let mut chr_ids = Vec::with_capacity(batch.lines.len());
                    let mut positions = Vec::with_capacity(batch.lines.len());
                    for min in mins {
                        if let Some(m) = min {
                            chr_ids.push(m.chr_id);
                            positions.push(m.pos);
                        } else {
                            chr_ids.push(u32::MAX);
                            positions.push(0);
                        }
                    }
                    let (offsets, counts) = gpu.lookup_pos_batch(&chr_ids, &positions, buffers)?;
                    let entry_indices = ani.pos_index().unwrap().entry_indices.as_slice();
                    build_bundles_from_pos_results(
                        ani,
                        &lines_view,
                        mins,
                        &offsets,
                        &counts,
                        entry_indices,
                        &field_meta,
                        need_info,
                        need_format,
                    )
                } else {
                    build_bundles_pos_index_batch(
                        ani,
                        &lines_view,
                        mins,
                        &field_meta,
                        need_info,
                        need_format,
                    )
                };
                bundles_per_batch.push(bundles);
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
                total_lines += bundles_per_batch
                    .iter()
                    .map(|b| b.iter().map(|v| v.len()).sum::<usize>())
                    .sum::<usize>();
                if last_report.elapsed().as_secs_f64() >= 2.0 {
                    eprintln!(
                        "[gpu] lines: {}, parse: {:.3}s, lookup: {:.3}s, bundle: {:.3}s, annotate: {:.3}s, send: {:.3}s",
                        total_lines,
                        parse_total,
                        lookup_total,
                        bundle_total,
                        annotate_total,
                        send_total
                    );
                    last_report = Instant::now();
                }
            }
            continue;
        }

        let keys_start = Instant::now();
        let mut keys: Vec<u64> = Vec::new();
        let mut key_batch_idx: Vec<usize> = Vec::new();
        let mut key_line_idx: Vec<usize> = Vec::new();
        let mut key_alt_idx: Vec<usize> = Vec::new();

        use crate::annotate::structs::ani::make_variant_key;
        for (batch_idx, batch) in batches.iter().enumerate() {
            let mins = &mins_per_batch[batch_idx];
            for (line_idx, line) in batch.lines.iter().enumerate() {
                let Some(ref min) = mins[line_idx] else {
                    continue;
                };
                let bytes = line.as_bytes();
                let ref_bytes = &bytes[min.ref_range.0..min.ref_range.1];

                let (alt_start, alt_end) = min.alt_range;
                let mut alt_idx = 0usize;
                let mut start = alt_start;
                for i in alt_start..alt_end {
                    if bytes[i] == b',' {
                        let key =
                            make_variant_key(min.chr_id, min.pos, ref_bytes, &bytes[start..i]);
                        keys.push(key);
                        key_batch_idx.push(batch_idx);
                        key_line_idx.push(line_idx);
                        key_alt_idx.push(alt_idx);
                        alt_idx += 1;
                        start = i + 1;
                    }
                }
                if start < alt_end {
                    let key = make_variant_key(
                        min.chr_id,
                        min.pos,
                        ref_bytes,
                        &bytes[start..alt_end],
                    );
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
                            if let Some(bundle) = ani.lookup_exact_by_chr_id_pos_index_opts(
                                min.chr_id,
                                min.pos,
                                ref_str,
                                ref_hash,
                                alt_str,
                                &field_meta,
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
                let bundle = ani.build_bundle_from_entry_idx_opts_with_meta(
                    *idx as usize,
                    &field_meta,
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
            bundle_read_total, bundle_info_total, bundle_optional_total, bundle_samples_total
        );
    }

    Ok(())
}

fn read_batches_gpu(
    reader: &mut StreamingVcfReader,
    read_tx: Sender<ReadBatch>,
    timing: bool,
) -> Result<()> {
    use crate::annotate::constants::{BATCH_MAX_LINES, BATCH_MIN_LINES, batch_target_bytes};

    let byte_target = batch_target_bytes();
    let bytes_cap = byte_target + 16 * 1024;
    let lines_cap = BATCH_MIN_LINES.max(1024);
    let mut batch = ReadBatch::with_capacity(bytes_cap, lines_cap);
    let mut total_lines = 0usize;
    let mut last_report = Instant::now();
    let mut blocked = false;

    if timing {
        eprintln!(
            "[gpu] reader: byte-bounded batches, target {} MB (min {} lines, max {} lines)",
            byte_target / (1024 * 1024),
            BATCH_MIN_LINES,
            BATCH_MAX_LINES
        );
    }

    while reader.read_line_into_batch(&mut batch)? {
        if batch.last_line_first_byte() == Some(b'#') {
            batch.pop_last_line();
            continue;
        }

        let oversized = batch.byte_len() >= byte_target && batch.len() >= BATCH_MIN_LINES;
        let too_many = batch.len() >= BATCH_MAX_LINES;
        if oversized || too_many {
            if timing {
                total_lines += batch.len();
                eprintln!("[gpu] read: {total_lines} lines");
                last_report = Instant::now();
            }
            if timing && !blocked && read_tx.is_full() {
                eprintln!("[gpu] read: channel full, waiting");
                blocked = true;
            }
            let next_lines_cap = batch.len().clamp(BATCH_MIN_LINES, BATCH_MAX_LINES);
            let next = ReadBatch::with_capacity(bytes_cap, next_lines_cap);
            if read_tx.send(std::mem::replace(&mut batch, next)).is_err() {
                break;
            }
            blocked = false;
        }
    }

    if !batch.is_empty() {
        if timing {
            total_lines += batch.len();
            eprintln!("[gpu] read: {total_lines} lines");
        }
        if timing && !blocked && read_tx.is_full() {
            eprintln!("[gpu] read: channel full, waiting");
        }
        let _ = read_tx.send(batch);
    }
    if timing && last_report.elapsed().as_secs_f64() >= 0.0 {
        eprintln!("[gpu] read done: {total_lines} lines");
    }
    Ok(())
}

/// RAII guard: unregisters a previously-pinned host memory region on drop.
struct PinnedHostRegion {
    ptr: *mut std::ffi::c_void,
}

unsafe impl Send for PinnedHostRegion {}
unsafe impl Sync for PinnedHostRegion {}

impl Drop for PinnedHostRegion {
    fn drop(&mut self) {
        unsafe {
            let _ = cust::sys::cuMemHostUnregister(self.ptr);
        }
    }
}

/// Pins `slice`'s host memory pages for DMA. Returns `None` on any error
/// (caller silently falls back to non-pinned copy).
fn pin_for_dma(slice: &[u8]) -> Option<PinnedHostRegion> {
    const MIN_PIN_BYTES: usize = 4 * 1024 * 1024;
    if slice.len() < MIN_PIN_BYTES {
        return None;
    }
    let ptr = slice.as_ptr() as *mut std::ffi::c_void;
    let rc = unsafe { cust::sys::cuMemHostRegister_v2(ptr, slice.len(), 0) };
    if rc != cust::sys::cudaError_enum::CUDA_SUCCESS {
        return None;
    }
    Some(PinnedHostRegion { ptr })
}

/// Build the per-entry MPH verification key array (parallel rayon scan).
fn build_entry_keys(ani: &AniIndex) -> Vec<u64> {
    use rayon::prelude::*;
    let n = ani.entries.len();
    let chunk = 16_384;
    let progress = std::sync::atomic::AtomicUsize::new(0);
    let report_every = chunk * 16;
    (0..n)
        .into_par_iter()
        .with_min_len(chunk)
        .map(|i| {
            let entry = &ani.entries[i];
            let ref_str = ani.read_cstring(entry.ref_ofs as usize);
            let alt_str = ani.read_cstring(entry.alt_ofs as usize);
            let key = make_key(entry.chr_id, entry.pos, ref_str.as_ref(), alt_str.as_ref());
            let prev = progress
                .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            if (prev + 1) % report_every == 0 {
                eprintln!(
                    "[gpu] init: entry_keys progress {}/{} ({:.0}%)",
                    prev + 1,
                    n,
                    100.0 * (prev + 1) as f64 / n as f64
                );
            }
            key
        })
        .collect()
}
