//! `.ktile` reader with adaptive strategy (mmap whole-file, sliding window, or compressed-chunk LRU).

use std::cell::RefCell;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;
use std::thread::JoinHandle;

use anyhow::{Context, Result};
use crossbeam_channel::{Receiver, Sender, bounded};
use flate2::{Decompress, FlushDecompress, Status};
use fxhash::FxHashMap;
use memmap2::{Mmap, MmapOptions};

use super::format::{CompressedChunkEntry, KtileError, KtileHeader};

const KTILE_MMAP_MAX_MB_DEFAULT: u64 = 16 * 1024;
const KTILE_CHUNK_MB_DEFAULT: u64 = 4 * 1024;

fn env_mb_or(default_mb: u64, key: &str) -> u64 {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(default_mb)
}

/// Opened `.ktile` file.
pub struct KtileReader {
    mmap: Arc<Mmap>,
    file: Arc<File>,
    header: KtileHeader,
    strategy: PoolStrategy,
}

enum PoolStrategy {
    Whole,
    Sliding {
        chunk_size: u64,
        cache: RefCell<Option<ChunkCache>>,
    },
    Compressed {
        chunks: Vec<CompressedChunkEntry>,
        lines_per_chunk: u32,
        cache: RefCell<CompressedCache>,
    },
}

struct CompressedCache {
    entries: Vec<CachedChunk>,
    decompressor: Decompress,
    prefetched: Receiver<(u32, Vec<u8>)>,
    request_tx: Option<Sender<u32>>,
    handles: Vec<JoinHandle<()>>,
    next_to_prefetch: u32,
    n_chunks: u32,
    pending: FxHashMap<u32, Vec<u8>>,
}

impl Drop for CompressedCache {
    fn drop(&mut self) {
        self.request_tx.take();
        for h in self.handles.drain(..) {
            let _ = h.join();
        }
    }
}

struct CachedChunk {
    chunk_idx: u32,
    data: Vec<u8>,
}

const CACHE_DEPTH: usize = 2;
const PENDING_MAX: usize = 16;

fn n_decompress_workers() -> usize {
    std::env::var("KIRA_BT_KTILE_DECOMPRESS_THREADS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&n| (1..=32).contains(&n))
        .unwrap_or(4)
}

fn spawn_decompress_pool(
    n_workers: usize,
    file: Arc<File>,
    line_pool_off: u64,
    chunks: Arc<Vec<CompressedChunkEntry>>,
    requests: Receiver<u32>,
    out: Sender<(u32, Vec<u8>)>,
) -> Vec<JoinHandle<()>> {
    (0..n_workers)
        .map(|wid| {
            let requests = requests.clone();
            let out = out.clone();
            let file = Arc::clone(&file);
            let chunks = Arc::clone(&chunks);
            std::thread::Builder::new()
                .name(format!("ktile-decompress-{wid}"))
                .spawn(move || {
                    let mut decompressor = Decompress::new(false);
                    while let Ok(chunk_idx) = requests.recv() {
                        let Some(chunk) = chunks.get(chunk_idx as usize) else {
                            continue;
                        };
                        let mut compressed = vec![0u8; chunk.compressed_size as usize];
                        let abs_off = line_pool_off + chunk.compressed_off;
                        if pread_exact(&file, &mut compressed, abs_off).is_err() {
                            continue;
                        }
                        let mut decompressed = vec![0u8; chunk.uncompressed_size as usize];
                        decompressor.reset(false);
                        let before_in = decompressor.total_in();
                        let before_out = decompressor.total_out();
                        let status = match decompressor.decompress(
                            &compressed,
                            &mut decompressed,
                            FlushDecompress::Finish,
                        ) {
                            Ok(s) => s,
                            Err(_) => continue,
                        };
                        if status != Status::StreamEnd
                            || (decompressor.total_in() - before_in) as u32 != chunk.compressed_size
                            || (decompressor.total_out() - before_out) as u32
                                != chunk.uncompressed_size
                        {
                            continue;
                        }
                        if out.send((chunk_idx, decompressed)).is_err() {
                            break;
                        }
                    }
                })
                .expect("spawn ktile-decompress thread")
        })
        .collect()
}

struct ChunkCache {
    data: Vec<u8>,
    pool_start: u64,
}

impl ChunkCache {
    fn covers(&self, pool_start: u64, pool_end: u64) -> bool {
        pool_start >= self.pool_start
            && pool_end <= self.pool_start + self.data.len() as u64
    }

    fn slice(&self, pool_start: u64, pool_end: u64) -> &[u8] {
        let local_start = (pool_start - self.pool_start) as usize;
        let local_end = (pool_end - self.pool_start) as usize;
        &self.data[local_start..local_end]
    }
}

#[cfg(unix)]
fn pread_exact(file: &File, buf: &mut [u8], offset: u64) -> std::io::Result<()> {
    use std::os::unix::fs::FileExt;
    file.read_exact_at(buf, offset)
}

#[cfg(windows)]
fn pread_exact(file: &File, buf: &mut [u8], offset: u64) -> std::io::Result<()> {
    use std::os::windows::fs::FileExt;
    let mut filled = 0usize;
    while filled < buf.len() {
        let n = file.seek_read(&mut buf[filled..], offset + filled as u64)?;
        if n == 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                "ktile chunk pread: unexpected EOF",
            ));
        }
        filled += n;
    }
    Ok(())
}

impl KtileReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path.as_ref())
            .with_context(|| format!("opening ktile {:?}", path.as_ref()))?;
        let file_size = file.metadata().context("stat ktile")?.len();
        // SAFETY: read-only mmap; concurrent writes unsupported.
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        #[cfg(unix)]
        let _ = mmap.advise(memmap2::Advice::Sequential);
        let mmap = Arc::new(mmap);
        let file = Arc::new(file);

        if mmap.len() < KtileHeader::SIZE {
            return Err(KtileError::Truncated { section: "header" }.into());
        }
        let header: KtileHeader = *bytemuck::from_bytes(&mmap[..KtileHeader::SIZE]);
        header.validate()?;

        let len = mmap.len() as u64;
        let check = |off: u64, byte_len: u64, name: &'static str| -> Result<()> {
            if off + byte_len > len {
                Err(KtileError::Truncated { section: name }.into())
            } else {
                Ok(())
            }
        };
        check(header.headers_off, header.headers_len, "headers")?;
        let n = header.n_records as usize;
        check(
            header.line_offsets_off,
            ((n + 1) * std::mem::size_of::<u64>()) as u64,
            "line_offsets",
        )?;
        check(
            header.chr_ids_off,
            (n * std::mem::size_of::<u32>()) as u64,
            "chr_ids",
        )?;
        check(
            header.positions_off,
            (n * std::mem::size_of::<u32>()) as u64,
            "positions",
        )?;
        check(header.line_pool_off, header.line_pool_len, "line_pool")?;
        // Phase 3 columns are optional — only validate when the flag is
        // set (writer pre-Phase 3 leaves the four fields zeroed).
        if header.has_ref_alt_columns() {
            let col_bytes = (n * std::mem::size_of::<u32>()) as u64;
            check(header.off_ref_offsets, col_bytes, "ref_offsets")?;
            check(header.off_ref_lens, col_bytes, "ref_lens")?;
            check(header.off_alt_offsets, col_bytes, "alt_offsets")?;
            check(header.off_alt_lens, col_bytes, "alt_lens")?;
        }

        // Choose strategy. Compressed mode is independent of the size
        // threshold — once a ktile was built compressed, we always need
        // the LRU-cached decompression path. Uncompressed ktiles go
        // through the file-size threshold to pick Whole vs Sliding.
        let mmap_max_bytes = env_mb_or(KTILE_MMAP_MAX_MB_DEFAULT, "KIRA_BT_KTILE_MMAP_MAX_MB")
            .saturating_mul(1024 * 1024);
        let strategy = if header.has_compressed_pool() {
            // Load the chunk index (small — 24 B per chunk).
            let n_chunks = header.n_chunks as usize;
            if n_chunks > 0 {
                let chunk_index_off = header.off_chunk_index as usize;
                let chunk_index_bytes = n_chunks * std::mem::size_of::<CompressedChunkEntry>();
                check(
                    header.off_chunk_index,
                    chunk_index_bytes as u64,
                    "chunk_index",
                )?;
                let chunks_vec: Vec<CompressedChunkEntry> = bytemuck::cast_slice(
                    &mmap[chunk_index_off..chunk_index_off + chunk_index_bytes],
                )
                .to_vec();
                // Spawn N decompress workers sharing a single request
                // channel. Channel capacity scales with worker count
                // so all workers can be busy with N requests in flight
                // plus a small slack for queueing.
                let n_workers = n_decompress_workers();
                let prefetch_cap = n_workers + 2;
                eprintln!(
                    "[ktile] compressed pool: {} chunks, {} lines/chunk, {} MB on disk → ~{} MB uncompressed, {} decompress workers",
                    n_chunks,
                    header.lines_per_chunk,
                    header.line_pool_len / (1024 * 1024),
                    chunks_vec.iter().map(|c| c.uncompressed_size as u64).sum::<u64>() / (1024 * 1024),
                    n_workers,
                );
                let chunks_arc = Arc::new(chunks_vec.clone());
                let (req_tx, req_rx) = bounded::<u32>(prefetch_cap);
                let (chunk_tx, chunk_rx) = bounded::<(u32, Vec<u8>)>(prefetch_cap);
                let handles = spawn_decompress_pool(
                    n_workers,
                    Arc::clone(&file),
                    header.line_pool_off,
                    Arc::clone(&chunks_arc),
                    req_rx,
                    chunk_tx,
                );
                // Pre-queue up to `prefetch_cap` chunks so every worker
                // has work immediately on startup. With N workers, the
                // first N requests start in parallel and the reader
                // sees most of them ready by the time it asks.
                let initial = (prefetch_cap as u32).min(n_chunks as u32);
                for i in 0..initial {
                    let _ = req_tx.send(i);
                }
                PoolStrategy::Compressed {
                    chunks: chunks_vec,
                    lines_per_chunk: header.lines_per_chunk,
                    cache: RefCell::new(CompressedCache {
                        entries: Vec::with_capacity(CACHE_DEPTH),
                        decompressor: Decompress::new(false),
                        prefetched: chunk_rx,
                        request_tx: Some(req_tx),
                        handles,
                        next_to_prefetch: initial,
                        n_chunks: n_chunks as u32,
                        pending: FxHashMap::default(),
                    }),
                }
            } else {
                PoolStrategy::Whole // empty ktile
            }
        } else if file_size <= mmap_max_bytes {
            PoolStrategy::Whole
        } else {
            let chunk_size = env_mb_or(KTILE_CHUNK_MB_DEFAULT, "KIRA_BT_KTILE_CHUNK_MB")
                .saturating_mul(1024 * 1024);
            eprintln!(
                "[ktile] file {} MB > {} MB threshold → sliding-window mode ({} MB chunks)",
                file_size / (1024 * 1024),
                mmap_max_bytes / (1024 * 1024),
                chunk_size / (1024 * 1024)
            );
            PoolStrategy::Sliding {
                chunk_size,
                cache: RefCell::new(None),
            }
        };

        Ok(Self {
            mmap,
            file,
            header,
            strategy,
        })
    }

    #[inline]
    pub fn n_records(&self) -> usize {
        self.header.n_records as usize
    }

    /// Read-only view of the parsed header. Used by freshness checks +
    /// debug introspection.
    #[inline]
    pub fn header(&self) -> &KtileHeader {
        &self.header
    }

    /// Raw VCF header block as it was captured during build (UTF-8,
    /// '\n'-separated, no trailing newline).
    pub fn headers_block(&self) -> &str {
        let start = self.header.headers_off as usize;
        let end = start + self.header.headers_len as usize;
        // SAFETY: bytes were written from a valid UTF-8 String during build.
        unsafe { std::str::from_utf8_unchecked(&self.mmap[start..end]) }
    }

    fn line_offsets(&self) -> &[u64] {
        let n = self.header.n_records as usize + 1;
        let start = self.header.line_offsets_off as usize;
        let end = start + n * std::mem::size_of::<u64>();
        bytemuck::cast_slice(&self.mmap[start..end])
    }

    fn chr_ids(&self) -> &[u32] {
        let n = self.header.n_records as usize;
        let start = self.header.chr_ids_off as usize;
        let end = start + n * std::mem::size_of::<u32>();
        bytemuck::cast_slice(&self.mmap[start..end])
    }

    fn positions(&self) -> &[u32] {
        let n = self.header.n_records as usize;
        let start = self.header.positions_off as usize;
        let end = start + n * std::mem::size_of::<u32>();
        bytemuck::cast_slice(&self.mmap[start..end])
    }

    /// Returns the i-th line as a String. Always allocates (uniformly
    /// across Whole and Sliding strategies) so the API is the same and
    /// the caller — `KtileSourceReader::read_line` — gets the owned
    /// String it needs to send through the bounded channel.
    pub fn line_owned(&self, idx: usize) -> String {
        let bytes = self.line_bytes_to_vec(idx);
        // SAFETY: ktile builder wrote bytes from a valid UTF-8 source line.
        unsafe { String::from_utf8_unchecked(bytes) }
    }

    /// Zero-copy fast path: copy line `idx` bytes directly into
    /// `batch` without going through a `String` intermediate. Saves
    /// one allocation + one memcpy per line vs `line_owned`.
    ///
    /// Used by the GPU/CPU annotate readers to feed ReadBatch
    /// without paying the `Vec<u8> → String → batch.bytes` chain. On
    /// a 6.5 M-record chr1 1000G this saves ~3-5 sec wall-clock.
    pub fn push_line_into_batch(
        &self,
        idx: usize,
        batch: &mut crate::annotate::cpu_v2::ReadBatch,
    ) {
        let offs = self.line_offsets();
        let pool_start = offs[idx];
        let pool_end = offs[idx + 1];
        match &self.strategy {
            PoolStrategy::Whole => {
                let abs_start = (self.header.line_pool_off + pool_start) as usize;
                let abs_end = (self.header.line_pool_off + pool_end) as usize;
                batch.push_line_bytes(&self.mmap[abs_start..abs_end]);
            }
            PoolStrategy::Sliding { chunk_size, cache } => {
                self.with_chunk(*chunk_size, cache, pool_start, pool_end, |bytes| {
                    batch.push_line_bytes(bytes);
                });
            }
            PoolStrategy::Compressed {
                chunks,
                lines_per_chunk,
                cache,
            } => {
                self.with_compressed_chunk(
                    idx,
                    chunks,
                    *lines_per_chunk,
                    cache,
                    |bytes| {
                        batch.push_line_bytes(bytes);
                    },
                );
            }
        }
    }

    /// Returns line `idx` bytes as a fresh Vec. Used internally by
    /// `line_owned`, `ref_slice`, `alt_slice`.
    fn line_bytes_to_vec(&self, idx: usize) -> Vec<u8> {
        let offs = self.line_offsets();
        let pool_start = offs[idx];
        let pool_end = offs[idx + 1];
        match &self.strategy {
            PoolStrategy::Whole => {
                let abs_start = (self.header.line_pool_off + pool_start) as usize;
                let abs_end = (self.header.line_pool_off + pool_end) as usize;
                self.mmap[abs_start..abs_end].to_vec()
            }
            PoolStrategy::Sliding { chunk_size, cache } => {
                self.with_chunk(*chunk_size, cache, pool_start, pool_end, |bytes| {
                    bytes.to_vec()
                })
            }
            PoolStrategy::Compressed {
                chunks,
                lines_per_chunk,
                cache,
            } => self.with_compressed_chunk(idx, chunks, *lines_per_chunk, cache, |bytes| {
                bytes.to_vec()
            }),
        }
    }

    /// Decompresses the chunk containing line `idx` (if not already
    /// cached) and runs `f` on the line's byte slice within the
    /// decompressed chunk buffer.
    ///
    /// Cache resolution order:
    ///   1. Already in LRU → hit (warm sequential / re-reads).
    ///   2. Already in the reorder buffer `pending` (workers got
    ///      ahead) → take it, promote to LRU.
    ///   3. Drain the worker output channel, store all in `pending`;
    ///      if our target chunk shows up, take it.
    ///   4. Fallback: synchronous decompress on the reader thread
    ///      (only when *no* worker has produced our chunk yet — rare
    ///      in steady state with N workers running ahead).
    ///
    /// Sequential access (the annotate read loop) hits path 2 or 3
    /// almost always, so the reader never waits on decompress — it
    /// runs in parallel with the worker pool.
    fn with_compressed_chunk<R>(
        &self,
        line_idx: usize,
        chunks: &[CompressedChunkEntry],
        lines_per_chunk: u32,
        cache: &RefCell<CompressedCache>,
        f: impl FnOnce(&[u8]) -> R,
    ) -> R {
        let chunk_idx = (line_idx / lines_per_chunk as usize) as u32;
        let chunk = chunks
            .get(chunk_idx as usize)
            .expect("line_idx must map to a valid chunk");

        let mut cache_borrow = cache.borrow_mut();

        // Always drain the worker output channel first. This:
        //   * Unblocks any worker stuck on `out.send()` (channel full),
        //     so the pool stays saturated.
        //   * Maximises the chance of finding `chunk_idx` in `pending`
        //     instead of synchronously decompressing.
        // Out-of-order arrivals go into `pending`; oldest gets evicted
        // when we hit PENDING_MAX (cap is generous, so this only
        // triggers under pathological out-of-order read patterns).
        while let Ok((idx, data)) = cache_borrow.prefetched.try_recv() {
            if cache_borrow.pending.len() >= PENDING_MAX
                && !cache_borrow.pending.contains_key(&idx)
            {
                if let Some(&smallest) = cache_borrow.pending.keys().min() {
                    cache_borrow.pending.remove(&smallest);
                }
            }
            cache_borrow.pending.insert(idx, data);
        }

        let mut entry_pos = cache_borrow
            .entries
            .iter()
            .position(|c| c.chunk_idx == chunk_idx);

        // Path 2: target chunk waiting in the reorder buffer.
        if entry_pos.is_none() {
            if let Some(data) = cache_borrow.pending.remove(&chunk_idx) {
                if cache_borrow.entries.len() >= CACHE_DEPTH {
                    cache_borrow.entries.remove(0);
                }
                cache_borrow.entries.push(CachedChunk {
                    chunk_idx,
                    data,
                });
                entry_pos = Some(cache_borrow.entries.len() - 1);
            }
        }

        // Path 4: synchronous fallback decompress.
        if entry_pos.is_none() {
            let decompressed = self.decompress_chunk(chunk, &mut cache_borrow.decompressor);
            if cache_borrow.entries.len() >= CACHE_DEPTH {
                cache_borrow.entries.remove(0);
            }
            cache_borrow.entries.push(CachedChunk {
                chunk_idx,
                data: decompressed,
            });
            entry_pos = Some(cache_borrow.entries.len() - 1);
            // Re-prime the worker pool to start from where we just
            // landed (in case workers fell behind / out-of-sync).
            cache_borrow.next_to_prefetch = chunk_idx + 1;
        }

        // Queue up the next batch of prefetch requests, capped at
        // `chunk_idx + MAX_AHEAD` so workers never race more than
        // ~N_workers chunks past the reader (prevents pending from
        // overflowing PENDING_MAX and evicting chunks the reader
        // will need next).
        Self::request_next_prefetch_compressed(&mut cache_borrow, chunk_idx);

        let entry_pos = entry_pos.expect("entry populated above");

        // Compute line-local offsets within the chunk's uncompressed
        // span. line_offsets[i] is uncompressed-pool-relative; the chunk
        // covers [chunk.uncompressed_off .. + uncompressed_size].
        let offs = self.line_offsets();
        let pool_start = offs[line_idx];
        let pool_end = offs[line_idx + 1];
        let local_start = (pool_start - chunk.uncompressed_off) as usize;
        let local_end = (pool_end - chunk.uncompressed_off) as usize;
        f(&cache_borrow.entries[entry_pos].data[local_start..local_end])
    }

    /// Pushes prefetch requests to the worker pool, capped at
    /// `reader_chunk_idx + max_ahead`. Workers never race more than
    /// ~max_ahead chunks past the reader cursor, so the reorder
    /// buffer `pending` stays well below PENDING_MAX even if workers
    /// finish far faster than the reader.
    ///
    /// `max_ahead` = `n_workers + 2`: enough to keep every worker
    /// busy with a different chunk plus a small slack for jitter,
    /// but not so deep that decompressed bytes pile up uselessly.
    fn request_next_prefetch_compressed(cache: &mut CompressedCache, reader_chunk_idx: u32) {
        let Some(tx) = cache.request_tx.as_ref() else {
            return;
        };
        let max_ahead = cache.handles.len() as u32 + 2;
        let cap = reader_chunk_idx.saturating_add(max_ahead + 1);
        while cache.next_to_prefetch < cache.n_chunks && cache.next_to_prefetch < cap {
            match tx.try_send(cache.next_to_prefetch) {
                Ok(()) => cache.next_to_prefetch += 1,
                Err(_) => break, // channel full or closed
            }
        }
    }

    /// Reads a chunk's compressed bytes via pread, decompresses with a
    /// reused decompressor. Returns the uncompressed bytes.
    fn decompress_chunk(
        &self,
        chunk: &CompressedChunkEntry,
        decompressor: &mut Decompress,
    ) -> Vec<u8> {
        let mut compressed = vec![0u8; chunk.compressed_size as usize];
        let abs_off = self.header.line_pool_off + chunk.compressed_off;
        pread_exact(&self.file, &mut compressed, abs_off).unwrap_or_else(|e| {
            panic!(
                "ktile compressed chunk pread failed at file offset {abs_off} length {}: {e}",
                chunk.compressed_size
            )
        });
        let mut decompressed = vec![0u8; chunk.uncompressed_size as usize];
        decompressor.reset(false);
        let before_in = decompressor.total_in();
        let before_out = decompressor.total_out();
        let status = decompressor
            .decompress(&compressed, &mut decompressed, FlushDecompress::Finish)
            .unwrap_or_else(|e| panic!("ktile chunk decompress: {e}"));
        if status != Status::StreamEnd
            || (decompressor.total_in() - before_in) as u32 != chunk.compressed_size
            || (decompressor.total_out() - before_out) as u32 != chunk.uncompressed_size
        {
            panic!(
                "ktile chunk decompress: unexpected status {:?} or size mismatch (in {} / {}, out {} / {})",
                status,
                decompressor.total_in() - before_in,
                chunk.compressed_size,
                decompressor.total_out() - before_out,
                chunk.uncompressed_size,
            );
        }
        decompressed
    }

    /// Ensures the sliding-window cache covers `[pool_start..pool_end]`,
    /// then runs `f` on the byte slice.
    ///
    /// Refill policy: if the current chunk doesn't cover the request,
    /// drop it and read a fresh `chunk_size`-byte chunk anchored at
    /// `pool_start`. Sequential access (the annotate loop) refills the
    /// chunk every ~chunk_size bytes; one extra refill is paid per line
    /// that straddles a chunk boundary (rare for VCF where most lines
    /// are 10s of KB).
    fn with_chunk<R>(
        &self,
        chunk_size: u64,
        cache: &RefCell<Option<ChunkCache>>,
        pool_start: u64,
        pool_end: u64,
        f: impl FnOnce(&[u8]) -> R,
    ) -> R {
        let line_len = pool_end - pool_start;
        let target_chunk = chunk_size.max(line_len);
        let pool_len = self.header.line_pool_len;

        let mut cache_borrow = cache.borrow_mut();
        let needs_refill = cache_borrow
            .as_ref()
            .map_or(true, |c| !c.covers(pool_start, pool_end));

        if needs_refill {
            let want = target_chunk.min(pool_len - pool_start) as usize;
            let mut data = vec![0u8; want];
            let abs_off = self.header.line_pool_off + pool_start;
            pread_exact(&self.file, &mut data, abs_off).unwrap_or_else(|e| {
                panic!("ktile pread failed at file offset {abs_off} length {want}: {e}")
            });
            *cache_borrow = Some(ChunkCache { data, pool_start });
        }

        let cache_ref = cache_borrow.as_ref().expect("chunk just refilled");
        f(cache_ref.slice(pool_start, pool_end))
    }

    /// Pre-parsed chromosome id for record `idx`.
    pub fn chr_id(&self, idx: usize) -> u32 {
        self.chr_ids()[idx]
    }

    /// Pre-parsed 1-based POS for record `idx`.
    pub fn position(&self, idx: usize) -> u32 {
        self.positions()[idx]
    }

    /// True iff the ktile was built with Phase 3 (per-record REF/ALT
    /// offset/length columns). Falls back to whole-line slicing via
    /// `line_owned()` when false; callers that depend on direct REF/ALT
    /// slices should check this first.
    pub fn has_ref_alt_columns(&self) -> bool {
        self.header.has_ref_alt_columns()
    }

    /// REF allele bytes for record `idx` (always owned — the underlying
    /// strategy may be sliding, in which case we can't borrow long-lived
    /// slices).
    pub fn ref_slice(&self, idx: usize) -> Option<Vec<u8>> {
        if !self.has_ref_alt_columns() {
            return None;
        }
        let line = self.line_bytes_to_vec(idx);
        let off = self.ref_offsets()[idx] as usize;
        let len = self.ref_lens()[idx] as usize;
        if off + len > line.len() {
            return None;
        }
        Some(line[off..off + len].to_vec())
    }

    /// ALT allele span bytes for record `idx`. Multi-allelic records
    /// return a comma-separated span — caller splits on `,`.
    pub fn alt_slice(&self, idx: usize) -> Option<Vec<u8>> {
        if !self.has_ref_alt_columns() {
            return None;
        }
        let line = self.line_bytes_to_vec(idx);
        let off = self.alt_offsets()[idx] as usize;
        let len = self.alt_lens()[idx] as usize;
        if off + len > line.len() {
            return None;
        }
        Some(line[off..off + len].to_vec())
    }

    fn ref_offsets(&self) -> &[u32] {
        let n = self.header.n_records as usize;
        let start = self.header.off_ref_offsets as usize;
        let end = start + n * std::mem::size_of::<u32>();
        bytemuck::cast_slice(&self.mmap[start..end])
    }
    fn ref_lens(&self) -> &[u32] {
        let n = self.header.n_records as usize;
        let start = self.header.off_ref_lens as usize;
        let end = start + n * std::mem::size_of::<u32>();
        bytemuck::cast_slice(&self.mmap[start..end])
    }
    fn alt_offsets(&self) -> &[u32] {
        let n = self.header.n_records as usize;
        let start = self.header.off_alt_offsets as usize;
        let end = start + n * std::mem::size_of::<u32>();
        bytemuck::cast_slice(&self.mmap[start..end])
    }
    fn alt_lens(&self) -> &[u32] {
        let n = self.header.n_records as usize;
        let start = self.header.off_alt_lens as usize;
        let end = start + n * std::mem::size_of::<u32>();
        bytemuck::cast_slice(&self.mmap[start..end])
    }

    /// True iff the reader is in sliding-window mode (huge uncompressed
    /// file). For tests / introspection / log messages.
    pub fn is_sliding(&self) -> bool {
        matches!(self.strategy, PoolStrategy::Sliding { .. })
    }

    /// True iff the reader is in compressed-pool mode (per-chunk
    /// deflate-on-disk). Decompression happens on demand into a small
    /// LRU cache.
    pub fn is_compressed(&self) -> bool {
        matches!(self.strategy, PoolStrategy::Compressed { .. })
    }
}

#[cfg(test)]
#[path = "../../../tests/unit/annotate_ktile_reader.rs"]
mod tests;
