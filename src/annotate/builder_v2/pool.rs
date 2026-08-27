//! Append-only compressed string pool for `.ani` building. Combines
//! short-string interning with deflate-compressed blocks.

use anyhow::Result;
use flate2::{Compress, Compression, Decompress, FlushCompress, FlushDecompress, Status};
use fxhash::FxHashMap;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

const DEFAULT_BLOCK_SIZE: usize = 1024 * 1024;

/// Intern window. Strings within `[INTERN_MIN_LEN..=INTERN_MAX_LEN]` bytes
/// participate in the dedup map; outside that range they bypass interning.
/// Set `KIRA_BT_DISABLE_INTERN=1` to disable interning entirely.
const INTERN_MIN_LEN: usize = 8;
const INTERN_MAX_LEN: usize = 64;

/// Deflate compression level for blob blocks.
const DEFLATE_LEVEL: u32 = 4;

/// Worst-case raw-deflate output bound (zlib's deflateBound formula).
#[inline]
fn deflate_compress_bound(input: usize) -> usize {
    input + (input >> 12) + (input >> 14) + (input >> 25) + 13
}

pub struct CompressedBlock {
    pub raw_start: u64,
    pub raw_len: u32,
    pub data: Vec<u8>,
}

pub struct StringPool {
    blocks: Vec<CompressedBlock>,
    current: Vec<u8>,
    current_start: u64,
    total_len: u64,
    block_size: usize,
    compressor: Compress,
    scratch: Vec<u8>,
    finished: bool,
    /// Content-keyed dedup for short strings. `None` when interning is disabled.
    intern: Option<FxHashMap<Box<[u8]>, u32>>,
    intern_hits: u64,
    intern_bytes_saved: u64,
}

impl StringPool {
    pub fn new() -> Self {
        let intern_enabled = std::env::var("KIRA_BT_DISABLE_INTERN")
            .ok()
            .is_none_or(|v| !matches!(v.as_str(), "1" | "true" | "yes" | "y"));
        Self {
            blocks: Vec::new(),
            current: Vec::with_capacity(DEFAULT_BLOCK_SIZE),
            current_start: 0,
            total_len: 0,
            block_size: DEFAULT_BLOCK_SIZE,
            compressor: Compress::new(Compression::new(DEFLATE_LEVEL), false),
            scratch: Vec::new(),
            finished: false,
            intern: intern_enabled.then(|| {
                FxHashMap::with_capacity_and_hasher(8192, Default::default())
            }),
            intern_hits: 0,
            intern_bytes_saved: 0,
        }
    }

    pub fn with_limit(_limit: Option<usize>, _spill_path: Option<PathBuf>) -> Self {
        Self::new()
    }

    pub fn len(&self) -> usize {
        self.total_len as usize
    }

    pub fn block_size(&self) -> usize {
        self.block_size
    }

    pub fn is_in_memory(&self) -> bool {
        true
    }

    pub fn intern_hits(&self) -> u64 {
        self.intern_hits
    }

    pub fn intern_bytes_saved(&self) -> u64 {
        self.intern_bytes_saved
    }

    pub fn append_cstr(&mut self, s: &str) -> usize {
        let bytes = s.as_bytes();
        let needed = bytes.len() + 1;
        let intern_eligible = (INTERN_MIN_LEN..=INTERN_MAX_LEN).contains(&bytes.len());

        if intern_eligible
            && let Some(intern) = self.intern.as_ref()
            && let Some(&existing_ofs) = intern.get(bytes)
        {
            self.intern_hits += 1;
            self.intern_bytes_saved += needed as u64;
            return existing_ofs as usize;
        }

        let ofs = self.total_len as usize;

        if !self.current.is_empty() && self.current.len() + needed > self.block_size {
            self.flush_block();
        }

        self.current.extend_from_slice(bytes);
        self.current.push(0);
        self.total_len += needed as u64;

        if intern_eligible
            && let Some(intern) = self.intern.as_mut()
            && let Ok(ofs_u32) = u32::try_from(ofs)
        {
            intern.insert(bytes.to_vec().into_boxed_slice(), ofs_u32);
        }

        if self.current.len() >= self.block_size {
            self.flush_block();
        }

        ofs
    }

    pub fn spilled(&self) -> bool {
        false
    }

    pub fn finish(&mut self) {
        if !self.finished {
            self.flush_block();
            self.finished = true;
        }
    }

    pub fn blocks(&mut self) -> &[CompressedBlock] {
        self.finish();
        &self.blocks
    }

    pub fn materialize(&mut self) -> Result<Vec<u8>> {
        self.finish();
        let mut out = Vec::with_capacity(self.total_len as usize);
        let mut decompressor = Decompress::new(false);

        for block in &self.blocks {
            let mut buf = vec![0u8; block.raw_len as usize];
            decompressor.reset(false);
            let before_in = decompressor.total_in();
            let before_out = decompressor.total_out();
            let status = decompressor
                .decompress(&block.data, &mut buf, FlushDecompress::Finish)
                .map_err(|_| anyhow::anyhow!("String pool decompression failed"))?;
            let consumed = (decompressor.total_in() - before_in) as usize;
            let produced = (decompressor.total_out() - before_out) as usize;
            if status != Status::StreamEnd
                || consumed != block.data.len()
                || produced != block.raw_len as usize
            {
                anyhow::bail!(
                    "String pool decompression incomplete: status={:?} consumed={}/{} produced={}/{}",
                    status,
                    consumed,
                    block.data.len(),
                    produced,
                    block.raw_len
                );
            }
            out.extend_from_slice(&buf);
        }

        Ok(out)
    }

    pub fn write_to(&mut self, out: &mut File) -> Result<()> {
        let raw = self.materialize()?;
        out.write_all(&raw)?;
        Ok(())
    }

    /// Drop the intern map once the pool is finalised.
    pub fn cleanup(&mut self) {
        self.intern = None;
    }

    fn flush_block(&mut self) {
        if self.current.is_empty() {
            return;
        }

        let raw_len = self.current.len();
        let bound = deflate_compress_bound(raw_len);
        if self.scratch.len() < bound {
            self.scratch.resize(bound, 0);
        }

        self.compressor.reset();
        let before_in = self.compressor.total_in();
        let before_out = self.compressor.total_out();
        let status = self
            .compressor
            .compress(&self.current, &mut self.scratch, FlushCompress::Finish)
            .expect("String pool compression failed");
        if status != Status::StreamEnd
            || (self.compressor.total_in() - before_in) as usize != raw_len
        {
            panic!(
                "String pool compression did not finish in one call (scratch too small?): status={:?}",
                status
            );
        }
        let comp_len = (self.compressor.total_out() - before_out) as usize;

        let mut data = vec![0u8; comp_len];
        data.copy_from_slice(&self.scratch[..comp_len]);

        self.blocks.push(CompressedBlock {
            raw_start: self.current_start,
            raw_len: raw_len as u32,
            data,
        });

        self.current_start += raw_len as u64;
        self.current.clear();
    }
}

#[cfg(test)]
#[path = "../../../tests/unit/annotate_builder_v2_pool.rs"]
mod tests;
