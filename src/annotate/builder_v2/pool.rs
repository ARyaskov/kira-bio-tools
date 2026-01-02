use anyhow::Result;
use libdeflater::{CompressionLvl, Compressor, Decompressor};
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

const DEFAULT_BLOCK_SIZE: usize = 1024 * 1024;

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
    compressor: Compressor,
    scratch: Vec<u8>,
    finished: bool,
}

impl StringPool {
    pub fn new() -> Self {
        let lvl = CompressionLvl::new(4).unwrap();
        Self {
            blocks: Vec::new(),
            current: Vec::with_capacity(DEFAULT_BLOCK_SIZE),
            current_start: 0,
            total_len: 0,
            block_size: DEFAULT_BLOCK_SIZE,
            compressor: Compressor::new(lvl),
            scratch: Vec::new(),
            finished: false,
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

    pub fn append_cstr(&mut self, s: &str) -> usize {
        let ofs = self.total_len as usize;
        let needed = s.len() + 1;

        if !self.current.is_empty() && self.current.len() + needed > self.block_size {
            self.flush_block();
        }

        self.current.extend_from_slice(s.as_bytes());
        self.current.push(0);
        self.total_len += needed as u64;

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
        let mut decompressor = Decompressor::new();

        for block in &self.blocks {
            let mut buf = vec![0u8; block.raw_len as usize];
            decompressor
                .deflate_decompress(&block.data, &mut buf)
                .map_err(|_| anyhow::anyhow!("String pool decompression failed"))?;
            out.extend_from_slice(&buf);
        }

        Ok(out)
    }

    pub fn write_to(&mut self, out: &mut File) -> Result<()> {
        let raw = self.materialize()?;
        out.write_all(&raw)?;
        Ok(())
    }

    pub fn cleanup(&mut self) {}

    fn flush_block(&mut self) {
        if self.current.is_empty() {
            return;
        }

        let raw_len = self.current.len();
        let bound = self.compressor.deflate_compress_bound(raw_len);
        if self.scratch.len() < bound {
            self.scratch.resize(bound, 0);
        }

        let comp_len = self
            .compressor
            .deflate_compress(&self.current, &mut self.scratch)
            .expect("String pool compression failed");

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
