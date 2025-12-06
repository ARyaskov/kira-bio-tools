use crossbeam_channel::{bounded, Receiver, Sender};
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::path::Path;
use std::sync::Arc;

use super::block::{BgzfBlock, BlockDecoder};
use super::{BATCH_SIZE, BGZF_BLOCK_SIZE};

pub struct ParallelBgzfReader {
    file: Arc<File>,
    file_size: u64,
    current_offset: u64,
    batch_size: usize,
}

impl ParallelBgzfReader {
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path)?;
        let file_size = file.metadata()?.len();

        Ok(Self {
            file: Arc::new(file),
            file_size,
            current_offset: 0,
            batch_size: BATCH_SIZE,
        })
    }

    pub fn read_batch(&mut self) -> io::Result<Vec<BgzfBlock>> {
        if self.current_offset >= self.file_size {
            return Ok(Vec::new());
        }

        let mut compressed_blocks = Vec::with_capacity(self.batch_size);
        let mut file_clone = (*self.file).try_clone()?;
        file_clone.seek(SeekFrom::Start(self.current_offset))?;

        // Read batch of compressed blocks
        for _ in 0..self.batch_size {
            if self.current_offset >= self.file_size {
                break;
            }

            let mut header = [0u8; 18];
            if file_clone.read_exact(&mut header).is_err() {
                break;
            }

            if BlockDecoder::is_eof_marker(&header) {
                break;
            }

            let bsize = u16::from_le_bytes([header[16], header[17]]) as usize + 1;

            let mut block = BgzfBlock::new(self.current_offset);
            block.compressed.extend_from_slice(&header);
            block.compressed.resize(bsize, 0);

            file_clone.read_exact(&mut block.compressed[18..])?;

            self.current_offset += bsize as u64;
            compressed_blocks.push(block);
        }

        if compressed_blocks.is_empty() {
            return Ok(Vec::new());
        }

        // Parallel decompression
        let decompressed: Vec<_> = compressed_blocks
            .into_par_iter()
            .map(|mut block| {
                BlockDecoder::decode(&mut block).ok();
                block
            })
            .collect();

        Ok(decompressed)
    }

    pub fn virtual_position(&self) -> u64 {
        self.current_offset << 16
    }

    pub fn seek(&mut self, vpos: u64) -> io::Result<()> {
        self.current_offset = vpos >> 16;
        Ok(())
    }
}

pub struct BatchedLineReader {
    reader: ParallelBgzfReader,
    buffer: Vec<u8>,
    buffer_pos: usize,
    batch_size: usize,
    vpos: u64,
}

impl BatchedLineReader {
    pub fn new(reader: ParallelBgzfReader, batch_size: usize) -> Self {
        Self {
            reader,
            buffer: Vec::with_capacity(batch_size * 1024),
            buffer_pos: 0,
            batch_size,
            vpos: 0,
        }
    }

    pub fn read_batch(&mut self) -> io::Result<Vec<(String, u64)>> {
        let blocks = self.reader.read_batch()?;
        if blocks.is_empty() {
            return Ok(Vec::new());
        }

        let mut lines = Vec::with_capacity(self.batch_size);

        for block in blocks {
            let block_vpos = block.virtual_offset();
            let mut line_start = 0;

            for (i, &byte) in block.uncompressed.iter().enumerate() {
                if byte == b'\n' {
                    let line =
                        String::from_utf8_lossy(&block.uncompressed[line_start..i]).to_string();
                    lines.push((line, block_vpos + line_start as u64));
                    line_start = i + 1;
                }
            }

            // Handle partial line at end of block
            if line_start < block.uncompressed.len() {
                self.buffer
                    .extend_from_slice(&block.uncompressed[line_start..]);
            }
        }

        Ok(lines)
    }

    pub fn virtual_position(&self) -> u64 {
        self.reader.virtual_position()
    }
}
