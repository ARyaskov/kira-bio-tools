use flate2::write::DeflateEncoder;
use flate2::Compression;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};

use super::BGZF_BLOCK_SIZE;

const BGZF_HEADER: [u8; 18] = [
    0x1f, 0x8b, 0x08, 0x04, // GZIP magic + compression method + flags
    0x00, 0x00, 0x00, 0x00, // MTIME
    0x00, 0xff, // XFL + OS
    0x06, 0x00, // XLEN = 6
    0x42, 0x43, // SI1='B', SI2='C'
    0x02, 0x00, // SLEN = 2
    0x00, 0x00, // BSIZE placeholder
];

const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

pub struct ParallelBgzfWriter {
    file: Arc<Mutex<File>>,
    buffer: Vec<u8>,
    compression_level: Compression,
}

impl ParallelBgzfWriter {
    pub fn create<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::create(path)?;

        Ok(Self {
            file: Arc::new(Mutex::new(file)),
            buffer: Vec::with_capacity(BGZF_BLOCK_SIZE),
            compression_level: Compression::fast(),
        })
    }

    pub fn write_all(&mut self, data: &[u8]) -> io::Result<()> {
        self.buffer.extend_from_slice(data);

        while self.buffer.len() >= BGZF_BLOCK_SIZE {
            self.flush_block()?;
        }

        Ok(())
    }

    fn flush_block(&mut self) -> io::Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }

        let chunk_size = self.buffer.len().min(BGZF_BLOCK_SIZE - 1024);
        let chunk = self.buffer.drain(..chunk_size).collect::<Vec<_>>();

        let compressed = Self::compress_block(&chunk, self.compression_level)?;

        let mut file = self.file.lock().unwrap();
        file.write_all(&compressed)?;

        Ok(())
    }

    fn compress_block(data: &[u8], level: Compression) -> io::Result<Vec<u8>> {
        let mut block = Vec::with_capacity(data.len() + 128);
        block.extend_from_slice(&BGZF_HEADER);

        let cdata_start = block.len();
        let mut encoder = DeflateEncoder::new(Vec::new(), level);
        encoder.write_all(data)?;
        let compressed = encoder.finish()?;

        block.extend_from_slice(&compressed);

        // CRC32
        let crc = crc32fast::hash(data);
        block.extend_from_slice(&crc.to_le_bytes());

        // ISIZE
        block.extend_from_slice(&(data.len() as u32).to_le_bytes());

        // Update BSIZE
        let bsize = (block.len() - 1) as u16;
        block[16..18].copy_from_slice(&bsize.to_le_bytes());

        Ok(block)
    }

    pub fn finish(mut self) -> io::Result<File> {
        // Flush remaining data
        if !self.buffer.is_empty() {
            self.flush_block()?;
        }

        // Write EOF marker
        let mut file = self.file.lock().unwrap();
        file.write_all(&BGZF_EOF)?;
        file.flush()?;

        drop(file);
        Arc::try_unwrap(self.file)
            .ok()
            .and_then(|m| m.into_inner().ok())
            .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "Failed to unwrap file"))
    }
}

impl Write for ParallelBgzfWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.write_all(buf)?;
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        self.flush_block()
    }
}
