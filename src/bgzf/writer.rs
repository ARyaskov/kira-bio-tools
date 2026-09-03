use crossbeam_channel::{Receiver, Sender, bounded};
use flate2::Compression;
use libdeflater::{CompressionLvl, Compressor};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::thread;

use super::simd::compute_crc32;
use super::structs::{BGZF_EOF, BGZF_HEADER, CHUNK_SIZE, CompressedBlock};

/// Bounded queue between the producer and the compressor pool: 512 blocks
/// (~32 MB of pre-compression data) keeps the pool fed without burning RAM.
const CHANNEL_DEPTH: usize = 512;
/// Output buffer for file targets. Larger values stall the writer thread on
/// periodic flushes and back-pressure the compressor pool.
pub const FILE_BUFFER_SIZE: usize = 128 * 1024 * 1024;
/// Output buffer for pipes: small so downstream consumers see data early.
pub const STREAM_BUFFER_SIZE: usize = 4 * 1024 * 1024;

/// Multithreaded BGZF writer (libdeflate). The output is finalized (remaining
/// blocks compressed, EOF marker written, buffers flushed) by
/// [`BgzfWriter::finish`]; dropping the writer without calling it finalizes as
/// well so an early return never leaves a truncated file, but errors are only
/// reported by `finish`. Worker count comes from [`crate::threads`].
pub struct BgzfWriter {
    tx: Option<Sender<(Vec<u8>, usize)>>,
    writer_thread: Option<thread::JoinHandle<io::Result<()>>>,
    sequence: usize,
    pending: Vec<u8>,
    finished: bool,
}

impl BgzfWriter {
    pub fn create<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::with_compression(path, Compression::new(1))
    }

    pub fn with_compression<P: AsRef<Path>>(path: P, compression: Compression) -> io::Result<Self> {
        let file = File::create(path)?;
        Self::from_writer_buffered(file, compression, FILE_BUFFER_SIZE)
    }

    /// BGZF stream to standard output.
    pub fn stdout(compression: Compression) -> io::Result<Self> {
        Self::from_writer_buffered(io::stdout(), compression, STREAM_BUFFER_SIZE)
    }

    pub fn from_writer<W: Write + Send + 'static>(writer: W, compression: Compression) -> io::Result<Self> {
        Self::from_writer_buffered(writer, compression, FILE_BUFFER_SIZE)
    }

    pub fn from_writer_buffered<W: Write + Send + 'static>(
        writer: W,
        compression: Compression,
        buffer_size: usize,
    ) -> io::Result<Self> {
        let num_workers = crate::threads::compress_workers();
        // zlib levels 0..=9 map one-to-one onto libdeflate levels (0 = stored).
        let level = compression.level().min(9) as i32;
        let lvl = CompressionLvl::new(level)
            .map_err(|e| io::Error::other(format!("BGZF compression level {level}: {e:?}")))?;
        let writer: Box<dyn Write + Send> = Box::new(writer);
        let writer = BufWriter::with_capacity(buffer_size.max(64 * 1024), writer);
        let (chunk_tx, chunk_rx) = bounded::<(Vec<u8>, usize)>(CHANNEL_DEPTH);
        let (block_tx, block_rx) = bounded::<CompressedBlock>(CHANNEL_DEPTH);

        let compression_workers: Vec<_> = (0..num_workers)
            .map(|_| {
                let rx = chunk_rx.clone();
                let tx = block_tx.clone();
                thread::spawn(move || {
                    let mut compressor = Compressor::new(lvl);
                    let mut scratch = Vec::new();
                    Self::compression_worker(rx, tx, &mut compressor, &mut scratch)
                })
            })
            .collect();

        drop(chunk_rx);
        drop(block_tx);

        let writer_thread =
            thread::spawn(move || Self::writer_worker(block_rx, writer, compression_workers));

        Ok(Self {
            tx: Some(chunk_tx),
            writer_thread: Some(writer_thread),
            sequence: 0,
            pending: Vec::with_capacity(CHUNK_SIZE),
            finished: false,
        })
    }

    fn compression_worker(
        chunk_rx: Receiver<(Vec<u8>, usize)>,
        block_tx: Sender<CompressedBlock>,
        compressor: &mut Compressor,
        scratch: &mut Vec<u8>,
    ) -> io::Result<()> {
        while let Ok((data, seq)) = chunk_rx.recv() {
            let compressed = Self::compress_block(&data, compressor, scratch)?;
            if block_tx.send(CompressedBlock { data: compressed, sequence: seq }).is_err() {
                break;
            }
        }
        Ok(())
    }

    fn writer_worker(
        block_rx: Receiver<CompressedBlock>,
        mut writer: BufWriter<Box<dyn Write + Send>>,
        workers: Vec<thread::JoinHandle<io::Result<()>>>,
    ) -> io::Result<()> {
        let mut pending = std::collections::BTreeMap::new();
        let mut next = 0usize;
        let mut mega_buffer = Vec::<u8>::with_capacity(8 * 1024 * 1024);

        for block in block_rx.iter() {
            pending.insert(block.sequence, block.data);
            while let Some(data) = pending.remove(&next) {
                mega_buffer.extend_from_slice(&data);
                if mega_buffer.len() > 4 * 1024 * 1024 {
                    writer.write_all(&mega_buffer)?;
                    mega_buffer.clear();
                }
                next += 1;
            }
        }

        if !mega_buffer.is_empty() {
            writer.write_all(&mega_buffer)?;
        }

        for handle in workers {
            handle
                .join()
                .map_err(|_| io::Error::other("BGZF compression worker panicked"))??;
        }
        if !pending.is_empty() {
            return Err(io::Error::other(format!(
                "BGZF block {next} was never compressed ({} blocks stranded)",
                pending.len()
            )));
        }

        writer.write_all(&BGZF_EOF)?;
        writer.flush()?;
        Ok(())
    }

    pub fn write_all(&mut self, data: &[u8]) -> io::Result<()> {
        let mut offset = 0;
        while offset < data.len() {
            let free = CHUNK_SIZE - self.pending.len();
            let take = free.min(data.len() - offset);
            self.pending.extend_from_slice(&data[offset..offset + take]);
            offset += take;
            if self.pending.len() == CHUNK_SIZE {
                self.flush_pending()?;
            }
        }
        Ok(())
    }

    fn flush_pending(&mut self) -> io::Result<()> {
        if self.pending.is_empty() {
            return Ok(());
        }
        let Some(tx) = self.tx.as_ref() else {
            return Err(io::Error::new(io::ErrorKind::BrokenPipe, "BGZF writer already finished"));
        };
        let seq = self.sequence;
        self.sequence += 1;
        let chunk = std::mem::replace(&mut self.pending, Vec::with_capacity(CHUNK_SIZE));
        if tx.send((chunk, seq)).is_err() {
            self.tx = None;
            return Err(self.thread_error("BGZF writer stopped early"));
        }
        Ok(())
    }

    /// The underlying cause when the pipeline stopped (disk full, closed
    /// pipe, compression failure), taken from the writer thread.
    fn thread_error(&mut self, what: &str) -> io::Error {
        match self.writer_thread.take().map(|h| h.join()) {
            Some(Ok(Err(e))) => io::Error::new(e.kind(), format!("{what}: {e}")),
            Some(Ok(Ok(()))) => io::Error::new(io::ErrorKind::BrokenPipe, what.to_string()),
            Some(Err(_)) => io::Error::other(format!("{what}: BGZF writer thread panicked")),
            None => io::Error::new(io::ErrorKind::BrokenPipe, what.to_string()),
        }
    }

    fn compress_block(data: &[u8], compressor: &mut Compressor, scratch: &mut Vec<u8>) -> io::Result<Vec<u8>> {
        let bound = compressor.deflate_compress_bound(data.len());
        if scratch.len() < bound {
            scratch.resize(bound, 0);
        }
        // Each BGZF block is an independent raw-deflate stream.
        let n = compressor
            .deflate_compress(data, scratch)
            .map_err(|e| io::Error::other(format!("BGZF compression failed: {e:?}")))?;

        let total = BGZF_HEADER.len() + n + 8;
        if total > 65536 {
            return Err(io::Error::other(format!("BGZF block of {total} bytes exceeds the 64 KiB limit")));
        }
        let mut block = Vec::with_capacity(total);
        block.extend_from_slice(&BGZF_HEADER);
        block.extend_from_slice(&scratch[..n]);
        block.extend_from_slice(&compute_crc32(data).to_le_bytes());
        block.extend_from_slice(&(data.len() as u32).to_le_bytes());
        let bsize = (block.len() - 1) as u16;
        block[16..18].copy_from_slice(&bsize.to_le_bytes());
        Ok(block)
    }

    fn finish_inner(&mut self) -> io::Result<()> {
        if self.finished {
            return Ok(());
        }
        self.finished = true;
        let flushed = self.flush_pending();
        drop(self.tx.take());
        let joined = match self.writer_thread.take() {
            Some(handle) => handle
                .join()
                .map_err(|_| io::Error::other("BGZF writer thread panicked"))
                .and_then(|r| r),
            None => Ok(()),
        };
        flushed?;
        joined
    }

    /// Compress the remaining data, write the EOF marker and flush.
    pub fn finish(mut self) -> io::Result<()> {
        self.finish_inner()
    }
}

impl Drop for BgzfWriter {
    fn drop(&mut self) {
        let _ = self.finish_inner();
    }
}

impl Write for BgzfWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.write_all(buf)?;
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        self.flush_pending()
    }
}

#[cfg(test)]
#[path = "../../tests/unit/bgzf_writer.rs"]
mod tests;
