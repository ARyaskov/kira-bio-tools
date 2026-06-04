use crossbeam_channel::{Receiver, Sender, bounded};
use flate2::{Compress, Compression, FlushCompress, Status};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

/// Worst-case raw-deflate output size for `input` bytes. Matches the historic
/// libdeflate bound (zlib's `deflateBound` formula) — safe for both miniz_oxide
/// and any C backend if we ever reintroduce one. Used to size scratch buffers.
#[inline]
fn deflate_compress_bound(input: usize) -> usize {
    input + (input >> 12) + (input >> 14) + (input >> 25) + 13
}

use super::simd::{compute_crc32, fast_copy_bgzf_header, fast_memcpy};
use super::structs::{BGZF_EOF, CHUNK_SIZE, CompressedBlock, WritePool};

/// Bounded queue between the writer pre-input stage and the
/// compressor pool. 512 deep — about 32 MB of in-flight pre-compression
/// blocks (CHUNK_SIZE × 512). Larger value gives the compressor pool
/// more parallelism opportunity but burns RAM; 256→512 keeps the
/// compressor pool fed without saturating RAM.
const CHANNEL_DEPTH: usize = 512;
/// Output BufWriter capacity. 128 MB sweet spot — larger values
/// (we tried 256 MB) cause transient stalls during periodic flushes
/// that block the writer thread for ~250 ms, which in turn fills the
/// upstream compressor channel and back-pressures the worker. The
/// observed regression was 16+ sec of `send` block time on GPU path
/// when the buffer was 256 MB.
const WRITER_BUFFER_SIZE: usize = 128 * 1024 * 1024;

pub struct BgzfWriter {
    tx: Sender<(Vec<u8>, usize)>,
    writer_thread: Option<thread::JoinHandle<io::Result<()>>>,
    sequence: AtomicUsize,
    pending: Vec<u8>,
    _num_workers: usize,
    _write_pool: WritePool,
}

impl BgzfWriter {
    pub fn create<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::with_compression(path, Compression::new(1))
    }

    pub fn with_compression<P: AsRef<Path>>(path: P, compression: Compression) -> io::Result<Self> {
        // Compressor pool size: up to 32 threads on big servers (was
        // capped at 16). For typical 8-16-core hosts this is min'd to
        // the actual CPU count. Doubles the parallel compression
        // throughput on machines with >16 logical cores.
        let num_workers = num_cpus::get().min(32).max(2);

        let file = File::create(path)?;
        let writer = BufWriter::with_capacity(WRITER_BUFFER_SIZE, file);
        let pool = WritePool::new(num_workers * 4, CHUNK_SIZE);
        let (chunk_tx, chunk_rx) = bounded::<(Vec<u8>, usize)>(CHANNEL_DEPTH);
        let (block_tx, block_rx) = bounded::<CompressedBlock>(CHANNEL_DEPTH);

        let compression_workers: Vec<_> = (0..num_workers)
            .map(|_| {
                let rx = chunk_rx.clone();
                let tx = block_tx.clone();
                let comp = compression;

                thread::spawn(move || {
                    // flate2::Compress with zlib_header=false → raw deflate, same
                    // wire format libdeflater produced. Level mapping is identical.
                    let mut compressor = Compress::new(comp, false);
                    let mut comp_buf = vec![0u8; 512 * 1024];

                    Self::compression_worker(rx, tx, &mut compressor, &mut comp_buf)
                })
            })
            .collect();

        drop(chunk_rx);
        drop(block_tx);

        let writer_thread =
            thread::spawn(move || Self::writer_worker(block_rx, writer, compression_workers));

        Ok(Self {
            tx: chunk_tx,
            writer_thread: Some(writer_thread),
            sequence: AtomicUsize::new(0),
            pending: Vec::with_capacity(CHUNK_SIZE),
            _num_workers: num_workers,
            _write_pool: pool,
        })
    }

    fn compression_worker(
        chunk_rx: Receiver<(Vec<u8>, usize)>,
        block_tx: Sender<CompressedBlock>,
        compressor: &mut Compress,
        scratch_buf: &mut Vec<u8>,
    ) -> io::Result<()> {
        while let Ok((data, seq)) = chunk_rx.recv() {
            let compressed = Self::compress_block(&data, compressor, scratch_buf)?;

            if block_tx
                .send(CompressedBlock {
                    data: compressed,
                    sequence: seq,
                })
                .is_err()
            {
                break;
            }
        }

        Ok(())
    }

    fn writer_worker(
        block_rx: Receiver<CompressedBlock>,
        mut writer: BufWriter<File>,
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
                .map_err(|_| io::Error::new(io::ErrorKind::Other, "Worker thread panicked"))??;
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
        let seq = self.sequence.fetch_add(1, Ordering::Relaxed);
        let chunk = std::mem::replace(&mut self.pending, Vec::with_capacity(CHUNK_SIZE));
        self.tx
            .send((chunk, seq))
            .map_err(|_| io::Error::new(io::ErrorKind::BrokenPipe, "Worker died"))
    }

    fn compress_block(
        data: &[u8],
        compressor: &mut Compress,
        scratch_buf: &mut Vec<u8>,
    ) -> io::Result<Vec<u8>> {
        let mut block = Vec::with_capacity(data.len() / 2 + 256);
        unsafe {
            block.resize(18, 0);
            fast_copy_bgzf_header(block.as_mut_ptr());
        }

        let bound = deflate_compress_bound(data.len());
        if scratch_buf.len() < bound {
            scratch_buf.resize(bound, 0);
        }

        // Reset the stateful Compress between blocks: each BGZF block is an
        // independent raw-deflate stream. Without reset() the second block
        // would inherit the LZ77 dictionary from the first one and the
        // FlushCompress::Finish from the prior call.
        compressor.reset();
        let before_out = compressor.total_out();
        let before_in = compressor.total_in();
        let status = compressor
            .compress(data, scratch_buf, FlushCompress::Finish)
            .map_err(|_| io::Error::new(io::ErrorKind::Other, "Compression failed"))?;
        // FlushCompress::Finish with a sufficient output buffer (deflate_compress_bound)
        // and a single-shot call must consume all input and reach StreamEnd. If it
        // didn't, the scratch buffer is undersized — a bug, not a runtime condition.
        if status != Status::StreamEnd
            || (compressor.total_in() - before_in) as usize != data.len()
        {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "Compression did not finish in one call (scratch too small?)",
            ));
        }
        let comp_len = (compressor.total_out() - before_out) as usize;

        unsafe {
            let old_len = block.len();
            block.resize(old_len + comp_len, 0);
            fast_memcpy(
                block.as_mut_ptr().add(old_len),
                scratch_buf.as_ptr(),
                comp_len,
            );
        }

        let crc = compute_crc32(data);
        block.extend_from_slice(&crc.to_le_bytes());
        block.extend_from_slice(&(data.len() as u32).to_le_bytes());

        let bsize = (block.len() - 1) as u16;
        block[16..18].copy_from_slice(&bsize.to_le_bytes());

        Ok(block)
    }

    pub fn finish(mut self) -> io::Result<()> {
        self.flush_pending()?;
        drop(self.tx);

        if let Some(handle) = self.writer_thread.take() {
            handle
                .join()
                .map_err(|_| io::Error::new(io::ErrorKind::Other, "Writer thread panicked"))??;
        }

        Ok(())
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
