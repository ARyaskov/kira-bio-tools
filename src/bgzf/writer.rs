use crossbeam_channel::{bounded, Receiver, Sender};
use flate2::Compression;
use libdeflater::Compressor;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

use super::simd::{compute_crc32, fast_copy_bgzf_header, fast_memcpy};
use super::structs::{CompressedBlock, WritePool, BGZF_EOF, CHUNK_SIZE};

const CHANNEL_DEPTH: usize = 256;
const WRITER_BUFFER_SIZE: usize = 128 * 1024 * 1024;

pub struct BgzfWriter {
    tx: Sender<(Vec<u8>, usize)>,
    writer_thread: Option<thread::JoinHandle<io::Result<()>>>,
    sequence: AtomicUsize,
    num_workers: usize,
    write_pool: WritePool,
}

impl BgzfWriter {
    pub fn create<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::with_compression(path, Compression::new(3))
    }

    pub fn with_compression<P: AsRef<Path>>(path: P, compression: Compression) -> io::Result<Self> {
        let num_workers = num_cpus::get().min(8).max(2);

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
                    let lvl = libdeflater::CompressionLvl::new(comp.level() as i32).unwrap();
                    let mut compressor = libdeflater::Compressor::new(lvl);
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
            num_workers,
            write_pool: pool,
        })
    }

    fn compression_worker(
        chunk_rx: Receiver<(Vec<u8>, usize)>,
        block_tx: Sender<CompressedBlock>,
        compressor: &mut Compressor,
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
            let end = (offset + CHUNK_SIZE).min(data.len());
            let seq = self.sequence.fetch_add(1, Ordering::Relaxed);

            let chunk = data[offset..end].to_vec();

            self.tx
                .send((chunk, seq))
                .map_err(|_| io::Error::new(io::ErrorKind::BrokenPipe, "Worker died"))?;

            offset = end;
        }

        Ok(())
    }

    fn compress_block(
        data: &[u8],
        compressor: &mut Compressor,
        scratch_buf: &mut Vec<u8>,
    ) -> io::Result<Vec<u8>> {
        let mut block = Vec::with_capacity(data.len() / 2 + 256);
        unsafe {
            block.resize(18, 0);
            fast_copy_bgzf_header(block.as_mut_ptr());
        }

        let bound = compressor.gzip_compress_bound(data.len());
        if scratch_buf.len() < bound {
            scratch_buf.resize(bound, 0);
        }

        let comp_len = compressor
            .gzip_compress(data, scratch_buf)
            .map_err(|_| io::Error::new(io::ErrorKind::Other, "Compression failed"))?;

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
        Ok(())
    }
}
