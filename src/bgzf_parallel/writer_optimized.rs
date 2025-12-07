use crossbeam_channel::{bounded, Receiver, Sender};
use flate2::write::DeflateEncoder;
use flate2::Compression;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::thread;

const BGZF_HEADER: [u8; 18] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x00, 0x00,
];

const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

const CHUNK_SIZE: usize = 56 * 1024;
const CHANNEL_DEPTH: usize = 256;
const WRITER_BUFFER_SIZE: usize = 128 * 1024 * 1024;

struct CompressedBlock {
    data: Vec<u8>,
    sequence: usize,
}

pub struct OptimizedBgzfWriter {
    tx: Sender<(Vec<u8>, usize)>,
    writer_thread: Option<thread::JoinHandle<io::Result<()>>>,
    sequence: AtomicUsize,
    num_workers: usize,
}

impl OptimizedBgzfWriter {
    pub fn create<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::with_compression(path, Compression::new(3))
    }

    pub fn with_compression<P: AsRef<Path>>(path: P, compression: Compression) -> io::Result<Self> {
        let num_workers = rayon::current_num_threads();

        let file = File::create(path)?;
        let writer = BufWriter::with_capacity(WRITER_BUFFER_SIZE, file);

        let (chunk_tx, chunk_rx) = bounded::<(Vec<u8>, usize)>(CHANNEL_DEPTH);
        let (block_tx, block_rx) = bounded::<CompressedBlock>(CHANNEL_DEPTH);

        let compression_workers: Vec<_> = (0..num_workers)
            .map(|_| {
                let rx = chunk_rx.clone();
                let tx = block_tx.clone();
                let comp = compression;

                thread::spawn(move || Self::compression_worker(rx, tx, comp))
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
        })
    }

    fn compression_worker(
        chunk_rx: Receiver<(Vec<u8>, usize)>,
        block_tx: Sender<CompressedBlock>,
        compression: Compression,
    ) -> io::Result<()> {
        while let Ok((data, seq)) = chunk_rx.recv() {
            let compressed = Self::compress_block(&data, compression)?;

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
        let mut pending_blocks = std::collections::BTreeMap::new();
        let mut next_expected = 0usize;

        for block in block_rx.iter() {
            pending_blocks.insert(block.sequence, block.data);

            while let Some(data) = pending_blocks.remove(&next_expected) {
                writer.write_all(&data)?;
                next_expected += 1;
            }
        }

        for data in pending_blocks.into_values() {
            writer.write_all(&data)?;
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
            let chunk_end = (offset + CHUNK_SIZE).min(data.len());
            let chunk = data[offset..chunk_end].to_vec();
            let seq = self.sequence.fetch_add(1, Ordering::Relaxed);

            self.tx
                .send((chunk, seq))
                .map_err(|_| io::Error::new(io::ErrorKind::BrokenPipe, "Worker died"))?;

            offset = chunk_end;
        }

        Ok(())
    }

    fn compress_block(data: &[u8], level: Compression) -> io::Result<Vec<u8>> {
        let mut block = Vec::with_capacity(data.len() / 2 + 256);
        block.extend_from_slice(&BGZF_HEADER);

        let mut encoder = DeflateEncoder::new(Vec::new(), level);
        encoder.write_all(data)?;
        let compressed = encoder.finish()?;

        block.extend_from_slice(&compressed);

        let crc = crc32fast::hash(data);
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

impl Write for OptimizedBgzfWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.write_all(buf)?;
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
