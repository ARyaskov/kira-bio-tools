use crossbeam_channel::{bounded, Receiver, Sender};
use flate2::write::DeflateEncoder;
use flate2::Compression;
use libdeflater::Compressor;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::thread;

#[cfg(target_arch = "x86_64")]
use core::arch::x86_64::_mm_crc32_u64;

#[target_feature(enable = "avx2")]
unsafe fn crc32_hw(mut crc: u32, data: &[u8]) -> u32 {
    use core::arch::x86_64::*;

    let mut i = 0;

    while i + 32 <= data.len() {
        let ptr = data.as_ptr().add(i) as *const __m256i;
        let chunk = _mm256_loadu_si256(ptr);

        let a = _mm256_extract_epi64(chunk, 0) as u64;
        let b = _mm256_extract_epi64(chunk, 1) as u64;
        let c = _mm256_extract_epi64(chunk, 2) as u64;
        let d = _mm256_extract_epi64(chunk, 3) as u64;

        crc = _mm_crc32_u64(crc as u64, a) as u32;
        crc = _mm_crc32_u64(crc as u64, b) as u32;
        crc = _mm_crc32_u64(crc as u64, c) as u32;
        crc = _mm_crc32_u64(crc as u64, d) as u32;

        i += 32;
    }

    while i + 8 <= data.len() {
        let chunk = *(data.as_ptr().add(i) as *const u64);
        crc = _mm_crc32_u64(crc as u64, chunk) as u32;
        i += 8;
    }

    while i < data.len() {
        crc = _mm_crc32_u8(crc, *data.get_unchecked(i));
        i += 1;
    }

    crc
}

#[inline(always)]
unsafe fn memcpy_avx(dst: *mut u8, src: *const u8, len: usize) {
    let mut i = 0;

    while i + 32 <= len {
        let v = core::arch::x86_64::_mm256_loadu_si256(src.add(i) as *const _);
        core::arch::x86_64::_mm256_storeu_si256(dst.add(i) as *mut _, v);
        i += 32;
    }

    while i < len {
        *dst.add(i) = *src.add(i);
        i += 1;
    }
}

#[inline(always)]
unsafe fn copy_bgzf_header(dst: *mut u8) {
    let mut tmp = [0u8; 32];
    tmp[..18].copy_from_slice(&BGZF_HEADER);
    let v = core::arch::x86_64::_mm256_loadu_si256(tmp.as_ptr() as *const _);
    core::arch::x86_64::_mm256_storeu_si256(dst as *mut _, v);
}

const BGZF_HEADER: [u8; 18] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x00, 0x00,
];

const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

const CHUNK_SIZE: usize = 256 * 1024;
const CHANNEL_DEPTH: usize = 256;
const WRITER_BUFFER_SIZE: usize = 128 * 1024 * 1024;

struct CompressedBlock {
    data: Vec<u8>,
    sequence: usize,
}

#[derive(Clone)]
pub struct WritePool {
    pool: std::sync::Arc<std::sync::Mutex<Vec<Vec<u8>>>>,
}

impl WritePool {
    pub fn new(n: usize) -> Self {
        let mut v = Vec::with_capacity(n);
        for _ in 0..n {
            let mut buf = Vec::with_capacity(CHUNK_SIZE);
            unsafe {
                buf.set_len(CHUNK_SIZE);
            }
            v.push(buf);
        }
        Self {
            pool: std::sync::Arc::new(std::sync::Mutex::new(v)),
        }
    }

    #[inline]
    pub fn get(&self) -> Vec<u8> {
        self.pool.lock().unwrap().pop().unwrap_or_else(|| {
            let mut b = Vec::with_capacity(CHUNK_SIZE);
            unsafe {
                b.set_len(CHUNK_SIZE);
            }
            b
        })
    }

    #[inline]
    pub fn put(&self, mut buf: Vec<u8>) {
        unsafe {
            buf.set_len(CHUNK_SIZE);
        }
        self.pool.lock().unwrap().push(buf);
    }
}

pub struct OptimizedBgzfWriter {
    tx: Sender<(Vec<u8>, usize)>,
    writer_thread: Option<thread::JoinHandle<io::Result<()>>>,
    sequence: AtomicUsize,
    num_workers: usize,
    write_pool: WritePool,
}

impl OptimizedBgzfWriter {
    pub fn create<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::with_compression(path, Compression::new(3))
    }

    pub fn with_compression<P: AsRef<Path>>(path: P, compression: Compression) -> io::Result<Self> {
        let num_workers = 4;

        let file = File::create(path)?;
        let writer = BufWriter::with_capacity(WRITER_BUFFER_SIZE, file);
        let pool = WritePool::new(num_workers * 4);
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

                    Self::compression_worker(rx, tx, &mut compressor, &mut comp_buf, comp)
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
        compression: Compression,
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

    #[inline(always)]
    fn find_next_newline_end(data: &[u8], mut i: usize) -> usize {
        while i < data.len() {
            if data[i] == b'\n' {
                return i + 1;
            }
            i += 1;
        }
        data.len()
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
            copy_bgzf_header(block.as_mut_ptr());
        }

        let bound = compressor.gzip_compress_bound(data.len());
        if scratch_buf.len() < bound {
            scratch_buf.resize(bound, 0);
        }

        let comp_len = compressor
            .gzip_compress(data, scratch_buf)
            .expect("gzip compression failed");

        unsafe {
            let old_len = block.len();
            block.resize(old_len + comp_len, 0);
            memcpy_avx(
                block.as_mut_ptr().add(old_len),
                scratch_buf.as_ptr(),
                comp_len,
            );
        }

        let crc = unsafe { crc32_hw(0, data) };
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
