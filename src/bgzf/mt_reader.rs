//! Multithreaded BGZF reader for streaming workloads.
//!
//! Thin wrapper over [`noodles_bgzf::io::MultithreadedReader`] that owns the
//! decompression worker pool. Worker count is configurable via
//! `KIRA_BT_BGZF_THREADS`.

use std::fs::File;
use std::io::{self, BufRead, Read};
use std::num::NonZeroUsize;
use std::path::Path;

use noodles_bgzf as bgzf;

use crate::bgzf::structs::{Result, VirtualPosition};

/// Default decompression worker count. Honours `KIRA_BT_BGZF_THREADS`;
/// otherwise a share of the process thread budget ([`crate::threads`]).
pub fn default_bgzf_workers() -> NonZeroUsize {
    if let Ok(v) = std::env::var("KIRA_BT_BGZF_THREADS")
        && let Ok(n) = v.parse::<usize>()
        && n > 0
    {
        return NonZeroUsize::new(n).unwrap();
    }
    NonZeroUsize::new(crate::threads::decompress_workers()).unwrap()
}

/// Multithreaded BGZF reader for streaming (annotate) workloads.
pub struct MtBgzfReader<R: Read + Send + 'static> {
    inner: bgzf::io::MultithreadedReader<R>,
}

impl MtBgzfReader<File> {
    /// Opens `path` with the default worker count.
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        Ok(Self::with_worker_count(file, default_bgzf_workers()))
    }
}

impl<R: Read + Send + 'static> MtBgzfReader<R> {
    pub fn new(reader: R) -> Self {
        Self::with_worker_count(reader, default_bgzf_workers())
    }

    pub fn with_worker_count(reader: R, workers: NonZeroUsize) -> Self {
        Self {
            inner: bgzf::io::MultithreadedReader::with_worker_count(workers, reader),
        }
    }

    /// Current virtual position (the currently-buffered block).
    pub fn virtual_position(&self) -> VirtualPosition {
        self.inner.virtual_position().into()
    }
}

impl<R: Read + Send + 'static> Read for MtBgzfReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.inner.read(buf)
    }
}

impl<R: Read + Send + 'static> BufRead for MtBgzfReader<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        self.inner.fill_buf()
    }

    fn consume(&mut self, amt: usize) {
        self.inner.consume(amt)
    }
}
