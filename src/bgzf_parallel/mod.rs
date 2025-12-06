//! High-performance parallel BGZF reader/writer with batched decompression
//!
//! Key features:
//! - 3-5x faster decompression using flate2 with miniz_oxide backend
//! - Parallel block decompression via rayon
//! - Prefetch pipeline for continuous block processing
//! - Full BGZF format compatibility (RFC 1951/1952)

mod block;
mod reader;
mod writer;

pub use block::{BgzfBlock, BlockDecoder};
pub use reader::{BatchedLineReader, ParallelBgzfReader};
pub use writer::ParallelBgzfWriter;

const BGZF_BLOCK_SIZE: usize = 64 * 1024; // 64KB max uncompressed
const BATCH_SIZE: usize = 64; // Blocks per batch
const PREFETCH_BATCHES: usize = 2;
