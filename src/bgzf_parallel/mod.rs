mod block;
mod reader;
mod writer;
mod writer_optimized;

pub use block::{BgzfBlock, BlockDecoder};
pub use reader::{BatchedLineReader, ParallelBgzfReader};
pub use writer::ParallelBgzfWriter;
pub use writer_optimized::OptimizedBgzfWriter;

const BGZF_BLOCK_SIZE: usize = 64 * 1024;
const BATCH_SIZE: usize = 64;
const PREFETCH_BATCHES: usize = 2;
