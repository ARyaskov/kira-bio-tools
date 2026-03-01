mod reader;
mod simd;
pub mod structs;
mod utils;
mod writer;

pub use reader::{BgzfLineReader, BgzfReader};
pub use structs::{
    BGZF_BLOCK_SIZE, BGZF_EOF, BGZF_HEADER, BgzfBlock, BgzfError, CHUNK_SIZE, CompressedBlock,
    Result, VirtualPosition, WritePool,
};
pub use utils::is_bgzf;
pub use writer::BgzfWriter;
