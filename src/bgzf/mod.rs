mod reader;
mod simd;
pub mod structs;
mod utils;
mod writer;

pub use reader::{BgzfLineReader, BgzfReader};
pub use structs::{
    BgzfBlock, BgzfError, CompressedBlock, Result, VirtualPosition, WritePool, BGZF_BLOCK_SIZE,
    BGZF_EOF, BGZF_HEADER, CHUNK_SIZE,
};
pub use utils::is_bgzf;
pub use writer::BgzfWriter;
