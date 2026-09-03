mod mt_reader;
mod reader;
pub mod seek;
mod simd;
pub mod structs;
mod utils;
mod writer;

pub use mt_reader::{MtBgzfReader, default_bgzf_workers};
pub use reader::{BgzfLineReader, BgzfReader};
pub use structs::{
    BGZF_BLOCK_SIZE, BGZF_EOF, BGZF_HEADER, BgzfBlock, BgzfError, CHUNK_SIZE, CompressedBlock,
    Result, VirtualPosition, is_bgzf_header, is_gzip_header,
};
pub use seek::BgzfSeekReader;
pub use utils::is_bgzf;
pub use writer::{BgzfWriter, FILE_BUFFER_SIZE, STREAM_BUFFER_SIZE};
