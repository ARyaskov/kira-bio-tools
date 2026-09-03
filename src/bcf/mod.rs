//! BCF 2.2 binary I/O. Supports -O u (uncompressed BCF) and -O b (BGZF-BCF).

pub mod typed;
pub mod header;
pub mod record;
pub mod reader;
pub mod writer;

pub use header::{BcfHeaderDict, HdrField, parse_header_to_dict};
pub use reader::{BcfInput, BcfReader};
pub use writer::BcfWriter;

pub const BCF_MAGIC: &[u8] = b"BCF\x02\x02";
