mod fetch;
mod mmap;
mod parser;
mod reader;
pub mod simd;
pub mod structs;
mod unified_reader;

pub use fetch::fetch_line;
pub use mmap::MmapVcfParser;
pub use parser::{extract_contig_id, parse_vcf_full_line, InfoParser, VcfParser};
pub use reader::VcfReader;
pub use structs::{Result, VcfError, VcfFields, VcfParsedRecord, VcfRecord};
pub use unified_reader::UnifiedVcfReader;
