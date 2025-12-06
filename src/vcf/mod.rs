mod fetch;
mod mmap;
mod parser;
mod reader;
pub mod structs;

pub use fetch::fetch_line;
pub use mmap::MmapVcfParser;
pub use parser::{extract_contig_id, parse_vcf_full_line};
pub use reader::{BgzfVcfReader, GzipVcfReader, PlainVcfReader, VcfReader};
pub use structs::{Result, VcfError, VcfParsedRecord, VcfRecord};
