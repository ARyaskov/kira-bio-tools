pub mod annotate;
pub mod bgzf;
pub mod cli;
pub mod csi;
pub mod filter;
pub mod kbi;
pub mod norm;
pub mod util;
pub mod vcf;

pub use bgzf::{BgzfReader, BgzfWriter, VirtualPosition};
pub use csi::{build_csi_index, read_csi_index, CsiQuery};
pub use kbi::{build_kbi_index, KbiBuilder, KbiIndex, KbiStats};
pub use util::{
    chr_id_to_name, chr_name_to_id, detect_format, ChrId, GenomicKey, Region, VcfFormat,
};
pub use vcf::{fetch_line, VcfReader, VcfRecord};

pub mod bgzf_parallel;
pub mod vcf_parser_fast;

// Re-export for convenience
pub use bgzf_parallel::{BatchedLineReader, ParallelBgzfReader, ParallelBgzfWriter};
pub use vcf_parser_fast::{FastVcfParser, InfoParser, VcfFields};
