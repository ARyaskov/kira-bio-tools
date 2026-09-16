mod binning;
mod builder;
mod fetch;
pub mod structs;
pub mod utils;

pub use binning::{BinEntry, BinIndex, IndexKind, RefIndex, RefMeta, TabixHeader};
pub use builder::{IndexBuilder, build_csi_index, build_index, build_index_in_memory, find_index_for, vcf_line_interval, TabixConf, build_tabix_index, tabix_line_interval, TBX_GENERIC, TBX_SAM, TBX_VCF, TBX_UCSC};
pub use fetch::IndexedVcfReader;
pub use structs::{CsiError, Result};
