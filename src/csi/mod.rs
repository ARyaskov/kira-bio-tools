mod binning;
mod builder;
mod fetch;
pub mod structs;
pub mod utils;

pub use binning::{BinEntry, BinIndex, IndexKind, RefIndex, RefMeta, TabixHeader};
pub use builder::{
    IndexBuilder, TBX_GENERIC, TBX_SAM, TBX_UCSC, TBX_VCF, TabixConf, build_csi_index, build_index,
    build_index_in_memory, build_tabix_index, find_index_for, tabix_line_interval,
    vcf_line_interval,
};
pub use fetch::IndexedVcfReader;
pub use structs::{CsiError, Result};
