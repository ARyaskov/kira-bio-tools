mod binning;
mod builder;
mod fetch;
pub mod structs;
pub mod utils;

pub use binning::{BinEntry, BinIndex, IndexKind, RefIndex, RefMeta, TabixHeader};
pub use builder::{
    IndexBuilder, build_csi_index, build_index, build_index_in_memory, find_index_for,
    vcf_line_interval,
};
pub use fetch::IndexedVcfReader;
pub use structs::{CsiError, Result};
