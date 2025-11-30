pub mod bgzf;
pub mod csi;
pub mod kbi;
pub mod util;
pub mod vcf;

pub use bgzf::{BgzfReader, BgzfWriter, VirtualPosition};
pub use csi::{build_csi_index, read_csi_index, CsiQuery};
pub use kbi::{build_kbi_index, KbiBuilder, KbiIndex, KbiStats};
pub use util::{chr_id_to_name, chr_name_to_id, detect_format, ChrId, GenomicKey, Region, VcfFormat};
pub use vcf::{fetch_line, VcfReader, VcfRecord};