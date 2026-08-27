#![allow(dead_code)] // WIP / alternate-path helpers (multiallelic LK, codon phasing, CNV spans, etc.) kept intentionally

pub mod annotate;
pub mod bam;
pub mod bcf;
pub mod bgzf;
pub mod call;
pub mod cli;
pub mod cnv;
pub mod csi;
pub mod filter;
pub mod filter_arch;
pub mod kbi;
pub mod norm;
pub mod roh;
pub mod util;
pub mod vcf;

pub use bgzf::{BgzfReader, BgzfWriter, VirtualPosition};
pub use csi::{CsiQuery, build_csi_index, read_csi_index};
pub use kbi::{KbiBuilder, KbiIndex, KbiStats, build_kbi_index};
pub use util::{
    ChrId, GenomicKey, Region, VcfFormat, chr_id_to_name, chr_name_to_id, detect_format,
};
pub use vcf::{InfoParser, VcfParser, VcfReader, VcfRecord, fetch_line};
