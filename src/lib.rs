// Structural lints this codebase deliberately does not follow: wide
// signatures mirror bcftools' C functions, index loops keep ported numeric
// code readable, and the k-way merges need their iterators dropped before
// the loop body mutates the readers. Everything else clippy reports is fixed.
#![allow(
    clippy::too_many_arguments,
    clippy::type_complexity,
    clippy::needless_range_loop,
    clippy::large_enum_variant,
    clippy::if_same_then_else,
    clippy::field_reassign_with_default,
    clippy::while_let_loop,
    clippy::single_match,
    clippy::explicit_counter_loop
)]

pub mod af_file;
pub mod align;
pub mod annotate;
pub mod bam;
pub mod bcf;
pub mod bgzf;
pub mod call;
pub mod cli;
pub mod cnv;
pub mod consensus;
pub mod csi;
pub mod csq;
pub mod fasta;
pub mod filter;
pub mod filter_arch;
pub mod fx;
pub mod hmm;
pub mod kbi;
pub mod norm;
pub mod regions;
pub mod roh;
pub mod sketch;
pub mod threads;
pub mod util;
pub mod vcf;
pub mod zarr;

pub use bgzf::{BgzfReader, BgzfWriter, VirtualPosition};
pub use csi::{
    BinIndex, IndexKind, IndexedVcfReader, build_csi_index, build_index, find_index_for,
};
pub use kbi::{KbiBuilder, KbiIndex, KbiStats, build_kbi_index};
pub use regions::RegionSet;
pub use util::{ChrId, GenomicKey, Region, VcfFormat, chr_name_to_id, detect_format};
pub use vcf::{
    ContigDict, HeaderInfo, InfoParser, LineFetcher, VcfParser, VcfReader, VcfRecord, VcfSink,
    fetch_line,
};

#[cfg(test)]
#[path = "../tests/unit/props.rs"]
mod props;
