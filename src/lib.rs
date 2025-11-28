//! # Kira Bio Tools
//!
//! High-performance VCF indexer using learned indexes (PGM + MPH hybrid).
//!
//! ## Features
//!
//! - O(1) lookup for genomic positions
//! - Memory-mapped persistence for instant loading
//! - Support for 100M+ variants
//!
//! ## Example
//!
//! ```no_run
//! use kira_bio_tools::{VcfIndex, VcfIndexBuilder, GenomicKey};
//!
//! // Build index from VCF
//! let mut builder = VcfIndexBuilder::new();
//! builder.add(GenomicKey::new(1, 12345), 1024)?;
//! builder.add(GenomicKey::new(1, 67890), 2048)?;
//! let index = builder.build()?;
//!
//! // Query
//! if let Some(offset) = index.get(GenomicKey::new(1, 12345)) {
//!     println!("Found at byte offset: {}", offset);
//! }
//!
//! // Save to disk
//! index.save("variants.kbi")?;
//!
//! // Load with mmap
//! let index = VcfIndex::load_mmap("variants.kbi")?;
//! ```

mod persistence;
mod vcf_index;

pub use persistence::{INDEX_MAGIC, INDEX_VERSION, IndexHeader};
pub use vcf_index::{ChrId, GenomicKey, VcfIndex, VcfIndexBuilder};

/// Standard chromosome name to ID mapping (1-24 for chr1-22, X, Y, MT=25)
pub fn chr_name_to_id(name: &str) -> Option<u8> {
    let normalized = name.strip_prefix("chr").unwrap_or(name).to_uppercase();

    match normalized.as_str() {
        "X" => Some(23),
        "Y" => Some(24),
        "M" | "MT" => Some(25),
        s => s.parse::<u8>().ok().filter(|&n| n >= 1 && n <= 22),
    }
}

/// ID to standard chromosome name
pub fn chr_id_to_name(id: u8) -> Option<&'static str> {
    match id {
        1 => Some("chr1"),
        2 => Some("chr2"),
        3 => Some("chr3"),
        4 => Some("chr4"),
        5 => Some("chr5"),
        6 => Some("chr6"),
        7 => Some("chr7"),
        8 => Some("chr8"),
        9 => Some("chr9"),
        10 => Some("chr10"),
        11 => Some("chr11"),
        12 => Some("chr12"),
        13 => Some("chr13"),
        14 => Some("chr14"),
        15 => Some("chr15"),
        16 => Some("chr16"),
        17 => Some("chr17"),
        18 => Some("chr18"),
        19 => Some("chr19"),
        20 => Some("chr20"),
        21 => Some("chr21"),
        22 => Some("chr22"),
        23 => Some("chrX"),
        24 => Some("chrY"),
        25 => Some("chrM"),
        _ => None,
    }
}
