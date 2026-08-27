//! `.ktile` — Kira columnar VCF sidecar format.
//!
//! See [`format`] for the on-disk binary layout, [`writer`] for the one-shot
//! VCF → `.ktile` builder, and [`reader`] for the mmap-backed reader.

pub mod format;
pub mod freshness;
pub mod reader;
pub mod writer;

pub use format::{KTILE_MAGIC, KTILE_VERSION, KtileHeader};
pub use freshness::{KtileFreshness, check_ktile_freshness, ktile_path_for};
pub use reader::KtileReader;
pub use writer::write_ktile_from_vcf;
