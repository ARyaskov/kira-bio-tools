pub mod bitmap;
mod build;
mod builder;
mod index;
pub mod structs;

pub use bitmap::{ContigBitmap, PosBitmap};
pub use build::{build_kbi_index, build_kbi_index_with_bitmap};
pub use builder::KbiBuilder;
pub use index::KbiIndex;
pub use structs::{KbiError, KbiStats, Result};
