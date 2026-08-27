//! `.ktile` on-disk binary layout. Little-endian; mmap-friendly.

use bytemuck::{Pod, Zeroable};

pub const KTILE_MAGIC: [u8; 8] = *b"KIRA_TIL";
pub const KTILE_VERSION: u32 = 1;

/// Fixed 128-byte header at file offset 0.
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct KtileHeader {
    pub magic: [u8; 8],
    pub version: u32,
    pub flags: u32,

    pub n_records: u64,

    pub headers_off: u64,
    pub headers_len: u64,

    pub line_offsets_off: u64,
    pub chr_ids_off: u64,
    pub positions_off: u64,

    pub line_pool_off: u64,
    pub line_pool_len: u64,

    /// Source file byte size at build time; 0 to skip freshness check.
    pub source_size: u64,
    /// Source mtime at build time (Unix seconds); 0 to skip freshness check.
    pub source_mtime_unix: u64,

    pub off_ref_offsets: u64,
    pub off_ref_lens: u64,
    pub off_alt_offsets: u64,
    pub off_alt_lens: u64,

    /// Lines per compressed chunk; 0 when uncompressed.
    pub lines_per_chunk: u32,
    pub n_chunks: u32,
    pub off_chunk_index: u64,
}

impl KtileHeader {
    pub const SIZE: usize = std::mem::size_of::<Self>();

    /// Returns `Ok(header)` if magic + version are recognized; otherwise
    /// returns a descriptive error. Does not validate offset bounds — the
    /// reader does that against the actual mmap length.
    pub fn validate(&self) -> Result<(), KtileError> {
        if self.magic != KTILE_MAGIC {
            return Err(KtileError::BadMagic);
        }
        if self.version != KTILE_VERSION {
            return Err(KtileError::UnsupportedVersion {
                found: self.version,
                expected: KTILE_VERSION,
            });
        }
        Ok(())
    }

    pub fn has_ref_alt_columns(&self) -> bool {
        (self.flags & flags::HAS_REF_ALT_COLUMNS) != 0
    }

    pub fn has_compressed_pool(&self) -> bool {
        (self.flags & flags::HAS_COMPRESSED_POOL) != 0
    }
}

pub mod flags {
    pub const HAS_REF_ALT_COLUMNS: u32 = 1 << 0;
    pub const HAS_COMPRESSED_POOL: u32 = 1 << 1;
}

#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct CompressedChunkEntry {
    pub compressed_off: u64,
    pub uncompressed_off: u64,
    pub compressed_size: u32,
    pub uncompressed_size: u32,
}

#[derive(Debug, thiserror::Error)]
pub enum KtileError {
    #[error("not a ktile file (bad magic)")]
    BadMagic,
    #[error("unsupported ktile version: found {found}, expected {expected}")]
    UnsupportedVersion { found: u32, expected: u32 },
    #[error("ktile file truncated: section {section} extends past file end")]
    Truncated { section: &'static str },
    #[error("i/o error: {0}")]
    Io(#[from] std::io::Error),
}

#[cfg(test)]
#[path = "../../../tests/unit/annotate_ktile_format.rs"]
mod tests;
