use std::io;
use std::mem;

use bytemuck::{Pod, Zeroable};
use kira_kv_engine::HybridError;
use thiserror::Error;

pub const KBI_MAGIC: [u8; 8] = *b"KBIV0003";
pub const KBI_VERSION: u32 = 3;
pub const ENDIAN_TAG: u32 = 0x01020304;

#[derive(Debug, Error)]
pub enum KbiError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
    #[error("KV error: {0}")]
    Kv(#[from] HybridError),
    #[error("VCF error: {0}")]
    Vcf(#[from] crate::vcf::VcfError),
    #[error("Empty dataset")]
    EmptyDataset,
    #[error("Invalid format: {0}")]
    InvalidFormat(String),
    #[error("Version mismatch: expected {expected}, got {got}")]
    VersionMismatch { expected: u32, got: u32 },
}

pub type Result<T> = std::result::Result<T, KbiError>;

#[repr(C, packed)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
pub struct KbiHeader {
    pub magic: [u8; 8],
    pub version: u32,
    pub endian: u32,
    pub n_entries: u64,
    pub index_len: u64,
    pub off_index: u64,
    pub off_keys: u64,
    pub off_offsets: u64,
    _reserved: [u8; 8],
}

const _: () = assert!(mem::size_of::<KbiHeader>() == 64);

impl KbiHeader {
    pub fn new(n_entries: usize, index_len: usize) -> Self {
        let header_size = mem::size_of::<Self>() as u64;
        let keys_size = (n_entries * mem::size_of::<u64>()) as u64;
        let index_size = index_len as u64;

        Self {
            magic: KBI_MAGIC,
            version: KBI_VERSION,
            endian: ENDIAN_TAG,
            n_entries: n_entries as u64,
            index_len: index_size,
            off_index: header_size,
            off_keys: header_size + index_size,
            off_offsets: header_size + index_size + keys_size,
            _reserved: [0; 8],
        }
    }

    pub fn validate(&self) -> Result<()> {
        if self.magic != KBI_MAGIC {
            return Err(KbiError::InvalidFormat("Invalid magic bytes".into()));
        }
        if self.version != KBI_VERSION {
            return Err(KbiError::VersionMismatch {
                expected: KBI_VERSION,
                got: self.version,
            });
        }
        if self.endian != ENDIAN_TAG {
            return Err(KbiError::InvalidFormat("Endianness mismatch".into()));
        }
        Ok(())
    }
}

pub struct KbiStats {
    pub entries: usize,
    pub memory_bytes: usize,
    pub bytes_per_key: f64,
    pub file_size: u64,
}
