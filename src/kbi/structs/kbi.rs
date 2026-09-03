use std::io;
use std::mem;

use bytemuck::{Pod, Zeroable};
use kira_kv_engine::IndexError;
use thiserror::Error;

pub const KBI_MAGIC: [u8; 8] = *b"KBIV0004";
pub const KBI_VERSION: u32 = 4;
pub const ENDIAN_TAG: u32 = 0x01020304;

#[derive(Debug, Error)]
pub enum KbiError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
    #[error("KV error: {0}")]
    Kv(#[from] IndexError),
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

/// On-disk layout: header, kv index, sorted keys, offsets, slot table, contig names.
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
    pub off_slots: u64,
    pub n_slots: u64,
    pub off_names: u64,
}

const _: () = assert!(mem::size_of::<KbiHeader>() == 80);

impl KbiHeader {
    pub fn new(n_entries: usize, index_len: usize, n_slots: usize) -> Self {
        let header_size = mem::size_of::<Self>() as u64;
        let keys_size = (n_entries * mem::size_of::<u64>()) as u64;
        let slots_size = (n_slots * mem::size_of::<u32>()) as u64;
        let index_size = index_len as u64;
        let off_index = header_size;
        let off_keys = off_index + index_size;
        let off_offsets = off_keys + keys_size;
        let off_slots = off_offsets + keys_size;
        let off_names = off_slots + slots_size;

        Self {
            magic: KBI_MAGIC,
            version: KBI_VERSION,
            endian: ENDIAN_TAG,
            n_entries: n_entries as u64,
            index_len: index_size,
            off_index,
            off_keys,
            off_offsets,
            off_slots,
            n_slots: n_slots as u64,
            off_names,
        }
    }

    pub fn validate(&self) -> Result<()> {
        if self.magic != KBI_MAGIC {
            return Err(KbiError::InvalidFormat(
                "Invalid magic bytes (rebuild the .kbi index with this version)".into(),
            ));
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
