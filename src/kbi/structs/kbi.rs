use std::io;
use std::mem;

use bytemuck::{Pod, Zeroable};
use kira_kv_engine::IndexError;
use thiserror::Error;

pub const KBI_MAGIC: [u8; 8] = *b"KBIV0005";
pub const KBI_VERSION: u32 = 5;
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
        let off_index = mem::size_of::<Self>() as u64;
        let values_size = (n_entries as u64).saturating_mul(8);
        let off_keys = off_index.saturating_add(index_len as u64);
        let off_offsets = off_keys.saturating_add(values_size);
        let off_slots = off_offsets.saturating_add(values_size);
        let off_names = off_slots.saturating_add((n_slots as u64).saturating_mul(4));
        Self {
            magic: *b"KBIV0004",
            version: 4,
            endian: ENDIAN_TAG,
            n_entries: n_entries as u64,
            index_len: index_len as u64,
            off_index,
            off_keys,
            off_offsets,
            off_slots,
            n_slots: n_slots as u64,
            off_names,
        }
    }

    pub(crate) fn new_compact(n_entries: usize, index_len: usize) -> Result<Self> {
        let n_entries = u64::try_from(n_entries)
            .map_err(|_| KbiError::InvalidFormat("Too many entries".into()))?;
        let index_len = u64::try_from(index_len)
            .map_err(|_| KbiError::InvalidFormat("Key index too large".into()))?;
        let off_index = mem::size_of::<Self>() as u64;
        let off_offsets = off_index
            .checked_add(index_len)
            .ok_or_else(|| KbiError::InvalidFormat("Key index size overflow".into()))?;
        let offsets_size = n_entries
            .checked_mul(8)
            .ok_or_else(|| KbiError::InvalidFormat("Offsets size overflow".into()))?;
        let off_names = off_offsets
            .checked_add(offsets_size)
            .ok_or_else(|| KbiError::InvalidFormat("Offsets size overflow".into()))?;
        let header = Self {
            magic: KBI_MAGIC,
            version: KBI_VERSION,
            endian: ENDIAN_TAG,
            n_entries,
            index_len,
            off_index,
            off_keys: off_offsets,
            off_offsets,
            off_slots: off_names,
            n_slots: 0,
            off_names,
        };
        header.validate()?;
        Ok(header)
    }

    pub fn validate(&self) -> Result<()> {
        let legacy = self.magic == *b"KBIV0004";
        if !legacy && self.magic != KBI_MAGIC {
            return Err(KbiError::InvalidFormat("Invalid magic bytes".into()));
        }
        let expected = if legacy { 4 } else { KBI_VERSION };
        if self.version != expected {
            return Err(KbiError::VersionMismatch {
                expected,
                got: self.version,
            });
        }
        if self.endian != ENDIAN_TAG {
            return Err(KbiError::InvalidFormat("Endianness mismatch".into()));
        }
        if self.n_entries == 0 || self.n_entries > u32::MAX as u64 || self.index_len == 0 {
            return Err(KbiError::InvalidFormat(
                "Invalid entry or index count".into(),
            ));
        }
        if self.off_index != mem::size_of::<Self>() as u64 {
            return Err(KbiError::InvalidFormat("Invalid index offset".into()));
        }
        let index_end = self
            .off_index
            .checked_add(self.index_len)
            .ok_or_else(|| KbiError::InvalidFormat("Index size overflow".into()))?;
        let values_size = self
            .n_entries
            .checked_mul(8)
            .ok_or_else(|| KbiError::InvalidFormat("Entry size overflow".into()))?;
        let expected_offsets = if legacy {
            index_end
                .checked_add(values_size)
                .ok_or_else(|| KbiError::InvalidFormat("Keys size overflow".into()))?
        } else {
            index_end
        };
        let expected_slots = expected_offsets
            .checked_add(values_size)
            .ok_or_else(|| KbiError::InvalidFormat("Offsets size overflow".into()))?;
        let slots_size = if legacy {
            if self.n_slots < self.n_entries || self.n_slots > u32::MAX as u64 {
                return Err(KbiError::InvalidFormat("Invalid slot count".into()));
            }
            self.n_slots
                .checked_mul(4)
                .ok_or_else(|| KbiError::InvalidFormat("Slots size overflow".into()))?
        } else {
            if self.n_slots != 0 {
                return Err(KbiError::InvalidFormat("Unexpected slot table".into()));
            }
            0
        };
        let expected_names = expected_slots
            .checked_add(slots_size)
            .ok_or_else(|| KbiError::InvalidFormat("Slots size overflow".into()))?;
        if self.off_keys != index_end
            || self.off_offsets != expected_offsets
            || self.off_slots != expected_slots
            || self.off_names != expected_names
        {
            return Err(KbiError::InvalidFormat(
                "Inconsistent section offsets".into(),
            ));
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
