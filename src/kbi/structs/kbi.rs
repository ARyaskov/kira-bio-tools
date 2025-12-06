use std::io;
use std::mem;

use bytemuck::{Pod, Zeroable};
use kira_kv_engine::{MphError, Mphf};
use thiserror::Error;

pub const KBI_MAGIC: [u8; 8] = *b"KBIV0002";
pub const KBI_VERSION: u32 = 2;
pub const ENDIAN_TAG: u32 = 0x01020304;

#[derive(Debug, Error)]
pub enum KbiError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
    #[error("MPH error: {0}")]
    Mph(#[from] MphError),
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
    pub mph_m: u32,
    pub mph_salt: u32,
    pub off_mph_g: u64,
    pub off_keys: u64,
    pub off_offsets: u64,
    _reserved: [u8; 8],
}

const _: () = assert!(mem::size_of::<KbiHeader>() == 64);

impl KbiHeader {
    pub fn new(n_entries: usize, mph: &Mphf) -> Self {
        let header_size = mem::size_of::<Self>() as u64;
        let mph_g_size = (mph.g.len() * mem::size_of::<u32>()) as u64;
        let keys_size = (n_entries * mem::size_of::<u64>()) as u64;

        Self {
            magic: KBI_MAGIC,
            version: KBI_VERSION,
            endian: ENDIAN_TAG,
            n_entries: n_entries as u64,
            mph_m: mph.m,
            mph_salt: mph.salt as u32,
            off_mph_g: header_size,
            off_keys: header_size + mph_g_size,
            off_offsets: header_size + mph_g_size + keys_size,
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
