use std::io;
use thiserror::Error;

use noodles_bgzf as bgzf;

#[derive(Debug, Error)]
pub enum BgzfError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
    #[error("Invalid BGZF format")]
    InvalidFormat,
}

pub type Result<T> = std::result::Result<T, BgzfError>;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct VirtualPosition(u64);

impl VirtualPosition {
    #[inline]
    pub fn new(compressed_offset: u64, uncompressed_offset: u16) -> Self {
        Self((compressed_offset << 16) | (uncompressed_offset as u64))
    }

    #[inline]
    pub fn compressed(&self) -> u64 {
        self.0 >> 16
    }

    #[inline]
    pub fn uncompressed(&self) -> u16 {
        (self.0 & 0xFFFF) as u16
    }

    #[inline]
    pub fn as_u64(&self) -> u64 {
        self.0
    }

    #[inline]
    pub fn from_u64(v: u64) -> Self {
        Self(v)
    }
}

impl From<bgzf::VirtualPosition> for VirtualPosition {
    fn from(vp: bgzf::VirtualPosition) -> Self {
        Self(u64::from(vp))
    }
}

impl From<VirtualPosition> for bgzf::VirtualPosition {
    fn from(vp: VirtualPosition) -> Self {
        bgzf::VirtualPosition::try_from(vp.0).unwrap_or_default()
    }
}
