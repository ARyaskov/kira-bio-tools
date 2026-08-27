use std::fmt;
use std::io;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct VirtualPosition {
    pub block_offset: u64,
    pub data_offset: u16,
}

impl VirtualPosition {
    pub fn new(block_offset: u64, data_offset: u16) -> Self {
        Self {
            block_offset,
            data_offset,
        }
    }

    pub fn from_raw(vpos: u64) -> Self {
        Self {
            block_offset: vpos >> 16,
            data_offset: (vpos & 0xFFFF) as u16,
        }
    }

    pub fn from_u64(vpos: u64) -> Self {
        Self::from_raw(vpos)
    }

    pub fn to_raw(&self) -> u64 {
        (self.block_offset << 16) | (self.data_offset as u64)
    }

    pub fn as_u64(&self) -> u64 {
        self.to_raw()
    }
}

impl From<noodles_bgzf::VirtualPosition> for VirtualPosition {
    fn from(pos: noodles_bgzf::VirtualPosition) -> Self {
        let raw: u64 = pos.into();
        Self::from_raw(raw)
    }
}

impl From<VirtualPosition> for noodles_bgzf::VirtualPosition {
    fn from(pos: VirtualPosition) -> Self {
        noodles_bgzf::VirtualPosition::from(pos.to_raw())
    }
}

#[derive(Debug)]
pub enum BgzfError {
    Io(io::Error),
    InvalidHeader,
    InvalidMagic,
    InvalidSignature,
    BlockTooSmall,
    IncompleteBlock,
    CompressionFailed,
}

impl fmt::Display for BgzfError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BgzfError::Io(e) => write!(f, "IO error: {}", e),
            BgzfError::InvalidHeader => write!(f, "Invalid BGZF header"),
            BgzfError::InvalidMagic => write!(f, "Invalid BGZF magic number"),
            BgzfError::InvalidSignature => write!(f, "Invalid BGZF signature"),
            BgzfError::BlockTooSmall => write!(f, "BGZF block too small"),
            BgzfError::IncompleteBlock => write!(f, "Incomplete BGZF block"),
            BgzfError::CompressionFailed => write!(f, "BGZF compression failed"),
        }
    }
}

impl std::error::Error for BgzfError {}

impl From<io::Error> for BgzfError {
    fn from(error: io::Error) -> Self {
        BgzfError::Io(error)
    }
}

pub type Result<T> = std::result::Result<T, BgzfError>;

#[derive(Clone)]
pub struct BgzfBlock {
    pub compressed: Vec<u8>,
    pub uncompressed: Vec<u8>,
    pub cdata_size: u16,
    pub file_offset: u64,
}

impl BgzfBlock {
    pub fn new(file_offset: u64) -> Self {
        const BGZF_BLOCK_SIZE: usize = 64 * 1024;
        Self {
            compressed: Vec::with_capacity(BGZF_BLOCK_SIZE),
            uncompressed: Vec::with_capacity(BGZF_BLOCK_SIZE),
            cdata_size: 0,
            file_offset,
        }
    }

    pub fn virtual_offset(&self) -> u64 {
        self.file_offset << 16
    }
}

pub struct WritePool {
    pool: std::sync::Arc<std::sync::Mutex<Vec<Vec<u8>>>>,
}

impl WritePool {
    pub fn new(n: usize, capacity: usize) -> Self {
        let mut v = Vec::with_capacity(n);
        for _ in 0..n {
            let mut buf = Vec::with_capacity(capacity);
            unsafe {
                buf.set_len(capacity);
            }
            v.push(buf);
        }
        Self {
            pool: std::sync::Arc::new(std::sync::Mutex::new(v)),
        }
    }

    #[inline]
    pub fn get(&self, capacity: usize) -> Vec<u8> {
        self.pool.lock().unwrap().pop().unwrap_or_else(|| {
            let mut b = Vec::with_capacity(capacity);
            unsafe {
                b.set_len(capacity);
            }
            b
        })
    }

    #[inline]
    pub fn put(&self, mut buf: Vec<u8>, capacity: usize) {
        unsafe {
            buf.set_len(capacity);
        }
        self.pool.lock().unwrap().push(buf);
    }
}

impl Clone for WritePool {
    fn clone(&self) -> Self {
        Self {
            pool: self.pool.clone(),
        }
    }
}

pub struct CompressedBlock {
    pub data: Vec<u8>,
    pub sequence: usize,
}

pub const BGZF_HEADER: [u8; 18] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x00, 0x00,
];

pub const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

pub const BGZF_MAGIC: [u8; 4] = [0x1f, 0x8b, 0x08, 0x04];
pub const BGZF_EXTRA_LEN: u16 = 6;
pub const BGZF_SI1: u8 = 0x42;
pub const BGZF_SI2: u8 = 0x43;
pub const BGZF_BLOCK_SIZE: usize = 64 * 1024;
pub const CHUNK_SIZE: usize = 64 * 1024 - 256;
