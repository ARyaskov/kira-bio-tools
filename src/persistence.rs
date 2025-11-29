//! Persistence layer for `VcfIndex`.
//!
//! - Cross-platform (Linux / macOS / Windows)
//! - No alignment assumptions for mmap on load
//! - On save uses BufWriter + bytemuck::cast_slice for fast bulk I/O
//! - On disk layout is:
//!   [ IndexHeader ][ MPH.g (u32[m]) ][ keys (u64[n]) ][ offsets (u64[n]) ]

use std::{
    fs::File,
    io::{BufWriter, Read, Seek, SeekFrom, Write},
    path::Path,
};

use kira_kv_engine::Mphf;
use memmap2::MmapOptions;

use crate::vcf_index::{IndexError, Result, VcfIndex};

/// Magic of `.kbi` files: 'KBI1'
pub const INDEX_MAGIC: u32 = 0x4B_42_49_31;

/// Index format version.
pub const INDEX_VERSION: u32 = 1;

/// Fixed-size on-disk header for `.kbi` files.
///
/// Layout (little-endian):
/// - magic:      u32  (4 bytes)
/// - version:    u32  (4 bytes)
/// - n_entries:  u64  (8 bytes)
/// - mph_n:      u64  (8 bytes)
/// - mph_m:      u32  (4 bytes)
/// - mph_salt:   u64  (8 bytes)
/// - off_mph_g:  u64  (8 bytes)
/// - off_keys:   u64  (8 bytes)
/// - off_offsets:u64  (8 bytes)
/// = 64 bytes total
#[derive(Clone, Copy, Debug)]
pub struct IndexHeader {
    pub magic: u32,
    pub version: u32,
    pub n_entries: u64,
    pub mph_n: u64,
    pub mph_m: u32,
    pub mph_salt: u64,
    pub off_mph_g: u64,
    pub off_keys: u64,
    pub off_offsets: u64,
}

impl IndexHeader {
    pub const ENCODED_LEN: usize = 64;

    pub fn new(
        n_entries: usize,
        mph_n: u64,
        mph_m: u32,
        mph_salt: u64,
        off_mph_g: usize,
        off_keys: usize,
        off_offsets: usize,
    ) -> Self {
        Self {
            magic: INDEX_MAGIC,
            version: INDEX_VERSION,
            n_entries: n_entries as u64,
            mph_n,
            mph_m,
            mph_salt,
            off_mph_g: off_mph_g as u64,
            off_keys: off_keys as u64,
            off_offsets: off_offsets as u64,
        }
    }

    pub fn validate(&self) -> Result<()> {
        if self.magic != INDEX_MAGIC {
            return Err(IndexError::InvalidFormat("Bad magic".to_string()));
        }
        if self.version != INDEX_VERSION {
            return Err(IndexError::VersionMismatch {
                expected: INDEX_VERSION,
                got: self.version,
            });
        }
        Ok(())
    }

    pub fn to_bytes(&self) -> [u8; Self::ENCODED_LEN] {
        let mut buf = [0u8; Self::ENCODED_LEN];
        let mut off = 0;

        buf[off..off + 4].copy_from_slice(&self.magic.to_le_bytes());
        off += 4;
        buf[off..off + 4].copy_from_slice(&self.version.to_le_bytes());
        off += 4;

        buf[off..off + 8].copy_from_slice(&self.n_entries.to_le_bytes());
        off += 8;
        buf[off..off + 8].copy_from_slice(&self.mph_n.to_le_bytes());
        off += 8;
        buf[off..off + 4].copy_from_slice(&self.mph_m.to_le_bytes());
        off += 4;
        buf[off..off + 8].copy_from_slice(&self.mph_salt.to_le_bytes());
        off += 8;
        buf[off..off + 8].copy_from_slice(&self.off_mph_g.to_le_bytes());
        off += 8;
        buf[off..off + 8].copy_from_slice(&self.off_keys.to_le_bytes());
        off += 8;
        buf[off..off + 8].copy_from_slice(&self.off_offsets.to_le_bytes());

        buf
    }

    pub fn from_bytes(buf: &[u8; Self::ENCODED_LEN]) -> Self {
        let mut off = 0;

        let magic = u32::from_le_bytes(buf[off..off + 4].try_into().unwrap());
        off += 4;
        let version = u32::from_le_bytes(buf[off..off + 4].try_into().unwrap());
        off += 4;

        let n_entries = u64::from_le_bytes(buf[off..off + 8].try_into().unwrap());
        off += 8;
        let mph_n = u64::from_le_bytes(buf[off..off + 8].try_into().unwrap());
        off += 8;
        let mph_m = u32::from_le_bytes(buf[off..off + 4].try_into().unwrap());
        off += 4;
        let mph_salt = u64::from_le_bytes(buf[off..off + 8].try_into().unwrap());
        off += 8;
        let off_mph_g = u64::from_le_bytes(buf[off..off + 8].try_into().unwrap());
        off += 8;
        let off_keys = u64::from_le_bytes(buf[off..off + 8].try_into().unwrap());
        off += 8;
        let off_offsets = u64::from_le_bytes(buf[off..off + 8].try_into().unwrap());

        Self {
            magic,
            version,
            n_entries,
            mph_n,
            mph_m,
            mph_salt,
            off_mph_g,
            off_keys,
            off_offsets,
        }
    }
}

/// Save index to `.kbi` file.
///
/// Layout:
/// - header (`IndexHeader`)
/// - MPH.g: `u32[m]` in little-endian
/// - keys: `u64[n]` in little-endian
/// - offsets: `u64[n]` in little-endian
pub fn save_index(path: &Path, keys: &[u64], offsets: &[u64], mph: &Mphf) -> Result<()> {
    use bytemuck::cast_slice;

    let n = keys.len();
    if n != offsets.len() {
        return Err(IndexError::InvalidFormat(
            "keys/offsets length mismatch".to_string(),
        ));
    }

    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    let hdr_size = IndexHeader::ENCODED_LEN;
    let off_mph_g = hdr_size;
    let off_keys = off_mph_g + (mph.m as usize) * std::mem::size_of::<u32>();
    let off_offsets = off_keys + n * std::mem::size_of::<u64>();

    // Header
    let header = IndexHeader::new(
        n,
        mph.n,
        mph.m,
        mph.salt,
        off_mph_g,
        off_keys,
        off_offsets,
    );
    writer.write_all(&header.to_bytes())?;

    // MPH.g (bulk write)
    writer.write_all(cast_slice::<u32, u8>(&mph.g))?;

    // Keys (bulk write)
    writer.write_all(cast_slice::<u64, u8>(keys))?;

    // Offsets (bulk write)
    writer.write_all(cast_slice::<u64, u8>(offsets))?;

    writer.flush()?;
    Ok(())
}

/// Load index using memory-mapped I/O.
///
/// Windows does not guarantee 8-byte alignment for mmap base address,
/// so we treat the mapped region as raw bytes and decode with
/// `from_le_bytes` instead of `cast_slice`.
pub fn load_index_mmap(path: &Path) -> Result<VcfIndex> {
    let file = File::open(path)?;
    let mmap = unsafe { MmapOptions::new().map(&file)? };

    if mmap.len() < IndexHeader::ENCODED_LEN {
        return Err(IndexError::InvalidFormat("File too small".to_string()));
    }

    let mut hdr_buf = [0u8; IndexHeader::ENCODED_LEN];
    hdr_buf.copy_from_slice(&mmap[..IndexHeader::ENCODED_LEN]);
    let header = IndexHeader::from_bytes(&hdr_buf);
    header.validate()?;

    let n = header.n_entries as usize;
    let m = header.mph_m as usize;

    // MPH.g
    let g_start = header.off_mph_g as usize;
    let g_end = g_start + m * std::mem::size_of::<u32>();
    if g_end > mmap.len() {
        return Err(IndexError::InvalidFormat("Truncated MPH.g".to_string()));
    }
    let mut g = Vec::with_capacity(m);
    for chunk in mmap[g_start..g_end].chunks_exact(4) {
        g.push(u32::from_le_bytes(chunk.try_into().unwrap()));
    }

    let mph = Mphf {
        n: header.mph_n,
        m: header.mph_m,
        salt: header.mph_salt,
        g,
    };

    // Keys
    let keys_start = header.off_keys as usize;
    let keys_end = keys_start + n * std::mem::size_of::<u64>();
    if keys_end > mmap.len() {
        return Err(IndexError::InvalidFormat("Truncated keys".to_string()));
    }
    let mut keys = Vec::with_capacity(n);
    for chunk in mmap[keys_start..keys_end].chunks_exact(8) {
        keys.push(u64::from_le_bytes(chunk.try_into().unwrap()));
    }

    // Offsets
    let offs_start = header.off_offsets as usize;
    let offs_end = offs_start + n * std::mem::size_of::<u64>();
    if offs_end > mmap.len() {
        return Err(IndexError::InvalidFormat("Truncated offsets".to_string()));
    }
    let mut offsets = Vec::with_capacity(n);
    for chunk in mmap[offs_start..offs_end].chunks_exact(8) {
        offsets.push(u64::from_le_bytes(chunk.try_into().unwrap()));
    }

    Ok(VcfIndex::from_parts(mph, keys, offsets))
}

/// Load index with regular buffered I/O (без mmap).
pub fn load_index(path: &Path) -> Result<VcfIndex> {
    let mut file = File::open(path)?;

    // Header
    let mut hdr_buf = [0u8; IndexHeader::ENCODED_LEN];
    file.read_exact(&mut hdr_buf)?;
    let header = IndexHeader::from_bytes(&hdr_buf);
    header.validate()?;

    let n = header.n_entries as usize;
    let m = header.mph_m as usize;

    // MPH.g
    file.seek(SeekFrom::Start(header.off_mph_g))?;
    let mut g = Vec::with_capacity(m);
    let mut buf4 = [0u8; 4];
    for _ in 0..m {
        file.read_exact(&mut buf4)?;
        g.push(u32::from_le_bytes(buf4));
    }

    let mph = Mphf {
        n: header.mph_n,
        m: header.mph_m,
        salt: header.mph_salt,
        g,
    };

    // Keys
    file.seek(SeekFrom::Start(header.off_keys))?;
    let mut keys = Vec::with_capacity(n);
    let mut buf8 = [0u8; 8];
    for _ in 0..n {
        file.read_exact(&mut buf8)?;
        keys.push(u64::from_le_bytes(buf8));
    }

    // Offsets
    file.seek(SeekFrom::Start(header.off_offsets))?;
    let mut offsets = Vec::with_capacity(n);
    for _ in 0..n {
        file.read_exact(&mut buf8)?;
        offsets.push(u64::from_le_bytes(buf8));
    }

    Ok(VcfIndex::from_parts(mph, keys, offsets))
}
