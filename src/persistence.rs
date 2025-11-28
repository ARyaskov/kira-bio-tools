//! Binary persistence for VcfIndex with mmap support.
//!
//! File format (.kbi - Kira Bio Index):
//! ```text
//! [IndexHeader]           - 64 bytes
//! [u32 mph_g[]; mph_m]    - MPH displacement table
//! [u64 keys[]; n]         - Sorted genomic keys
//! [u64 offsets[]; n]      - VCF byte offsets
//! ```

use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::mem;
use std::path::Path;
use std::slice;

use bytemuck::{Pod, Zeroable};
use kira_kv_engine::Mphf;
use memmap2::{Mmap, MmapOptions};

use crate::vcf_index::{IndexError, Result, VcfIndex};

/// Magic bytes identifying KBI format
pub const INDEX_MAGIC: [u8; 8] = *b"KBIV0001";

/// Current format version
pub const INDEX_VERSION: u32 = 1;

/// Endianness marker
const ENDIAN_TAG: u32 = 0x01020304;

/// Index file header (64 bytes, cache-line aligned)
///
/// Layout (no padding):
///   magic      : [u8; 8]   ( 8)
///   n_entries  : u64       ( 8) = 16
///   mph_salt   : u64       ( 8) = 24
///   off_mph_g  : u64       ( 8) = 32
///   off_keys   : u64       ( 8) = 40
///   off_offsets: u64       ( 8) = 48
///   version    : u32       ( 4) = 52
///   endian     : u32       ( 4) = 56
///   mph_m      : u32       ( 4) = 60
///   _reserved  : [u8; 4]   ( 4) = 64
#[repr(C, align(64))]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
pub struct IndexHeader {
    pub magic: [u8; 8],

    /// Number of indexed entries (also used as MPH `n`)
    pub n_entries: u64,

    /// MPH salt
    pub mph_salt: u64,

    /// Section offsets from file start
    pub off_mph_g: u64,
    pub off_keys: u64,
    pub off_offsets: u64,

    /// File format version
    pub version: u32,

    /// Endianness marker
    pub endian: u32,

    /// MPH `m` parameter (length of displacement table `g`)
    pub mph_m: u32,

    /// Reserved for future use
    pub _reserved: [u8; 4],
}

impl IndexHeader {
    fn new(n_entries: usize, mph: &Mphf) -> Self {
        let header_size = mem::size_of::<IndexHeader>() as u64;
        let mph_g_size = (mph.g.len() * mem::size_of::<u32>()) as u64;
        let keys_size = (n_entries * mem::size_of::<u64>()) as u64;

        Self {
            magic: INDEX_MAGIC,
            n_entries: n_entries as u64,
            mph_salt: mph.salt,
            off_mph_g: header_size,
            off_keys: header_size + mph_g_size,
            off_offsets: header_size + mph_g_size + keys_size,
            version: INDEX_VERSION,
            endian: ENDIAN_TAG,
            mph_m: mph.m,
            _reserved: [0; 4],
        }
    }

    fn validate(&self) -> Result<()> {
        if self.magic != INDEX_MAGIC {
            return Err(IndexError::InvalidFormat("Bad magic bytes".into()));
        }
        if self.version != INDEX_VERSION {
            return Err(IndexError::VersionMismatch {
                expected: INDEX_VERSION,
                got: self.version,
            });
        }
        if self.endian != ENDIAN_TAG {
            return Err(IndexError::InvalidFormat("Endianness mismatch".into()));
        }
        Ok(())
    }
}

/// Save index to file
pub fn save_index(path: &Path, keys: &[u64], offsets: &[u64], mph: &Mphf) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::with_capacity(8 * 1024 * 1024, file);

    let header = IndexHeader::new(keys.len(), mph);

    // Write header
    writer.write_all(bytemuck::bytes_of(&header))?;

    // Write MPH g array
    let g_bytes = unsafe {
        slice::from_raw_parts(
            mph.g.as_ptr() as *const u8,
            mph.g.len() * mem::size_of::<u32>(),
        )
    };
    writer.write_all(g_bytes)?;

    // Write keys
    let keys_bytes = unsafe {
        slice::from_raw_parts(
            keys.as_ptr() as *const u8,
            keys.len() * mem::size_of::<u64>(),
        )
    };
    writer.write_all(keys_bytes)?;

    // Write offsets
    let offsets_bytes = unsafe {
        slice::from_raw_parts(
            offsets.as_ptr() as *const u8,
            offsets.len() * mem::size_of::<u64>(),
        )
    };
    writer.write_all(offsets_bytes)?;

    writer.flush()?;
    Ok(())
}

/// Load index using mmap (zero-copy, instant load)
pub fn load_index_mmap(path: &Path) -> Result<VcfIndex> {
    let file = File::open(path)?;
    let mmap = unsafe { MmapOptions::new().map(&file)? };

    if mmap.len() < mem::size_of::<IndexHeader>() {
        return Err(IndexError::InvalidFormat("File too small".into()));
    }

    // Parse header
    let header: &IndexHeader = bytemuck::from_bytes(&mmap[..mem::size_of::<IndexHeader>()]);
    header.validate()?;

    let n = header.n_entries as usize;
    let m = header.mph_m as usize;

    // Extract MPH g array (must copy for Mphf ownership)
    let g_start = header.off_mph_g as usize;
    let g_end = g_start + m * mem::size_of::<u32>();
    let g: Vec<u32> = bytemuck::cast_slice(&mmap[g_start..g_end]).to_vec();

    let mph = Mphf {
        n: header.n_entries,
        m: header.mph_m,
        salt: header.mph_salt,
        g,
    };

    // Extract keys (copy)
    let keys_start = header.off_keys as usize;
    let keys_end = keys_start + n * mem::size_of::<u64>();
    let keys: Vec<u64> = bytemuck::cast_slice(&mmap[keys_start..keys_end]).to_vec();

    // Extract offsets (copy)
    let offsets_start = header.off_offsets as usize;
    let offsets_end = offsets_start + n * mem::size_of::<u64>();
    let offsets: Vec<u64> = bytemuck::cast_slice(&mmap[offsets_start..offsets_end]).to_vec();

    Ok(VcfIndex::from_parts(mph, keys, offsets))
}

/// Load index with full read (no mmap)
pub fn load_index(path: &Path) -> Result<VcfIndex> {
    let mut file = File::open(path)?;

    // Read header
    let mut header_bytes = [0u8; mem::size_of::<IndexHeader>()];
    file.read_exact(&mut header_bytes)?;
    let header: IndexHeader = *bytemuck::from_bytes(&header_bytes);
    header.validate()?;

    let n = header.n_entries as usize;
    let m = header.mph_m as usize;

    // Read MPH g
    file.seek(SeekFrom::Start(header.off_mph_g))?;
    let mut g = vec![0u32; m];
    let g_bytes =
        unsafe { slice::from_raw_parts_mut(g.as_mut_ptr() as *mut u8, m * mem::size_of::<u32>()) };
    file.read_exact(g_bytes)?;

    let mph = Mphf {
        n: header.n_entries,
        m: header.mph_m,
        salt: header.mph_salt,
        g,
    };

    // Read keys
    file.seek(SeekFrom::Start(header.off_keys))?;
    let mut keys = vec![0u64; n];
    let keys_bytes = unsafe {
        slice::from_raw_parts_mut(keys.as_mut_ptr() as *mut u8, n * mem::size_of::<u64>())
    };
    file.read_exact(keys_bytes)?;

    // Read offsets
    file.seek(SeekFrom::Start(header.off_offsets))?;
    let mut offsets = vec![0u64; n];
    let offsets_bytes = unsafe {
        slice::from_raw_parts_mut(offsets.as_mut_ptr() as *mut u8, n * mem::size_of::<u64>())
    };
    file.read_exact(offsets_bytes)?;

    Ok(VcfIndex::from_parts(mph, keys, offsets))
}

/// Calculate file size for given parameters
pub fn calculate_file_size(n_entries: usize, mph_m: usize) -> u64 {
    let header = mem::size_of::<IndexHeader>() as u64;
    let mph_g = (mph_m * mem::size_of::<u32>()) as u64;
    let keys = (n_entries * mem::size_of::<u64>()) as u64;
    let offsets = (n_entries * mem::size_of::<u64>()) as u64;
    header + mph_g + keys + offsets
}
