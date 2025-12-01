use std::fs::File;
use std::io::{BufWriter, Write};
use std::mem;
use std::path::Path;
use std::slice;

use bytemuck::{Pod, Zeroable};
use kira_kv_engine::{BuildConfig, Builder, MphError, Mphf};
use memmap2::{Mmap, MmapOptions};
use rayon::prelude::*;
use thiserror::Error;

use crate::util::{ChrId, GenomicKey};
use crate::vcf::{VcfReader, VcfRecord};

pub const KBI_MAGIC: [u8; 8] = *b"KBIV0002";
pub const KBI_VERSION: u32 = 2;
const ENDIAN_TAG: u32 = 0x01020304;

#[derive(Debug, Error)]
pub enum KbiError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
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
    fn new(n_entries: usize, mph: &Mphf) -> Self {
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

    fn validate(&self) -> Result<()> {
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

pub struct KbiIndex {
    mph: Mphf,
    keys: Vec<u64>,
    offsets: Vec<u64>,
}

impl KbiIndex {
    #[inline]
    pub fn get(&self, key: GenomicKey) -> Option<u64> {
        let key_bytes = key.as_u64().to_le_bytes();
        let idx = self.mph.index(&key_bytes) as usize;

        if idx < self.keys.len() && self.keys[idx] == key.as_u64() {
            Some(self.offsets[idx])
        } else {
            None
        }
    }

    #[inline]
    pub fn contains(&self, key: GenomicKey) -> bool {
        self.get(key).is_some()
    }

    pub fn len(&self) -> usize {
        self.keys.len()
    }

    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    pub fn get_batch(&self, keys: &[GenomicKey]) -> Vec<Option<u64>> {
        keys.par_iter().map(|&k| self.get(k)).collect()
    }

    pub fn range(&self, chr: ChrId, start_pos: u32, end_pos: u32) -> Vec<(u32, u64)> {
        let start_key = GenomicKey::new(chr, start_pos).as_u64();
        let end_key = GenomicKey::new(chr, end_pos).as_u64();

        let start_idx = self.keys.partition_point(|&k| k < start_key);
        let end_idx = self.keys.partition_point(|&k| k <= end_key);

        self.keys[start_idx..end_idx]
            .iter()
            .zip(&self.offsets[start_idx..end_idx])
            .map(|(&k, &off)| (GenomicKey::from_u64(k).position(), off))
            .collect()
    }

    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::with_capacity(8 * 1024 * 1024, file);

        let header = KbiHeader::new(self.keys.len(), &self.mph);

        writer.write_all(bytemuck::bytes_of(&header))?;

        let g_bytes = unsafe {
            slice::from_raw_parts(
                self.mph.g.as_ptr() as *const u8,
                self.mph.g.len() * mem::size_of::<u32>(),
            )
        };
        writer.write_all(g_bytes)?;

        let keys_bytes = unsafe {
            slice::from_raw_parts(
                self.keys.as_ptr() as *const u8,
                self.keys.len() * mem::size_of::<u64>(),
            )
        };
        writer.write_all(keys_bytes)?;

        let offsets_bytes = unsafe {
            slice::from_raw_parts(
                self.offsets.as_ptr() as *const u8,
                self.offsets.len() * mem::size_of::<u64>(),
            )
        };
        writer.write_all(offsets_bytes)?;

        writer.flush()?;
        Ok(())
    }

    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let mmap = unsafe {
            let file = File::open(path)?;
            MmapOptions::new().map(&file)?
        };

        Self::from_mmap(&mmap)
    }

    pub fn from_mmap(mmap: &Mmap) -> Result<Self> {
        if mmap.len() < mem::size_of::<KbiHeader>() {
            return Err(KbiError::InvalidFormat("File too small".into()));
        }

        let header: &KbiHeader = bytemuck::from_bytes(&mmap[..mem::size_of::<KbiHeader>()]);
        header.validate()?;

        let n = header.n_entries as usize;
        let m = header.mph_m as usize;

        let g_start = header.off_mph_g as usize;
        let g_end = g_start + m * mem::size_of::<u32>();
        let g: Vec<u32> = bytemuck::cast_slice(&mmap[g_start..g_end]).to_vec();

        let mph = Mphf {
            n: header.n_entries,
            m: header.mph_m,
            salt: header.mph_salt as u64,
            g,
        };

        let keys_start = header.off_keys as usize;
        let keys_end = keys_start + n * mem::size_of::<u64>();
        let keys: Vec<u64> = bytemuck::cast_slice(&mmap[keys_start..keys_end]).to_vec();

        let offsets_start = header.off_offsets as usize;
        let offsets_end = offsets_start + n * mem::size_of::<u64>();
        let offsets: Vec<u64> = bytemuck::cast_slice(&mmap[offsets_start..offsets_end]).to_vec();

        Ok(Self { mph, keys, offsets })
    }

    pub fn memory_usage(&self) -> usize {
        mem::size_of::<Self>()
            + self.keys.len() * mem::size_of::<u64>()
            + self.offsets.len() * mem::size_of::<u64>()
            + self.mph.g.len() * mem::size_of::<u32>()
    }

    pub fn bytes_per_key(&self) -> f64 {
        self.memory_usage() as f64 / self.len().max(1) as f64
    }

    fn from_parts(mph: Mphf, keys: Vec<u64>, offsets: Vec<u64>) -> Self {
        Self { mph, keys, offsets }
    }
}

pub struct KbiBuilder {
    entries: Vec<(u64, u64)>,
    config: BuildConfig,
}

impl KbiBuilder {
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
            config: BuildConfig {
                gamma: 1.27,
                rehash_limit: 16,
                salt: 0xC0FF_EE00_D15E_A5E,
            },
        }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            entries: Vec::with_capacity(capacity),
            ..Self::new()
        }
    }

    pub fn gamma(mut self, gamma: f64) -> Self {
        self.config.gamma = gamma;
        self
    }

    #[inline]
    pub fn add(&mut self, key: GenomicKey, offset: u64) {
        self.entries.push((key.as_u64(), offset));
    }

    #[inline]
    pub fn add_record(&mut self, record: &VcfRecord) {
        self.add(record.key(), record.offset);
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn build(mut self) -> Result<KbiIndex> {
        if self.entries.is_empty() {
            return Err(KbiError::EmptyDataset);
        }

        self.entries.par_sort_unstable_by_key(|(k, _)| *k);

        let mut deduped = Vec::with_capacity(self.entries.len());
        let mut last_key: Option<u64> = None;
        let mut dup_count = 0usize;

        for (key, offset) in self.entries {
            if Some(key) == last_key {
                dup_count += 1;
                continue;
            }
            deduped.push((key, offset));
            last_key = Some(key);
        }

        if dup_count > 0 {
            eprintln!("Deduplicated {} entries", dup_count);
        }

        let keys: Vec<u64> = deduped.iter().map(|(k, _)| *k).collect();
        let offsets: Vec<u64> = deduped.iter().map(|(_, o)| *o).collect();

        let key_bytes: Vec<[u8; 8]> = keys.iter().map(|k| k.to_le_bytes()).collect();

        let mph = Builder::new()
            .with_config(self.config)
            .build(key_bytes.iter().map(|b| b.as_slice()))?;

        let n = keys.len();
        let sorted_keys = vec![0u64; n];
        let sorted_offsets = vec![0u64; n];

        keys.par_iter()
            .zip(offsets.par_iter())
            .for_each(|(&key, &offset)| {
                let idx = mph.index(&key.to_le_bytes()) as usize;
                unsafe {
                    let keys_ptr = sorted_keys.as_ptr() as *mut u64;
                    let offsets_ptr = sorted_offsets.as_ptr() as *mut u64;
                    *keys_ptr.add(idx) = key;
                    *offsets_ptr.add(idx) = offset;
                }
            });

        Ok(KbiIndex::from_parts(mph, sorted_keys, sorted_offsets))
    }
}

impl Default for KbiBuilder {
    fn default() -> Self {
        Self::new()
    }
}

pub fn build_kbi_index<P: AsRef<Path>>(vcf_path: P, output_path: P) -> Result<KbiIndex> {
    let mut reader = VcfReader::open(vcf_path)?;
    reader.header()?;

    let estimated_capacity = 10_000_000;
    let mut builder = KbiBuilder::with_capacity(estimated_capacity);

    let mut count = 0usize;
    for record in reader.records() {
        let record = record?;
        builder.add_record(&record);
        count += 1;

        if count % 1_000_000 == 0 {
            eprintln!("Processed {} records...", count);
        }
    }

    eprintln!("Building MPH for {} entries...", builder.len());
    let index = builder.build()?;

    eprintln!("Saving index to {:?}...", output_path.as_ref());
    index.save(&output_path)?;

    Ok(index)
}

pub struct KbiStats {
    pub entries: usize,
    pub memory_bytes: usize,
    pub bytes_per_key: f64,
    pub file_size: u64,
}

impl KbiStats {
    pub fn from_index(index: &KbiIndex, file_size: u64) -> Self {
        Self {
            entries: index.len(),
            memory_bytes: index.memory_usage(),
            bytes_per_key: index.bytes_per_key(),
            file_size,
        }
    }
}
