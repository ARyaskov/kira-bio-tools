//! VCF Index implementation using kira_kv_engine hybrid index.

use kira_kv_engine::Builder as KvBuilder;
use kira_kv_engine::{BuildConfig, Builder, MphError, Mphf};
use std::path::Path;
use thiserror::Error;

use crate::persistence::{self, INDEX_VERSION, IndexHeader};

/// Chromosome identifier (1-25 for standard chromosomes)
pub type ChrId = u8;

/// Genomic coordinate encoded as single u64 for efficient indexing.
/// Format: (chr_id << 32) | position
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GenomicKey(u64);

impl GenomicKey {
    /// Create key from chromosome ID and position.
    #[inline]
    pub fn new(chr: ChrId, position: u32) -> Self {
        Self(((chr as u64) << 32) | (position as u64))
    }

    /// Extract chromosome ID.
    #[inline]
    pub fn chr(&self) -> ChrId {
        (self.0 >> 32) as ChrId
    }

    /// Extract position within chromosome.
    #[inline]
    pub fn position(&self) -> u32 {
        self.0 as u32
    }

    /// Raw u64 value for indexing.
    #[inline]
    pub fn as_u64(&self) -> u64 {
        self.0
    }

    /// Create from raw u64.
    #[inline]
    pub fn from_u64(v: u64) -> Self {
        Self(v)
    }
}

#[derive(Debug, Error)]
pub enum IndexError {
    #[error("Empty dataset")]
    EmptyDataset,

    #[error("Duplicate key: chr{0}:{1}")]
    DuplicateKey(ChrId, u32),

    #[error("MPH construction failed: {0}")]
    MphError(#[from] MphError),

    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),

    #[error("Invalid index format: {0}")]
    InvalidFormat(String),

    #[error("Version mismatch: expected {expected}, got {got}")]
    VersionMismatch { expected: u32, got: u32 },
}

pub type Result<T> = std::result::Result<T, IndexError>;

/// High-performance VCF index using minimal perfect hashing.
///
/// Provides O(1) lookup from genomic coordinates to VCF file byte offsets.
pub struct VcfIndex {
    /// MPH for key lookup
    mph: Mphf,
    /// Sorted keys (for verification and range queries)
    keys: Vec<u64>,
    /// VCF byte offsets corresponding to keys
    offsets: Vec<u64>,
}

impl VcfIndex {
    /// Look up byte offset for genomic position.
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

    /// Check if key exists in index.
    #[inline]
    pub fn contains(&self, key: GenomicKey) -> bool {
        self.get(key).is_some()
    }

    /// Number of indexed positions.
    #[inline]
    pub fn len(&self) -> usize {
        self.keys.len()
    }

    /// Check if index is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    /// Batch lookup for multiple keys.
    pub fn get_batch(&self, keys: &[GenomicKey]) -> Vec<Option<u64>> {
        keys.iter().map(|&k| self.get(k)).collect()
    }

    /// Range query: find all offsets for positions in [start, end] on given chromosome.
    pub fn range(&self, chr: ChrId, start_pos: u32, end_pos: u32) -> Vec<(u32, u64)> {
        let start_key = GenomicKey::new(chr, start_pos).as_u64();
        let end_key = GenomicKey::new(chr, end_pos).as_u64();

        let start_idx = self.keys.partition_point(|&k| k < start_key);
        let end_idx = self.keys.partition_point(|&k| k <= end_key);

        self.keys[start_idx..end_idx]
            .iter()
            .zip(&self.offsets[start_idx..end_idx])
            .map(|(&k, &off)| (k as u32, off))
            .collect()
    }

    /// Save index to file.
    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        persistence::save_index(path.as_ref(), &self.keys, &self.offsets, &self.mph)
    }

    /// Load index from file using memory mapping.
    pub fn load_mmap<P: AsRef<Path>>(path: P) -> Result<Self> {
        persistence::load_index_mmap(path.as_ref())
    }

    /// Load index from file (full read into memory).
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        persistence::load_index(path.as_ref())
    }

    /// Memory usage estimate in bytes.
    pub fn memory_usage(&self) -> usize {
        std::mem::size_of::<Self>()
            + self.keys.len() * std::mem::size_of::<u64>()
            + self.offsets.len() * std::mem::size_of::<u64>()
            + self.mph.g.len() * std::mem::size_of::<u32>()
    }

    /// Bytes per key (efficiency metric).
    pub fn bytes_per_key(&self) -> f64 {
        self.memory_usage() as f64 / self.len().max(1) as f64
    }

    /// Internal constructor from components.
    pub(crate) fn from_parts(mph: Mphf, keys: Vec<u64>, offsets: Vec<u64>) -> Self {
        Self { mph, keys, offsets }
    }

    /// Get reference to internal MPH (for persistence).
    pub(crate) fn mph(&self) -> &Mphf {
        &self.mph
    }

    /// Get reference to keys (for persistence).
    pub(crate) fn keys(&self) -> &[u64] {
        &self.keys
    }

    /// Get reference to offsets (for persistence).
    pub(crate) fn offsets(&self) -> &[u64] {
        &self.offsets
    }
}

/// Builder for constructing VcfIndex incrementally.
pub struct VcfIndexBuilder {
    entries: Vec<(u64, u64)>, // (key, offset)
    config: BuildConfig,
}

impl VcfIndexBuilder {
    /// Create new builder with default configuration.
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

    /// Create builder with pre-allocated capacity.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            entries: Vec::with_capacity(capacity),
            config: BuildConfig::default(),
        }
    }

    /// Set MPH gamma parameter (vertex ratio m/n).
    pub fn gamma(mut self, gamma: f64) -> Self {
        self.config.gamma = gamma;
        self
    }

    /// Add genomic position with its VCF byte offset.
    #[inline]
    pub fn add(&mut self, key: GenomicKey, offset: u64) -> Result<()> {
        let packed = key.0;
        self.entries.push((packed, offset));
        Ok(())
    }

    /// Add raw key-offset pair.
    #[inline]
    pub fn add_raw(&mut self, key: u64, offset: u64) -> Result<()> {
        self.entries.push((key, offset));
        Ok(())
    }

    /// Number of entries added.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Check if builder is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn build(mut self) -> Result<VcfIndex> {
        if self.entries.is_empty() {
            return Err(IndexError::EmptyDataset);
        }

        self.entries.sort_unstable_by_key(|(k, _)| *k);

        self.entries.dedup_by(|a, b| a.0 == b.0);

        let keys: Vec<u64> = self.entries.iter().map(|(k, _)| *k).collect();
        let offsets: Vec<u64> = self.entries.iter().map(|(_, o)| *o).collect();

        let key_bytes: Vec<[u8; 8]> = keys.iter().map(|k| k.to_le_bytes()).collect();
        let mph = KvBuilder::new()
            .with_config(self.config)
            .build(key_bytes.iter().map(|b| b.as_slice()))?;

        let n = keys.len();
        let mut sorted_keys = vec![0u64; n];
        let mut sorted_offsets = vec![0u64; n];

        for (i, &key) in keys.iter().enumerate() {
            let idx = mph.index(&key.to_le_bytes()) as usize;
            sorted_keys[idx] = key;
            sorted_offsets[idx] = offsets[i];
        }

        Ok(VcfIndex::from_parts(mph, sorted_keys, sorted_offsets))
    }
}

impl Default for VcfIndexBuilder {
    fn default() -> Self {
        Self::new()
    }
}
