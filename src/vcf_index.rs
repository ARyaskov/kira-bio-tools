//! VCF index implementation using a hybrid MPH + PGM index.
//!
//! Design:
//! - Minimal perfect hash (MPH, `kira_kv_engine::Mphf`) for O(1) point lookups.
//! - PGM-Index (`pgm_index::PGMIndex`) built over sorted genomic keys for
//!   fast range scans and introspection.
//!
//! On disk we persist only MPH + raw key/offset arrays; PGM is rebuilt
//! on load from the sorted keys.

use std::path::Path;

use kira_kv_engine::{BuildConfig, Builder, MphError, Mphf};
use pgm_index::{PGMIndex, PGMStats};
use thiserror::Error;

use crate::persistence;

/// Chromosome identifier (1–25 for standard chromosomes)
pub type ChrId = u8;

/// Genomic coordinate encoded as a single `u64` for efficient indexing.
///
/// Layout:
/// ```text
/// [ 8 bits chr_id ][ 32 bits position ][ 24 bits reserved ]
/// ```
///
/// In practice we only need chr_id + 32-bit POS; the high bits are kept
/// for potential future extensions and for better separation of chromosomes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GenomicKey(u64);

impl GenomicKey {
    /// Create a new genomic key from chromosome id and 1-based position.
    #[inline]
    pub fn new(chr: ChrId, pos: u32) -> Self {
        let v = ((chr as u64) << 32) | (pos as u64);
        Self(v)
    }

    /// Chromosome id encoded in this key.
    #[inline]
    pub fn chr(&self) -> ChrId {
        (self.0 >> 32) as u8
    }

    /// 1-based position encoded in this key.
    #[inline]
    pub fn position(&self) -> u32 {
        (self.0 & 0xFFFF_FFFF) as u32
    }

    /// Raw `u64` value used for indexing.
    #[inline]
    pub fn as_u64(&self) -> u64 {
        self.0
    }

    /// Construct from raw `u64` (inverse of `as_u64`).
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

/// Default epsilon for PGM index.
///
/// Trade-offs:
/// - Smaller epsilon -> more segments, higher memory, slightly faster.
/// - Larger epsilon  -> fewer segments, lower memory, slightly slower.
///
/// For 10^5–10^8 keys, 32 is a good general-purpose value.
const PGM_EPSILON: usize = 32;

/// High-performance VCF index using a minimal perfect hash + PGM-Index.
///
/// Responsibilities:
/// - Point lookup `GenomicKey -> offset` via MPH (O(1) average).
/// - Range lookup `[chr, start..=end]` via PGM-guided bounds + sequential scan.
/// - Persistence via `persistence::{save_index, load_index, load_index_mmap}`.
#[derive(Debug)]
pub struct VcfIndex {
    /// Minimal perfect hash for key → index mapping.
    mph: Mphf,
    /// PGM-Index built over sorted genomic keys (same logical keys as persisted).
    pgm: PGMIndex<u64>,
    /// VCF byte offsets corresponding to keys in the same order as PGM data.
    offsets: Vec<u64>,
}

impl VcfIndex {
    /// Load index from file using mmap (recommended for large indexes).
    pub fn load_mmap<P: AsRef<Path>>(path: P) -> Result<Self> {
        persistence::load_index_mmap(path.as_ref())
    }

    /// Load index from file with full read into memory.
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        persistence::load_index(path.as_ref())
    }

    /// Get raw slice of all keys (sorted by genomic coordinate).
    #[inline]
    pub fn keys(&self) -> &[u64] {
        self.pgm.data.as_ref()
    }

    /// Get raw slice of all offsets.
    #[inline]
    pub fn offsets(&self) -> &[u64] {
        &self.offsets
    }

    /// Get reference to the underlying MPH structure (used by persistence).
    #[inline]
    pub(crate) fn mph(&self) -> &Mphf {
        &self.mph
    }

    /// Number of indexed positions.
    #[inline]
    pub fn len(&self) -> usize {
        self.offsets.len()
    }

    /// Check if index is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.offsets.is_empty()
    }

    /// Point lookup: get byte offset for a genomic key.
    ///
    /// Uses MPH for O(1) index, with a verification against the key array
    /// to protect against corrupted data.
    #[inline]
    pub fn get(&self, key: GenomicKey) -> Option<u64> {
        if self.is_empty() {
            return None;
        }

        let key_u64 = key.as_u64();
        let idx = self.mph.index(&key_u64.to_le_bytes()) as usize;

        let keys = self.keys();
        if idx < keys.len() && keys[idx] == key_u64 {
            Some(self.offsets[idx])
        } else {
            None
        }
    }

    /// Check if a genomic key exists in the index.
    #[inline]
    pub fn contains(&self, key: GenomicKey) -> bool {
        self.get(key).is_some()
    }

    /// Batch lookup for multiple keys.
    pub fn get_batch(&self, keys: &[GenomicKey]) -> Vec<Option<u64>> {
        keys.iter().copied().map(|k| self.get(k)).collect()
    }

    /// Range query: find all positions in `[start_pos, end_pos]` on a given chromosome.
    ///
    /// Implementation:
    /// - Map chr,start/end to packed u64 keys.
    /// - Use PGM to get a tight bound when the boundary key exists.
    /// - Fallback to binary search (`partition_point`) otherwise.
    /// - Sequentially scan the resulting slice and decode positions back for output.
    pub fn range(&self, chr: ChrId, start_pos: u32, end_pos: u32) -> Vec<(u32, u64)> {
        if self.is_empty() || start_pos > end_pos {
            return Vec::new();
        }

        let start_key = GenomicKey::new(chr, start_pos).as_u64();
        let end_key = GenomicKey::new(chr, end_pos).as_u64();

        let keys = self.keys();
        let n = keys.len();
        if n == 0 {
            return Vec::new();
        }

        let start_idx = self.lower_bound(start_key);
        let end_idx = self.upper_bound(end_key);

        if start_idx >= end_idx {
            return Vec::new();
        }

        let mut result = Vec::with_capacity(end_idx.saturating_sub(start_idx));
        for (k, &off) in keys[start_idx..end_idx].iter().zip(&self.offsets[start_idx..end_idx]) {
            let key = GenomicKey::from_u64(*k);
            if key.chr() == chr {
                result.push((key.position(), off));
            }
        }
        result
    }

    /// Memory usage estimate in bytes (index + keys + offsets).
    ///
    /// Includes:
    /// - `self` struct
    /// - PGM-Index internal storage (keys + segments)
    /// - MPH displacement table `g`
    /// - offsets array
    pub fn memory_usage(&self) -> usize {
        std::mem::size_of::<Self>()
            + self.pgm.memory_usage()
            + self.offsets.len() * std::mem::size_of::<u64>()
            + self.mph.g.len() * std::mem::size_of::<u32>()
    }

    /// Bytes per key (efficiency metric).
    pub fn bytes_per_key(&self) -> f64 {
        if self.len() == 0 {
            return 0.0;
        }
        self.memory_usage() as f64 / self.len() as f64
    }

    /// Expose PGM stats for diagnostics / tuning.
    #[inline]
    pub fn pgm_stats(&self) -> PGMStats {
        self.pgm.stats()
    }

    /// Internal constructor from raw components (used by builder and persistence).
    ///
    /// Invariants:
    /// - `keys` must be sorted ascending.
    /// - `offsets.len() == keys.len()`.
    /// - `mph` must be built over the same set of keys.
    pub(crate) fn from_parts(mph: Mphf, keys: Vec<u64>, offsets: Vec<u64>) -> Self {
        debug_assert_eq!(keys.len(), offsets.len());

        let pgm = PGMIndex::new(keys, PGM_EPSILON);

        Self { mph, pgm, offsets }
    }

    /// Lower bound helper: first index `i` such that `keys[i] >= key`.
    ///
    /// Strategy:
    /// - If PGM knows the position of `key` exactly (it exists), use that
    ///   and scan left within equal-key run.
    /// - Otherwise, fall back to binary search via `partition_point`.
    fn lower_bound(&self, key: u64) -> usize {
        let keys = self.keys();
        if keys.is_empty() {
            return 0;
        }

        if let Some(mut idx) = self.pgm.get(key) {
            // PGM guarantees `keys[idx] == key` when it returns Some.
            while idx > 0 && keys[idx - 1] >= key {
                idx -= 1;
            }
            idx
        } else {
            // Key is not present; use binary search.
            keys.partition_point(|&k| k < key)
        }
    }

    /// Upper bound helper: first index `i` such that `keys[i] > key`.
    ///
    /// Strategy:
    /// - If PGM knows the position of `key` exactly, start from that index
    ///   and scan right within equal-key run.
    /// - Otherwise, fall back to binary search via `partition_point`.
    fn upper_bound(&self, key: u64) -> usize {
        let keys = self.keys();
        if keys.is_empty() {
            return 0;
        }

        let n = keys.len();

        if let Some(mut idx) = self.pgm.get(key) {
            while idx + 1 < n && keys[idx + 1] <= key {
                idx += 1;
            }
            idx + 1
        } else {
            keys.partition_point(|&k| k <= key)
        }
    }

    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let keys = self.pgm.data.as_ref();
        let offsets = &self.offsets;
        persistence::save_index(path.as_ref(), keys, offsets, &self.mph)
    }
}

/// Builder for constructing `VcfIndex` incrementally from VCF scan.
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
        self.entries.push((key.as_u64(), offset));
        Ok(())
    }

    /// Add raw key–offset pair.
    #[inline]
    pub fn add_raw(&mut self, key: u64, offset: u64) -> Result<()> {
        self.entries.push((key, offset));
        Ok(())
    }

    /// Number of entries added.
    #[inline]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Check if builder is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Build the hybrid MPH + PGM index.
    pub fn build(mut self) -> Result<VcfIndex> {
        if self.entries.is_empty() {
            return Err(IndexError::EmptyDataset);
        }

        // Sort by key for deterministic layout and PGM construction.
        self.entries.sort_unstable_by_key(|(k, _)| *k);

        // Check for duplicate keys – we require a functional mapping.
        for window in self.entries.windows(2) {
            if window[0].0 == window[1].0 {
                let key = GenomicKey::from_u64(window[0].0);
                return Err(IndexError::DuplicateKey(key.chr(), key.position()));
            }
        }

        let keys: Vec<u64> = self.entries.iter().map(|(k, _)| *k).collect();
        let offsets: Vec<u64> = self.entries.iter().map(|(_, o)| *o).collect();

        // Build MPH over raw key bytes.
        let key_bytes: Vec<[u8; 8]> = keys.iter().map(|k| k.to_le_bytes()).collect();
        let mph = Builder::new()
            .with_config(self.config)
            .build(key_bytes.iter().map(|b| b.as_slice()))?;

        // Reorder keys and offsets according to MPH to get true O(1) layout.
        let n = keys.len();
        let mut mph_keys = vec![0u64; n];
        let mut mph_offsets = vec![0u64; n];

        for (i, &key) in keys.iter().enumerate() {
            let idx = mph.index(&key.to_le_bytes()) as usize;
            mph_keys[idx] = key;
            mph_offsets[idx] = offsets[i];
        }

        Ok(VcfIndex::from_parts(mph, mph_keys, mph_offsets))
    }
}

impl Default for VcfIndexBuilder {
    fn default() -> Self {
        Self::new()
    }
}
