use crate::pgm_index::Key;
use serde::{Deserialize, Serialize};

/// Key-Value entry
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct KVEntry<K: Key, V> {
    pub key: K,
    pub value: V,
}

impl<K: Key, V> KVEntry<K, V> {
    pub fn new(key: K, value: V) -> Self {
        Self { key, value }
    }
}

/// Error types for all index structures
#[derive(Debug, thiserror::Error)]
pub enum KVError {
    #[error("Key not found: {0}")]
    KeyNotFound(String),

    #[error("Invalid epsilon value: {0}")]
    InvalidEpsilon(usize),

    #[error("Empty dataset provided")]
    EmptyDataset,

    #[error("Segment size too large: {0}")]
    SegmentTooLarge(usize),

    #[error("Perfect hash construction failed: {0}")]
    PerfectHashFailed(String),

    #[error("Hybrid index construction failed: {0}")]
    HybridIndexFailed(String),

    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),

    #[error("Serialization error: {0}")]
    SerializationError(String),

    #[error("Configuration error: {0}")]
    ConfigError(String),
}

pub type KVResult<T> = Result<T, KVError>;

/// Performance metrics for any index structure
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerformanceMetrics {
    pub structure_type: IndexType,
    pub total_keys: usize,
    pub memory_usage_bytes: usize,
    pub build_time_ms: u64,
    pub total_queries: u64,
    pub cache_hits: u64,
    pub cache_hit_rate: f64,
    pub avg_query_time_ns: f64,
    pub throughput_queries_per_sec: f64,
}

/// Type of index structure
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum IndexType {
    PGMIndex,
    PerfectHash,
    HybridIndex,
    KeyValueStore,
    CompressedStore,
}

/// Configuration for building optimal index structures
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IndexConfig {
    /// PGM epsilon parameter
    pub pgm_epsilon: usize,
    /// Maximum segment size for perfect hashing
    pub max_segment_size: usize,
    /// Target number of segments
    pub target_segments: Option<usize>,
    /// Enable parallel construction
    pub parallel_build: bool,
    /// Number of threads for parallel operations
    pub num_threads: Option<usize>,
    /// Enable compression
    pub compression: bool,
    /// Data pattern hint for optimization
    pub data_pattern: Option<DataPattern>,
}

impl Default for IndexConfig {
    fn default() -> Self {
        Self {
            pgm_epsilon: 64,
            max_segment_size: 2000,
            target_segments: None,
            parallel_build: true,
            num_threads: None,
            compression: false,
            data_pattern: None,
        }
    }
}

/// Data pattern types for optimization
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum DataPattern {
    Sequential,  // 0, 1, 2, 3, ...
    SmallGaps,   // 0, 2, 4, 6, ...
    LargeGaps,   // 0, 100, 200, 300, ...
    Exponential, // 1, 2, 4, 8, 16, ...
    Random,      // Random but sorted
    Clustered,   // Data grouped in clusters
    Genomic,     // Genomic positions (sparse with clusters)
}

/// Query statistics for performance analysis
#[derive(Debug, Clone, Default)]
pub struct QueryStats {
    pub sequential_queries: u64,
    pub random_queries: u64,
    pub range_queries: u64,
    pub batch_queries: u64,
    pub total_time_ns: u64,
    pub fastest_query_ns: u64,
    pub slowest_query_ns: u64,
}

/// Range query specification
#[derive(Debug, Clone, Copy)]
pub struct Range<K: Key> {
    pub start: K,
    pub end: K,
    pub inclusive: bool,
}

impl<K: Key> Range<K> {
    pub fn new(start: K, end: K) -> Self {
        Self {
            start,
            end,
            inclusive: true,
        }
    }

    pub fn exclusive(start: K, end: K) -> Self {
        Self {
            start,
            end,
            inclusive: false,
        }
    }

    pub fn contains(&self, key: K) -> bool {
        if self.inclusive {
            key >= self.start && key <= self.end
        } else {
            key >= self.start && key < self.end
        }
    }
}

/// Batch operation specification
#[derive(Debug, Clone)]
pub struct BatchOperation<K: Key, V> {
    pub operation_type: OperationType,
    pub keys: Vec<K>,
    pub values: Option<Vec<V>>,
    pub ranges: Option<Vec<Range<K>>>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum OperationType {
    Get,
    Put,
    Delete,
    RangeQuery,
    Contains,
}

/// Configuration builder for easy setup
pub struct IndexConfigBuilder {
    config: IndexConfig,
}

impl IndexConfigBuilder {
    pub fn new() -> Self {
        Self {
            config: IndexConfig::default(),
        }
    }

    pub fn pgm_epsilon(mut self, epsilon: usize) -> Self {
        self.config.pgm_epsilon = epsilon;
        self
    }

    pub fn max_segment_size(mut self, size: usize) -> Self {
        self.config.max_segment_size = size;
        self
    }

    pub fn target_segments(mut self, segments: usize) -> Self {
        self.config.target_segments = Some(segments);
        self
    }

    pub fn parallel_build(mut self, enable: bool) -> Self {
        self.config.parallel_build = enable;
        self
    }

    pub fn num_threads(mut self, threads: usize) -> Self {
        self.config.num_threads = Some(threads);
        self
    }

    pub fn compression(mut self, enable: bool) -> Self {
        self.config.compression = enable;
        self
    }

    pub fn data_pattern(mut self, pattern: DataPattern) -> Self {
        self.config.data_pattern = Some(pattern);
        self
    }

    /// Auto-configure based on data analysis
    pub fn auto_configure<T: Key>(mut self, data: &[T]) -> Self {
        let pattern = crate::types::key_utils::detect_pattern(data);
        let epsilon = match pattern {
            DataPattern::Sequential => (data.len() / 1000).max(64).min(512),
            DataPattern::SmallGaps => (data.len() / 500).max(32).min(256),
            DataPattern::LargeGaps => (data.len() / 200).max(16).min(128),
            DataPattern::Exponential => (data.len() / 100).max(8).min(64),
            DataPattern::Random => (data.len() / 50).max(16).min(64),
            DataPattern::Clustered => (data.len() / 300).max(32).min(128),
            DataPattern::Genomic => (data.len() / 400).max(16).min(64),
        };
        let segment_size = (data.len() / self.config.target_segments.unwrap_or(data.len() / 1000))
            .max(100)
            .min(10_000);

        self.config.data_pattern = Some(pattern);
        self.config.pgm_epsilon = epsilon;
        self.config.max_segment_size = segment_size;
        self
    }

    pub fn build(self) -> IndexConfig {
        self.config
    }
}

impl Default for IndexConfigBuilder {
    fn default() -> Self {
        Self::new()
    }
}

/// Comprehensive benchmark results
#[derive(Debug, Clone)]
pub struct BenchmarkResults {
    pub index_type: IndexType,
    pub data_size: usize,
    pub data_pattern: DataPattern,
    pub config: IndexConfig,
    pub metrics: PerformanceMetrics,
    pub query_stats: QueryStats,
    pub comparison_baselines: Vec<(String, PerformanceMetrics)>,
}

/// Memory usage breakdown
#[derive(Debug, Clone)]
pub struct MemoryBreakdown {
    pub index_structure: usize,
    pub data_storage: usize,
    pub metadata: usize,
    pub cache: usize,
    pub total: usize,
}

impl MemoryBreakdown {
    pub fn bytes_per_key(&self, num_keys: usize) -> f64 {
        self.total as f64 / num_keys.max(1) as f64
    }

    pub fn efficiency_ratio(&self) -> f64 {
        self.data_storage as f64 / self.total as f64
    }
}

/// Utility functions for working with different key types
pub mod key_utils {
    use super::*;

    /// Convert any key to u64 for internal processing
    pub fn key_to_u64<K: Key>(key: K) -> u64 {
        key.to_u64().unwrap_or(0)
    }

    /// Estimate key range for optimization
    pub fn estimate_key_range<K: Key>(keys: &[K]) -> Option<(K, K)> {
        if keys.is_empty() {
            None
        } else {
            Some((keys[0], keys[keys.len() - 1]))
        }
    }

    /// Calculate key density (useful for optimization)
    pub fn calculate_key_density<K: Key>(keys: &[K]) -> f64 {
        if keys.len() < 2 {
            return 1.0;
        }

        let start = key_to_u64(keys[0]) as f64;
        let end = key_to_u64(keys[keys.len() - 1]) as f64;
        let range = end - start;

        if range > 0.0 {
            keys.len() as f64 / range
        } else {
            1.0
        }
    }

    /// Detect data pattern automatically
    pub fn detect_pattern<K: Key>(keys: &[K]) -> DataPattern {
        if keys.len() < 10 {
            return DataPattern::Sequential;
        }

        let sample_size = keys.len().min(1000);
        let sample = &keys[0..sample_size];

        // Calculate gaps
        let mut gaps = Vec::new();
        for i in 1..sample.len() {
            let curr = key_to_u64(sample[i]) as f64;
            let prev = key_to_u64(sample[i - 1]) as f64;
            gaps.push(curr - prev);
        }

        if gaps.is_empty() {
            return DataPattern::Sequential;
        }

        let avg_gap = gaps.iter().sum::<f64>() / gaps.len() as f64;
        let variance = gaps.iter().map(|&g| (g - avg_gap).powi(2)).sum::<f64>() / gaps.len() as f64;

        let coefficient_of_variation = if avg_gap.abs() > f64::EPSILON {
            variance.sqrt() / avg_gap.abs()
        } else {
            0.0
        };

        // Check for sequential pattern (gaps mostly = 1)
        let unit_gaps = gaps.iter().filter(|&&g| (g - 1.0).abs() < 0.1).count();
        if unit_gaps as f64 / gaps.len() as f64 > 0.8 {
            return DataPattern::Sequential;
        }

        // Check for clustering
        let large_gaps = gaps.iter().filter(|&&g| g > avg_gap * 10.0).count();
        if large_gaps as f64 / gaps.len() as f64 > 0.1 {
            return DataPattern::Clustered;
        }

        // Classify by coefficient of variation
        match coefficient_of_variation {
            cv if cv < 0.1 => DataPattern::SmallGaps,
            cv if cv < 1.0 => DataPattern::LargeGaps,
            cv if cv < 10.0 => DataPattern::Exponential,
            _ => DataPattern::Random,
        }
    }
}

/// Validation utilities
pub mod validation {
    use super::*;

    /// Validate index configuration
    pub fn validate_config(config: &IndexConfig) -> KVResult<()> {
        if config.pgm_epsilon == 0 {
            return Err(KVError::InvalidEpsilon(config.pgm_epsilon));
        }

        if config.max_segment_size == 0 {
            return Err(KVError::SegmentTooLarge(config.max_segment_size));
        }

        if config.max_segment_size > 100_000 {
            return Err(KVError::SegmentTooLarge(config.max_segment_size));
        }

        if let Some(threads) = config.num_threads {
            if threads == 0 {
                return Err(KVError::ConfigError(
                    "Number of threads cannot be zero".to_string(),
                ));
            }
        }

        Ok(())
    }

    /// Validate key-value data
    pub fn validate_data<K: Key, V>(data: &[KVEntry<K, V>]) -> KVResult<()> {
        if data.is_empty() {
            return Err(KVError::EmptyDataset);
        }

        // Check for key ordering (should be sorted)
        for i in 1..data.len() {
            if data[i].key < data[i - 1].key {
                return Err(KVError::ConfigError(
                    "Data must be sorted by keys".to_string(),
                ));
            }
        }

        Ok(())
    }

    /// Validate keys only
    pub fn validate_keys<K: Key>(keys: &[K]) -> KVResult<()> {
        if keys.is_empty() {
            return Err(KVError::EmptyDataset);
        }

        // Check for key ordering
        for i in 1..keys.len() {
            if keys[i] < keys[i - 1] {
                return Err(KVError::ConfigError("Keys must be sorted".to_string()));
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_config_builder() {
        let config = IndexConfigBuilder::new()
            .pgm_epsilon(32)
            .max_segment_size(1000)
            .parallel_build(true)
            .compression(false)
            .build();

        assert_eq!(config.pgm_epsilon, 32);
        assert_eq!(config.max_segment_size, 1000);
        assert_eq!(config.parallel_build, true);
        assert_eq!(config.compression, false);
    }

    #[test]
    fn test_auto_configure() {
        let keys: Vec<u64> = (0..10000).collect();
        let config = IndexConfigBuilder::new().auto_configure(&keys).build();

        assert!(config.pgm_epsilon > 0);
        assert!(config.max_segment_size > 0);
        assert!(config.data_pattern.is_some());
    }

    #[test]
    fn test_range() {
        let range = Range::new(10u64, 20u64);
        assert!(range.contains(15));
        assert!(range.contains(10));
        assert!(range.contains(20));
        assert!(!range.contains(5));
        assert!(!range.contains(25));

        let exclusive_range = Range::exclusive(10u64, 20u64);
        assert!(exclusive_range.contains(15));
        assert!(exclusive_range.contains(10));
        assert!(!exclusive_range.contains(20));
    }

    #[test]
    fn test_key_utils() {
        use key_utils::*;

        let keys = vec![1u64, 5, 10, 15, 20];
        let (min, max) = estimate_key_range(&keys).unwrap();
        assert_eq!(min, 1);
        assert_eq!(max, 20);

        let density = calculate_key_density(&keys);
        assert!(density > 0.0);

        let pattern = detect_pattern(&keys);
        assert_eq!(pattern, DataPattern::LargeGaps);
    }

    #[test]
    fn test_validation() {
        use validation::*;

        // Valid configuration
        let config = IndexConfig::default();
        assert!(validate_config(&config).is_ok());

        // Invalid epsilon
        let mut bad_config = IndexConfig::default();
        bad_config.pgm_epsilon = 0;
        assert!(validate_config(&bad_config).is_err());

        // Valid keys
        let keys = vec![1u64, 2, 3, 4, 5];
        assert!(validate_keys(&keys).is_ok());

        // Invalid keys (not sorted)
        let bad_keys = vec![5u64, 1, 3, 2, 4];
        assert!(validate_keys(&bad_keys).is_err());
    }

    #[test]
    fn test_memory_breakdown() {
        let breakdown = MemoryBreakdown {
            index_structure: 1000,
            data_storage: 8000,
            metadata: 500,
            cache: 500,
            total: 10000,
        };

        assert_eq!(breakdown.bytes_per_key(1000), 10.0);
        assert_eq!(breakdown.efficiency_ratio(), 0.8);
    }

    #[test]
    fn test_kv_entry() {
        let entry = KVEntry::new(42u64, "hello".to_string());
        assert_eq!(entry.key, 42);
        assert_eq!(entry.value, "hello");

        // Test serialization
        let serialized = serde_json::to_string(&entry).unwrap();
        let deserialized: KVEntry<u64, String> = serde_json::from_str(&serialized).unwrap();
        assert_eq!(entry, deserialized);
    }
}
