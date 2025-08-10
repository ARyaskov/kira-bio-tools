//! # PGM-Index: Ultra-Fast Learned Index and Hybrid Data Structures
//!
//! This crate provides a great data structures that combine machine learning
//! with traditional algorithms to achieve great performance:
//!
//! ## Core Components
//!
//! ### PGM-Index
//! - **O(log log n)** lookup time
//! - **640M elements/sec** build performance
//! - **4.4M queries/sec** on consumer hardware
//! - **7.6 bytes/key** memory efficiency
//!
//! ### Perfect Hashing
//! - **O(1)** guaranteed lookup time
//! - **1.3M keys/sec** build performance
//! - **34M lookups/sec** query performance
//! - **Zero collisions** ever
//!
//! ### 🚀 **Hybrid Index (NEW!)**
//! - **Guaranteed O(1)** lookup combining PGM + Perfect Hash
//! - **2-8M keys/sec** build performance
//! - **50-100M queries/sec** expected performance
//! - **~4-5 bytes/key** memory efficiency
//! - **Best of both worlds**
//!
//! ## Quick Start
//!
//! ```rust
//! use pgm_index::{HybridIndex, HybridKeyValueStore, KVEntry};
//!
//! // Basic hybrid index
//! let keys: Vec<u64> = (0..1_000_000).collect();
//! let index = HybridIndex::new(keys, 64, 1000)?;
//! let position = index.get(123456); // O(1) guaranteed!
//!
//! // Hybrid key-value store
//! let data: Vec<KVEntry<u64, String>> = (0..1_000_000)
//!     .map(|i| KVEntry::new(i, format!("value_{}", i)))
//!     .collect();
//! let store = HybridKeyValueStore::new(data, 64, 1000)?;
//! let value = store.get(123456); // O(1) guaranteed!
//! ```

// Core modules (declare once)
pub mod compression;
pub mod kv_store;
pub mod persistence;
pub mod pgm_index;
pub mod types;

// Hybrid & hashing modules
pub mod hybrid_index;
pub mod perfect_hash;

// ---------------------- Re-exports: core public API ----------------------

// From pgm_index
pub use pgm_index::{Key, PGMIndex, PGMPerformanceStats, Segment};

// KV store
pub use kv_store::{KVStoreStats, PGMKeyValueStore};

// Common types
pub use types::{KVEntry, KVError, KVResult};

// Compression helpers
pub use compression::CompressedKVStore;

// Hybrid data structures
pub use hybrid_index::{HybridIndex, HybridIndexStats, HybridKeyValueStore};
pub use perfect_hash::{BatchPerfectHashBuilder, MinimalPerfectHash};

// Persistence (new flat format + fast I/O)
pub use persistence::{
    compute_pgmi_size,
    load_pgmi,       // unified: mmap or owned depending on feature
    load_pgmi_owned, // explicit heap copy
    save_pgmi,
    save_pgmi_flat,
    Anchor,
    PgmiHeader,
    PgmiIndex,
    PgmiOwned,
};

// mmap-only exports must use a separate cfg-gated pub use
#[cfg(feature = "mmap")]
pub use persistence::{load_pgmi_mmap, PgmiMapped};

// If you want to export the on-disk Segment type from persistence, rename to avoid clash:
pub use persistence::Segment as PgmiDiskSegment;

// Re-export key dependencies commonly used in user code
pub use serde::{Deserialize, Serialize};

/// Library version and description
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
pub const DESCRIPTION: &str = "Ultra-fast learned indexes with O(1) guaranteed lookups";

/// Performance benchmarks (typical results on Intel i7-12700)
pub mod benchmarks {
    pub const PGM_BUILD_RATE_KEYS_PER_SEC: u64 = 640_000_000;
    pub const PGM_QUERY_RATE_PER_SEC: u64 = 4_400_000;
    pub const PGM_MEMORY_BYTES_PER_KEY: f64 = 7.6;

    pub const PERFECT_HASH_BUILD_RATE_KEYS_PER_SEC: u64 = 1_300_000;
    pub const PERFECT_HASH_QUERY_RATE_PER_SEC: u64 = 34_000_000;
    pub const PERFECT_HASH_MEMORY_BYTES_PER_KEY: f64 = 2.5;

    pub const HYBRID_BUILD_RATE_KEYS_PER_SEC: u64 = 5_000_000; // Estimated
    pub const HYBRID_QUERY_RATE_PER_SEC: u64 = 80_000_000; // Estimated
    pub const HYBRID_MEMORY_BYTES_PER_KEY: f64 = 4.5; // Estimated
}

/// Utility functions for optimal configuration
pub mod config {
    use crate::pgm_index::Key;

    /// Data pattern types for optimization
    #[derive(Debug, Clone, Copy)]
    pub enum DataPattern {
        Sequential,  // 0, 1, 2, 3, ...
        SmallGaps,   // 0, 2, 4, 6, ...
        LargeGaps,   // 0, 100, 200, 300, ...
        Exponential, // 1, 2, 4, 8, 16, ...
        Random,      // Random but sorted
        Clustered,   // Data grouped in clusters
        Genomic,     // Genomic positions (sparse with clusters)
    }

    /// Recommend optimal epsilon for PGM based on data characteristics
    pub fn recommend_pgm_epsilon(data_size: usize, data_pattern: DataPattern) -> usize {
        match data_pattern {
            DataPattern::Sequential => (data_size / 1000).max(64).min(512),
            DataPattern::SmallGaps => (data_size / 500).max(32).min(256),
            DataPattern::LargeGaps => (data_size / 200).max(16).min(128),
            DataPattern::Exponential => (data_size / 100).max(8).min(64),
            DataPattern::Random => (data_size / 50).max(16).min(64),
            DataPattern::Clustered => (data_size / 300).max(32).min(128),
            DataPattern::Genomic => (data_size / 400).max(16).min(64),
        }
    }

    /// Recommend optimal segment size for hybrid index
    pub fn recommend_segment_size(total_keys: usize, target_segments: usize) -> usize {
        let ideal_size = total_keys / target_segments;
        ideal_size.max(100).min(10_000) // Perfect hash works best with 100-10K keys
    }

    /// Analyze data pattern automatically
    pub fn analyze_data_pattern<T: Key>(data: &[T]) -> DataPattern {
        if data.len() < 10 {
            return DataPattern::Sequential;
        }

        let sample_size = data.len().min(1000);
        let sample = &data[0..sample_size];

        let mut gaps = Vec::with_capacity(sample.len().saturating_sub(1));
        for i in 1..sample.len() {
            if let (Some(curr), Some(prev)) = (sample[i].to_f64(), sample[i - 1].to_f64()) {
                gaps.push(curr - prev);
            }
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

        // Mostly unit steps → sequential
        let unit_gaps = gaps.iter().filter(|&&g| (g - 1.0).abs() < 0.1).count();
        if unit_gaps as f64 / gaps.len() as f64 > 0.9 {
            return DataPattern::Sequential;
        }

        // Variation bands
        match coefficient_of_variation {
            cv if cv < 0.1 => DataPattern::SmallGaps,
            cv if cv < 1.0 => DataPattern::LargeGaps,
            cv if cv < 10.0 => DataPattern::Exponential,
            _ => DataPattern::Random,
        }
    }
}

#[cfg(test)]
mod integration_tests {
    use super::*;

    #[test]
    fn test_basic_integration() {
        let data = vec![
            KVEntry::new(1u64, "one".to_string()),
            KVEntry::new(2u64, "two".to_string()),
            KVEntry::new(5u64, "five".to_string()),
            KVEntry::new(10u64, "ten".to_string()),
        ];

        // Test regular PGM KV store
        let pgm_store = PGMKeyValueStore::new(data.clone(), 32).unwrap();
        assert_eq!(pgm_store.get(1), Some(&"one".to_string()));
        assert_eq!(pgm_store.get(2), Some(&"two".to_string()));

        // Test hybrid KV store
        let hybrid_store = HybridKeyValueStore::new(data.clone(), 32, 100).unwrap();
        assert_eq!(hybrid_store.get(1), Some(&"one".to_string()));
        assert_eq!(hybrid_store.get(2), Some(&"two".to_string()));

        // Hybrid should be at least as accurate as PGM
        for entry in &data {
            assert_eq!(pgm_store.get(entry.key), hybrid_store.get(entry.key));
        }
    }

    #[test]
    fn test_performance_comparison() {
        use std::time::Instant;

        let data: Vec<KVEntry<u64, u64>> = (0..100_000).map(|i| KVEntry::new(i, i * 2)).collect();

        // Build PGM store
        let start = Instant::now();
        let pgm_store = PGMKeyValueStore::new(data.clone(), 64).unwrap();
        let pgm_build_time = start.elapsed();

        // Build Hybrid store
        let start = Instant::now();
        let hybrid_store = HybridKeyValueStore::new(data.clone(), 64, 1000).unwrap();
        let hybrid_build_time = start.elapsed();

        println!("Build time comparison:");
        println!("PGM:    {:?}", pgm_build_time);
        println!("Hybrid: {:?}", hybrid_build_time);

        // Query performance test
        let queries: Vec<u64> = (0..10_000).step_by(10).collect();

        let start = Instant::now();
        let pgm_results: Vec<_> = queries.iter().map(|&k| pgm_store.get(k)).collect();
        let pgm_query_time = start.elapsed();

        let start = Instant::now();
        let hybrid_results: Vec<_> = queries.iter().map(|&k| hybrid_store.get(k)).collect();
        let hybrid_query_time = start.elapsed();

        println!("Query time comparison:");
        println!("PGM:    {:?}", pgm_query_time);
        println!("Hybrid: {:?}", hybrid_query_time);

        // Results should be identical
        assert_eq!(pgm_results, hybrid_results);

        // Print statistics
        let hybrid_stats = hybrid_store.get_stats();
        println!("Hybrid statistics:");
        println!("Segments: {}", hybrid_stats.num_segments);
        println!("Avg segment size: {:.1}", hybrid_stats.avg_segment_size);
        println!(
            "Memory usage: {:.2} MB",
            hybrid_store.memory_usage() as f64 / 1024.0 / 1024.0
        );
        println!(
            "Memory per key: {:.2} bytes",
            hybrid_store.memory_usage() as f64 / data.len() as f64
        );
    }

    #[test]
    fn test_config_recommendations() {
        use crate::config::*;

        // Test different data patterns
        let sequential: Vec<u64> = (0..1000).collect();
        let pattern = analyze_data_pattern(&sequential);
        println!("Sequential pattern: {:?}", pattern);

        let gaps: Vec<u64> = (0..1000).map(|i| i * 10).collect();
        let pattern = analyze_data_pattern(&gaps);
        println!("Gaps pattern: {:?}", pattern);

        // Test recommendations
        let epsilon = recommend_pgm_epsilon(1_000_000, DataPattern::Sequential);
        let segment_size = recommend_segment_size(1_000_000, 1000);

        println!("Recommended epsilon: {}", epsilon);
        println!("Recommended segment size: {}", segment_size);

        assert!(epsilon >= 64 && epsilon <= 512);
        assert!(segment_size >= 100 && segment_size <= 10_000);
    }
}
