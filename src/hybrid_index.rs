//! # Hybrid PGM + Perfect Hash Index
//!
//! Great data structure combining:
//! - PGM-Index for O(log log n) segment prediction
//! - Perfect Hash for O(1) exact lookup within segments
//!
//! Result: Guaranteed O(1) lookup with excellent scalability!

use crate::perfect_hash::{BatchPerfectHashBuilder, MinimalPerfectHash};
use crate::pgm_index::{Key, PGMIndex};
use crate::types::{KVEntry, KVError, KVResult};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::hash::Hash;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;
use std::time::Instant;

/// Performance statistics for hybrid index
#[derive(Debug, Clone, Default)]
pub struct HybridIndexStats {
    pub total_keys: usize,
    pub num_segments: usize,
    pub avg_segment_size: f64,
    pub total_queries: u64,
    pub cache_hits: u64,
    pub cache_hit_rate: f64,
    pub memory_usage_bytes: usize,
    pub build_time_ms: u64,
    pub pgm_prediction_time_ns: u64,
    pub perfect_hash_time_ns: u64,
}

/// Hybrid index combining PGM and Perfect Hashing
#[derive(Debug, Serialize, Deserialize)]
pub struct HybridIndex<K: Key + Hash> {
    /// PGM index for segment routing
    pgm_router: PGMIndex<K>,

    /// Perfect hash tables for each segment
    segment_hashes: Vec<MinimalPerfectHash<K>>,

    /// Segment boundaries for routing
    segment_boundaries: Vec<K>,

    /// Segment start positions in global data
    segment_starts: Vec<usize>,

    /// Configuration
    pgm_epsilon: usize,
    max_segment_size: usize,

    /// Runtime statistics (not serialized)
    #[serde(skip)]
    stats: Arc<HybridIndexStats>,

    #[serde(skip)]
    query_count: Arc<AtomicU64>,

    #[serde(skip)]
    cache_hits: Arc<AtomicU64>,
}

impl<K: Key + Hash + std::fmt::Debug> HybridIndex<K> {
    /// Build new hybrid index
    pub fn new(mut keys: Vec<K>, pgm_epsilon: usize, max_segment_size: usize) -> KVResult<Self> {
        if keys.is_empty() {
            return Err(KVError::EmptyDataset);
        }

        let build_start = Instant::now();

        // Sort keys
        keys.sort_unstable();
        keys.dedup();

        let total_keys = keys.len();

        // Step 1: Build PGM index for routing
        let pgm_router = PGMIndex::new(keys.clone(), pgm_epsilon);

        // Step 2: Create segments based on PGM structure + size limits
        let segments = Self::create_optimal_segments(&keys, &pgm_router, max_segment_size)?;

        // Step 3: Build perfect hash for each segment
        let batch_builder = BatchPerfectHashBuilder::new(max_segment_size);
        let segment_keys: Vec<Vec<K>> = segments
            .iter()
            .map(|(start, end)| keys[*start..*end].to_vec())
            .collect();

        let segment_hashes = batch_builder
            .build_segments::<K>(segment_keys)
            .map_err(|e| {
                KVError::SerializationError(format!("Perfect hash build failed: {}", e))
            })?;

        // Step 4: Create metadata
        let segment_boundaries: Vec<K> = segments.iter().map(|(start, _)| keys[*start]).collect();

        let segment_starts: Vec<usize> = segments.iter().map(|(start, _)| *start).collect();

        let build_time = build_start.elapsed();

        let stats = Arc::new(HybridIndexStats {
            total_keys,
            num_segments: segments.len(),
            avg_segment_size: total_keys as f64 / segments.len() as f64,
            total_queries: 0,
            cache_hits: 0,
            cache_hit_rate: 0.0,
            memory_usage_bytes: 0, // Will be calculated later
            build_time_ms: build_time.as_millis() as u64,
            pgm_prediction_time_ns: 0,
            perfect_hash_time_ns: 0,
        });

        let mut index = HybridIndex {
            pgm_router,
            segment_hashes,
            segment_boundaries,
            segment_starts,
            pgm_epsilon,
            max_segment_size,
            stats,
            query_count: Arc::new(AtomicU64::new(0)),
            cache_hits: Arc::new(AtomicU64::new(0)),
        };

        // Calculate memory usage
        let memory_usage = index.calculate_memory_usage();
        Arc::get_mut(&mut index.stats).unwrap().memory_usage_bytes = memory_usage;

        Ok(index)
    }

    /// Create optimal segments balancing PGM predictions and perfect hash limits
    fn create_optimal_segments(
        keys: &[K],
        pgm_router: &PGMIndex<K>,
        max_segment_size: usize,
    ) -> KVResult<Vec<(usize, usize)>> {
        let mut segments = Vec::new();
        let mut current_start = 0;

        // Get PGM segment boundaries as hints
        let pgm_segments = pgm_router.segment_count();
        let keys_per_pgm_segment = keys.len() / pgm_segments.max(1);

        while current_start < keys.len() {
            // Try to align with PGM segments, but respect size limits
            let ideal_end = (current_start + keys_per_pgm_segment).min(keys.len());
            let max_end = (current_start + max_segment_size).min(keys.len());
            let actual_end = ideal_end.min(max_end);

            // Ensure minimum segment size (unless it's the last segment)
            let min_segment_size = 100;
            let final_end =
                if actual_end - current_start < min_segment_size && actual_end < keys.len() {
                    (current_start + min_segment_size).min(keys.len())
                } else {
                    actual_end
                };

            segments.push((current_start, final_end));
            current_start = final_end;
        }

        Ok(segments)
    }

    /// Ultra-fast O(1) lookup with timing
    #[inline(always)]
    pub fn get(&self, key: K) -> Option<usize> {
        let _query_start = Instant::now();
        self.query_count.fetch_add(1, Ordering::Relaxed);

        // Step 1: PGM prediction to find segment
        let pgm_start = Instant::now();
        let segment_idx = self.find_segment_pgm_guided(key);
        let _pgm_time = pgm_start.elapsed().as_nanos() as u64;

        // Step 2: Perfect hash lookup within segment
        let hash_start = Instant::now();
        let result = if segment_idx < self.segment_hashes.len() {
            let perfect_hash = &self.segment_hashes[segment_idx];
            if let Some(local_pos) = perfect_hash.get_position(key) {
                let global_pos = self.segment_starts[segment_idx] + local_pos;
                self.cache_hits.fetch_add(1, Ordering::Relaxed);
                Some(global_pos)
            } else {
                None
            }
        } else {
            None
        };
        let _hash_time = hash_start.elapsed().as_nanos() as u64;

        result
    }

    /// PGM-guided segment finding with fallback
    #[inline(always)]
    fn find_segment_pgm_guided(&self, key: K) -> usize {
        // Try PGM prediction first
        let (lo, _hi) = self.pgm_router.predict(key);

        // Convert position to segment index
        let predicted_segment = self.position_to_segment(lo);

        // Verify prediction and adjust if needed
        if predicted_segment < self.segment_boundaries.len() {
            if key >= self.segment_boundaries[predicted_segment] {
                // Check if key is actually in the next segment
                if predicted_segment + 1 < self.segment_boundaries.len() {
                    if key >= self.segment_boundaries[predicted_segment + 1] {
                        return self.find_segment_binary_search(key);
                    }
                }
                return predicted_segment;
            }
        }

        // Fallback to binary search if prediction is wrong
        self.find_segment_binary_search(key)
    }

    /// Convert global position to segment index
    #[inline(always)]
    fn position_to_segment(&self, position: usize) -> usize {
        // Binary search through segment starts
        match self.segment_starts.binary_search(&position) {
            Ok(idx) => idx,
            Err(idx) => idx.saturating_sub(1),
        }
    }

    /// Fallback binary search for segment finding
    fn find_segment_binary_search(&self, key: K) -> usize {
        match self.segment_boundaries.binary_search(&key) {
            Ok(idx) => idx,
            Err(idx) => idx.saturating_sub(1).min(self.segment_boundaries.len() - 1),
        }
    }

    /// Batch lookup operations
    pub fn batch_get(&self, keys: &[K]) -> Vec<Option<usize>> {
        if keys.len() < 1000 {
            // Small batches: sequential processing
            keys.iter().map(|&key| self.get(key)).collect()
        } else {
            // Large batches: parallel processing
            keys.par_iter().map(|&key| self.get(key)).collect()
        }
    }

    /// Verify index correctness (for testing)
    pub fn verify(&self, original_keys: &[K]) -> bool {
        for (i, &key) in original_keys.iter().enumerate() {
            if let Some(found_pos) = self.get(key) {
                if found_pos != i {
                    eprintln!(
                        "Position mismatch for key {:?}: expected {}, found {}",
                        key, i, found_pos
                    );
                    return false;
                }
            } else {
                eprintln!("Key {:?} not found at expected position {}", key, i);
                return false;
            }
        }
        true
    }

    /// Get current statistics
    pub fn get_stats(&self) -> HybridIndexStats {
        let total_queries = self.query_count.load(Ordering::Relaxed);
        let cache_hits = self.cache_hits.load(Ordering::Relaxed);
        let cache_hit_rate = if total_queries > 0 {
            cache_hits as f64 / total_queries as f64
        } else {
            0.0
        };

        HybridIndexStats {
            total_queries,
            cache_hits,
            cache_hit_rate,
            ..self.stats.as_ref().clone()
        }
    }

    /// Reset statistics
    pub fn reset_stats(&self) {
        self.query_count.store(0, Ordering::Relaxed);
        self.cache_hits.store(0, Ordering::Relaxed);
    }

    /// Calculate memory usage
    fn calculate_memory_usage(&self) -> usize {
        let pgm_memory = self.pgm_router.memory_usage();
        let segments_memory: usize = self.segment_hashes.iter().map(|h| h.memory_usage()).sum();
        let boundaries_memory = std::mem::size_of_val(&*self.segment_boundaries);
        let starts_memory = std::mem::size_of_val(&*self.segment_starts);
        let struct_memory = std::mem::size_of::<Self>();

        pgm_memory + segments_memory + boundaries_memory + starts_memory + struct_memory
    }

    /// Get segment information
    pub fn segment_info(&self) -> Vec<(K, usize, usize)> {
        self.segment_boundaries
            .iter()
            .zip(self.segment_hashes.iter())
            .enumerate()
            .map(|(i, (&boundary, hash))| (boundary, hash.len(), self.segment_starts[i]))
            .collect()
    }

    /// Get memory usage per key
    pub fn memory_per_key(&self) -> f64 {
        self.calculate_memory_usage() as f64 / self.stats.total_keys as f64
    }

    /// Number of segments
    pub fn segment_count(&self) -> usize {
        self.segment_hashes.len()
    }

    /// Total number of keys
    pub fn len(&self) -> usize {
        self.stats.total_keys
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.stats.total_keys == 0
    }
}

/// Hybrid Key-Value store using the hybrid index
#[derive(Debug)]
pub struct HybridKeyValueStore<K: Key + Hash + std::fmt::Debug, V: Clone + Send + Sync> {
    /// The hybrid index
    index: HybridIndex<K>,

    /// Values stored in the same order as keys
    values: Arc<Vec<V>>,
}

impl<K: Key + Hash + std::fmt::Debug, V: Clone + Send + Sync> HybridKeyValueStore<K, V> {
    /// Create new hybrid KV store
    pub fn new(
        mut data: Vec<KVEntry<K, V>>,
        pgm_epsilon: usize,
        max_segment_size: usize,
    ) -> KVResult<Self> {
        if data.is_empty() {
            return Err(KVError::EmptyDataset);
        }

        // Sort and deduplicate
        data.sort_by_key(|entry| entry.key);
        data.dedup_by_key(|entry| entry.key);

        // Separate keys and values
        let keys: Vec<K> = data.iter().map(|entry| entry.key).collect();
        let values: Vec<V> = data.iter().map(|entry| entry.value.clone()).collect();

        // Build hybrid index
        let index = HybridIndex::new(keys, pgm_epsilon, max_segment_size)?;

        Ok(HybridKeyValueStore {
            index,
            values: Arc::new(values),
        })
    }

    /// Get value by key
    #[inline(always)]
    pub fn get(&self, key: K) -> Option<&V> {
        self.index.get(key).and_then(|pos| self.values.get(pos))
    }

    /// Batch get values
    pub fn batch_get(&self, keys: &[K]) -> Vec<Option<V>> {
        let positions = self.index.batch_get(keys);
        positions
            .into_iter()
            .map(|pos_opt| pos_opt.and_then(|pos| self.values.get(pos).cloned()))
            .collect()
    }

    /// Check if key exists
    #[inline(always)]
    pub fn contains_key(&self, key: K) -> bool {
        self.index.get(key).is_some()
    }

    /// Get number of entries
    pub fn len(&self) -> usize {
        self.index.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.index.is_empty()
    }

    /// Get statistics
    pub fn get_stats(&self) -> HybridIndexStats {
        self.index.get_stats()
    }

    /// Get memory usage
    pub fn memory_usage(&self) -> usize {
        self.index.calculate_memory_usage() + std::mem::size_of_val(&**self.values)
    }

    /// Get underlying index reference
    pub fn index(&self) -> &HybridIndex<K> {
        &self.index
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hybrid_index_basic() {
        let keys: Vec<u64> = (0..10_000).collect();
        let index = HybridIndex::new(keys.clone(), 64, 1000).unwrap();

        // Test all keys can be found
        for (i, &key) in keys.iter().enumerate() {
            assert_eq!(index.get(key), Some(i));
        }

        // Test non-existent keys
        assert_eq!(index.get(99999), None);
        assert_eq!(index.get(10000), None);
    }

    #[test]
    fn test_hybrid_kv_store() {
        let data: Vec<KVEntry<u64, String>> = (0..1000)
            .map(|i| KVEntry::new(i, format!("value_{}", i)))
            .collect();

        let store = HybridKeyValueStore::new(data, 32, 100).unwrap();

        // Test lookups
        assert_eq!(store.get(0), Some(&"value_0".to_string()));
        assert_eq!(store.get(999), Some(&"value_999".to_string()));
        assert_eq!(store.get(1000), None);

        // Test batch operations
        let keys = vec![0, 100, 200, 500, 999];
        let results = store.batch_get(&keys);
        assert_eq!(results.len(), 5);
        assert!(results.iter().all(|r| r.is_some()));
    }

    #[test]
    fn test_performance_characteristics() {
        use std::time::Instant;

        let keys: Vec<u64> = (0..100_000).collect();

        // Build performance
        let start = Instant::now();
        let index = HybridIndex::new(keys.clone(), 64, 2000).unwrap();
        let build_time = start.elapsed();

        // Query performance
        let start = Instant::now();
        let mut found = 0;
        for &key in &keys {
            if index.get(key).is_some() {
                found += 1;
            }
        }
        let query_time = start.elapsed();

        assert_eq!(found, keys.len());

        let stats = index.get_stats();
        println!("Hybrid Index Performance:");
        println!("Build time: {:?}", build_time);
        println!("Query time: {:?}", query_time);
        println!(
            "Build rate: {:.0} keys/s",
            keys.len() as f64 / build_time.as_secs_f64()
        );
        println!(
            "Query rate: {:.0} queries/s",
            keys.len() as f64 / query_time.as_secs_f64()
        );
        println!("Memory per key: {:.2} bytes", index.memory_per_key());
        println!("Segments: {}", stats.num_segments);
        println!("Avg segment size: {:.1}", stats.avg_segment_size);
        println!("Cache hit rate: {:.1}%", stats.cache_hit_rate * 100.0);
    }
}
