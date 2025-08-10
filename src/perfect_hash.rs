//! # Minimal Perfect Hash Implementation
//!
//! Ultra-fast minimal perfect hashing using CHD (Compress, Hash, Displace) algorithm
//! optimized for small to medium sized key sets (up to ~100K keys per hash table).

use crate::pgm_index::Key;
use serde::{Deserialize, Serialize};
use std::hash::Hash;
use std::marker::PhantomData;

/// Fast hash function for internal use
#[inline(always)]
fn fast_hash(key: u64, seed: u32) -> u32 {
    // MurmurHash3-inspired fast hash
    let mut h = key.wrapping_mul(0xc6a4a7935bd1e995u64);
    h ^= h >> 47;
    h = h.wrapping_mul(0xc6a4a7935bd1e995u64);
    h ^= seed as u64;
    h ^= h >> 47;
    (h as u32).wrapping_mul(0x5bd1e995)
}

/// Minimal Perfect Hash Table for small key sets
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MinimalPerfectHash<K: Key> {
    /// Displacement table for resolving conflicts
    displacements: Vec<u32>,
    /// Hash seeds for different levels
    seeds: [u32; 3],
    /// Original keys for verification
    keys: Vec<K>,
    /// Bucket size (usually 4-8)
    bucket_size: usize,
    /// Number of buckets
    num_buckets: usize,
}

impl<K: Key + Hash> MinimalPerfectHash<K> {
    /// Build minimal perfect hash for given keys
    pub fn new(keys: Vec<K>) -> Result<Self, &'static str> {
        if keys.is_empty() {
            return Err("Cannot build perfect hash for empty key set");
        }

        if keys.len() > 100_000 {
            return Err("Key set too large for minimal perfect hash");
        }

        let bucket_size = 6; // Optimal for most cases
        let num_buckets = (keys.len() + bucket_size - 1) / bucket_size;

        // Try different seed combinations
        for seed1 in 1..=100 {
            for seed2 in 1..=100 {
                for seed3 in 1..=100 {
                    let seeds = [seed1, seed2, seed3];

                    if let Ok(displacements) =
                        Self::try_build(&keys, bucket_size, num_buckets, seeds)
                    {
                        return Ok(MinimalPerfectHash {
                            displacements,
                            seeds,
                            keys: keys.clone(),
                            bucket_size,
                            num_buckets,
                        });
                    }
                }
            }
        }

        Err("Failed to build perfect hash after trying all seed combinations")
    }

    /// Try to build with specific seeds
    fn try_build(
        keys: &[K],
        _bucket_size: usize,
        num_buckets: usize,
        seeds: [u32; 3],
    ) -> Result<Vec<u32>, ()> {
        let mut buckets: Vec<Vec<(K, usize)>> = vec![Vec::new(); num_buckets];

        // Distribute keys into buckets
        for (idx, &key) in keys.iter().enumerate() {
            let bucket_id = Self::hash_to_bucket(key, seeds[0], num_buckets);
            buckets[bucket_id].push((key, idx));
        }

        // Sort buckets by size (largest first) for better success rate
        let mut bucket_indices: Vec<usize> = (0..num_buckets).collect();
        bucket_indices.sort_by_key(|&i| std::cmp::Reverse(buckets[i].len()));

        let mut displacements = vec![0u32; num_buckets];
        let mut used_positions = vec![false; keys.len()];

        // Process each bucket
        for &bucket_idx in &bucket_indices {
            let bucket = &buckets[bucket_idx];
            if bucket.is_empty() {
                continue;
            }

            // Find displacement that avoids conflicts
            let mut displacement = 0u32;
            'displacement_loop: loop {
                let mut positions = Vec::new();

                // Check if this displacement works
                for &(key, _) in bucket {
                    let pos = Self::hash_with_displacement(
                        key,
                        seeds[1],
                        seeds[2],
                        displacement,
                        keys.len(),
                    );
                    if used_positions[pos] {
                        displacement += 1;
                        if displacement > 1000 {
                            return Err(()); // Too many attempts
                        }
                        continue 'displacement_loop;
                    }
                    positions.push(pos);
                }

                // Mark positions as used
                for pos in positions {
                    used_positions[pos] = true;
                }

                displacements[bucket_idx] = displacement;
                break;
            }
        }

        Ok(displacements)
    }

    /// Hash key to bucket
    #[inline(always)]
    fn hash_to_bucket(key: K, seed: u32, num_buckets: usize) -> usize {
        let key_u64 = key.to_u64().unwrap_or(0);
        let hash = fast_hash(key_u64, seed);
        (hash as usize) % num_buckets
    }

    /// Hash with displacement
    #[inline(always)]
    fn hash_with_displacement(
        key: K,
        seed1: u32,
        seed2: u32,
        displacement: u32,
        table_size: usize,
    ) -> usize {
        let key_u64 = key.to_u64().unwrap_or(0);
        let hash1 = fast_hash(key_u64, seed1.wrapping_add(displacement));
        let hash2 = fast_hash(key_u64, seed2.wrapping_add(displacement));
        let combined = hash1.wrapping_add(hash2);
        (combined as usize) % table_size
    }

    /// Get position of key (returns None if key not found)
    #[inline(always)]
    pub fn get_position(&self, key: K) -> Option<usize> {
        let bucket_id = Self::hash_to_bucket(key, self.seeds[0], self.num_buckets);
        let displacement = self.displacements[bucket_id];
        let position = Self::hash_with_displacement(
            key,
            self.seeds[1],
            self.seeds[2],
            displacement,
            self.keys.len(),
        );

        // Verify the key matches (perfect hash should guarantee this)
        if position < self.keys.len() && self.keys[position] == key {
            Some(position)
        } else {
            None
        }
    }

    /// Check if key exists
    #[inline(always)]
    pub fn contains(&self, key: K) -> bool {
        self.get_position(key).is_some()
    }

    /// Get number of keys
    pub fn len(&self) -> usize {
        self.keys.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    /// Get memory usage in bytes
    pub fn memory_usage(&self) -> usize {
        std::mem::size_of_val(&*self.displacements)
            + std::mem::size_of_val(&*self.keys)
            + std::mem::size_of::<Self>()
    }

    /// Get keys slice
    pub fn keys(&self) -> &[K] {
        &self.keys
    }
}

/// Batch perfect hash builder for multiple small key sets
pub struct BatchPerfectHashBuilder {
    max_segment_size: usize,
    _phantom: PhantomData<u8>, // Placeholder since we removed the K parameter
}

impl BatchPerfectHashBuilder {
    pub fn new(max_segment_size: usize) -> Self {
        Self {
            max_segment_size: max_segment_size.min(50_000), // Reasonable limit
            _phantom: PhantomData,
        }
    }

    /// Build multiple perfect hashes for segments
    pub fn build_segments<K: Key + Hash>(
        &self,
        segments: Vec<Vec<K>>,
    ) -> Result<Vec<MinimalPerfectHash<K>>, String> {
        use rayon::prelude::*;

        // Parallel construction of perfect hashes
        segments
            .into_par_iter()
            .map(|keys| {
                if keys.len() > self.max_segment_size {
                    return Err(format!(
                        "Segment too large: {} > {}",
                        keys.len(),
                        self.max_segment_size
                    ));
                }

                MinimalPerfectHash::new(keys)
                    .map_err(|e| format!("Perfect hash construction failed: {}", e))
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_minimal_perfect_hash() {
        let keys: Vec<u64> = vec![1, 5, 10, 15, 20, 25, 30];
        let mph = MinimalPerfectHash::new(keys.clone()).unwrap();

        // Test all keys can be found
        for &key in &keys {
            assert!(mph.contains(key));
            let pos = mph.get_position(key).unwrap();
            assert_eq!(mph.keys()[pos], key);
        }

        // Test non-existent keys
        assert!(!mph.contains(999));
        assert!(mph.get_position(999).is_none());
    }

    #[test]
    fn test_perfect_hash_performance() {
        use std::time::Instant;

        // Test with 10K keys
        let keys: Vec<u64> = (0..10_000).collect();

        let start = Instant::now();
        let mph = MinimalPerfectHash::new(keys.clone()).unwrap();
        let build_time = start.elapsed();

        println!("Perfect hash build time for 10K keys: {:?}", build_time);
        println!(
            "Memory usage: {} bytes ({:.2} bytes/key)",
            mph.memory_usage(),
            mph.memory_usage() as f64 / keys.len() as f64
        );

        // Test lookup performance
        let start = Instant::now();
        let mut found = 0;
        for &key in &keys {
            if mph.contains(key) {
                found += 1;
            }
        }
        let lookup_time = start.elapsed();

        assert_eq!(found, keys.len());
        println!("Lookup time for 10K queries: {:?}", lookup_time);
        println!(
            "Lookup rate: {:.0} queries/sec",
            keys.len() as f64 / lookup_time.as_secs_f64()
        );
    }

    #[test]
    fn test_batch_builder() {
        let builder = BatchPerfectHashBuilder::new(1000);

        let segments = vec![
            vec![1u64, 2, 3, 4, 5],
            vec![10, 20, 30, 40, 50],
            vec![100, 200, 300],
        ];

        let perfect_hashes = builder.build_segments(segments).unwrap();
        assert_eq!(perfect_hashes.len(), 3);

        // Test first segment
        let mph = &perfect_hashes[0];
        assert!(mph.contains(1));
        assert!(mph.contains(5));
        assert!(!mph.contains(10));
    }
}
