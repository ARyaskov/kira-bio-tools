//! # Compressed storage implementations

use crate::pgm_index::{Key, PGMIndex};
use crate::types::{KVError, KVResult};
use serde::{de::DeserializeOwned, Serialize};
use std::sync::Arc;

/// KV store with compressed values
pub struct CompressedKVStore<K: Key> {
    index: PGMIndex<K>,
    compressed_values: Arc<Vec<Vec<u8>>>,
}

impl<K: Key> CompressedKVStore<K> {
    /// Create compressed KV store
    pub fn new<V: Serialize>(data: Vec<(K, V)>, epsilon: usize) -> KVResult<Self> {
        if data.is_empty() {
            return Err(KVError::EmptyDataset);
        }

        let mut sorted_data = data;
        sorted_data.sort_by_key(|(k, _)| *k);

        let keys: Vec<K> = sorted_data.iter().map(|(k, _)| *k).collect();
        let compressed_values: Result<Vec<Vec<u8>>, KVError> =
            sorted_data.iter().map(|(_, v)| compress_value(v)).collect();

        let compressed_values = compressed_values?;
        let index = PGMIndex::new(keys, epsilon);

        Ok(Self {
            index,
            compressed_values: Arc::new(compressed_values),
        })
    }

    /// Get decompressed value
    pub fn get_decompressed<V: DeserializeOwned>(&self, key: K) -> Option<V> {
        self.index
            .get(key)
            .and_then(|pos| self.compressed_values.get(pos))
            .and_then(|compressed| decompress_value(compressed).ok())
    }

    /// Batch get decompressed values
    pub fn batch_get_decompressed<V: DeserializeOwned>(&self, keys: &[K]) -> Vec<Option<V>> {
        let positions = self.index.batch_get(keys);
        positions
            .into_iter()
            .map(|pos_opt| {
                pos_opt
                    .and_then(|pos| self.compressed_values.get(pos))
                    .and_then(|compressed| decompress_value(compressed).ok())
            })
            .collect()
    }

    /// Get compressed data size
    pub fn compressed_size(&self) -> usize {
        self.compressed_values.iter().map(|v| v.len()).sum()
    }

    /// Get compression ratio
    pub fn compression_ratio(&self, original_size: usize) -> f64 {
        if original_size == 0 {
            0.0
        } else {
            self.compressed_size() as f64 / original_size as f64
        }
    }

    /// Get number of entries
    pub fn len(&self) -> usize {
        self.compressed_values.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.compressed_values.is_empty()
    }

    /// Get memory usage
    pub fn memory_usage(&self) -> usize {
        self.index.memory_usage()
            + std::mem::size_of_val(&**self.compressed_values)
            + self.compressed_size()
            + std::mem::size_of::<Self>()
    }
}

// Compression functions with conditional compilation
fn compress_value<V: Serialize>(value: &V) -> KVResult<Vec<u8>> {
    use flate2::{write::GzEncoder, Compression};
    use std::io::Write;

    let serialized =
        bincode::serialize(value).map_err(|e| KVError::SerializationError(e.to_string()))?;

    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder
        .write_all(&serialized)
        .map_err(|e| KVError::IoError(e))?;

    encoder.finish().map_err(|e| KVError::IoError(e))
}

fn decompress_value<V: DeserializeOwned>(compressed: &[u8]) -> KVResult<V> {
    use flate2::read::GzDecoder;
    use std::io::Read;

    let mut decoder = GzDecoder::new(compressed);
    let mut decompressed = Vec::new();
    decoder
        .read_to_end(&mut decompressed)
        .map_err(|e| KVError::IoError(e))?;

    bincode::deserialize(&decompressed).map_err(|e| KVError::SerializationError(e.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compressed_store() {
        let data = vec![
            (1u64, "hello world".repeat(100)),
            (2u64, "test data".repeat(100)),
        ];

        let store = CompressedKVStore::new(data, 32).unwrap();

        let value: Option<String> = store.get_decompressed(1);
        assert!(value.is_some());
        assert_eq!(value.unwrap(), "hello world".repeat(100));
    }
}
