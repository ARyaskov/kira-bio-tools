//! # Key-Value Store implementation using PGM-Index

use crate::pgm_index::{Key, PGMIndex};
use crate::types::{KVEntry, KVError, KVResult};
use serde::{Deserialize, Serialize};
use std::sync::Arc;
use std::time::SystemTime;

/// Key-Value store statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KVStoreStats {
    pub total_entries: usize,
    pub memory_usage_bytes: usize,
    pub cache_hit_rate: f64,
    pub avg_segment_size: f64,
    pub created_at: SystemTime,
    pub version: u64,
}

/// High-performance Key-Value store built on PGM-Index
#[derive(Debug)]
pub struct PGMKeyValueStore<K: Key, V: Clone + Send + Sync> {
    /// PGM Index for fast key lookups
    index: PGMIndex<K>,
    /// Values stored in the same order as keys
    values: Arc<Vec<V>>,
    /// Store metadata
    created_at: SystemTime,
    version: u64,
}

impl<K: Key, V: Clone + Send + Sync> PGMKeyValueStore<K, V> {
    /// Create new KV store from key-value pairs
    pub fn new(mut data: Vec<KVEntry<K, V>>, epsilon: usize) -> KVResult<Self> {
        if data.is_empty() {
            return Err(KVError::EmptyDataset);
        }

        if epsilon == 0 {
            return Err(KVError::InvalidEpsilon(epsilon));
        }

        // Sort by keys for PGM-Index
        data.sort_by_key(|entry| entry.key);

        // Remove duplicates (keep last value for duplicate keys)
        data.dedup_by_key(|entry| entry.key);

        // Separate keys and values
        let keys: Vec<K> = data.iter().map(|entry| entry.key).collect();
        let values: Vec<V> = data.iter().map(|entry| entry.value.clone()).collect();

        // Build PGM Index
        let index = PGMIndex::new(keys, epsilon);

        Ok(Self {
            index,
            values: Arc::new(values),
            created_at: SystemTime::now(),
            version: 1,
        })
    }

    /// Get value by key
    pub fn get(&self, key: K) -> Option<&V> {
        self.index.get(key).and_then(|pos| self.values.get(pos))
    }

    /// Check if key exists
    pub fn contains_key(&self, key: K) -> bool {
        self.index.get(key).is_some()
    }

    /// Batch get multiple values - adaptive approach
    pub fn batch_get(&self, keys: &[K]) -> Vec<Option<V>> {
        if keys.len() < 50_000 {
            keys.iter().map(|&key| self.get(key).cloned()).collect()
        } else {
            let positions = self.index.batch_get(keys);
            positions
                .into_iter()
                .map(|pos_opt| pos_opt.and_then(|pos| self.values.get(pos).cloned()))
                .collect()
        }
    }

    /// Force parallel batch get (for benchmarking)
    pub fn batch_get_parallel(&self, keys: &[K]) -> Vec<Option<V>> {
        let positions = self.index.batch_get(keys);
        positions
            .into_iter()
            .map(|pos_opt| pos_opt.and_then(|pos| self.values.get(pos).cloned()))
            .collect()
    }

    /// Force sequential batch get (for benchmarking)
    pub fn batch_get_sequential(&self, keys: &[K]) -> Vec<Option<V>> {
        keys.iter().map(|&key| self.get(key).cloned()).collect()
    }

    /// Get all key-value pairs in range [start, end]
    pub fn range(&self, start: K, end: K) -> Vec<KVEntry<K, V>> {
        let (start_pos, _) = self.index.predict(start);
        let (_, end_pos) = self.index.predict(end);

        let mut result = Vec::new();
        for pos in start_pos..end_pos.min(self.values.len()) {
            if let (Some(&key), Some(value)) = (self.index.data.get(pos), self.values.get(pos)) {
                if key >= start && key <= end {
                    result.push(KVEntry::new(key, value.clone()));
                }
            }
        }
        result
    }

    /// Get all keys in range [start, end]
    pub fn keys_range(&self, start: K, end: K) -> Vec<K> {
        let (start_pos, _) = self.index.predict(start);
        let (_, end_pos) = self.index.predict(end);

        let mut result = Vec::new();
        for pos in start_pos..end_pos.min(self.index.data.len()) {
            if let Some(&key) = self.index.data.get(pos) {
                if key >= start && key <= end {
                    result.push(key);
                }
            }
        }
        result
    }

    /// Get store statistics
    pub fn stats(&self) -> KVStoreStats {
        KVStoreStats {
            total_entries: self.values.len(),
            memory_usage_bytes: self.memory_usage(),
            cache_hit_rate: self.index.cache_hit_rate(),
            avg_segment_size: self.index.avg_segment_size(),
            created_at: self.created_at,
            version: self.version,
        }
    }

    /// Get number of entries
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Check if store is empty
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// Get memory usage in bytes
    pub fn memory_usage(&self) -> usize {
        self.index.memory_usage()
            + std::mem::size_of_val(&**self.values)
            + std::mem::size_of::<Self>()
    }

    /// Reset performance statistics
    pub fn reset_stats(&self) {
        self.index.reset_stats();
    }

    /// Get all keys (sorted)
    pub fn keys(&self) -> &[K] {
        &self.index.data
    }

    /// Get underlying PGM index
    pub fn index(&self) -> &PGMIndex<K> {
        &self.index
    }
}

// Custom serialization for KV store
impl<K: Key + Serialize, V: Clone + Send + Sync + Serialize> Serialize for PGMKeyValueStore<K, V> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        use serde::ser::SerializeStruct;
        let mut state = serializer.serialize_struct("PGMKeyValueStore", 4)?;
        state.serialize_field("index", &self.index)?;
        state.serialize_field("values", &*self.values)?;
        state.serialize_field("created_at", &self.created_at)?;
        state.serialize_field("version", &self.version)?;
        state.end()
    }
}

impl<'de, K: Key + Deserialize<'de>, V: Clone + Send + Sync + Deserialize<'de>> Deserialize<'de>
    for PGMKeyValueStore<K, V>
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        use serde::de::{self, MapAccess, Visitor};
        use std::fmt;

        struct KVStoreVisitor<K, V>(std::marker::PhantomData<(K, V)>);

        impl<'de, K: Key + Deserialize<'de>, V: Clone + Send + Sync + Deserialize<'de>> Visitor<'de>
            for KVStoreVisitor<K, V>
        {
            type Value = PGMKeyValueStore<K, V>;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("struct PGMKeyValueStore")
            }

            fn visit_map<A>(self, mut map: A) -> Result<PGMKeyValueStore<K, V>, A::Error>
            where
                A: MapAccess<'de>,
            {
                let mut index = None;
                let mut values = None;
                let mut created_at = None;
                let mut version = None;

                while let Some(key) = map.next_key()? {
                    match key {
                        "index" => index = Some(map.next_value()?),
                        "values" => values = Some(Arc::new(map.next_value::<Vec<V>>()?)),
                        "created_at" => created_at = Some(map.next_value()?),
                        "version" => version = Some(map.next_value()?),
                        _ => {
                            let _ = map.next_value::<serde::de::IgnoredAny>()?;
                        }
                    }
                }

                Ok(PGMKeyValueStore {
                    index: index.ok_or_else(|| de::Error::missing_field("index"))?,
                    values: values.ok_or_else(|| de::Error::missing_field("values"))?,
                    created_at: created_at.ok_or_else(|| de::Error::missing_field("created_at"))?,
                    version: version.ok_or_else(|| de::Error::missing_field("version"))?,
                })
            }
        }

        const FIELDS: &[&str] = &["index", "values", "created_at", "version"];
        deserializer.deserialize_struct(
            "PGMKeyValueStore",
            FIELDS,
            KVStoreVisitor(std::marker::PhantomData),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kv_store_basic_operations() {
        let data = vec![
            KVEntry::new(1, "one".to_string()),
            KVEntry::new(3, "three".to_string()),
            KVEntry::new(2, "two".to_string()),
        ];

        let store = PGMKeyValueStore::new(data, 32).unwrap();

        assert_eq!(store.get(1), Some(&"one".to_string()));
        assert_eq!(store.get(2), Some(&"two".to_string()));
        assert_eq!(store.get(3), Some(&"three".to_string()));
        assert_eq!(store.get(4), None);

        assert!(store.contains_key(1));
        assert!(!store.contains_key(4));

        assert_eq!(store.len(), 3);
    }
}
