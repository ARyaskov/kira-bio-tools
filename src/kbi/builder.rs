use kira_kv_engine::{IndexBuilder, IndexConfig};
use rayon::prelude::*;

use crate::kbi::index::KbiIndex;
use crate::kbi::structs::{KbiError, Result};
use crate::util::GenomicKey;
use crate::vcf::VcfRecord;

pub struct KbiBuilder {
    entries: Vec<(u64, u64)>,
    gamma: f64,
    // Renamed in kira_kv_engine 0.6.0: BuildConfig::rehash_limit -> max_rehash.
    // Kept locally as max_rehash so we map 1:1 to the underlying field.
    max_rehash: u32,
    salt: u64,
}

impl KbiBuilder {
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
            gamma: 1.27,
            max_rehash: 16,
            salt: 0xC0FF_EE00_D15E_A5E,
        }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            entries: Vec::with_capacity(capacity),
            ..Self::new()
        }
    }

    pub fn gamma(mut self, gamma: f64) -> Self {
        self.gamma = gamma;
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

        let mut config = IndexConfig::default();
        config.mph_config.gamma = self.gamma;
        // kira_kv_engine 0.6.0: BuildConfig field renamed rehash_limit -> max_rehash.
        config.mph_config.max_rehash = self.max_rehash;
        config.mph_config.seed = self.salt;
        config.auto_detect_numeric = true;

        let index = IndexBuilder::new()
            .with_config(config)
            .build_index(key_bytes.clone())?;

        let n = keys.len();
        let sorted_keys = vec![0u64; n];
        let sorted_offsets = vec![0u64; n];

        keys.par_iter()
            .zip(offsets.par_iter())
            .for_each(|(&key, &offset)| {
                let idx = index.lookup_u64(key).unwrap();
                unsafe {
                    let keys_ptr = sorted_keys.as_ptr() as *mut u64;
                    let offsets_ptr = sorted_offsets.as_ptr() as *mut u64;
                    *keys_ptr.add(idx) = key;
                    *offsets_ptr.add(idx) = offset;
                }
            });

        Ok(KbiIndex::from_parts(index, sorted_keys, sorted_offsets))
    }
}

impl Default for KbiBuilder {
    fn default() -> Self {
        Self::new()
    }
}
