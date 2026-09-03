use kira_kv_engine::{IndexBuilder, IndexConfig};
use rayon::prelude::*;

use crate::kbi::index::KbiIndex;
use crate::kbi::structs::{KbiError, Result};
use crate::util::GenomicKey;
use crate::vcf::VcfRecord;
use crate::vcf::header::ContigDict;

pub struct KbiBuilder {
    entries: Vec<(u64, u64)>,
    gamma: f64,
    max_rehash: u32,
    salt: u64,
    contigs: ContigDict,
}

impl KbiBuilder {
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
            gamma: 1.27,
            max_rehash: 16,
            salt: 0xC0FF_EE00_D15E_A5E,
            contigs: ContigDict::new(),
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

    /// Contig dictionary the keys' chr ids refer to; stored in the index so
    /// queries can be made by name.
    pub fn set_contigs(&mut self, contigs: ContigDict) {
        self.contigs = contigs;
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

        // Keep the first record at each (contig, pos); later ones share the key.
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
            eprintln!("KBI: {} records share a position with an earlier record and are not point-addressable", dup_count);
        }

        let keys: Vec<u64> = deduped.iter().map(|(k, _)| *k).collect();
        let offsets: Vec<u64> = deduped.iter().map(|(_, o)| *o).collect();

        let key_bytes: Vec<[u8; 8]> = keys.iter().map(|k| k.to_le_bytes()).collect();

        let mut config = IndexConfig::default();
        config.mph_config.gamma = self.gamma;
        config.mph_config.max_rehash = self.max_rehash;
        config.mph_config.seed = self.salt;
        config.auto_detect_numeric = true;

        let index = IndexBuilder::new()
            .with_config(config)
            .build_index(key_bytes)?;

        // Keys stay sorted; the MPH slot of each key is recorded so point
        // lookups can jump straight to the sorted position.
        let n = keys.len();
        let slots: Vec<usize> = keys
            .par_iter()
            .map(|&key| index.lookup_u64(key).unwrap_or(usize::MAX))
            .collect();
        let n_slots = slots.iter().copied().filter(|&s| s != usize::MAX).max().map(|m| m + 1).unwrap_or(0);
        let mut slot_to_idx = vec![u32::MAX; n_slots];
        for (i, &slot) in slots.iter().enumerate() {
            if slot == usize::MAX || slot >= n_slots {
                return Err(KbiError::InvalidFormat("MPH lookup failed during build".into()));
            }
            slot_to_idx[slot] = i as u32;
        }
        let _ = n;

        Ok(KbiIndex::from_parts(index, keys, offsets, slot_to_idx, self.contigs))
    }
}

impl Default for KbiBuilder {
    fn default() -> Self {
        Self::new()
    }
}
