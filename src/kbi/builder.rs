use rayon::prelude::*;

use crate::kbi::index::KbiIndex;
use crate::kbi::structs::{KbiError, Result};
use crate::util::GenomicKey;
use crate::vcf::VcfRecord;
use crate::vcf::header::ContigDict;

pub struct KbiBuilder {
    entries: Vec<(u64, u64)>,
    contigs: ContigDict,
}

impl KbiBuilder {
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
            contigs: ContigDict::new(),
        }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            entries: Vec::with_capacity(capacity),
            ..Self::new()
        }
    }

    #[deprecated(
        since = "0.8.1",
        note = "no effect: the MPH output range is always minimal"
    )]
    pub fn gamma(self, _gamma: f64) -> Self {
        self
    }

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
        self.entries.par_sort_by_key(|(key, _)| *key);
        let original_len = self.entries.len();
        self.entries.dedup_by_key(|entry| entry.0);
        let dup_count = original_len - self.entries.len();
        if dup_count > 0 {
            eprintln!(
                "KBI: {} records share a position with an earlier record and are not point-addressable",
                dup_count
            );
        }
        if self.entries.len() > u32::MAX as usize {
            return Err(KbiError::InvalidFormat("Too many entries".into()));
        }
        let (keys, offsets) = self.entries.into_iter().unzip();
        Ok(KbiIndex::from_parts(keys, offsets, self.contigs))
    }
}

impl Default for KbiBuilder {
    fn default() -> Self {
        Self::new()
    }
}
