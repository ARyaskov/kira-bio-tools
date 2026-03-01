use std::fs::File;
use std::io::{BufWriter, Write};
use std::mem;
use std::path::Path;
use std::slice;

use kira_kv_engine::Index;
use memmap2::{Mmap, MmapOptions};
use rayon::prelude::*;

use crate::kbi::structs::{KbiHeader, KbiStats, Result};
use crate::util::{ChrId, GenomicKey};

pub struct KbiIndex {
    index: Index,
    keys: Vec<u64>,
    offsets: Vec<u64>,
}

impl KbiIndex {
    #[inline]
    pub fn get(&self, key: GenomicKey) -> Option<u64> {
        let idx = match self.index.lookup_u64(key.as_u64()) {
            Ok(v) => v,
            Err(_) => return None,
        };

        if idx < self.keys.len() && self.keys[idx] == key.as_u64() {
            Some(self.offsets[idx])
        } else {
            None
        }
    }

    #[inline]
    pub fn contains(&self, key: GenomicKey) -> bool {
        self.get(key).is_some()
    }

    pub fn len(&self) -> usize {
        self.keys.len()
    }

    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    pub fn get_batch(&self, keys: &[GenomicKey]) -> Vec<Option<u64>> {
        keys.par_iter().map(|&k| self.get(k)).collect()
    }

    pub fn range(&self, chr: ChrId, start_pos: u32, end_pos: u32) -> Vec<(u32, u64)> {
        let start_key = GenomicKey::new(chr, start_pos).as_u64();
        let end_key = GenomicKey::new(chr, end_pos).as_u64();

        let start_idx = self.keys.partition_point(|&k| k < start_key);
        let end_idx = self.keys.partition_point(|&k| k <= end_key);

        self.keys[start_idx..end_idx]
            .iter()
            .zip(&self.offsets[start_idx..end_idx])
            .map(|(&k, &off)| (GenomicKey::from_u64(k).position(), off))
            .collect()
    }

    pub fn has_range(&self, chr: ChrId, start_pos: u32, end_pos: u32) -> bool {
        let start_key = GenomicKey::new(chr, start_pos).as_u64();
        let end_key = GenomicKey::new(chr, end_pos).as_u64();

        let start_idx = self.keys.partition_point(|&k| k < start_key);
        let end_idx = self.keys.partition_point(|&k| k <= end_key);
        start_idx < end_idx
    }

    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::with_capacity(8 * 1024 * 1024, file);

        let index_bytes = self.index.serialize()?;
        let header = KbiHeader::new(self.keys.len(), index_bytes.len());

        writer.write_all(bytemuck::bytes_of(&header))?;

        writer.write_all(&index_bytes)?;

        let keys_bytes = unsafe {
            slice::from_raw_parts(
                self.keys.as_ptr() as *const u8,
                self.keys.len() * mem::size_of::<u64>(),
            )
        };
        writer.write_all(keys_bytes)?;

        let offsets_bytes = unsafe {
            slice::from_raw_parts(
                self.offsets.as_ptr() as *const u8,
                self.offsets.len() * mem::size_of::<u64>(),
            )
        };
        writer.write_all(offsets_bytes)?;

        writer.flush()?;
        Ok(())
    }

    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let mmap = unsafe {
            let file = File::open(path)?;
            MmapOptions::new().map(&file)?
        };

        Self::from_mmap(&mmap)
    }

    pub fn from_mmap(mmap: &Mmap) -> Result<Self> {
        if mmap.len() < mem::size_of::<KbiHeader>() {
            return Err(crate::kbi::structs::KbiError::InvalidFormat(
                "File too small".into(),
            ));
        }

        let header: &KbiHeader = bytemuck::from_bytes(&mmap[..mem::size_of::<KbiHeader>()]);
        header.validate()?;

        let n = header.n_entries as usize;
        let index_start = header.off_index as usize;
        let index_end = index_start + header.index_len as usize;
        if index_end > mmap.len() {
            return Err(crate::kbi::structs::KbiError::InvalidFormat(
                "Index out of range".into(),
            ));
        }
        let index = Index::deserialize(&mmap[index_start..index_end])?;

        let keys_start = header.off_keys as usize;
        let keys_end = keys_start + n * mem::size_of::<u64>();
        let keys: Vec<u64> = bytemuck::cast_slice(&mmap[keys_start..keys_end]).to_vec();

        let offsets_start = header.off_offsets as usize;
        let offsets_end = offsets_start + n * mem::size_of::<u64>();
        let offsets: Vec<u64> = bytemuck::cast_slice(&mmap[offsets_start..offsets_end]).to_vec();

        Ok(Self {
            index,
            keys,
            offsets,
        })
    }

    pub fn memory_usage(&self) -> usize {
        mem::size_of::<Self>()
            + self.keys.len() * mem::size_of::<u64>()
            + self.offsets.len() * mem::size_of::<u64>()
            + self.index.stats().total_memory
    }

    pub fn bytes_per_key(&self) -> f64 {
        self.memory_usage() as f64 / self.len().max(1) as f64
    }

    pub(crate) fn from_parts(index: Index, keys: Vec<u64>, offsets: Vec<u64>) -> Self {
        Self {
            index,
            keys,
            offsets,
        }
    }

    pub fn stats(&self, file_size: u64) -> KbiStats {
        KbiStats {
            entries: self.len(),
            memory_bytes: self.memory_usage(),
            bytes_per_key: self.bytes_per_key(),
            file_size,
        }
    }
}
