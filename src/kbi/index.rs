use std::fs::File;
use std::io::{BufWriter, Write};
use std::mem;
use std::path::Path;

use kira_kv_engine::Index;
use memmap2::{Mmap, MmapOptions};
use rayon::prelude::*;

use crate::kbi::structs::{KbiError, KbiHeader, KbiStats, Result};
use crate::util::{ChrId, GenomicKey};
use crate::vcf::header::ContigDict;

/// Point index over `(contig id, position)`. Keys and offsets are kept in
/// sorted order (so range scans are binary searches); the minimal perfect
/// hash maps a key to a slot, and `slot_to_idx` maps that slot back into the
/// sorted arrays.
pub struct KbiIndex {
    index: Index,
    keys: Vec<u64>,
    offsets: Vec<u64>,
    slot_to_idx: Vec<u32>,
    contigs: ContigDict,
}

impl KbiIndex {
    #[inline]
    pub fn get(&self, key: GenomicKey) -> Option<u64> {
        let slot = match self.index.lookup_u64(key.as_u64()) {
            Ok(v) => v,
            Err(_) => return None,
        };
        let idx = *self.slot_to_idx.get(slot)? as usize;
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

    pub fn contigs(&self) -> &ContigDict {
        &self.contigs
    }

    pub fn get_batch(&self, keys: &[GenomicKey]) -> Vec<Option<u64>> {
        keys.par_iter().map(|&k| self.get(k)).collect()
    }

    /// Records with `start_pos <= pos <= end_pos` on contig `chr`, in position
    /// order, as `(pos, file offset)` pairs.
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

    pub fn range_by_name(&self, chrom: &str, start_pos: u32, end_pos: u32) -> Vec<(u32, u64)> {
        match self.contigs.id(chrom) {
            Some(id) => self.range(id, start_pos, end_pos),
            None => Vec::new(),
        }
    }

    pub fn has_range(&self, chr: ChrId, start_pos: u32, end_pos: u32) -> bool {
        let start_key = GenomicKey::new(chr, start_pos).as_u64();
        let end_key = GenomicKey::new(chr, end_pos).as_u64();
        let start_idx = self.keys.partition_point(|&k| k < start_key);
        let end_idx = self.keys.partition_point(|&k| k <= end_key);
        start_idx < end_idx
    }

    /// Contig names that have at least one indexed record, in dictionary order.
    pub fn contigs_with_records(&self) -> Vec<String> {
        let mut out = Vec::new();
        for (id, name, _) in self.contigs.iter() {
            if self.has_range(id, 0, u32::MAX) {
                out.push(name.to_string());
            }
        }
        out
    }

    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::with_capacity(8 * 1024 * 1024, file);

        let index_bytes = self.index.serialize()?;
        let header = KbiHeader::new(self.keys.len(), index_bytes.len(), self.slot_to_idx.len());

        writer.write_all(bytemuck::bytes_of(&header))?;
        writer.write_all(&index_bytes)?;
        writer.write_all(bytemuck::cast_slice(&self.keys))?;
        writer.write_all(bytemuck::cast_slice(&self.offsets))?;
        writer.write_all(bytemuck::cast_slice(&self.slot_to_idx))?;

        let names = self.contigs.names();
        writer.write_all(&(names.len() as u32).to_le_bytes())?;
        for n in names {
            let b = n.as_bytes();
            let len = u16::try_from(b.len())
                .map_err(|_| KbiError::InvalidFormat(format!("contig name too long: {n}")))?;
            writer.write_all(&len.to_le_bytes())?;
            writer.write_all(b)?;
        }

        writer.flush()?;
        Ok(())
    }

    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        // SAFETY: read-only mapping of an index file that is only ever
        // replaced atomically, never modified in place while mapped.
        let mmap = unsafe { MmapOptions::new().map(&file)? };

        Self::from_mmap(&mmap)
    }

    pub fn from_mmap(mmap: &Mmap) -> Result<Self> {
        if mmap.len() < mem::size_of::<KbiHeader>() {
            return Err(KbiError::InvalidFormat("File too small".into()));
        }

        let header: KbiHeader = *bytemuck::from_bytes(&mmap[..mem::size_of::<KbiHeader>()]);
        header.validate()?;

        let n = header.n_entries as usize;
        let n_slots = header.n_slots as usize;
        let index_start = header.off_index as usize;
        let index_end = index_start + header.index_len as usize;
        if index_end > mmap.len() {
            return Err(KbiError::InvalidFormat("Index out of range".into()));
        }
        let index = Index::deserialize(&mmap[index_start..index_end])?;

        let keys_start = header.off_keys as usize;
        let keys_end = keys_start + n * mem::size_of::<u64>();
        let offsets_start = header.off_offsets as usize;
        let offsets_end = offsets_start + n * mem::size_of::<u64>();
        let slots_start = header.off_slots as usize;
        let slots_end = slots_start + n_slots * mem::size_of::<u32>();
        let names_start = header.off_names as usize;
        if keys_end > mmap.len() || offsets_end > mmap.len() || slots_end > mmap.len() || names_start + 4 > mmap.len() {
            return Err(KbiError::InvalidFormat("Section out of range".into()));
        }
        // Sections are not necessarily 8-aligned inside the file: decode by value.
        let keys: Vec<u64> = mmap[keys_start..keys_end]
            .chunks_exact(8)
            .map(|c| u64::from_le_bytes(c.try_into().unwrap()))
            .collect();
        let offsets: Vec<u64> = mmap[offsets_start..offsets_end]
            .chunks_exact(8)
            .map(|c| u64::from_le_bytes(c.try_into().unwrap()))
            .collect();
        let slot_to_idx: Vec<u32> = mmap[slots_start..slots_end]
            .chunks_exact(4)
            .map(|c| u32::from_le_bytes(c.try_into().unwrap()))
            .collect();

        let mut contigs = ContigDict::new();
        let mut p = names_start;
        let n_names = u32::from_le_bytes(mmap[p..p + 4].try_into().unwrap()) as usize;
        p += 4;
        for _ in 0..n_names {
            if p + 2 > mmap.len() {
                return Err(KbiError::InvalidFormat("Contig table truncated".into()));
            }
            let len = u16::from_le_bytes([mmap[p], mmap[p + 1]]) as usize;
            p += 2;
            if p + len > mmap.len() {
                return Err(KbiError::InvalidFormat("Contig table truncated".into()));
            }
            let name = std::str::from_utf8(&mmap[p..p + len])
                .map_err(|_| KbiError::InvalidFormat("Contig name not UTF-8".into()))?;
            contigs.insert(name);
            p += len;
        }

        Ok(Self { index, keys, offsets, slot_to_idx, contigs })
    }

    pub fn memory_usage(&self) -> usize {
        mem::size_of::<Self>()
            + self.keys.len() * mem::size_of::<u64>()
            + self.offsets.len() * mem::size_of::<u64>()
            + self.slot_to_idx.len() * mem::size_of::<u32>()
            + self.index.stats().total_memory
    }

    pub fn bytes_per_key(&self) -> f64 {
        self.memory_usage() as f64 / self.len().max(1) as f64
    }

    pub(crate) fn from_parts(
        index: Index,
        keys: Vec<u64>,
        offsets: Vec<u64>,
        slot_to_idx: Vec<u32>,
        contigs: ContigDict,
    ) -> Self {
        Self { index, keys, offsets, slot_to_idx, contigs }
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
