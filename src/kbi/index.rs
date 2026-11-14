use std::fs::File;
use std::io::{BufWriter, Write};
use std::mem;
use std::path::Path;

use kira_kv_engine::{EliasFano, Index};
use memmap2::{Mmap, MmapOptions};
use rayon::prelude::*;

use crate::kbi::structs::{KbiError, KbiHeader, KbiStats, Result};
use crate::util::{ChrId, GenomicKey};
use crate::vcf::header::ContigDict;

enum KbiKeys {
    Plain(Vec<u64>),
    Compact(EliasFano),
}

impl KbiKeys {
    fn from_sorted(keys: Vec<u64>) -> Self {
        if let Some(compact) = EliasFano::from_sorted(&keys)
            && compact.memory_usage() < keys.len() * mem::size_of::<u64>()
        {
            return Self::Compact(compact);
        }
        Self::Plain(keys)
    }

    fn len(&self) -> usize {
        match self {
            Self::Plain(keys) => keys.len(),
            Self::Compact(keys) => keys.len(),
        }
    }

    fn get(&self, index: usize) -> u64 {
        match self {
            Self::Plain(keys) => keys[index],
            Self::Compact(keys) => keys.get(index),
        }
    }

    fn lower_bound(&self, target: u64) -> usize {
        if let Self::Plain(keys) = self {
            return keys.partition_point(|&key| key < target);
        }
        let mut start = 0;
        let mut end = self.len();
        while start < end {
            let middle = start + (end - start) / 2;
            if self.get(middle) < target {
                start = middle + 1;
            } else {
                end = middle;
            }
        }
        start
    }

    fn bounds(&self, start: u64, end: u64) -> std::ops::Range<usize> {
        if start > end {
            return 0..0;
        }
        let from = self.lower_bound(start);
        let to = end
            .checked_add(1)
            .map_or(self.len(), |next| self.lower_bound(next));
        from..to
    }

    fn memory_usage(&self) -> usize {
        match self {
            Self::Plain(keys) => keys.capacity() * mem::size_of::<u64>(),
            Self::Compact(keys) => keys.memory_usage() - mem::size_of::<EliasFano>(),
        }
    }

    fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();
        match self {
            Self::Plain(keys) => {
                bytes.push(0);
                for key in keys {
                    bytes.extend_from_slice(&key.to_le_bytes());
                }
            }
            Self::Compact(keys) => {
                bytes.push(1);
                keys.write_to(&mut bytes);
            }
        }
        bytes
    }

    fn from_bytes(bytes: &[u8], n: usize) -> Result<Self> {
        match bytes.first() {
            Some(0) if bytes.len().checked_sub(1) == n.checked_mul(8) => {
                let keys: Vec<_> = bytes[1..]
                    .as_chunks::<8>()
                    .0
                    .iter()
                    .map(|chunk| u64::from_le_bytes(*chunk))
                    .collect();
                if keys.windows(2).any(|pair| pair[0] >= pair[1]) {
                    return Err(KbiError::InvalidFormat(
                        "Keys are not strictly ordered".into(),
                    ));
                }
                Ok(Self::from_sorted(keys))
            }
            Some(1) => {
                let mut position = 1;
                let compact = EliasFano::read_from(bytes, &mut position)
                    .ok_or_else(|| KbiError::InvalidFormat("Invalid compact keys".into()))?;
                if position != bytes.len() || compact.len() != n {
                    return Err(KbiError::InvalidFormat(
                        "Compact key count or size mismatch".into(),
                    ));
                }
                if compact.memory_usage() < n * mem::size_of::<u64>() {
                    Ok(Self::Compact(compact))
                } else {
                    let mut keys = Vec::with_capacity(n);
                    compact.materialize_range(0, n, &mut keys);
                    Ok(Self::Plain(keys))
                }
            }
            _ => Err(KbiError::InvalidFormat(
                "Invalid key encoding or size".into(),
            )),
        }
    }
}

pub struct KbiIndex {
    keys: KbiKeys,
    offsets: Vec<u64>,
    contigs: ContigDict,
}

impl KbiIndex {
    #[inline]
    pub fn get(&self, key: GenomicKey) -> Option<u64> {
        let key = key.as_u64();
        let index = self.keys.lower_bound(key);
        if index < self.len() && self.keys.get(index) == key {
            Some(self.offsets[index])
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
        self.len() == 0
    }

    pub fn contigs(&self) -> &ContigDict {
        &self.contigs
    }

    pub fn get_batch(&self, keys: &[GenomicKey]) -> Vec<Option<u64>> {
        keys.par_iter().map(|&key| self.get(key)).collect()
    }

    pub fn range(&self, chr: ChrId, start_pos: u32, end_pos: u32) -> Vec<(u32, u64)> {
        let bounds = self.keys.bounds(
            GenomicKey::new(chr, start_pos).as_u64(),
            GenomicKey::new(chr, end_pos).as_u64(),
        );
        match &self.keys {
            KbiKeys::Plain(keys) => keys[bounds.clone()]
                .iter()
                .zip(&self.offsets[bounds])
                .map(|(&key, &offset)| (GenomicKey::from_u64(key).position(), offset))
                .collect(),
            KbiKeys::Compact(keys) => {
                let mut out = Vec::with_capacity(bounds.len());
                let mut buffer = Vec::with_capacity(bounds.len().min(1024));
                for from in bounds.clone().step_by(1024) {
                    let count = (bounds.end - from).min(1024);
                    keys.materialize_range(from, count, &mut buffer);
                    out.extend(
                        buffer
                            .iter()
                            .zip(&self.offsets[from..from + count])
                            .map(|(&key, &offset)| (GenomicKey::from_u64(key).position(), offset)),
                    );
                }
                out
            }
        }
    }

    pub fn range_by_name(&self, chrom: &str, start_pos: u32, end_pos: u32) -> Vec<(u32, u64)> {
        match self.contigs.id(chrom) {
            Some(id) => self.range(id, start_pos, end_pos),
            None => Vec::new(),
        }
    }

    pub fn has_range(&self, chr: ChrId, start_pos: u32, end_pos: u32) -> bool {
        !self
            .keys
            .bounds(
                GenomicKey::new(chr, start_pos).as_u64(),
                GenomicKey::new(chr, end_pos).as_u64(),
            )
            .is_empty()
    }

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
        let key_bytes = self.keys.to_bytes();
        let header = KbiHeader::new_compact(self.len(), key_bytes.len())?;
        let names = self.contigs.names();
        let count = u32::try_from(names.len())
            .map_err(|_| KbiError::InvalidFormat("Too many contigs".into()))?;
        let mut name_bytes = count.to_le_bytes().to_vec();
        for name in names {
            let len = u16::try_from(name.len())
                .map_err(|_| KbiError::InvalidFormat(format!("Contig name too long: {name}")))?;
            name_bytes.extend_from_slice(&len.to_le_bytes());
            name_bytes.extend_from_slice(name.as_bytes());
        }
        header
            .off_names
            .checked_add(
                u64::try_from(name_bytes.len())
                    .map_err(|_| KbiError::InvalidFormat("Contig table too large".into()))?,
            )
            .ok_or_else(|| KbiError::InvalidFormat("File size overflow".into()))?;
        let file = File::create(path)?;
        let mut writer = BufWriter::with_capacity(8 * 1024 * 1024, file);
        writer.write_all(bytemuck::bytes_of(&header))?;
        writer.write_all(&key_bytes)?;
        for offset in &self.offsets {
            writer.write_all(&offset.to_le_bytes())?;
        }
        writer.write_all(&name_bytes)?;
        writer.flush()?;
        Ok(())
    }

    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        Self::from_mmap(&mmap)
    }

    pub fn from_mmap(mmap: &Mmap) -> Result<Self> {
        let header_bytes = mmap
            .get(..mem::size_of::<KbiHeader>())
            .ok_or_else(|| KbiError::InvalidFormat("File too small".into()))?;
        let header: KbiHeader = bytemuck::pod_read_unaligned(header_bytes);
        header.validate()?;
        let n = usize::try_from(header.n_entries)
            .map_err(|_| KbiError::InvalidFormat("Entry count exceeds address space".into()))?;
        let values_size = header
            .n_entries
            .checked_mul(8)
            .ok_or_else(|| KbiError::InvalidFormat("Entry size overflow".into()))?;
        let index_bytes = section(mmap, header.off_index, header.index_len)?;
        let offset_bytes = section(mmap, header.off_offsets, values_size)?;
        let names_start = usize::try_from(header.off_names)
            .map_err(|_| KbiError::InvalidFormat("Contig offset exceeds address space".into()))?;
        let name_bytes = mmap
            .get(names_start..)
            .ok_or_else(|| KbiError::InvalidFormat("Contig table out of range".into()))?;
        let contigs = read_contigs(name_bytes)?;
        let keys = if header.version == 4 {
            let key_bytes = section(mmap, header.off_keys, values_size)?;
            let keys: Vec<_> = key_bytes
                .as_chunks::<8>()
                .0
                .iter()
                .map(|chunk| u64::from_le_bytes(*chunk))
                .collect();
            if keys.windows(2).any(|pair| pair[0] >= pair[1]) {
                return Err(KbiError::InvalidFormat(
                    "Keys are not strictly ordered".into(),
                ));
            }
            let slots_size = header
                .n_slots
                .checked_mul(4)
                .ok_or_else(|| KbiError::InvalidFormat("Slot size overflow".into()))?;
            let slot_bytes = section(mmap, header.off_slots, slots_size)?;
            let slots: Vec<_> = slot_bytes
                .as_chunks::<4>()
                .0
                .iter()
                .map(|chunk| u32::from_le_bytes(*chunk))
                .collect();
            let index = Index::from_bytes(index_bytes)?;
            if index.len() != n
                || slots.len() > index.slot_capacity()
                || slots.last() == Some(&u32::MAX)
                || slots.iter().filter(|&&rank| rank != u32::MAX).count() != n
                || slots
                    .iter()
                    .any(|&rank| rank != u32::MAX && rank as usize >= n)
            {
                return Err(KbiError::InvalidFormat(
                    "Invalid legacy slot metadata".into(),
                ));
            }
            for (rank, &key) in keys.iter().enumerate() {
                let slot = index
                    .lookup_u64(key)
                    .map_err(|_| KbiError::InvalidFormat("Legacy key lookup failed".into()))?;
                if slots.get(slot).copied() != Some(rank as u32) {
                    return Err(KbiError::InvalidFormat(
                        "Legacy slot mapping mismatch".into(),
                    ));
                }
            }
            KbiKeys::from_sorted(keys)
        } else {
            KbiKeys::from_bytes(index_bytes, n)?
        };
        let offsets = offset_bytes
            .as_chunks::<8>()
            .0
            .iter()
            .map(|chunk| u64::from_le_bytes(*chunk))
            .collect();
        Ok(Self {
            keys,
            offsets,
            contigs,
        })
    }

    pub fn memory_usage(&self) -> usize {
        mem::size_of::<Self>()
            + self.keys.memory_usage()
            + self.offsets.capacity() * mem::size_of::<u64>()
    }

    pub fn bytes_per_key(&self) -> f64 {
        self.memory_usage() as f64 / self.len().max(1) as f64
    }

    pub(crate) fn from_parts(keys: Vec<u64>, offsets: Vec<u64>, contigs: ContigDict) -> Self {
        Self {
            keys: KbiKeys::from_sorted(keys),
            offsets,
            contigs,
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

fn section(bytes: &[u8], start: u64, len: u64) -> Result<&[u8]> {
    let end = start
        .checked_add(len)
        .ok_or_else(|| KbiError::InvalidFormat("Section size overflow".into()))?;
    let start = usize::try_from(start)
        .map_err(|_| KbiError::InvalidFormat("Section offset exceeds address space".into()))?;
    let end = usize::try_from(end)
        .map_err(|_| KbiError::InvalidFormat("Section end exceeds address space".into()))?;
    bytes
        .get(start..end)
        .ok_or_else(|| KbiError::InvalidFormat("Section out of range".into()))
}

fn read_contigs(bytes: &[u8]) -> Result<ContigDict> {
    let count_bytes = bytes
        .get(..4)
        .ok_or_else(|| KbiError::InvalidFormat("Contig table truncated".into()))?;
    let count = u32::from_le_bytes(count_bytes.try_into().unwrap()) as usize;
    if count > (bytes.len() - 4) / 2 {
        return Err(KbiError::InvalidFormat("Invalid contig count".into()));
    }
    let mut position = 4usize;
    let mut contigs = ContigDict::new();
    for id in 0..count {
        let length_end = position
            .checked_add(2)
            .ok_or_else(|| KbiError::InvalidFormat("Contig length overflow".into()))?;
        let length_bytes = bytes
            .get(position..length_end)
            .ok_or_else(|| KbiError::InvalidFormat("Contig table truncated".into()))?;
        let len = u16::from_le_bytes(length_bytes.try_into().unwrap()) as usize;
        let end = length_end
            .checked_add(len)
            .ok_or_else(|| KbiError::InvalidFormat("Contig name size overflow".into()))?;
        let name_bytes = bytes
            .get(length_end..end)
            .ok_or_else(|| KbiError::InvalidFormat("Contig table truncated".into()))?;
        let name = std::str::from_utf8(name_bytes)
            .map_err(|_| KbiError::InvalidFormat("Contig name not UTF-8".into()))?;
        if contigs.insert(name) as usize != id {
            return Err(KbiError::InvalidFormat("Duplicate contig name".into()));
        }
        position = end;
    }
    if position != bytes.len() {
        return Err(KbiError::InvalidFormat("Trailing contig data".into()));
    }
    Ok(contigs)
}

#[cfg(test)]
#[path = "../../tests/unit/kbi_index.rs"]
mod tests;
