use anyhow::{Result, anyhow};
use bytemuck::{self, Pod, Zeroable};
#[cfg(feature = "gpu")]
use cust::memory::DeviceCopy;
use kira_kv_engine::Index;
use libdeflater::Decompressor;
use memchr::memchr;
use memmap2::{Mmap, MmapOptions};
use std::cell::RefCell;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::File;
use std::mem;
use std::ops::Deref;
use std::path::Path;
use std::sync::{Arc, Mutex};

use super::header::{
    ANI_HEADER_END, ANI_MAGIC, ANI_STR_NONE, ANI_VERSION, AniBlockEntry, AniEntry, AniEntryV2,
    AniHeader, AniHeaderV3, AniHeaderV6,
};
use crate::annotate::structs::bundle::{
    AnnotationBundle, FieldNumber, FieldType, StructuredInfoField,
};

const BLOCK_CACHE_CAP: usize = 16;
const CSTR_CACHE_CAP: usize = 256;

pub enum StringStorage {
    Owned(Vec<u8>),
    Mmap { offset: usize, len: usize },
}

impl StringStorage {
    pub fn as_slice<'a>(&'a self, mmap: &'a Mmap) -> &'a [u8] {
        match self {
            StringStorage::Owned(v) => v.as_slice(),
            StringStorage::Mmap { offset, len } => &mmap[*offset..(*offset + *len)],
        }
    }

    pub fn len(&self) -> usize {
        match self {
            StringStorage::Owned(v) => v.len(),
            StringStorage::Mmap { len, .. } => *len,
        }
    }
}

struct BlockCacheEntry {
    idx: usize,
    data: Arc<Vec<u8>>,
}

struct BlockCache {
    entries: Vec<BlockCacheEntry>,
}

impl BlockCache {
    fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    fn get(&mut self, idx: usize) -> Option<Arc<Vec<u8>>> {
        if let Some(pos) = self.entries.iter().position(|e| e.idx == idx) {
            let entry = self.entries.remove(pos);
            let data = entry.data.clone();
            self.entries.push(entry);
            return Some(data);
        }
        None
    }

    fn insert(&mut self, idx: usize, data: Arc<Vec<u8>>) {
        self.entries.push(BlockCacheEntry { idx, data });
        if self.entries.len() > BLOCK_CACHE_CAP {
            self.entries.remove(0);
        }
    }
}

enum StringSource {
    Raw(StringStorage),
    Compressed {
        blocks: Vec<AniBlockEntry>,
        cache: Mutex<BlockCache>,
        total_len: usize,
    },
}

#[derive(Clone, Copy)]
struct CStrCacheEntry {
    offset: u64,
    end: usize,
    block_idx: Option<usize>,
}

struct CStrCache {
    entries: Vec<CStrCacheEntry>,
}

impl CStrCache {
    fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    fn get(&mut self, offset: u64) -> Option<CStrCacheEntry> {
        if let Some(pos) = self.entries.iter().position(|e| e.offset == offset) {
            let entry = self.entries.remove(pos);
            self.entries.push(entry);
            Some(entry)
        } else {
            None
        }
    }

    fn insert(&mut self, entry: CStrCacheEntry) {
        if self.entries.len() >= CSTR_CACHE_CAP {
            self.entries.remove(0);
        }
        self.entries.push(entry);
    }
}

thread_local! {
    static CSTR_CACHE: RefCell<CStrCache> = RefCell::new(CStrCache::new());
}

pub struct CStrRef<'a> {
    data: CStrData<'a>,
    start: usize,
    end: usize,
}

enum CStrData<'a> {
    Borrowed(&'a [u8]),
    Shared(Arc<Vec<u8>>),
}

impl<'a> CStrRef<'a> {
    pub fn empty() -> CStrRef<'a> {
        CStrRef {
            data: CStrData::Borrowed(&[]),
            start: 0,
            end: 0,
        }
    }

    pub fn as_str(&self) -> &str {
        match &self.data {
            CStrData::Borrowed(bytes) => unsafe {
                std::str::from_utf8_unchecked(&bytes[self.start..self.end])
            },
            CStrData::Shared(buf) => unsafe {
                std::str::from_utf8_unchecked(&buf[self.start..self.end])
            },
        }
    }
}

impl<'a> AsRef<str> for CStrRef<'a> {
    fn as_ref(&self) -> &str {
        self.as_str()
    }
}

impl<'a> Deref for CStrRef<'a> {
    type Target = str;

    fn deref(&self) -> &Self::Target {
        self.as_str()
    }
}

impl std::fmt::Display for CStrRef<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

#[derive(Clone, Copy)]
pub struct IntervalEntry {
    pub start: u32,
    pub end: u32,
    pub entry_idx: usize,
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct AniPosIndexHeader {
    pub contig_count: u32,
    pub block_count: u32,
    pub pos_count: u32,
    pub entry_index_count: u32,
    pub off_contigs: u32,
    pub off_blocks: u32,
    pub off_pos_offsets: u32,
    pub off_pos_counts: u32,
    pub off_entry_indices: u32,
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct AniPosContig {
    pub chr_id: u16,
    pub _pad: u16,
    pub min_pos: u32,
    pub max_pos: u32,
    pub block_start: u32,
    pub block_count: u32,
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct AniPosBlock {
    pub base_pos: u32,
    pub _pad: u32,
    pub masks: [u64; 8],
    pub offsets_start: u32,
    pub _pad2: u32,
}

pub struct AniPosIndex {
    pub contigs: Vec<AniPosContig>,
    pub blocks: Vec<AniPosBlock>,
    pub pos_offsets: Vec<u32>,
    pub pos_counts: Vec<u16>,
    pub entry_indices: Vec<u32>,
}

#[cfg(feature = "gpu")]
unsafe impl DeviceCopy for AniPosContig {}

#[cfg(feature = "gpu")]
unsafe impl DeviceCopy for AniPosBlock {}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct AniInfoBlobHeader {
    pub n_entries: u64,
    pub dict_count: u32,
    pub pair_count: u32,
    pub off_dict_offsets: u64,
    pub off_dict_data: u64,
    pub off_entry_offsets: u64,
    pub off_entry_counts: u64,
    pub off_pairs: u64,
    pub off_values: u64,
}

pub struct AniInfoBlob {
    pub dict_offsets: Vec<u32>,
    pub dict_data: Vec<u8>,
    pub dict_strings: Vec<String>,
    pub entry_offsets: Vec<u32>,
    pub entry_counts: Vec<u16>,
    pub pairs: Vec<AniInfoPair>,
    pub values: Vec<u8>,
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct AniInfoPair {
    pub tag_id: u32,
    pub value_off: u32,
    pub value_len: u32,
}

pub struct AniInfoCache {
    pub tag_types: Vec<FieldType>,
    pub entry_offsets: Vec<u32>,
    pub entry_counts: Vec<u16>,
    pub pair_tag_ids: Vec<u32>,
    pub pair_offsets: Vec<u32>,
    pub pair_counts: Vec<u16>,
    pub int_values: Vec<i32>,
    pub float_values: Vec<f32>,
    pub str_offsets: Vec<u32>,
    pub str_lens: Vec<u32>,
    pub str_data: Vec<u8>,
}

pub struct AniIndex {
    pub header: AniHeader,
    pub index: Index,
    pub entries: Vec<AniEntry>,
    intervals: Vec<Vec<IntervalEntry>>,
    strings: StringSource,
    mmap: Mmap,
    pos_index: Option<AniPosIndex>,
    info_blob: Option<AniInfoBlob>,
    info_cache: Option<AniInfoCache>,
}

impl AniIndex {
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };

        if mmap.len() < mem::size_of::<AniHeaderV3>() {
            return Err(anyhow!("ANI file too small"));
        }

        let header = read_header(&mmap)?;
        header.validate()?;

        let index = load_index(&mmap, &header)?;
        let entries = load_entries(&mmap, &header)?;
        let strings = load_strings(&mmap, &header)?;

        let (pos_index, info_blob) = if header.version == ANI_VERSION {
            (
                load_pos_index(&mmap, &header)?,
                load_info_blob(&mmap, &header)?,
            )
        } else {
            (None, None)
        };

        let mut ani = Self {
            header,
            index,
            entries,
            intervals: Vec::new(),
            strings,
            mmap,
            pos_index,
            info_blob,
            info_cache: None,
        };
        ani.build_interval_index();
        if let Some(blob) = &ani.info_blob {
            ani.info_cache = Some(build_info_cache(&ani, blob));
        }

        Ok(ani)
    }

    pub fn strings_len(&self) -> usize {
        match &self.strings {
            StringSource::Raw(storage) => storage.len(),
            StringSource::Compressed { total_len, .. } => *total_len,
        }
    }

    pub fn strings_slice(&self) -> &[u8] {
        match &self.strings {
            StringSource::Raw(storage) => storage.as_slice(&self.mmap),
            StringSource::Compressed { .. } => &[],
        }
    }

    pub fn strings_owned(&self) -> Vec<u8> {
        match &self.strings {
            StringSource::Raw(storage) => storage.as_slice(&self.mmap).to_vec(),
            StringSource::Compressed {
                blocks, total_len, ..
            } => {
                let mut out = Vec::with_capacity(*total_len);
                for (idx, block) in blocks.iter().enumerate() {
                    if let Ok(data) = self.decompress_block(idx, block) {
                        out.extend_from_slice(&data);
                    }
                }
                out
            }
        }
    }

    pub fn lookup_pos_index(&self, chr_id: u8, pos: u32) -> Option<&[u32]> {
        let pos_index = self.pos_index.as_ref()?;
        let contig = pos_index
            .contigs
            .iter()
            .find(|c| c.chr_id == chr_id as u16)?;
        if pos < contig.min_pos || pos > contig.max_pos {
            return None;
        }
        let start = contig.block_start as usize;
        let end = start + contig.block_count as usize;
        let blocks = &pos_index.blocks[start..end];
        let base = (pos / 512) * 512;
        let block_idx = match blocks.binary_search_by_key(&base, |b| b.base_pos) {
            Ok(v) => v,
            Err(_) => return None,
        };
        let block = &blocks[block_idx];
        let bit = (pos - base) as usize;
        let word = bit / 64;
        let bit_in_word = bit % 64;
        let mask = block.masks[word];
        if (mask >> bit_in_word) & 1 == 0 {
            return None;
        }
        let mut rank = 0usize;
        for w in 0..word {
            rank += block.masks[w].count_ones() as usize;
        }
        let lower_mask = if bit_in_word == 0 {
            0
        } else {
            (1u64 << bit_in_word) - 1
        };
        rank += (mask & lower_mask).count_ones() as usize;
        let pos_idx = block.offsets_start as usize + rank;
        let count = *pos_index.pos_counts.get(pos_idx)? as usize;
        let offset = *pos_index.pos_offsets.get(pos_idx)? as usize;
        let end = offset + count;
        pos_index.entry_indices.get(offset..end)
    }

    pub fn has_pos_index(&self) -> bool {
        self.pos_index.is_some()
    }

    pub fn has_info_blob(&self) -> bool {
        self.info_blob.is_some()
    }

    pub fn info_cache(&self) -> Option<&AniInfoCache> {
        self.info_cache.as_ref()
    }

    pub fn info_blob(&self) -> Option<&AniInfoBlob> {
        self.info_blob.as_ref()
    }

    pub fn pos_index(&self) -> Option<&AniPosIndex> {
        self.pos_index.as_ref()
    }

    pub fn read_cstring(&self, offset: usize) -> CStrRef<'_> {
        match &self.strings {
            StringSource::Raw(storage) => {
                let bytes = storage.as_slice(&self.mmap);
                if offset >= bytes.len() {
                    return CStrRef {
                        data: CStrData::Borrowed(&[]),
                        start: 0,
                        end: 0,
                    };
                }
                let offset_u64 = offset as u64;
                let cached = CSTR_CACHE.with(|c| c.borrow_mut().get(offset_u64));
                let end = if let Some(entry) = cached {
                    entry.end
                } else {
                    let end = find_cstr_end(bytes, offset);
                    CSTR_CACHE.with(|c| {
                        c.borrow_mut().insert(CStrCacheEntry {
                            offset: offset_u64,
                            end,
                            block_idx: None,
                        })
                    });
                    end
                };
                CStrRef {
                    data: CStrData::Borrowed(bytes),
                    start: offset,
                    end,
                }
            }
            StringSource::Compressed { blocks, cache, .. } => {
                let offset = offset as u64;
                let cached = CSTR_CACHE.with(|c| c.borrow_mut().get(offset));
                let cached_block_idx = cached.and_then(|e| e.block_idx);
                let idx = if let Some(v) = cached_block_idx {
                    v
                } else {
                    match find_block_index(blocks, offset) {
                        Some(v) => v,
                        None => {
                            return CStrRef {
                                data: CStrData::Borrowed(&[]),
                                start: 0,
                                end: 0,
                            };
                        }
                    }
                };
                let data = {
                    let mut guard = cache.lock().unwrap();
                    if let Some(v) = guard.get(idx) {
                        v
                    } else {
                        drop(guard);
                        let block = &blocks[idx];
                        let decompressed = match self.decompress_block(idx, block) {
                            Ok(v) => Arc::new(v),
                            Err(_) => {
                                return CStrRef {
                                    data: CStrData::Borrowed(&[]),
                                    start: 0,
                                    end: 0,
                                };
                            }
                        };
                        let mut guard = cache.lock().unwrap();
                        guard.insert(idx, decompressed.clone());
                        decompressed
                    }
                };

                let block = &blocks[idx];
                let start = (offset - block.raw_start) as usize;
                let end = if let Some(entry) = cached {
                    if entry.block_idx == Some(idx) {
                        entry.end
                    } else {
                        find_cstr_end(data.as_slice(), start)
                    }
                } else {
                    find_cstr_end(data.as_slice(), start)
                };
                if cached.is_none() || cached.map(|e| e.block_idx != Some(idx)).unwrap_or(false) {
                    CSTR_CACHE.with(|c| {
                        c.borrow_mut().insert(CStrCacheEntry {
                            offset,
                            end,
                            block_idx: Some(idx),
                        })
                    });
                }

                CStrRef {
                    data: CStrData::Shared(data),
                    start,
                    end,
                }
            }
        }
    }

    pub fn find_interval_entry(&self, chr_id: u8, pos: u32) -> Option<usize> {
        if self.intervals.is_empty() {
            return None;
        }
        let list = match self.intervals.get(chr_id as usize) {
            Some(v) => v,
            None => return None,
        };
        if list.is_empty() {
            return None;
        }

        let mut lo = 0usize;
        let mut hi = list.len();
        while lo < hi {
            let mid = (lo + hi) / 2;
            if list[mid].start > pos {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        if lo == 0 {
            return None;
        }

        let mut idx = lo - 1;
        loop {
            let ent = &list[idx];
            if pos < ent.start {
                break;
            }
            if pos <= ent.end {
                return Some(ent.entry_idx);
            }
            if idx == 0 {
                break;
            }
            idx -= 1;
        }

        None
    }

    fn decompress_block(&self, idx: usize, block: &AniBlockEntry) -> Result<Vec<u8>> {
        let start = block.data_off as usize;
        let end = start + block.data_len as usize;
        if end > self.mmap.len() {
            return Err(anyhow!("ANI block {} out of range", idx));
        }

        let compressed = &self.mmap[start..end];
        let mut out = vec![0u8; block.raw_len as usize];
        let mut decompressor = Decompressor::new();
        decompressor
            .deflate_decompress(compressed, &mut out)
            .map_err(|_| anyhow!("ANI block {} decompression failed", idx))?;
        Ok(out)
    }

    fn build_interval_index(&mut self) {
        if !self.has_interval_headers() {
            return;
        }

        let mut per_chr: Vec<Vec<IntervalEntry>> = vec![Vec::new(); 256];
        for (idx, entry) in self.entries.iter().enumerate() {
            let rf = self.read_cstring(entry.ref_ofs as usize);
            if rf.as_ref() != "." {
                continue;
            }
            let alt = self.read_cstring(entry.alt_ofs as usize);
            let end = match alt.as_ref().parse::<u32>() {
                Ok(v) => v,
                Err(_) => continue,
            };
            per_chr[entry.chr_id as usize].push(IntervalEntry {
                start: entry.pos,
                end,
                entry_idx: idx,
            });
        }

        for list in &mut per_chr {
            if list.len() > 1 {
                list.sort_by_key(|e| e.start);
            }
        }

        self.intervals = per_chr;
    }

    fn has_interval_headers(&self) -> bool {
        let mut idx = 0usize;
        while idx < self.strings_len() {
            let line = self.read_cstring(idx);
            let s = line.as_ref();
            if s.is_empty() {
                idx += 1;
                continue;
            }
            if s == "##KIRA_BT_ANI_INTERVALS" {
                return true;
            }
            if s == "##KIRA_BT_ANI_HEADER_END" {
                break;
            }
            idx += s.len() + 1;
        }
        false
    }

    pub fn info_fields_from_blob(
        &self,
        entry_idx: usize,
        field_meta: &HashMap<String, FieldNumber>,
    ) -> Vec<StructuredInfoField> {
        let Some(blob) = &self.info_blob else {
            return Vec::new();
        };
        if entry_idx >= blob.entry_offsets.len() {
            return Vec::new();
        }
        let offset = blob.entry_offsets[entry_idx] as usize;
        let count = blob.entry_counts[entry_idx] as usize;
        if offset + count > blob.pairs.len() {
            return Vec::new();
        }
        let mut fields = Vec::with_capacity(count);
        for pair in &blob.pairs[offset..offset + count] {
            let key = blob
                .dict_strings
                .get(pair.tag_id as usize)
                .cloned()
                .unwrap_or_default();
            if key.is_empty() {
                continue;
            }
            if pair.value_len == 0 {
                let number = field_meta.get(&key).copied().unwrap_or(FieldNumber::Zero);
                fields.push(StructuredInfoField {
                    key,
                    number,
                    ty: FieldType::Flag,
                    values: Vec::new(),
                });
                continue;
            }
            let start = pair.value_off as usize;
            let end = start + pair.value_len as usize;
            if end > blob.values.len() {
                continue;
            }
            let raw = std::str::from_utf8(&blob.values[start..end]).unwrap_or("");
            let values: Vec<String> = raw.split(',').map(|s| s.to_string()).collect();
            let number = field_meta.get(&key).copied().unwrap_or_else(|| {
                if values.len() == 1 {
                    FieldNumber::One
                } else {
                    FieldNumber::Many
                }
            });
            let ty = infer_info_type(&key);
            fields.push(StructuredInfoField {
                key,
                number,
                ty,
                values,
            });
        }
        fields
    }

    pub fn build_bundle_from_entry_idx_opts_with_meta(
        &self,
        entry_idx: usize,
        field_meta: &HashMap<String, FieldNumber>,
        need_info: bool,
        need_format: bool,
    ) -> AnnotationBundle {
        let e = &self.entries[entry_idx];
        let alt_str = self.read_cstring(e.alt_ofs as usize);
        let id_str = self.read_cstring(e.id_ofs as usize);
        let qual_str = self.read_cstring(e.qual_ofs as usize);
        let filter_str = self.read_cstring(e.filter_ofs as usize);
        let info = if need_info {
            if self.info_blob.is_some() {
                self.info_fields_from_blob(entry_idx, field_meta)
            } else {
                let info_str = self.read_cstring(e.info_ofs as usize);
                crate::annotate::structs::bundle::infer_structured_info_fields(
                    info_str.as_ref(),
                    field_meta,
                )
            }
        } else {
            Vec::new()
        };

        let (format_opt, samples) = if need_format && e.format_ofs != ANI_STR_NONE {
            let format_str = self.read_cstring(e.format_ofs as usize);
            let samples_str = if e.samples_ofs != ANI_STR_NONE {
                self.read_cstring(e.samples_ofs as usize)
            } else {
                CStrRef::empty()
            };
            let format_opt = parse_optional(format_str.as_ref());
            let samples = if format_opt.is_some() && !samples_str.as_ref().is_empty() {
                let s = samples_str.as_ref();
                let mut count = 1usize;
                for b in s.as_bytes() {
                    if *b == b'\t' {
                        count += 1;
                    }
                }
                let mut out = Vec::with_capacity(count);
                for v in s.split('\t') {
                    out.push(v.to_string());
                }
                out
            } else {
                Vec::new()
            };
            (format_opt, samples)
        } else {
            (None, Vec::new())
        };

        AnnotationBundle {
            alt: alt_str.to_string(),
            id: parse_optional(id_str.as_ref()),
            qual: parse_optional(qual_str.as_ref()),
            filter: parse_optional(filter_str.as_ref()),
            info,
            format_str: format_opt,
            format_samples: samples,
        }
    }
}

fn infer_info_type(key: &str) -> FieldType {
    if key.starts_with('I') && key.chars().nth(1).map_or(false, |c| c.is_uppercase()) {
        FieldType::Integer
    } else if key.starts_with('F') && key.chars().nth(1).map_or(false, |c| c.is_uppercase()) {
        FieldType::Float
    } else {
        FieldType::String
    }
}

fn build_info_cache(ani: &AniIndex, blob: &AniInfoBlob) -> AniInfoCache {
    let tag_types = build_tag_types(ani, blob);
    let entry_offsets = blob.entry_offsets.clone();
    let entry_counts = blob.entry_counts.clone();
    let pair_tag_ids: Vec<u32> = blob.pairs.iter().map(|p| p.tag_id).collect();
    let mut pair_offsets = Vec::with_capacity(blob.pairs.len());
    let mut pair_counts = Vec::with_capacity(blob.pairs.len());
    let mut int_values = Vec::new();
    let mut float_values = Vec::new();
    let mut str_offsets = Vec::new();
    let mut str_lens = Vec::new();
    let mut str_data = Vec::new();

    for pair in &blob.pairs {
        let tag_type = tag_types
            .get(pair.tag_id as usize)
            .copied()
            .unwrap_or(FieldType::String);
        if pair.value_len == 0 {
            pair_offsets.push(0);
            pair_counts.push(0);
            continue;
        }
        let start = pair.value_off as usize;
        let end = start + pair.value_len as usize;
        if end > blob.values.len() {
            pair_offsets.push(0);
            pair_counts.push(0);
            continue;
        }
        let raw = std::str::from_utf8(&blob.values[start..end]).unwrap_or("");
        match tag_type {
            FieldType::Flag => {
                pair_offsets.push(0);
                pair_counts.push(0);
            }
            FieldType::Integer => {
                let off = int_values.len() as u32;
                let mut count = 0u16;
                for token in raw.split(',') {
                    let v = if is_missing_value(token) {
                        i32::MIN
                    } else {
                        token.parse::<i32>().unwrap_or(i32::MIN)
                    };
                    int_values.push(v);
                    count = count.wrapping_add(1);
                }
                pair_offsets.push(off);
                pair_counts.push(count);
            }
            FieldType::Float => {
                let off = float_values.len() as u32;
                let mut count = 0u16;
                for token in raw.split(',') {
                    let v = if is_missing_value(token) {
                        f32::NAN
                    } else {
                        token.parse::<f32>().unwrap_or(f32::NAN)
                    };
                    float_values.push(v);
                    count = count.wrapping_add(1);
                }
                pair_offsets.push(off);
                pair_counts.push(count);
            }
            FieldType::String => {
                let off = str_offsets.len() as u32;
                let mut count = 0u16;
                for token in raw.split(',') {
                    if is_missing_value(token) {
                        str_offsets.push(str_data.len() as u32);
                        str_lens.push(0);
                    } else {
                        str_offsets.push(str_data.len() as u32);
                        str_lens.push(token.len() as u32);
                        str_data.extend_from_slice(token.as_bytes());
                    }
                    count = count.wrapping_add(1);
                }
                pair_offsets.push(off);
                pair_counts.push(count);
            }
        }
    }

    AniInfoCache {
        tag_types,
        entry_offsets,
        entry_counts,
        pair_tag_ids,
        pair_offsets,
        pair_counts,
        int_values,
        float_values,
        str_offsets,
        str_lens,
        str_data,
    }
}

fn build_tag_types(ani: &AniIndex, blob: &AniInfoBlob) -> Vec<FieldType> {
    let header_types = load_info_types_from_headers(ani);
    blob.dict_strings
        .iter()
        .map(|k| {
            header_types
                .get(k)
                .copied()
                .unwrap_or_else(|| infer_info_type(k))
        })
        .collect()
}

fn load_info_types_from_headers(ani: &AniIndex) -> HashMap<String, FieldType> {
    let mut out = HashMap::new();
    for line in iter_ani_header_lines_local(ani) {
        if !line.starts_with("##INFO=") {
            continue;
        }
        let key = extract_info_key_local(&line);
        let ty = extract_info_type_local(&line);
        if let (Some(k), Some(t)) = (key, ty) {
            out.insert(k, t);
        }
    }
    out
}

fn extract_info_key_local(line: &str) -> Option<String> {
    if let Some(start) = line.find("ID=") {
        let rest = &line[start + 3..];
        if let Some(end) = rest.find(',') {
            return Some(rest[..end].to_string());
        }
    }
    None
}

fn extract_info_type_local(line: &str) -> Option<FieldType> {
    if let Some(start) = line.find("Type=") {
        let rest = &line[start + 5..];
        let end = rest.find(',').unwrap_or(rest.len());
        let ty = &rest[..end];
        return match ty {
            "Integer" => Some(FieldType::Integer),
            "Float" => Some(FieldType::Float),
            "String" => Some(FieldType::String),
            "Flag" => Some(FieldType::Flag),
            _ => None,
        };
    }
    None
}

fn iter_ani_header_lines_local(ani: &AniIndex) -> Vec<String> {
    let mut headers = Vec::new();
    let mut saw_header = false;
    let mut idx = 0usize;
    while idx < ani.strings_len() {
        let line_ref = ani.read_cstring(idx);
        let line = line_ref.as_ref();
        if line.is_empty() {
            idx += 1;
            continue;
        }
        if line == ANI_HEADER_END {
            break;
        }
        let is_header = line.starts_with("##INFO=")
            || line.starts_with("##FORMAT=")
            || line.starts_with("##FILTER=")
            || line.starts_with("#CHROM");
        if is_header {
            headers.push(line.to_string());
            saw_header = true;
        } else if saw_header {
            break;
        }
        idx += line.len() + 1;
    }
    headers
}

fn is_missing_value(val: &str) -> bool {
    val.is_empty() || val == "."
}

fn parse_optional(s: &str) -> Option<String> {
    if s.is_empty() || s == "." {
        None
    } else {
        Some(s.to_string())
    }
}

fn read_pod_slice<T: Pod>(bytes: &[u8]) -> Result<Vec<T>> {
    let size = mem::size_of::<T>();
    if size == 0 {
        return Ok(Vec::new());
    }
    if bytes.len() % size != 0 {
        return Err(anyhow!("ANI blob size misaligned"));
    }
    let mut out = Vec::with_capacity(bytes.len() / size);
    let mut i = 0usize;
    while i < bytes.len() {
        let v = bytemuck::pod_read_unaligned(&bytes[i..i + size]);
        out.push(v);
        i += size;
    }
    Ok(out)
}

fn read_header(mmap: &Mmap) -> Result<AniHeader> {
    let h3: AniHeaderV3 = unsafe { *(mmap.as_ptr() as *const AniHeaderV3) };
    if h3.magic != ANI_MAGIC {
        return Err(anyhow!("Bad ANI magic"));
    }

    match h3.version {
        ANI_VERSION => {
            if mmap.len() < mem::size_of::<AniHeaderV6>() {
                return Err(anyhow!("ANI file too small for v6 header"));
            }
            let h6: AniHeaderV6 = unsafe { *(mmap.as_ptr() as *const AniHeaderV6) };
            Ok(AniHeader {
                magic: h6.magic,
                version: h6.version,
                n_entries: h6.n_entries,
                index_len: h6.index_len,
                off_index: h6.off_index,
                off_entries: h6.off_entries,
                off_strings: h6.off_strings,
                off_block_index: h6.off_block_index,
                n_blocks: h6.n_blocks,
                block_size: h6.block_size,
                off_pos_index: h6.off_pos_index,
                pos_index_len: h6.pos_index_len,
                off_blob: h6.off_blob,
                blob_len: h6.blob_len,
            })
        }
        2 | 3 | 4 | 5 => Err(anyhow!(
            "Unsupported ANI version {} (rebuild index)",
            h3.version
        )),
        _ => Err(anyhow!("Unsupported ANI version {}", h3.version)),
    }
}

fn load_index(mmap: &Mmap, header: &AniHeader) -> Result<Index> {
    let start = header.off_index as usize;
    let end = start + header.index_len as usize;
    if end > mmap.len() {
        return Err(anyhow!("ANI index out of range"));
    }
    Index::deserialize(&mmap[start..end]).map_err(|e| anyhow!(e))
}

fn load_entries(mmap: &Mmap, header: &AniHeader) -> Result<Vec<AniEntry>> {
    let ent_start = header.off_entries as usize;
    match header.version {
        2 => load_entries_v2(mmap, ent_start, header.n_entries as usize),
        3 | ANI_VERSION => load_entries_v3(mmap, ent_start, header.n_entries as usize),
        _ => Err(anyhow!("Unsupported ANI version {}", header.version)),
    }
}

fn load_entries_v2(mmap: &Mmap, ent_start: usize, n_entries: usize) -> Result<Vec<AniEntry>> {
    let ent_size = mem::size_of::<AniEntryV2>();
    let ent_end = ent_start + n_entries * ent_size;
    let mut entries = Vec::with_capacity(n_entries);

    for chunk in mmap[ent_start..ent_end].chunks_exact(ent_size) {
        let e: AniEntryV2 = unsafe { *(chunk.as_ptr() as *const AniEntryV2) };
        entries.push(AniEntry {
            chr_id: e.chr_id,
            pos: e.pos,
            ref_ofs: e.ref_ofs,
            alt_ofs: e.alt_ofs,
            id_ofs: e.id_ofs,
            qual_ofs: e.qual_ofs,
            filter_ofs: e.filter_ofs,
            info_ofs: e.info_ofs,
            info_len: e.info_len,
            format_ofs: ANI_STR_NONE,
            samples_ofs: ANI_STR_NONE,
        });
    }

    Ok(entries)
}

fn load_entries_v3(mmap: &Mmap, ent_start: usize, n_entries: usize) -> Result<Vec<AniEntry>> {
    let ent_size = mem::size_of::<AniEntry>();
    let ent_end = ent_start + n_entries * ent_size;
    let mut entries = Vec::with_capacity(n_entries);

    for chunk in mmap[ent_start..ent_end].chunks_exact(ent_size) {
        let e: AniEntry = unsafe { *(chunk.as_ptr() as *const AniEntry) };
        entries.push(e);
    }

    Ok(entries)
}

fn load_strings(mmap: &Mmap, header: &AniHeader) -> Result<StringSource> {
    if header.version == ANI_VERSION {
        let start = header.off_block_index as usize;
        let count = header.n_blocks as usize;
        let ent_size = mem::size_of::<AniBlockEntry>();
        let end = start + count * ent_size;
        if end > mmap.len() {
            return Err(anyhow!("ANI block index out of range"));
        }

        let mut blocks = Vec::with_capacity(count);
        for chunk in mmap[start..end].chunks_exact(ent_size) {
            let e: AniBlockEntry = unsafe { *(chunk.as_ptr() as *const AniBlockEntry) };
            blocks.push(e);
        }

        let total_len = blocks
            .last()
            .map(|b| (b.raw_start + b.raw_len as u64) as usize)
            .unwrap_or(0);

        Ok(StringSource::Compressed {
            blocks,
            cache: Mutex::new(BlockCache::new()),
            total_len,
        })
    } else {
        let str_start = header.off_strings as usize;
        let len = mmap.len().saturating_sub(str_start);
        Ok(StringSource::Raw(StringStorage::Mmap {
            offset: str_start,
            len,
        }))
    }
}

fn load_pos_index(mmap: &Mmap, header: &AniHeader) -> Result<Option<AniPosIndex>> {
    if header.pos_index_len == 0 {
        return Ok(None);
    }
    let start = header.off_pos_index as usize;
    let end = start + header.pos_index_len as usize;
    if end > mmap.len() {
        return Err(anyhow!("ANI pos index out of range"));
    }
    let bytes = &mmap[start..end];
    if bytes.len() < mem::size_of::<AniPosIndexHeader>() {
        return Err(anyhow!("ANI pos index too small"));
    }
    let h = unsafe { *(bytes.as_ptr() as *const AniPosIndexHeader) };
    let contig_start = h.off_contigs as usize;
    let block_start = h.off_blocks as usize;
    let pos_offsets_start = h.off_pos_offsets as usize;
    let pos_counts_start = h.off_pos_counts as usize;
    let entry_indices_start = h.off_entry_indices as usize;

    let contig_len = h.contig_count as usize * mem::size_of::<AniPosContig>();
    let block_len = h.block_count as usize * mem::size_of::<AniPosBlock>();
    let pos_offsets_len = h.pos_count as usize * mem::size_of::<u32>();
    let pos_counts_len = h.pos_count as usize * mem::size_of::<u16>();
    let entry_indices_len = h.entry_index_count as usize * mem::size_of::<u32>();

    let contig_end = contig_start + contig_len;
    let block_end = block_start + block_len;
    let pos_offsets_end = pos_offsets_start + pos_offsets_len;
    let pos_counts_end = pos_counts_start + pos_counts_len;
    let entry_indices_end = entry_indices_start + entry_indices_len;

    if contig_end > bytes.len()
        || block_end > bytes.len()
        || pos_offsets_end > bytes.len()
        || pos_counts_end > bytes.len()
        || entry_indices_end > bytes.len()
    {
        return Err(anyhow!("ANI pos index offsets invalid"));
    }

    let contig_bytes = &bytes[contig_start..contig_end];
    let block_bytes = &bytes[block_start..block_end];
    let pos_offsets_bytes = &bytes[pos_offsets_start..pos_offsets_end];
    let pos_counts_bytes = &bytes[pos_counts_start..pos_counts_end];
    let entry_indices_bytes = &bytes[entry_indices_start..entry_indices_end];

    let contigs: Vec<AniPosContig> = read_pod_slice(contig_bytes)?;
    let blocks: Vec<AniPosBlock> = read_pod_slice(block_bytes)?;
    let pos_offsets: Vec<u32> = read_pod_slice(pos_offsets_bytes)?;
    let pos_counts: Vec<u16> = read_pod_slice(pos_counts_bytes)?;
    let entry_indices: Vec<u32> = read_pod_slice(entry_indices_bytes)?;

    if contigs.len() != h.contig_count as usize
        || blocks.len() != h.block_count as usize
        || pos_offsets.len() != h.pos_count as usize
        || pos_counts.len() != h.pos_count as usize
        || entry_indices.len() != h.entry_index_count as usize
    {
        return Err(anyhow!("ANI pos index count mismatch"));
    }

    Ok(Some(AniPosIndex {
        contigs,
        blocks,
        pos_offsets,
        pos_counts,
        entry_indices,
    }))
}

fn load_info_blob(mmap: &Mmap, header: &AniHeader) -> Result<Option<AniInfoBlob>> {
    if header.blob_len == 0 {
        return Ok(None);
    }
    let start = header.off_blob as usize;
    let end = start + header.blob_len as usize;
    if end > mmap.len() {
        return Err(anyhow!("ANI blob out of range"));
    }
    let bytes = &mmap[start..end];
    if bytes.len() < mem::size_of::<AniInfoBlobHeader>() {
        return Err(anyhow!("ANI blob too small"));
    }
    let h = unsafe { *(bytes.as_ptr() as *const AniInfoBlobHeader) };
    let dict_offsets_start = h.off_dict_offsets as usize;
    let dict_data_start = h.off_dict_data as usize;
    let entry_offsets_start = h.off_entry_offsets as usize;
    let entry_counts_start = h.off_entry_counts as usize;
    let pairs_start = h.off_pairs as usize;
    let values_start = h.off_values as usize;
    if dict_offsets_start > bytes.len()
        || dict_data_start > bytes.len()
        || entry_offsets_start > bytes.len()
        || entry_counts_start > bytes.len()
        || pairs_start > bytes.len()
        || values_start > bytes.len()
    {
        return Err(anyhow!("ANI blob offsets invalid"));
    }
    let dict_offsets_len = h.dict_count as usize * mem::size_of::<u32>();
    let entry_offsets_len = h.n_entries as usize * mem::size_of::<u32>();
    let entry_counts_len = h.n_entries as usize * mem::size_of::<u16>();
    let pairs_len = h.pair_count as usize * mem::size_of::<AniInfoPair>();

    let dict_offsets_end = dict_offsets_start + dict_offsets_len;
    let entry_offsets_end = entry_offsets_start + entry_offsets_len;
    let entry_counts_end = entry_counts_start + entry_counts_len;
    let pairs_end = pairs_start + pairs_len;

    if dict_offsets_end > bytes.len()
        || entry_offsets_end > bytes.len()
        || entry_counts_end > bytes.len()
        || pairs_end > bytes.len()
    {
        return Err(anyhow!("ANI blob offsets invalid"));
    }

    let dict_offsets_bytes = &bytes[dict_offsets_start..dict_offsets_end];
    let dict_data = bytes[dict_data_start..entry_offsets_start].to_vec();
    let entry_offsets_bytes = &bytes[entry_offsets_start..entry_offsets_end];
    let entry_counts_bytes = &bytes[entry_counts_start..entry_counts_end];
    let pairs_bytes = &bytes[pairs_start..pairs_end];
    let values = bytes[values_start..].to_vec();

    let dict_offsets: Vec<u32> = read_pod_slice(dict_offsets_bytes)?;
    let entry_offsets: Vec<u32> = read_pod_slice(entry_offsets_bytes)?;
    let entry_counts: Vec<u16> = read_pod_slice(entry_counts_bytes)?;
    let pairs: Vec<AniInfoPair> = read_pod_slice(pairs_bytes)?;

    if dict_offsets.len() != h.dict_count as usize
        || entry_offsets.len() != h.n_entries as usize
        || entry_counts.len() != h.n_entries as usize
        || pairs.len() != h.pair_count as usize
    {
        return Err(anyhow!("ANI blob entry count mismatch"));
    }

    let mut dict_strings = Vec::with_capacity(dict_offsets.len());
    for &ofs in &dict_offsets {
        let start = ofs as usize;
        if start >= dict_data.len() {
            dict_strings.push(String::new());
            continue;
        }
        let end = match memchr(0, &dict_data[start..]) {
            Some(pos) => start + pos,
            None => dict_data.len(),
        };
        let s = std::str::from_utf8(&dict_data[start..end])
            .unwrap_or("")
            .to_string();
        dict_strings.push(s);
    }

    Ok(Some(AniInfoBlob {
        dict_offsets,
        dict_data,
        dict_strings,
        entry_offsets,
        entry_counts,
        pairs,
        values,
    }))
}

fn find_cstr_end(bytes: &[u8], start: usize) -> usize {
    if start >= bytes.len() {
        return bytes.len();
    }
    match memchr(0, &bytes[start..]) {
        Some(pos) => start + pos,
        None => bytes.len(),
    }
}

fn find_block_index(blocks: &[AniBlockEntry], offset: u64) -> Option<usize> {
    if blocks.is_empty() {
        return None;
    }

    let idx = match blocks.binary_search_by(|b| {
        if offset < b.raw_start {
            Ordering::Greater
        } else if offset >= b.raw_start + b.raw_len as u64 {
            Ordering::Less
        } else {
            Ordering::Equal
        }
    }) {
        Ok(i) => i,
        Err(i) => i.saturating_sub(1),
    };

    let b = &blocks[idx];
    if offset >= b.raw_start && offset < b.raw_start + b.raw_len as u64 {
        Some(idx)
    } else {
        None
    }
}
