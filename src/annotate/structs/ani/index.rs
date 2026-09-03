use anyhow::{Result, anyhow};
use bytemuck::{self, Pod, Zeroable};
#[cfg(feature = "gpu")]
use cust::memory::DeviceCopy;
use kira_kv_engine::Index;
use flate2::{Decompress, FlushDecompress, Status};
use memchr::memchr;
use memmap2::{Mmap, MmapOptions};
use std::cell::RefCell;
use std::cmp::Ordering;
use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::mem;
use std::ops::Deref;
use std::path::Path;
use std::sync::{Arc, Mutex, OnceLock};

use super::contig_dict::ContigDict;
use super::header::{
    ANI_HEADER_END, ANI_MAGIC, ANI_SENTINEL_CHR_ID, ANI_STR_NONE, ANI_VERSION, AniBlockEntry,
    AniEntry, AniEntryV2, AniHeader, AniHeaderV3, AniHeaderV6,
};
use crate::annotate::structs::bundle::{
    AnnotationBundle, FieldNumber, FieldType, StructuredInfoField,
};

const BLOCK_CACHE_CAP: usize = 16;
const CSTR_CACHE_CAP: usize = 256;
/// Max number of dense `Vec`-indexed contigs in the pos-index lookup table.
/// Above this, falls back to `HashMap<u32, usize>`.
const POS_INDEX_DENSE_CONTIG_CAP: usize = 1024;

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

/// LRU cache of decompressed string-pool blocks.
struct BlockCacheEntry {
    idx: usize,
    data: Arc<Vec<u8>>,
}

struct BlockCache {
    entries: VecDeque<BlockCacheEntry>,
}

impl BlockCache {
    fn new() -> Self {
        Self {
            entries: VecDeque::with_capacity(BLOCK_CACHE_CAP),
        }
    }

    fn get(&mut self, idx: usize) -> Option<Arc<Vec<u8>>> {
        let pos = self.entries.iter().position(|e| e.idx == idx)?;
        let entry = self.entries.remove(pos)?;
        let data = entry.data.clone();
        self.entries.push_back(entry);
        Some(data)
    }

    fn insert(&mut self, idx: usize, data: Arc<Vec<u8>>) {
        if self.entries.len() == BLOCK_CACHE_CAP {
            self.entries.pop_front();
        }
        self.entries.push_back(BlockCacheEntry { idx, data });
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
        // SAFETY: both buffers hold ANI string sections, which the loader
        // validates as UTF-8 once when it maps or decompresses them, and
        // `start..end` are the bounds of one string within that section.
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
    pub chr_id: u32,
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

// SAFETY: plain `#[repr(C)]` integer structs (bytemuck `Pod`), copied to the
// device byte for byte.
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
    /// O(1) chr_id → index-into-`pos_index.contigs`. Dense `Vec` lookup table
    /// when there are ≤1024 contigs; `HashMap` otherwise.
    pos_contig_lut: PosContigLut,
    info_blob: Option<AniInfoBlob>,
    /// Lazily computed on first access.
    info_cache: OnceLock<AniInfoCache>,
    /// Header-derived `name → id` map for variant lookups.
    contigs: ContigDict,
}

enum PosContigLut {
    Empty,
    Dense(Vec<Option<u32>>),
    Sparse(HashMap<u32, u32>),
}

impl PosContigLut {
    fn build(contigs: &[AniPosContig]) -> Self {
        if contigs.is_empty() {
            return Self::Empty;
        }
        let max_id = contigs.iter().map(|c| c.chr_id).max().unwrap_or(0);
        if (max_id as usize) < POS_INDEX_DENSE_CONTIG_CAP {
            let mut lut = vec![None; max_id as usize + 1];
            for (idx, c) in contigs.iter().enumerate() {
                lut[c.chr_id as usize] = Some(idx as u32);
            }
            Self::Dense(lut)
        } else {
            let mut lut = HashMap::with_capacity(contigs.len());
            for (idx, c) in contigs.iter().enumerate() {
                lut.insert(c.chr_id, idx as u32);
            }
            Self::Sparse(lut)
        }
    }

    #[inline]
    fn get(&self, chr_id: u32) -> Option<usize> {
        match self {
            Self::Empty => None,
            Self::Dense(v) => v.get(chr_id as usize).and_then(|o| o.map(|x| x as usize)),
            Self::Sparse(m) => m.get(&chr_id).map(|&x| x as usize),
        }
    }
}

impl AniIndex {
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        // SAFETY: read-only mapping of an index that is replaced atomically on
        // rebuild, never modified in place while mapped.
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

        let contigs = load_contig_dict(&mmap, &header)?.unwrap_or_default();

        let pos_contig_lut = match pos_index.as_ref() {
            Some(pi) => PosContigLut::build(&pi.contigs),
            None => PosContigLut::Empty,
        };

        let mut ani = Self {
            header,
            index,
            entries,
            intervals: Vec::new(),
            strings,
            mmap,
            pos_index,
            pos_contig_lut,
            info_blob,
            info_cache: OnceLock::new(),
            contigs,
        };
        if ani.header.has_intervals() {
            ani.build_interval_index();
        }

        Ok(ani)
    }

    /// O(1) header-derived contig name → id. Returns `None` if the dict
    /// wasn't stored or the name isn't in the dict.
    #[inline]
    pub fn contig_id(&self, name: &str) -> Option<u32> {
        self.contigs.id(name)
    }

    /// Number of contigs in the dict. `0` means legacy `.ani` without a dict.
    #[inline]
    pub fn contig_count(&self) -> usize {
        self.contigs.len()
    }

    /// O(1): does the index hold any entry for the given chromosome?
    #[inline]
    pub fn chrom_has_entries(&self, chr_id: u32) -> bool {
        if self.pos_index.is_none() {
            return true;
        }
        self.pos_contig_lut.get(chr_id).is_some()
    }

    /// Pre-computed per-entry MPH verification keys (`u64[n_entries]`).
    /// Returns `None` for legacy `.ani` files without the cached section.
    pub fn cached_entry_keys(&self) -> Option<&[u64]> {
        if !self.header.has_entry_keys() {
            return None;
        }
        let start = self.header.off_entry_keys as usize;
        let len_bytes = self.header.entry_keys_len as usize;
        let end = start.checked_add(len_bytes)?;
        if end > self.mmap.len() {
            return None;
        }
        if (start % std::mem::align_of::<u64>()) != 0 {
            return None;
        }
        let n_expected = self.entries.len();
        if len_bytes != n_expected * std::mem::size_of::<u64>() {
            return None;
        }
        Some(bytemuck::cast_slice(&self.mmap[start..end]))
    }

    /// Reference to the embedded contig dict.
    #[inline]
    pub fn contig_dict(&self) -> &ContigDict {
        &self.contigs
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

    pub fn lookup_pos_index(&self, chr_id: u32, pos: u32) -> Option<&[u32]> {
        let pos_index = self.pos_index.as_ref()?;
        let contig_idx = self.pos_contig_lut.get(chr_id)?;
        let contig = &pos_index.contigs[contig_idx];
        if pos < contig.min_pos || pos > contig.max_pos {
            return None;
        }
        let start = contig.block_start as usize;
        let end = start + contig.block_count as usize;
        let blocks = &pos_index.blocks[start..end];
        let base = (pos / 512) * 512;
        let Ok(block_idx) = blocks.binary_search_by_key(&base, |b| b.base_pos) else {
            return None;
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

    /// Lazy info cache. First access pays the build cost; subsequent threads
    /// see the cached `&AniInfoCache` via `OnceLock`.
    pub fn info_cache(&self) -> Option<&AniInfoCache> {
        let blob = self.info_blob.as_ref()?;
        Some(self.info_cache.get_or_init(|| build_info_cache(self, blob)))
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
                            // Validated once per block so `CStrRef::as_str` can skip the check.
                            Ok(v) if std::str::from_utf8(&v).is_ok() => Arc::new(v),
                            _ => {
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

    pub fn find_interval_entry(&self, chr_id: u32, pos: u32) -> Option<usize> {
        if self.intervals.is_empty() {
            return None;
        }
        let list = self.intervals.get(chr_id as usize)?;
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
        let mut decompressor = Decompress::new(false);
        let status = decompressor
            .decompress(compressed, &mut out, FlushDecompress::Finish)
            .map_err(|_| anyhow!("ANI block {} decompression failed", idx))?;
        if status != Status::StreamEnd
            || decompressor.total_in() as usize != compressed.len()
            || decompressor.total_out() as usize != block.raw_len as usize
        {
            return Err(anyhow!(
                "ANI block {} decompression incomplete: status={:?} in={}/{} out={}/{}",
                idx,
                status,
                decompressor.total_in(),
                compressed.len(),
                decompressor.total_out(),
                block.raw_len
            ));
        }
        Ok(out)
    }

    fn build_interval_index(&mut self) {
        let n_chr = self.contigs.len().max(256);
        let mut per_chr: Vec<Vec<IntervalEntry>> = vec![Vec::new(); n_chr];
        for (idx, entry) in self.entries.iter().enumerate() {
            if entry.chr_id == ANI_SENTINEL_CHR_ID {
                continue;
            }
            let rf = self.read_cstring(entry.ref_ofs as usize);
            if rf.as_ref() != "." {
                continue;
            }
            let alt = self.read_cstring(entry.alt_ofs as usize);
            let Ok(end) = alt.as_ref().parse::<u32>() else {
                continue;
            };
            let bucket = entry.chr_id as usize;
            if bucket >= per_chr.len() {
                per_chr.resize_with(bucket + 1, Vec::new);
            }
            per_chr[bucket].push(IntervalEntry {
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

    pub fn info_fields_from_blob(
        &self,
        entry_idx: usize,
        field_meta: &HashMap<String, FieldNumber>,
    ) -> Vec<StructuredInfoField> {
        self.info_fields_from_blob_filtered(entry_idx, field_meta, None)
    }

    /// Selective variant of [`Self::info_fields_from_blob`]. When `filter` is
    /// `Some(keys)`, only fields whose tag string is in `keys` get parsed.
    pub fn info_fields_from_blob_filtered(
        &self,
        entry_idx: usize,
        field_meta: &HashMap<String, FieldNumber>,
        filter: Option<&std::collections::HashSet<&str>>,
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
        let cap = filter.map_or(count, |s| s.len().min(count));
        let mut fields = Vec::with_capacity(cap);
        for pair in &blob.pairs[offset..offset + count] {
            let Some(key_ref) = blob.dict_strings.get(pair.tag_id as usize) else {
                continue;
            };
            if key_ref.is_empty() {
                continue;
            }
            if let Some(set) = filter
                && !set.contains(key_ref.as_str())
            {
                continue;
            }
            let key = key_ref.clone();
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

    pub fn entry_info_string(&self, entry_idx: usize) -> Option<String> {
        if let Some(blob) = &self.info_blob {
            if entry_idx >= blob.entry_offsets.len() {
                return None;
            }
            let offset = blob.entry_offsets[entry_idx] as usize;
            let count = blob.entry_counts[entry_idx] as usize;
            if count == 0 || offset + count > blob.pairs.len() {
                return None;
            }
            let pairs = &blob.pairs[offset..offset + count];
            let mut len = pairs.len().saturating_sub(1);
            for pair in pairs {
                if let Some(key) = blob.dict_strings.get(pair.tag_id as usize) {
                    len += key.len();
                }
                if pair.value_len > 0 {
                    len += 1 + pair.value_len as usize;
                }
            }
            let mut out = String::with_capacity(len);
            for pair in pairs {
                let Some(key) = blob.dict_strings.get(pair.tag_id as usize) else {
                    continue;
                };
                if key.is_empty() {
                    continue;
                }
                if !out.is_empty() {
                    out.push(';');
                }
                out.push_str(key);
                if pair.value_len > 0 {
                    let start = pair.value_off as usize;
                    let end = start + pair.value_len as usize;
                    if end > blob.values.len() {
                        continue;
                    }
                    let value = std::str::from_utf8(&blob.values[start..end]).unwrap_or("");
                    out.push('=');
                    out.push_str(value);
                }
            }
            if out.is_empty() {
                return None;
            }
            return Some(out);
        }

        let e = self.entries.get(entry_idx)?;
        if e.info_ofs == ANI_STR_NONE || e.info_len == 0 {
            return None;
        }
        let info = self.read_cstring(e.info_ofs as usize);
        let info = info.as_ref();
        if info.is_empty() || info == "." {
            None
        } else {
            Some(info.to_string())
        }
    }

    pub fn build_bundle_from_entry_idx_opts_with_meta(
        &self,
        entry_idx: usize,
        field_meta: &HashMap<String, FieldNumber>,
        need_info: bool,
        need_format: bool,
    ) -> AnnotationBundle {
        self.build_bundle_from_entry_idx_opts_with_meta_filtered(
            entry_idx,
            field_meta,
            need_info,
            need_format,
            None,
        )
    }

    /// Selective bundle builder. When `info_filter` is `Some(keys)`, only
    /// matching INFO tags are populated.
    pub fn build_bundle_from_entry_idx_opts_with_meta_filtered(
        &self,
        entry_idx: usize,
        field_meta: &HashMap<String, FieldNumber>,
        need_info: bool,
        need_format: bool,
        info_filter: Option<&std::collections::HashSet<&str>>,
    ) -> AnnotationBundle {
        let e = &self.entries[entry_idx];
        let ref_str = self.read_cstring(e.ref_ofs as usize);
        let alt_str = self.read_cstring(e.alt_ofs as usize);
        let id_str = self.read_cstring(e.id_ofs as usize);
        let qual_str = self.read_cstring(e.qual_ofs as usize);
        let filter_str = self.read_cstring(e.filter_ofs as usize);
        let info = if need_info {
            if self.info_blob.is_some() {
                self.info_fields_from_blob_filtered(entry_idx, field_meta, info_filter)
            } else {
                let info_str = self.read_cstring(e.info_ofs as usize);
                let all = crate::annotate::structs::bundle::infer_structured_info_fields(
                    info_str.as_ref(),
                    field_meta,
                );
                if let Some(set) = info_filter {
                    all.into_iter().filter(|f| set.contains(f.key.as_str())).collect()
                } else {
                    all
                }
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
            db_ref: ref_str.to_string(),
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

/// Build the typed-value cache from a raw `AniInfoBlob`. Chunked-parallel
/// via rayon: per-chunk caches built independently, merged with offset rebasing.
fn build_info_cache(ani: &AniIndex, blob: &AniInfoBlob) -> AniInfoCache {
    use rayon::prelude::*;
    let tag_types = build_tag_types(ani, blob);
    let entry_offsets = blob.entry_offsets.clone();
    let entry_counts = blob.entry_counts.clone();
    let pair_tag_ids: Vec<u32> = blob.pairs.iter().map(|p| p.tag_id).collect();

    const CHUNK: usize = 32_768;
    let pairs = blob.pairs.as_slice();

    let chunks: Vec<ChunkCache> = pairs
        .par_chunks(CHUNK)
        .map(|chunk| build_chunk(chunk, blob, &tag_types))
        .collect();

    let mut int_base: u32 = 0;
    let mut float_base: u32 = 0;
    let mut str_base: u32 = 0;
    let mut str_data_base: u32 = 0;
    let mut total_pairs = 0;
    let mut total_ints = 0;
    let mut total_floats = 0;
    let mut total_strs = 0;
    let mut total_str_data = 0;
    for c in &chunks {
        total_pairs += c.pair_offsets.len();
        total_ints += c.int_values.len();
        total_floats += c.float_values.len();
        total_strs += c.str_offsets.len();
        total_str_data += c.str_data.len();
    }

    let mut pair_offsets = Vec::with_capacity(total_pairs);
    let mut pair_counts = Vec::with_capacity(total_pairs);
    let mut int_values = Vec::with_capacity(total_ints);
    let mut float_values = Vec::with_capacity(total_floats);
    let mut str_offsets = Vec::with_capacity(total_strs);
    let mut str_lens = Vec::with_capacity(total_strs);
    let mut str_data = Vec::with_capacity(total_str_data);

    for c in chunks {
        for (off, count) in c.pair_offsets.iter().zip(c.pair_counts.iter()) {
            let off = *off;
            let count = *count;
            let rebased = if count == 0 { 0u32 } else { off };
            pair_offsets.push(rebased);
            pair_counts.push(count);
        }
        let pairs_start = pair_offsets.len() - c.pair_offsets.len();
        for (i, kind) in c.pair_kinds.iter().enumerate() {
            let p = pairs_start + i;
            match kind {
                ChunkPairKind::None => {}
                ChunkPairKind::Int => {
                    pair_offsets[p] = pair_offsets[p].wrapping_add(int_base);
                }
                ChunkPairKind::Float => {
                    pair_offsets[p] = pair_offsets[p].wrapping_add(float_base);
                }
                ChunkPairKind::Str => {
                    pair_offsets[p] = pair_offsets[p].wrapping_add(str_base);
                }
            }
        }
        int_values.extend_from_slice(&c.int_values);
        float_values.extend_from_slice(&c.float_values);
        for s_off in c.str_offsets.iter() {
            str_offsets.push(s_off.wrapping_add(str_data_base));
        }
        str_lens.extend_from_slice(&c.str_lens);
        str_data.extend_from_slice(&c.str_data);

        int_base = int_base.wrapping_add(c.int_values.len() as u32);
        float_base = float_base.wrapping_add(c.float_values.len() as u32);
        str_base = str_base.wrapping_add(c.str_offsets.len() as u32);
        str_data_base = str_data_base.wrapping_add(c.str_data.len() as u32);
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

/// Per-pair type marker used when rebasing `pair_offsets[i]` during merge.
#[derive(Clone, Copy)]
enum ChunkPairKind {
    None,
    Int,
    Float,
    Str,
}

struct ChunkCache {
    pair_offsets: Vec<u32>,
    pair_counts: Vec<u16>,
    pair_kinds: Vec<ChunkPairKind>,
    int_values: Vec<i32>,
    float_values: Vec<f32>,
    str_offsets: Vec<u32>,
    str_lens: Vec<u32>,
    str_data: Vec<u8>,
}

fn build_chunk(
    pairs: &[AniInfoPair],
    blob: &AniInfoBlob,
    tag_types: &[FieldType],
) -> ChunkCache {
    let mut pair_offsets = Vec::with_capacity(pairs.len());
    let mut pair_counts = Vec::with_capacity(pairs.len());
    let mut pair_kinds = Vec::with_capacity(pairs.len());
    let mut int_values = Vec::new();
    let mut float_values = Vec::new();
    let mut str_offsets = Vec::new();
    let mut str_lens = Vec::new();
    let mut str_data = Vec::new();

    for pair in pairs {
        let tag_type = tag_types
            .get(pair.tag_id as usize)
            .copied()
            .unwrap_or(FieldType::String);
        if pair.value_len == 0 {
            pair_offsets.push(0);
            pair_counts.push(0);
            pair_kinds.push(ChunkPairKind::None);
            continue;
        }
        let start = pair.value_off as usize;
        let end = start + pair.value_len as usize;
        if end > blob.values.len() {
            pair_offsets.push(0);
            pair_counts.push(0);
            pair_kinds.push(ChunkPairKind::None);
            continue;
        }
        let raw = std::str::from_utf8(&blob.values[start..end]).unwrap_or("");
        match tag_type {
            FieldType::Flag => {
                pair_offsets.push(0);
                pair_counts.push(0);
                pair_kinds.push(ChunkPairKind::None);
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
                pair_kinds.push(ChunkPairKind::Int);
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
                pair_kinds.push(ChunkPairKind::Float);
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
                pair_kinds.push(ChunkPairKind::Str);
            }
        }
    }

    ChunkCache {
        pair_offsets,
        pair_counts,
        pair_kinds,
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
    if mmap.len() < mem::size_of::<AniHeaderV3>() {
        return Err(anyhow!("ANI file too small for a header"));
    }
    let h3: AniHeaderV3 = bytemuck::pod_read_unaligned(&mmap[..mem::size_of::<AniHeaderV3>()]);
    if h3.magic != ANI_MAGIC {
        return Err(anyhow!("Bad ANI magic"));
    }

    match h3.version {
        ANI_VERSION => {
            if mmap.len() < mem::size_of::<AniHeaderV6>() {
                return Err(anyhow!("ANI file too small for v6 header (rebuild)"));
            }
            let h6: AniHeaderV6 = bytemuck::pod_read_unaligned(&mmap[..mem::size_of::<AniHeaderV6>()]);
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
                flags: h6.flags,
                off_pos_index: h6.off_pos_index,
                pos_index_len: h6.pos_index_len,
                off_blob: h6.off_blob,
                blob_len: h6.blob_len,
                off_contigs: h6.off_contigs,
                contigs_len: h6.contigs_len,
                off_entry_keys: h6.off_entry_keys,
                entry_keys_len: h6.entry_keys_len,
            })
        }
        2 | 3 | 4 | 5 => Err(anyhow!(
            "Unsupported ANI version {} (rebuild index)",
            h3.version
        )),
        _ => Err(anyhow!("Unsupported ANI version {}", h3.version)),
    }
}

fn load_contig_dict(mmap: &Mmap, header: &AniHeader) -> Result<Option<ContigDict>> {
    if header.contigs_len == 0 {
        return Ok(None);
    }
    let start = header.off_contigs as usize;
    let end = start
        .checked_add(header.contigs_len as usize)
        .ok_or_else(|| anyhow!("ANI contig dict offset overflow"))?;
    if end > mmap.len() {
        return Err(anyhow!("ANI contig dict out of range"));
    }
    ContigDict::parse_bytes(&mmap[start..end])
        .map(Some)
        .map_err(|e| anyhow!("ANI contig dict: {e}"))
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
        // SAFETY: chunk is exactly size_of::<AniEntryV2>() bytes of a plain-data struct.
        let e: AniEntryV2 = unsafe { std::ptr::read_unaligned(chunk.as_ptr() as *const AniEntryV2) };
        entries.push(AniEntry {
            chr_id: e.chr_id as u32,
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
        // SAFETY: chunk is exactly size_of::<AniEntry>() bytes of a plain-data struct.
        let e: AniEntry = unsafe { std::ptr::read_unaligned(chunk.as_ptr() as *const AniEntry) };
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
            // SAFETY: chunk is exactly size_of::<AniBlockEntry>() bytes of a plain-data struct.
            let e: AniBlockEntry = unsafe { std::ptr::read_unaligned(chunk.as_ptr() as *const AniBlockEntry) };
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
        // Validated once here so `CStrRef::as_str` can skip the check.
        if std::str::from_utf8(&mmap[str_start.min(mmap.len())..]).is_err() {
            return Err(anyhow!("ANI string section is not valid UTF-8 (rebuild the index)"));
        }
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
    // SAFETY: length checked above; the header is a plain-data struct.
    let h = unsafe { std::ptr::read_unaligned(bytes.as_ptr() as *const AniPosIndexHeader) };
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
    // SAFETY: length checked above; the header is a plain-data struct.
    let h = unsafe { std::ptr::read_unaligned(bytes.as_ptr() as *const AniInfoBlobHeader) };
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

#[cfg(test)]
#[path = "../../../../tests/unit/annotate_structs_ani_index.rs"]
mod tests;
