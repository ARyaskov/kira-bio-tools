use anyhow::{anyhow, Result};
use kira_kv_engine::Mphf;
use libdeflater::Decompressor;
use memchr::memchr;
use memmap2::{Mmap, MmapOptions};
use std::cmp::Ordering;
use std::fs::File;
use std::mem;
use std::ops::Deref;
use std::path::Path;
use std::sync::{Arc, Mutex};

use super::header::{
    AniBlockEntry, AniEntry, AniEntryV2, AniHeader, AniHeaderV3, AniHeaderV4, ANI_MAGIC,
    ANI_STR_NONE, ANI_VERSION,
};

const BLOCK_CACHE_CAP: usize = 16;

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

pub struct AniIndex {
    pub header: AniHeader,
    pub mph: Mphf,
    pub entries: Vec<AniEntry>,
    intervals: Vec<Vec<IntervalEntry>>,
    strings: StringSource,
    mmap: Mmap,
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

        let mph = load_mph(&mmap, &header)?;
        let entries = load_entries(&mmap, &header)?;
        let strings = load_strings(&mmap, &header)?;

        let mut ani = Self {
            header,
            mph,
            entries,
            intervals: Vec::new(),
            strings,
            mmap,
        };
        ani.build_interval_index();

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
                let end = find_cstr_end(bytes, offset);
                CStrRef {
                    data: CStrData::Borrowed(bytes),
                    start: offset,
                    end,
                }
            }
            StringSource::Compressed { blocks, cache, .. } => {
                let offset = offset as u64;
                let idx = match find_block_index(blocks, offset) {
                    Some(v) => v,
                    None => {
                        return CStrRef {
                            data: CStrData::Borrowed(&[]),
                            start: 0,
                            end: 0,
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
                                }
                            }
                        };
                        let mut guard = cache.lock().unwrap();
                        guard.insert(idx, decompressed.clone());
                        decompressed
                    }
                };

                let block = &blocks[idx];
                let start = (offset - block.raw_start) as usize;
                let end = find_cstr_end(data.as_slice(), start);

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
}

fn read_header(mmap: &Mmap) -> Result<AniHeader> {
    let h3: AniHeaderV3 = unsafe { *(mmap.as_ptr() as *const AniHeaderV3) };
    if h3.magic != ANI_MAGIC {
        return Err(anyhow!("Bad ANI magic"));
    }

    match h3.version {
        2 | 3 => Ok(AniHeader {
            magic: h3.magic,
            version: h3.version,
            n_entries: h3.n_entries,
            mph_m: h3.mph_m,
            mph_salt: h3.mph_salt,
            off_mph_g: h3.off_mph_g,
            off_entries: h3.off_entries,
            off_strings: h3.off_strings,
            off_block_index: 0,
            n_blocks: 0,
            block_size: 0,
        }),
        ANI_VERSION => {
            if mmap.len() < mem::size_of::<AniHeaderV4>() {
                return Err(anyhow!("ANI file too small for v4 header"));
            }
            let h4: AniHeaderV4 = unsafe { *(mmap.as_ptr() as *const AniHeaderV4) };
            Ok(AniHeader {
                magic: h4.magic,
                version: h4.version,
                n_entries: h4.n_entries,
                mph_m: h4.mph_m,
                mph_salt: h4.mph_salt,
                off_mph_g: h4.off_mph_g,
                off_entries: h4.off_entries,
                off_strings: h4.off_strings,
                off_block_index: h4.off_block_index,
                n_blocks: h4.n_blocks,
                block_size: h4.block_size,
            })
        }
        _ => Err(anyhow!("Unsupported ANI version {}", h3.version)),
    }
}

fn load_mph(mmap: &Mmap, header: &AniHeader) -> Result<Mphf> {
    let g_start = header.off_mph_g as usize;
    let g_end = g_start + header.mph_m as usize * 4;
    let mut g = Vec::with_capacity(header.mph_m as usize);

    for chunk in mmap[g_start..g_end].chunks_exact(4) {
        g.push(u32::from_le_bytes(chunk.try_into().unwrap()));
    }

    Ok(Mphf {
        n: header.n_entries,
        m: header.mph_m as u32,
        salt: header.mph_salt,
        g,
    })
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
