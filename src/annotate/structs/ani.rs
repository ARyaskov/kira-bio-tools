use anyhow::{anyhow, Result};
use kira_kv_engine::Mphf;
use memmap2::Mmap;
use memmap2::MmapOptions;
use std::fs::File;
use std::mem;
use std::path::Path;

use super::bundle::{parse_info_field, AnnotationBundle};

pub const ANI_MAGIC: &[u8; 8] = b"ANI00002";
pub const ANI_VERSION: u32 = 2;

#[repr(C)]
#[derive(Clone, Copy)]
pub struct AniHeader {
    pub magic: [u8; 8],
    pub version: u32,
    pub n_entries: u64,
    pub mph_m: u32,
    pub mph_salt: u64,
    pub off_mph_g: u64,
    pub off_entries: u64,
    pub off_strings: u64,
}

impl AniHeader {
    pub fn new(n: usize, mph: &Mphf, g_size: usize, entries_size: usize) -> Self {
        let head = mem::size_of::<Self>() as u64;
        Self {
            magic: *ANI_MAGIC,
            version: ANI_VERSION,
            n_entries: n as u64,
            mph_m: mph.m,
            mph_salt: mph.salt,
            off_mph_g: head,
            off_entries: head + g_size as u64,
            off_strings: head + g_size as u64 + entries_size as u64,
        }
    }

    pub fn validate(&self) -> Result<()> {
        if self.magic != *ANI_MAGIC {
            return Err(anyhow!("Bad ANI magic"));
        }
        if self.version != ANI_VERSION {
            return Err(anyhow!("ANI version mismatch"));
        }
        Ok(())
    }
}

#[repr(C)]
#[derive(Clone, Copy, Default)]
pub struct AniEntry {
    pub chr_id: u8,
    pub pos: u32,
    pub ref_ofs: u32,
    pub alt_ofs: u32,
    pub id_ofs: u32,
    pub qual_ofs: u32,
    pub filter_ofs: u32,
    pub info_ofs: u32,
    pub info_len: u32,
}

pub struct AniIndex {
    pub header: AniHeader,
    pub mph: Mphf,
    pub entries: Vec<AniEntry>,
    pub strings: Vec<u8>,
    pub mmap: Mmap,
}

impl AniIndex {
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };

        if mmap.len() < mem::size_of::<AniHeader>() {
            return Err(anyhow!("ANI file too small"));
        }

        let header: AniHeader = unsafe { *(mmap.as_ptr() as *const AniHeader) };
        header.validate()?;

        let g_start = header.off_mph_g as usize;
        let g_end = g_start + header.mph_m as usize * 4;
        let mut g = Vec::with_capacity(header.mph_m as usize);
        for chunk in mmap[g_start..g_end].chunks_exact(4) {
            g.push(u32::from_le_bytes(chunk.try_into().unwrap()));
        }

        let mph = Mphf {
            n: header.n_entries,
            m: header.mph_m,
            salt: header.mph_salt,
            g,
        };

        let ent_start = header.off_entries as usize;
        let ent_end = ent_start + header.n_entries as usize * mem::size_of::<AniEntry>();

        let mut entries = Vec::new();
        for chunk in mmap[ent_start..ent_end].chunks_exact(mem::size_of::<AniEntry>()) {
            let e: AniEntry = unsafe { *(chunk.as_ptr() as *const AniEntry) };
            entries.push(e);
        }

        let str_start = header.off_strings as usize;
        let strings = mmap[str_start..].to_vec();

        Ok(Self {
            header,
            mph,
            entries,
            strings,
            mmap,
        })
    }

    /// Lookup annotation for specific chr:pos:ref:alt
    /// Returns annotation bundle for EACH matching alt allele
    ///
    /// Note: alt_list is derived from input alt_raw, so it has input lifetime
    pub fn lookup_full<'a>(
        &'a self,
        chr: &str,
        pos: u32,
        rf: &str,
        alt_raw: &'a str,
    ) -> Option<(AnnotationBundle<'a>, Vec<&'a str>)> {
        use crate::chr_name_to_id;
        use fxhash::hash64;

        let chr_id = chr_name_to_id(chr)? as u8;

        // alt_list lifetime is tied to alt_raw input
        let alt_list: Vec<&'a str> = alt_raw.split(',').collect();

        // Try to find annotation for EACH alt allele
        for alt in &alt_list {
            let mut h = hash64(&[chr_id]);
            h ^= hash64(pos.to_le_bytes().as_ref());
            h ^= hash64(rf.as_bytes());
            h ^= hash64(alt.as_bytes());

            let idx = self.mph.index(&h.to_le_bytes()) as usize;
            if idx >= self.entries.len() {
                continue;
            }

            let e = &self.entries[idx];

            let rf_str = read_cstring(&self.strings, e.ref_ofs as usize);
            let alt_str = read_cstring(&self.strings, e.alt_ofs as usize);

            // Verify exact match
            if rf_str != rf || alt_str != *alt {
                continue;
            }

            // Found match - extract all fields (all have lifetime 'a from &'a self)
            let id_str = read_cstring(&self.strings, e.id_ofs as usize);
            let qual_str = read_cstring(&self.strings, e.qual_ofs as usize);
            let filter_str = read_cstring(&self.strings, e.filter_ofs as usize);
            let info_str = read_cstring(&self.strings, e.info_ofs as usize);

            let info_fields = parse_info_field(info_str);

            let bundle = AnnotationBundle {
                id: if id_str == "." || id_str.is_empty() {
                    None
                } else {
                    Some(id_str)
                },
                qual: if qual_str == "." || qual_str.is_empty() {
                    None
                } else {
                    Some(qual_str)
                },
                filter: if filter_str == "." || filter_str.is_empty() {
                    None
                } else {
                    Some(filter_str)
                },
                info: info_fields,
            };

            // bundle has lifetime 'a (from self.strings)
            // alt_list has lifetime 'a (from input alt_raw)
            return Some((bundle, alt_list));
        }

        None
    }
}

pub fn read_cstring<'a>(data: &'a [u8], mut pos: usize) -> &'a str {
    let start = pos;
    while pos < data.len() && data[pos] != 0 {
        pos += 1;
    }
    unsafe { std::str::from_utf8_unchecked(&data[start..pos]) }
}
