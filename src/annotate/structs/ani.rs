use anyhow::{anyhow, Result};
use kira_kv_engine::Mphf;
use memmap2::{Mmap, MmapOptions};
use std::fs::File;
use std::mem;
use std::path::Path;

use super::bundle::{parse_info_field, AnnotationBundle};

pub const ANI_MAGIC: u64 = 0x494E4149524B4256;
pub const ANI_VERSION: u64 = 2;

#[repr(C)]
#[derive(Clone, Copy)]
pub struct AniHeader {
    pub magic: u64,
    pub version: u64,
    pub n_entries: u64,
    pub mph_m: u64,
    pub mph_salt: u64,
    pub off_mph_g: u64,
    pub off_entries: u64,
    pub off_strings: u64,
}

impl AniHeader {
    pub fn validate(&self) -> Result<()> {
        if self.magic != ANI_MAGIC {
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
            m: header.mph_m as u32,
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

    pub fn lookup(&self, chr: &str, pos: u32, rf: &str, alt: &str) -> Option<AnnotationBundle> {
        use crate::chr_name_to_id;
        use fxhash::hash64;

        let debug = std::env::var("KIRA_BT_DEBUG").is_ok();
        let chr_id = chr_name_to_id(chr)? as u8;

        let mut h = (chr_id as u64) << 32 | (pos as u64);
        h ^= hash64(rf.as_bytes());
        h ^= hash64(alt.as_bytes());

        let idx = self.mph.index(&h.to_le_bytes()) as usize;
        if idx >= self.entries.len() {
            if debug {
                eprintln!(
                    "[DEBUG-LOOKUP] MPH returned idx {} >= entries.len() {}",
                    idx,
                    self.entries.len()
                );
            }
            return None;
        }

        let e = &self.entries[idx];

        if e.chr_id != chr_id || e.pos != pos {
            return None;
        }

        let rf_str = read_cstring(&self.strings, e.ref_ofs as usize);
        if rf_str != rf {
            if debug {
                eprintln!(
                    "[DEBUG-LOOKUP] REF mismatch: expected {}, got {}",
                    rf, rf_str
                );
            }
            return None;
        }

        let alt_str = read_cstring(&self.strings, e.alt_ofs as usize);
        if alt_str != alt {
            if debug {
                eprintln!(
                    "[DEBUG-LOOKUP] ALT mismatch: expected {}, got {}",
                    alt, alt_str
                );
            }
            return None;
        }

        let id_str = read_cstring(&self.strings, e.id_ofs as usize);
        let qual_str = read_cstring(&self.strings, e.qual_ofs as usize);
        let filter_str = read_cstring(&self.strings, e.filter_ofs as usize);
        let info_str = read_cstring(&self.strings, e.info_ofs as usize);

        let info_fields = parse_info_field(info_str);

        let bundle = AnnotationBundle {
            alt: alt_str.to_string(),
            id: if id_str == "." || id_str.is_empty() {
                None
            } else {
                Some(id_str.to_string())
            },
            qual: if qual_str == "." || qual_str.is_empty() {
                None
            } else {
                Some(qual_str.to_string())
            },
            filter: if filter_str == "." || filter_str.is_empty() {
                None
            } else {
                Some(filter_str.to_string())
            },
            info: info_fields,
        };

        Some(bundle)
    }

    pub fn lookup_any_alt(&self, chr: &str, pos: u32, rf: &str) -> Option<AnnotationBundle> {
        use crate::chr_name_to_id;

        let chr_id = chr_name_to_id(chr)? as u8;

        let mut found: Option<&AniEntry> = None;

        for e in &self.entries {
            if e.chr_id != chr_id || e.pos != pos {
                continue;
            }

            let rf_str = read_cstring(&self.strings, e.ref_ofs as usize);
            if rf_str != rf {
                continue;
            }

            if found.is_some() {
                return None;
            }
            found = Some(e);
        }

        let e = found?;

        let alt_str = read_cstring(&self.strings, e.alt_ofs as usize);
        let id_str = read_cstring(&self.strings, e.id_ofs as usize);
        let qual_str = read_cstring(&self.strings, e.qual_ofs as usize);
        let filter_str = read_cstring(&self.strings, e.filter_ofs as usize);
        let info_str = read_cstring(&self.strings, e.info_ofs as usize);

        let info_fields = parse_info_field(info_str);

        Some(AnnotationBundle {
            alt: alt_str.to_string(),
            id: if id_str == "." || id_str.is_empty() {
                None
            } else {
                Some(id_str.to_string())
            },
            qual: if qual_str == "." || qual_str.is_empty() {
                None
            } else {
                Some(qual_str.to_string())
            },
            filter: if filter_str == "." || filter_str.is_empty() {
                None
            } else {
                Some(filter_str.to_string())
            },
            info: info_fields,
        })
    }
}

pub fn read_cstring<'a>(data: &'a [u8], mut pos: usize) -> &'a str {
    let start = pos;
    while pos < data.len() && data[pos] != 0 {
        pos += 1;
    }
    std::str::from_utf8(&data[start..pos]).unwrap_or("")
}
