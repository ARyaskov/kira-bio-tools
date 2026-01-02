use anyhow::{anyhow, Result};
use kira_kv_engine::Mphf;
use memmap2::{Mmap, MmapOptions};
use std::fs::File;
use std::mem;
use std::path::Path;

use super::header::{AniEntry, AniEntryV2, AniHeader, ANI_STR_NONE};

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

        let mph = load_mph(&mmap, &header)?;
        let entries = load_entries(&mmap, &header)?;
        let strings = load_strings(&mmap, &header);

        Ok(Self {
            header,
            mph,
            entries,
            strings,
            mmap,
        })
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
        3 => load_entries_v3(mmap, ent_start, header.n_entries as usize),
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

fn load_strings(mmap: &Mmap, header: &AniHeader) -> Vec<u8> {
    let str_start = header.off_strings as usize;
    mmap[str_start..].to_vec()
}
