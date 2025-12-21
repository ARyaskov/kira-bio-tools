use anyhow::{anyhow, Result};
use kira_kv_engine::Mphf;
use memmap2::{Mmap, MmapOptions};
use std::fs::File;
use std::mem;
use std::path::Path;

use super::header::{AniEntry, AniHeader};

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
    let ent_end = ent_start + header.n_entries as usize * mem::size_of::<AniEntry>();

    let mut entries = Vec::new();
    for chunk in mmap[ent_start..ent_end].chunks_exact(mem::size_of::<AniEntry>()) {
        let e: AniEntry = unsafe { *(chunk.as_ptr() as *const AniEntry) };
        entries.push(e);
    }

    Ok(entries)
}

fn load_strings(mmap: &Mmap, header: &AniHeader) -> Vec<u8> {
    let str_start = header.off_strings as usize;
    mmap[str_start..].to_vec()
}
