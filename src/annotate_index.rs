use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::mem;
use std::path::Path;
use std::slice;

use anyhow::Result;
use kira_kv_engine::{BuildConfig, Builder, Mphf};
use memmap2::{Mmap, MmapOptions};

use crate::chr_name_to_id;

const ANI_MAGIC: &[u8; 8] = b"ANI00001";
const VERSION: u32 = 1;

#[repr(C)]
#[derive(Clone, Copy)]
pub struct AniHeader {
    magic: [u8; 8],
    version: u32,
    n_entries: u64,
    mph_m: u32,
    mph_salt: u32,
    off_mph_g: u64,
    off_entries: u64,
    off_strings: u64,
}

impl AniHeader {
    fn new(n: usize, mph: &Mphf, g_size: usize, entries_size: usize) -> Self {
        let header_size = mem::size_of::<Self>() as u64;

        Self {
            magic: *ANI_MAGIC,
            version: VERSION,
            n_entries: n as u64,
            mph_m: mph.m,
            mph_salt: mph.salt as u32,
            off_mph_g: header_size,
            off_entries: header_size + g_size as u64,
            off_strings: header_size + g_size as u64 + entries_size as u64,
        }
    }

    pub fn magic(&self) -> &[u8; 8] {
        &self.magic
    }
    pub fn version(&self) -> u32 {
        self.version
    }
    pub fn n_entries(&self) -> u64 {
        self.n_entries
    }
}

#[repr(C)]
#[derive(Clone, Copy)]
pub struct AniEntry {
    pub chr_id: u8,
    pub pos: u32,
    pub ref_ofs: u32,
    pub alt_ofs: u32,
    pub info_ofs: u32,
}

pub fn build_ani_index(input_vcf: &Path, output_ani: &Path) -> Result<()> {
    let f = File::open(input_vcf)?;
    let rdr = BufReader::new(f);

    let mut entries: Vec<(u64, AniEntry)> = Vec::new();
    let mut string_pool = String::new();

    eprintln!("[ani] reading db VCF...");

    for line in rdr.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }

        let mut cols = line.split('\t');

        let chrom = cols.next().unwrap();
        let pos = cols.next().unwrap().parse::<u32>().unwrap();
        let _id = cols.next().unwrap();
        let ref_ = cols.next().unwrap();
        let alt = cols.next().unwrap();
        let _qual = cols.next().unwrap();
        let _filter = cols.next().unwrap();
        let info = cols.next().unwrap_or("");

        let chr_id = match chr_name_to_id(chrom) {
            Some(v) => v,
            None => continue,
        };

        let ref_ofs = string_pool.len() as u32;
        string_pool.push_str(ref_);
        string_pool.push('\0');

        let alt_ofs = string_pool.len() as u32;
        string_pool.push_str(alt);
        string_pool.push('\0');

        let info_ofs = string_pool.len() as u32;
        string_pool.push_str(info);
        string_pool.push('\0');

        let mut h = fxhash::hash64(&[chr_id]);
        h ^= fxhash::hash64(pos.to_le_bytes().as_ref());
        h ^= fxhash::hash64(ref_.as_bytes());
        h ^= fxhash::hash64(alt.as_bytes());

        entries.push((
            h,
            AniEntry {
                chr_id,
                pos,
                ref_ofs,
                alt_ofs,
                info_ofs,
            },
        ));
    }

    eprintln!("[ani] building MPH...");

    let keys_bytes: Vec<[u8; 8]> = entries.iter().map(|(k, _)| k.to_le_bytes()).collect();

    let mph = Builder::new()
        .with_config(BuildConfig {
            gamma: 1.27,
            rehash_limit: 16,
            salt: 0xabcddccd00112233,
        })
        .build(keys_bytes.iter().map(|b| b.as_slice()))?;

    let n = entries.len();
    let mut arr_entries: Vec<AniEntry> = vec![
        AniEntry {
            chr_id: 0,
            pos: 0,
            ref_ofs: 0,
            alt_ofs: 0,
            info_ofs: 0,
        };
        n
    ];

    for (k, e) in &entries {
        let idx = mph.index(&k.to_le_bytes()) as usize;
        arr_entries[idx] = *e;
    }

    let g_size = mph.g.len() * mem::size_of::<u32>();
    let entries_size = arr_entries.len() * mem::size_of::<AniEntry>();

    let header = AniHeader::new(n, &mph, g_size, entries_size);

    let out = File::create(output_ani)?;
    let mut bw = BufWriter::new(out);

    unsafe {
        bw.write_all(slice::from_raw_parts(
            (&header as *const AniHeader) as *const u8,
            mem::size_of::<AniHeader>(),
        ))?;
    }

    let g_bytes = unsafe {
        slice::from_raw_parts(
            mph.g.as_ptr() as *const u8,
            mph.g.len() * mem::size_of::<u32>(),
        )
    };
    bw.write_all(g_bytes)?;

    let entries_bytes = unsafe {
        slice::from_raw_parts(
            arr_entries.as_ptr() as *const u8,
            arr_entries.len() * mem::size_of::<AniEntry>(),
        )
    };
    bw.write_all(entries_bytes)?;

    bw.write_all(string_pool.as_bytes())?;
    bw.flush()?;

    eprintln!("[ani] DONE: {} variants", n);

    Ok(())
}

pub struct AniIndex {
    pub header: AniHeader,
    pub mph: Mphf,
    pub entries: Vec<AniEntry>,
    pub string_block: Vec<u8>,
    pub mmap: Mmap,
}

impl AniIndex {
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };

        let header: AniHeader = unsafe { *(mmap.as_ptr() as *const AniHeader) };

        if &header.magic != ANI_MAGIC {
            anyhow::bail!("ANI: bad magic");
        }
        if header.version != VERSION {
            anyhow::bail!("ANI: bad version {}", header.version);
        }

        let g_start = header.off_mph_g as usize;
        let g_end = g_start + header.mph_m as usize * 4;

        let g_slice = &mmap[g_start..g_end];
        let mut g = Vec::with_capacity(header.mph_m as usize);
        for chunk in g_slice.chunks_exact(4) {
            g.push(u32::from_le_bytes(chunk.try_into().unwrap()));
        }

        let mph = Mphf {
            n: header.n_entries,
            m: header.mph_m,
            salt: header.mph_salt as u64,
            g,
        };

        let ent_start = header.off_entries as usize;
        let ent_end = ent_start + header.n_entries as usize * mem::size_of::<AniEntry>();
        let ent_bytes = &mmap[ent_start..ent_end];

        let mut entries = Vec::with_capacity(header.n_entries as usize);
        for chunk in ent_bytes.chunks_exact(mem::size_of::<AniEntry>()) {
            let e = unsafe { *(chunk.as_ptr() as *const AniEntry) };
            entries.push(e);
        }

        let str_start = header.off_strings as usize;
        let string_block = mmap[str_start..].to_vec();

        Ok(Self {
            header,
            mph,
            entries,
            string_block,
            mmap,
        })
    }

    pub fn lookup(&self, chr: &str, pos: u32, ref_: &str, alt: &str) -> Option<&str> {
        let chr_id = chr_name_to_id(chr)?;

        let mut h = fxhash::hash64(&[chr_id]);
        h ^= fxhash::hash64(pos.to_le_bytes().as_ref());
        h ^= fxhash::hash64(ref_.as_bytes());
        h ^= fxhash::hash64(alt.as_bytes());

        let idx = self.mph.index(&h.to_le_bytes()) as usize;
        if idx >= self.entries.len() {
            return None;
        }

        let e = &self.entries[idx];

        if e.chr_id != chr_id || e.pos != pos {
            return None;
        }

        let ref_s = read_cstring(&self.string_block, e.ref_ofs as usize);
        let alt_s = read_cstring(&self.string_block, e.alt_ofs as usize);

        if ref_s != ref_ || alt_s != alt {
            return None;
        }

        let info = read_cstring(&self.string_block, e.info_ofs as usize);
        Some(info)
    }
}

fn read_cstring<'a>(data: &'a [u8], mut pos: usize) -> &'a str {
    let start = pos;
    while pos < data.len() && data[pos] != 0 {
        pos += 1;
    }
    unsafe { std::str::from_utf8_unchecked(&data[start..pos]) }
}
