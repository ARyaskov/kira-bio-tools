use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::mem;
use std::path::Path;
use std::slice;

use anyhow::{anyhow, Result};
use kira_kv_engine::Mphf;
use memmap2::{Mmap, MmapOptions};

use crate::chr_name_to_id;

const ANI_MAGIC: &[u8; 8] = b"ANI00001";
const VERSION: u32 = 1;

#[repr(C)]
#[derive(Clone, Copy)]
struct AniHeader {
    magic: [u8; 8],
    version: u32,
    n_entries: u64,
    mph_m: u32,
    mph_salt: u32,
    off_mph_g: u64,
    off_entries: u64,
    off_strings: u64,
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

pub struct AniIndex {
    mmap: Mmap,
    pub mph: Mphf,
    pub entries: &'static [AniEntry],
    pub string_block: &'static [u8],
}

impl AniIndex {
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };

        let data: &'static [u8] = unsafe { std::mem::transmute::<&[u8], &'static [u8]>(&mmap[..]) };

        if data.len() < mem::size_of::<AniHeader>() {
            return Err(anyhow!("ANI file too small"));
        }

        let header: &AniHeader = unsafe { &*(data.as_ptr() as *const AniHeader) };

        if &header.magic != ANI_MAGIC {
            return Err(anyhow!("Invalid ANI magic"));
        }
        if header.version != VERSION {
            return Err(anyhow!("ANI version mismatch"));
        }

        let n = header.n_entries as usize;

        let g_start = header.off_mph_g as usize;
        let g_end = g_start + header.mph_m as usize * 4;

        let g_slice = &data[g_start..g_end];
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
        let ent_end = ent_start + n * mem::size_of::<AniEntry>();

        let entries: &'static [AniEntry] = unsafe {
            slice::from_raw_parts(data[ent_start..ent_end].as_ptr() as *const AniEntry, n)
        };

        let string_block: &'static [u8] = &data[header.off_strings as usize..];

        Ok(Self {
            mmap,
            mph,
            entries,
            string_block,
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

        let r = read_cstring(self.string_block, e.ref_ofs as usize);
        let a = read_cstring(self.string_block, e.alt_ofs as usize);
        let i = read_cstring(self.string_block, e.info_ofs as usize);

        if r == ref_ && a == alt {
            Some(i)
        } else {
            None
        }
    }
}

fn read_cstring<'a>(data: &'a [u8], mut pos: usize) -> &'a str {
    let start = pos;
    while pos < data.len() && data[pos] != 0 {
        pos += 1;
    }
    unsafe { std::str::from_utf8_unchecked(&data[start..pos]) }
}

fn extract_info(line: &str) -> &str {
    let mut tabs = 0;
    let mut start = 0;
    for (i, c) in line.char_indices() {
        if c == '\t' {
            tabs += 1;
            if tabs == 8 {
                start = i + 1;
            }
            if tabs == 9 {
                return &line[start..i];
            }
        }
    }
    ""
}

fn merge_info(base: &str, add: &str) -> String {
    if add.is_empty() {
        return base.to_string();
    }
    if base.is_empty() {
        return add.to_string();
    }

    let mut out = String::with_capacity(base.len() + add.len() + 2);
    out.push_str(base);
    out.push(';');
    out.push_str(add);
    out
}

fn parse_fields(line: &str) -> Option<(&str, u32, &str, &str)> {
    let mut c = line.split('\t');
    let chrom = c.next()?;
    let pos = c.next()?.parse::<u32>().ok()?;
    let _id = c.next()?;
    let ref_ = c.next()?;
    let alt = c.next()?;
    Some((chrom, pos, ref_, alt))
}

pub fn annotate_vcf_ani(db_ani: &Path, input_vcf: &Path, output_vcf: &Path) -> Result<()> {
    eprintln!("[annotate] Loading ANI index...");
    let ani = AniIndex::open(db_ani)?;

    eprintln!("[annotate] Annotating...");

    let fin = File::open(input_vcf)?;
    let rdr = BufReader::new(fin);

    let fout = File::create(output_vcf)?;
    let mut bw = BufWriter::new(fout);

    let mut processed = 0usize;
    let mut annotated = 0usize;
    let start = std::time::Instant::now();

    for line in rdr.lines() {
        let line = line?;

        if line.starts_with('#') {
            bw.write_all(line.as_bytes())?;
            bw.write_all(b"\n")?;
            continue;
        }

        processed += 1;

        if processed % 100_000 == 0 {
            eprintln!(
                "[annotate] processed={} annotated={} {:.3}s",
                processed,
                annotated,
                start.elapsed().as_secs_f64()
            );
        }

        if let Some((chr, pos, ref_, alt)) = parse_fields(&line) {
            if let Some(info2) = ani.lookup(chr, pos, ref_, alt) {
                let base_info = extract_info(&line);
                let merged = merge_info(base_info, info2);

                let mut cols: Vec<&str> = line.split('\t').collect();
                if cols.len() >= 8 {
                    cols[7] = &merged;
                    let new_line = cols.join("\t");
                    bw.write_all(new_line.as_bytes())?;
                    bw.write_all(b"\n")?;
                    annotated += 1;
                    continue;
                }
            }
        }

        // fallback
        bw.write_all(line.as_bytes())?;
        bw.write_all(b"\n")?;
    }

    eprintln!(
        "[annotate] DONE: processed={} annotated={} total={:.3}s",
        processed,
        annotated,
        start.elapsed().as_secs_f64()
    );

    Ok(())
}
