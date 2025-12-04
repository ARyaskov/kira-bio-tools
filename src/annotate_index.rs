// annotate_index.rs — Full ANI v2 Builder (bcftools‑compatible)
// High‑performance annotation index used by annotate.rs

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::mem;
use std::path::Path;
use std::slice;

use anyhow::{anyhow, Result};
use fxhash::hash64;
use kira_kv_engine::{BuildConfig, Builder, Mphf};
use memmap2::{Mmap, MmapOptions};

use crate::chr_name_to_id;

// ============================================================
// ANI v2 HEADER
// ============================================================

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
    fn new(n: usize, mph: &Mphf, g_size: usize, entries_size: usize) -> Self {
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
}

// ============================================================
// FIELD descriptors
// ============================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FieldNumber {
    Zero,
    One,
    Many,
    A,
    R,
    G,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FieldType {
    Int,
    Float,
    Str,
    Flag,
}

#[derive(Debug, Clone)]
pub struct StructuredInfoField<'a> {
    pub key: &'a str,
    pub number: FieldNumber,
    pub ty: FieldType,
    pub values: Vec<&'a str>,
}

#[derive(Debug)]
pub struct AnnotationBundle<'a> {
    pub id: Option<&'a str>,
    pub qual: Option<&'a str>,
    pub filter: Option<&'a str>,
    pub info: Vec<StructuredInfoField<'a>>,
}

// ============================================================
// ANI ENTRY
// ============================================================

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

// ============================================================
// CSTRING reader
// ============================================================

pub fn read_cstring<'a>(data: &'a [u8], mut pos: usize) -> &'a str {
    let start = pos;
    while pos < data.len() && data[pos] != 0 {
        pos += 1;
    }
    unsafe { std::str::from_utf8_unchecked(&data[start..pos]) }
}

// ============================================================
// INFO PARSER
// ============================================================

fn parse_info_field<'a>(raw: &'a str) -> Vec<StructuredInfoField<'a>> {
    if raw == "." || raw.is_empty() {
        return vec![];
    }

    let mut out = Vec::new();

    for part in raw.split(';') {
        if part.is_empty() {
            continue;
        }

        let mut kv = part.splitn(2, '=');
        let key = kv.next().unwrap();
        let vals = kv.next().unwrap_or("");

        let number = if vals.is_empty() {
            FieldNumber::Zero
        } else if !vals.contains(',') {
            FieldNumber::One
        } else {
            FieldNumber::Many
        };

        let ty = if number == FieldNumber::Zero {
            FieldType::Flag
        } else if vals.parse::<i64>().is_ok() {
            FieldType::Int
        } else if vals.parse::<f64>().is_ok() {
            FieldType::Float
        } else {
            FieldType::Str
        };

        let values: Vec<&str> = if number == FieldNumber::Zero {
            vec![]
        } else {
            vals.split(',').collect()
        };

        out.push(StructuredInfoField {
            key,
            number,
            ty,
            values,
        });
    }

    out
}

// ============================================================
// KEY BUILDER
// ============================================================

fn make_key(chr_id: u8, pos: u32, rf: &str, alt: &str) -> u64 {
    let mut h = hash64(&[chr_id]);
    h ^= hash64(pos.to_le_bytes().as_ref());
    h ^= hash64(rf.as_bytes());
    h ^= hash64(alt.as_bytes());
    h
}

/// Build ANI index from a simplified annotation table (TSV)
/// Columns: CHR POS ID REF ALT QUAL FILTER INFO
pub fn build_ani_index_from_tab(input_tab: &Path, output_ani: &Path) -> Result<()> {
    let f = File::open(input_tab)?;
    let rdr = BufReader::new(f);

    let mut rows: Vec<(u64, AniEntry)> = Vec::new();
    let mut pool = Vec::<u8>::new();

    for line in rdr.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }

        let mut c = line.split('\t');

        let chr = c.next().unwrap_or("");
        let pos = c.next().unwrap_or("0").parse::<u32>().unwrap();
        let id = c.next().unwrap_or("");
        let rf = c.next().unwrap_or("");
        let alt_raw = c.next().unwrap_or("");
        let qual = c.next().unwrap_or("");
        let filt = c.next().unwrap_or("");
        let info_raw = c.next().unwrap_or("");

        let chr_id = match chr_name_to_id(chr) {
            Some(v) => v,
            None => continue,
        };

        for a in alt_raw.split(',') {
            let ref_ofs = append_cstr(&mut pool, rf);
            let alt_ofs = append_cstr(&mut pool, a);
            let id_ofs = append_cstr(&mut pool, id);
            let qual_ofs = append_cstr(&mut pool, qual);
            let filter_ofs = append_cstr(&mut pool, filt);

            let info_start = pool.len();
            let fields = parse_info_field(info_raw);
            let encoded = encode_structured_info(&fields);

            pool.extend_from_slice(encoded.as_bytes());
            pool.push(0);

            let info_ofs = info_start as u32;
            let info_len = encoded.len() as u32;

            let key = {
                let mut h = hash64(&[chr_id]);
                h ^= hash64(pos.to_le_bytes().as_ref());
                h ^= hash64(rf.as_bytes());
                h ^= hash64(a.as_bytes());
                h
            };

            rows.push((
                key,
                AniEntry {
                    chr_id,
                    pos,
                    ref_ofs,
                    alt_ofs,
                    id_ofs,
                    qual_ofs,
                    filter_ofs,
                    info_ofs,
                    info_len,
                },
            ));
        }
    }

    // --- Build MPH identical to standard build_ani_index() ---
    let keys_bytes: Vec<[u8; 8]> = rows.iter().map(|(k, _)| k.to_le_bytes()).collect();

    let mph = Builder::new()
        .with_config(BuildConfig {
            gamma: 1.2,
            rehash_limit: 16,
            salt: 0x9E3779B185EBCA87,
        })
        .build(keys_bytes.iter().map(|b| b.as_slice()))?;

    let n = rows.len();
    let mut arr = vec![AniEntry::default(); n];

    for (k, e) in &rows {
        let idx = mph.index(&k.to_le_bytes()) as usize;
        arr[idx] = *e;
    }

    // --- Write ANI ---
    let g_size = mph.g.len() * mem::size_of::<u32>();
    let ent_size = arr.len() * mem::size_of::<AniEntry>();

    let header = AniHeader::new(n, &mph, g_size, ent_size);

    let out = File::create(output_ani)?;
    let mut bw = BufWriter::new(out);

    unsafe {
        bw.write_all(slice::from_raw_parts(
            (&header as *const _) as *const u8,
            mem::size_of::<AniHeader>(),
        ))?;
    }

    let g_bytes = unsafe { slice::from_raw_parts(mph.g.as_ptr() as *const u8, g_size) };
    bw.write_all(g_bytes)?;

    let ent_bytes = unsafe { slice::from_raw_parts(arr.as_ptr() as *const u8, ent_size) };
    bw.write_all(ent_bytes)?;

    bw.write_all(&pool)?;
    bw.flush()?;

    Ok(())
}

// ============================================================
// BUILD ANI INDEX
// ============================================================

pub fn build_ani_index(input_vcf: &Path, output_ani: &Path) -> Result<()> {
    let f = File::open(input_vcf)?;
    let rdr = BufReader::new(f);

    let mut rows: Vec<(u64, AniEntry)> = Vec::new();
    let mut pool = Vec::<u8>::new();

    for line in rdr.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let mut c = line.split('\t');
        let chr = c.next().unwrap();
        let pos = c.next().unwrap().parse::<u32>().unwrap();
        let id = c.next().unwrap_or("");
        let rf = c.next().unwrap();
        let alt = c.next().unwrap();
        let qual = c.next().unwrap_or("");
        let filt = c.next().unwrap_or("");
        let info_raw = c.next().unwrap_or("");

        let chr_id = match chr_name_to_id(chr) {
            Some(v) => v,
            None => continue,
        };

        for a in alt.split(',') {
            let ref_ofs = append_cstr(&mut pool, rf);
            let alt_ofs = append_cstr(&mut pool, a);
            let id_ofs = append_cstr(&mut pool, id);
            let qual_ofs = append_cstr(&mut pool, qual);
            let filter_ofs = append_cstr(&mut pool, filt);

            let info_start = pool.len();
            let bundle = parse_info_field(info_raw);
            let encoded_info = encode_structured_info(&bundle);
            pool.extend_from_slice(encoded_info.as_bytes());
            pool.push(0);
            let info_ofs = info_start as u32;
            let info_len = encoded_info.len() as u32;

            let key = make_key(chr_id, pos, rf, a);

            rows.push((
                key,
                AniEntry {
                    chr_id,
                    pos,
                    ref_ofs,
                    alt_ofs,
                    id_ofs,
                    qual_ofs,
                    filter_ofs,
                    info_ofs,
                    info_len,
                },
            ));
        }
    }

    let keys_bytes: Vec<[u8; 8]> = rows.iter().map(|(k, _)| k.to_le_bytes()).collect();

    let mph = Builder::new()
        .with_config(BuildConfig {
            gamma: 1.2,
            rehash_limit: 16,
            salt: 0x9E3779B185EBCA87,
        })
        .build(keys_bytes.iter().map(|b| b.as_slice()))?;

    let n = rows.len();
    let mut arr = vec![AniEntry::default(); n];

    for (k, e) in &rows {
        let idx = mph.index(&k.to_le_bytes()) as usize;
        arr[idx] = *e;
    }

    // -------------------------
    // WRITE ANI FILE
    // -------------------------

    let g_size = mph.g.len() * mem::size_of::<u32>();
    let ent_size = arr.len() * mem::size_of::<AniEntry>();

    let header = AniHeader::new(n, &mph, g_size, ent_size);

    let out = File::create(output_ani)?;
    let mut bw = BufWriter::new(out);

    unsafe {
        bw.write_all(slice::from_raw_parts(
            (&header as *const _) as *const u8,
            mem::size_of::<AniHeader>(),
        ))?;
    }

    let g_bytes = unsafe { slice::from_raw_parts(mph.g.as_ptr() as *const u8, g_size) };
    bw.write_all(g_bytes)?;

    let ent_bytes = unsafe { slice::from_raw_parts(arr.as_ptr() as *const u8, ent_size) };
    bw.write_all(ent_bytes)?;

    bw.write_all(&pool)?;
    bw.flush()?;

    Ok(())
}

pub struct AniIndex {
    pub header: AniHeader,
    pub mph: Mphf,
    pub entries: Vec<AniEntry>,
    pub strings: Vec<u8>,
    mmap: Mmap,
}

impl AniIndex {
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };

        if mmap.len() < mem::size_of::<AniHeader>() {
            return Err(anyhow!("ANI file too small"));
        }

        let header: AniHeader = unsafe { *(mmap.as_ptr() as *const AniHeader) };

        if &header.magic != ANI_MAGIC {
            return Err(anyhow!("Bad ANI magic"));
        }
        if header.version != ANI_VERSION {
            return Err(anyhow!("ANI version mismatch"));
        }

        // MPH g[]
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

        // entries[]
        let ent_start = header.off_entries as usize;
        let ent_end = ent_start + header.n_entries as usize * mem::size_of::<AniEntry>();

        let mut entries = Vec::new();
        for chunk in mmap[ent_start..ent_end].chunks_exact(mem::size_of::<AniEntry>()) {
            let e: AniEntry = unsafe { *(chunk.as_ptr() as *const AniEntry) };
            entries.push(e);
        }

        // strings
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

    pub fn lookup_full(
        &self,
        chr: &str,
        pos: u32,
        rf: &str,
        alt_raw: &str,
    ) -> Option<(AnnotationBundle<'_>, Vec<&'_ str>)> {
        let chr_id = chr_name_to_id(chr)? as u8;

        let mut h = hash64(&[chr_id]);
        h ^= hash64(pos.to_le_bytes().as_ref());
        h ^= hash64(rf.as_bytes());
        h ^= hash64(alt_raw.as_bytes());

        let idx = self.mph.index(&h.to_le_bytes()) as usize;
        if idx >= self.entries.len() {
            return None;
        }

        let e = &self.entries[idx];

        let id_str = read_cstring(&self.strings, e.id_ofs as usize);
        let qual_str = read_cstring(&self.strings, e.qual_ofs as usize);
        let filter_str = read_cstring(&self.strings, e.filter_ofs as usize);
        let info_str = read_cstring(&self.strings, e.info_ofs as usize);
        let alt_str = read_cstring(&self.strings, e.alt_ofs as usize);
        let rf_str = read_cstring(&self.strings, e.ref_ofs as usize);

        if alt_str != alt_raw {
            return None;
        }
        if rf_str != rf {
            return None;
        }

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

        let ann_alt_list: Vec<&str> = alt_str.split(',').collect();

        Some((bundle, ann_alt_list))
    }
}

pub fn encode_structured_info(f: &[StructuredInfoField]) -> String {
    let mut out = String::new();

    for (i, fld) in f.iter().enumerate() {
        if i > 0 {
            out.push(';');
        }

        match fld.number {
            FieldNumber::Zero => {
                out.push_str(fld.key);
            }

            FieldNumber::One => {
                out.push_str(fld.key);
                out.push('=');
                if let Some(v) = fld.values.first() {
                    out.push_str(v);
                }
            }

            FieldNumber::Many | FieldNumber::A | FieldNumber::R | FieldNumber::G => {
                out.push_str(fld.key);
                out.push('=');
                out.push_str(&fld.values.join(","));
            }
        }
    }

    out
}

fn append_cstr(pool: &mut Vec<u8>, s: &str) -> u32 {
    let ofs = pool.len() as u32;
    pool.extend_from_slice(s.as_bytes());
    pool.push(0);
    ofs
}
