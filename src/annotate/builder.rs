use anyhow::Result;
use fxhash::hash64;
use kira_kv_engine::{BuildConfig, Builder};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::mem;
use std::path::Path;
use std::slice;

use super::structs::*;
use crate::chr_name_to_id;

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

fn make_key(chr_id: u8, pos: u32, rf: &str, alt: &str) -> u64 {
    let mut h = hash64(&[chr_id]);
    h ^= hash64(pos.to_le_bytes().as_ref());
    h ^= hash64(rf.as_bytes());
    h ^= hash64(alt.as_bytes());
    h
}

fn append_cstr(pool: &mut Vec<u8>, s: &str) -> u32 {
    let ofs = pool.len() as u32;
    pool.extend_from_slice(s.as_bytes());
    pool.push(0);
    ofs
}
