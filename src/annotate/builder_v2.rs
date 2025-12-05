use anyhow::Result;
use fxhash::hash64;
use kira_kv_engine::{BuildConfig, Builder};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::mem;
use std::path::Path;
use std::slice;

use super::reader::VcfAnnotationReader;
use super::structs::*;
use crate::chr_name_to_id;

pub fn build_ani_index_auto_v2(input: &Path, output: &Path) -> Result<()> {
    let timing = std::env::var("KIRA_BT_TIMING").is_ok();
    let start = std::time::Instant::now();

    let mut reader = VcfAnnotationReader::open(input)?;

    let mut rows: Vec<(u64, AniEntry)> = Vec::new();
    let mut pool = Vec::<u8>::new();

    while let Some(line) = reader.read_line()? {
        if line.starts_with("#CHROM") {
            break;
        }
        if !line.starts_with('#') {
            process_vcf_line(&line, &mut rows, &mut pool)?;
        }
    }

    let mut count = 0usize;
    while let Some(line) = reader.read_line()? {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        process_vcf_line(&line, &mut rows, &mut pool)?;
        count += 1;

        if timing && count % 100_000 == 0 {
            eprintln!("[ani-build] Processed {} variants...", count);
        }
    }

    if timing {
        eprintln!(
            "[ani-build] Parse complete: {} variants in {:.3}s",
            count,
            start.elapsed().as_secs_f64()
        );
    }

    finalize_ani_index(rows, pool, output, timing)?;

    Ok(())
}

fn process_vcf_line(line: &str, rows: &mut Vec<(u64, AniEntry)>, pool: &mut Vec<u8>) -> Result<()> {
    let mut cols = line.split('\t');
    let chr = cols.next().unwrap();
    let pos = cols.next().unwrap().parse::<u32>().unwrap();
    let id = cols.next().unwrap_or("");
    let rf = cols.next().unwrap();
    let alt = cols.next().unwrap();
    let qual = cols.next().unwrap_or("");
    let filt = cols.next().unwrap_or("");
    let info_raw = cols.next().unwrap_or("");

    let chr_id = match chr_name_to_id(chr) {
        Some(v) => v,
        None => return Ok(()),
    };

    for a in alt.split(',') {
        let ref_ofs = append_cstr(pool, rf);
        let alt_ofs = append_cstr(pool, a);
        let id_ofs = append_cstr(pool, id);
        let qual_ofs = append_cstr(pool, qual);
        let filter_ofs = append_cstr(pool, filt);

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

    Ok(())
}

fn finalize_ani_index(
    rows: Vec<(u64, AniEntry)>,
    pool: Vec<u8>,
    output_ani: &Path,
    timing: bool,
) -> Result<()> {
    let mph_start = std::time::Instant::now();

    let keys_bytes: Vec<[u8; 8]> = rows.iter().map(|(k, _)| k.to_le_bytes()).collect();

    let mph = Builder::new()
        .with_config(BuildConfig {
            gamma: 1.2,
            rehash_limit: 16,
            salt: 0x9E3779B185EBCA87,
        })
        .build(keys_bytes.iter().map(|b| b.as_slice()))?;

    if timing {
        eprintln!(
            "[ani-build] MPH build: {:.3}s",
            mph_start.elapsed().as_secs_f64()
        );
    }

    let n = rows.len();
    let mut arr = vec![AniEntry::default(); n];

    for (k, e) in &rows {
        let idx = mph.index(&k.to_le_bytes()) as usize;
        arr[idx] = *e;
    }

    let g_size = mph.g.len() * mem::size_of::<u32>();
    let ent_size = arr.len() * mem::size_of::<AniEntry>();

    let header = AniHeader::new(n, &mph, g_size, ent_size);

    let write_start = std::time::Instant::now();

    let out = File::create(output_ani)?;
    let mut bw = BufWriter::with_capacity(64 * 1024 * 1024, out);

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

    if timing {
        eprintln!(
            "[ani-build] Write: {:.3}s",
            write_start.elapsed().as_secs_f64()
        );
        eprintln!("[ani-build] DONE: {} variants indexed", n);
    }

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
