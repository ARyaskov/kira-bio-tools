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

/// Detect format: VCF starts with ##, TAB is plain columns
fn detect_annotation_format(path: &Path) -> Result<AnnotationFormat> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut first_line = String::new();
    reader.read_line(&mut first_line)?;

    if first_line.starts_with("##") || first_line.starts_with("#CHROM") {
        Ok(AnnotationFormat::Vcf)
    } else {
        Ok(AnnotationFormat::Tab)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnnotationFormat {
    Vcf,
    Tab,
}

/// Unified entry point - auto-detects format
pub fn build_ani_index_auto(input: &Path, output: &Path) -> Result<()> {
    let format = detect_annotation_format(input)?;

    match format {
        AnnotationFormat::Vcf => build_ani_index_from_vcf(input, output),
        AnnotationFormat::Tab => build_ani_index_from_tab(input, output),
    }
}

/// Build ANI index from VCF format
fn build_ani_index_from_vcf(input_vcf: &Path, output_ani: &Path) -> Result<()> {
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

    finalize_ani_index(rows, pool, output_ani)
}

/// Build ANI index from TAB format (bcftools-compatible)
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

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            eprintln!(
                "Warning: skipping malformed TAB line (expected 9 cols): {}",
                line
            );
            continue;
        }

        // TAB format: CHROM POS REF ALT ID QUAL T_INT T_FLOAT T_STR/FLAG
        let chr = cols[0];
        let pos = cols[1].parse::<u32>().unwrap_or(0);
        let rf = cols[2];
        let alt_raw = cols[3];
        let id = cols[4];
        let qual = cols[5];
        let t_int = cols[6];
        let t_float = cols[7];
        let t_str_or_flag = cols[8];

        let chr_id = match chr_name_to_id(chr) {
            Some(v) => v,
            None => {
                eprintln!("Warning: unknown chromosome: {}", chr);
                continue;
            }
        };

        // Build INFO field from TAB columns
        let info_raw = build_info_from_tab_columns(t_int, t_float, t_str_or_flag);

        // Process each ALT allele separately
        for a in alt_raw.split(',') {
            let ref_ofs = append_cstr(&mut pool, rf);
            let alt_ofs = append_cstr(&mut pool, a);
            let id_ofs = append_cstr(&mut pool, id);
            let qual_ofs = append_cstr(&mut pool, qual);
            let filter_ofs = append_cstr(&mut pool, ".");

            let info_start = pool.len();
            let fields = parse_info_field(&info_raw);
            let encoded = encode_structured_info(&fields);

            pool.extend_from_slice(encoded.as_bytes());
            pool.push(0);

            let info_ofs = info_start as u32;
            let info_len = encoded.len() as u32;

            // KEY: hash(chr_id, pos, ref, SINGLE_ALT)
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

    finalize_ani_index(rows, pool, output_ani)
}

/// Convert TAB columns 7-9 into INFO string
/// Column 7: T_INT (. or values)
/// Column 8: T_FLOAT (. or values)
/// Column 9: T_STR value OR flag (0/1)
fn build_info_from_tab_columns(t_int: &str, t_float: &str, t_str_or_flag: &str) -> String {
    let mut info_parts = Vec::new();

    // T_INT field
    if t_int != "." && !t_int.is_empty() {
        info_parts.push(format!("T_INT={}", t_int));
    }

    // T_FLOAT field
    if t_float != "." && !t_float.is_empty() {
        info_parts.push(format!("T_FLOAT={}", t_float));
    }

    // T_STR field or flag
    // If value is "1" or "0" → this is NOT T_STR, it's just a flag indicator
    // Only add T_STR if it's an actual string value
    if t_str_or_flag != "."
        && t_str_or_flag != "0"
        && t_str_or_flag != "1"
        && !t_str_or_flag.is_empty()
    {
        info_parts.push(format!("T_STR={}", t_str_or_flag));
    }

    if info_parts.is_empty() {
        ".".to_string()
    } else {
        info_parts.join(";")
    }
}

fn finalize_ani_index(rows: Vec<(u64, AniEntry)>, pool: Vec<u8>, output_ani: &Path) -> Result<()> {
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

    eprintln!("[ani] DONE: {} variants indexed", n);

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
