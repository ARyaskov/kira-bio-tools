use anyhow::Result;
use fxhash::{hash64, FxHashMap};
use kira_kv_engine::{BuildConfig, Builder};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::mem;
use std::path::Path;
use std::slice;
use std::sync::atomic::{AtomicUsize, Ordering};

use super::structs::*;
use crate::chr_name_to_id;
use crate::vcf::UnifiedVcfReader;
use crate::vcf_parser_fast::FastVcfParser;

pub fn build_ani_index_auto_v2(input: &Path, output: &Path) -> Result<()> {
    let timing = std::env::var("KIRA_BT_TIMING").is_ok();
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok();
    let start = std::time::Instant::now();

    let mut reader = UnifiedVcfReader::open(input)?;
    let _headers = reader.header()?;

    let mut entries_map: FxHashMap<u64, (AniEntry, usize)> = FxHashMap::default();
    let mut pool = Vec::<u8>::new();

    let total_variants = AtomicUsize::new(0);
    let duplicates_skipped = AtomicUsize::new(0);
    let multiallelic_count = AtomicUsize::new(0);

    let mut insertion_order = 0usize;
    let mut count = 0usize;
    let report_interval = if timing { 100_000 } else { 1_000_000 };

    while let Some(line) = reader.read_line()? {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        total_variants.fetch_add(1, Ordering::Relaxed);

        let processed = process_vcf_line_dedup(
            &line,
            &mut entries_map,
            &mut pool,
            &mut insertion_order,
            &duplicates_skipped,
            &multiallelic_count,
            debug,
        )?;

        count += processed;

        if count > 0 && count % report_interval == 0 {
            let total = total_variants.load(Ordering::Relaxed);
            let unique = entries_map.len();
            let dups = duplicates_skipped.load(Ordering::Relaxed);
            let dup_rate = (dups as f64 / total as f64) * 100.0;

            eprintln!(
                "[ani-build] Progress: {} variants → {} unique entries ({} dups, {:.1}%)",
                total, unique, dups, dup_rate
            );
        }
    }

    let total = total_variants.load(Ordering::Relaxed);
    let dups = duplicates_skipped.load(Ordering::Relaxed);
    let multi = multiallelic_count.load(Ordering::Relaxed);

    if timing {
        eprintln!("[ani-build] Parse complete:");
        eprintln!("  Total variants:     {}", total);
        eprintln!("  Unique entries:     {}", entries_map.len());
        eprintln!(
            "  Duplicates skipped: {} ({:.1}%)",
            dups,
            (dups as f64 / total.max(1) as f64) * 100.0
        );
        eprintln!(
            "  Multiallelic:       {} ({:.1}%)",
            multi,
            (multi as f64 / total.max(1) as f64) * 100.0
        );
        eprintln!("  Time: {:.3}s", start.elapsed().as_secs_f64());
    }

    let mut rows: Vec<(u64, AniEntry, usize)> = entries_map
        .into_iter()
        .map(|(key, (entry, order))| (key, entry, order))
        .collect();

    rows.sort_by_key(|(_, _, order)| *order);

    let rows: Vec<(u64, AniEntry)> = rows.into_iter().map(|(k, e, _)| (k, e)).collect();

    finalize_ani_index(rows, pool, output, timing)?;

    Ok(())
}

fn process_vcf_line_dedup(
    line: &str,
    entries_map: &mut FxHashMap<u64, (AniEntry, usize)>,
    pool: &mut Vec<u8>,
    insertion_order: &mut usize,
    duplicates_skipped: &AtomicUsize,
    multiallelic_count: &AtomicUsize,
    debug: bool,
) -> Result<usize> {
    let mut parser = FastVcfParser::new(line);

    let fields = match parser.parse_standard_fields() {
        Some(f) => f,
        None => return Ok(0),
    };

    let chr = fields.chrom;
    let pos = fields.pos.parse::<u32>().unwrap_or(0);
    let id = fields.id;
    let rf = fields.ref_allele;
    let alt = fields.alt;
    let qual = fields.qual;
    let filt = fields.filter;
    let info_raw = fields.info;

    let chr_id = match chr_name_to_id(chr) {
        Some(v) => v,
        None => return Ok(0),
    };

    let alt_alleles: Vec<&str> = alt.split(',').collect();

    if alt_alleles.len() > 1 {
        multiallelic_count.fetch_add(1, Ordering::Relaxed);
    }

    let mut processed = 0usize;

    for a in alt_alleles {
        let key = make_key(chr_id, pos, rf, a);

        if entries_map.contains_key(&key) {
            duplicates_skipped.fetch_add(1, Ordering::Relaxed);

            if debug {
                eprintln!(
                    "[ani-build] Skipping duplicate: {}:{} {}→{}",
                    chr, pos, rf, a
                );
            }
            continue;
        }

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

        let entry = AniEntry {
            chr_id,
            pos,
            ref_ofs,
            alt_ofs,
            id_ofs,
            qual_ofs,
            filter_ofs,
            info_ofs,
            info_len,
        };

        entries_map.insert(key, (entry, *insertion_order));
        *insertion_order += 1;
        processed += 1;
    }

    Ok(processed)
}

fn finalize_ani_index(
    rows: Vec<(u64, AniEntry)>,
    pool: Vec<u8>,
    output_ani: &Path,
    timing: bool,
) -> Result<()> {
    if rows.is_empty() {
        anyhow::bail!("No valid entries to index");
    }

    let mph_start = std::time::Instant::now();
    let n = rows.len();

    let (gamma, rehash_limit) = match n {
        0..=100_000 => (1.2, 16),
        100_001..=1_000_000 => (1.5, 32),
        1_000_001..=10_000_000 => (2.0, 64),
        _ => (2.5, 100),
    };

    if timing {
        eprintln!(
            "[ani-build] MPHF config: gamma={:.1}, rehash_limit={}, entries={}",
            gamma, rehash_limit, n
        );
    }

    let keys_bytes: Vec<[u8; 8]> = rows.iter().map(|(k, _)| k.to_le_bytes()).collect();

    let mph = Builder::new()
        .with_config(BuildConfig {
            gamma,
            rehash_limit,
            salt: 0x9E3779B185EBCA87,
        })
        .build(keys_bytes.iter().map(|b| b.as_slice()))?;

    if timing {
        eprintln!(
            "[ani-build] MPH build: {:.3}s",
            mph_start.elapsed().as_secs_f64()
        );

        let mph_size = mph.g.len() * mem::size_of::<u32>();
        let entries_size = n * mem::size_of::<AniEntry>();
        let pool_size = pool.len();
        let total_size = mph_size + entries_size + pool_size;

        eprintln!(
            "[ani-build] Index size: {:.2} MB (MPHF: {:.2} MB, entries: {:.2} MB, pool: {:.2} MB)",
            total_size as f64 / (1024.0 * 1024.0),
            mph_size as f64 / (1024.0 * 1024.0),
            entries_size as f64 / (1024.0 * 1024.0),
            pool_size as f64 / (1024.0 * 1024.0)
        );
    }

    let n = rows.len();
    let mut arr = vec![AniEntry::default(); n];

    for (k, e) in &rows {
        let idx = mph.index(&k.to_le_bytes()) as usize;
        if idx >= n {
            anyhow::bail!("MPHF index out of bounds: {} >= {}", idx, n);
        }
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
        eprintln!("[ani-build] DONE: {} unique variants indexed", n);
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
