use anyhow::Result;
use fxhash::FxHashMap;
use kira_kv_engine::Mphf;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::mem;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::annotate::structs::{AniEntry, AniHeader, TabSchema, ANI_MAGIC, ANI_VERSION};
use crate::util::{append_cstr, chr_name_to_id, url_encode_info_value};
use crate::vcf::UnifiedVcfReader;

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

        let processed = process_vcf_line_multiallelic(
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

pub fn build_ani_index_from_tab(input: &Path, output: &Path, columns: Option<&str>) -> Result<()> {
    let timing = std::env::var("KIRA_BT_TIMING").is_ok();
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok();
    let start = std::time::Instant::now();

    let schema = TabSchema::parse(input, columns)?;

    if timing {
        eprintln!("[ani-build-tab] Schema: {:?}", schema);
    }

    let file = File::open(input)?;
    let reader = BufReader::new(file);

    let mut entries_map: FxHashMap<u64, (AniEntry, usize)> = FxHashMap::default();
    let mut pool = Vec::<u8>::new();

    let total_variants = AtomicUsize::new(0);
    let duplicates_skipped = AtomicUsize::new(0);
    let multiallelic_count = AtomicUsize::new(0);

    let mut insertion_order = 0usize;
    let mut count = 0usize;
    let report_interval = if timing { 100_000 } else { 1_000_000 };

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        total_variants.fetch_add(1, Ordering::Relaxed);

        let processed = process_tab_line_multiallelic(
            &line,
            &schema,
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
                "[ani-build-tab] Progress: {} variants → {} unique entries ({} dups, {:.1}%)",
                total, unique, dups, dup_rate
            );
        }
    }

    let total = total_variants.load(Ordering::Relaxed);
    let dups = duplicates_skipped.load(Ordering::Relaxed);
    let multi = multiallelic_count.load(Ordering::Relaxed);

    if timing {
        eprintln!("[ani-build-tab] Parse complete:");
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

fn process_vcf_line_multiallelic(
    line: &str,
    entries_map: &mut FxHashMap<u64, (AniEntry, usize)>,
    pool: &mut Vec<u8>,
    insertion_order: &mut usize,
    duplicates_skipped: &AtomicUsize,
    multiallelic_count: &AtomicUsize,
    debug: bool,
) -> Result<usize> {
    let parts: Vec<&str> = line.split('\t').collect();
    if parts.len() < 8 {
        return Ok(0);
    }

    let chr = parts[0];
    let pos_str = parts[1];
    let rf = parts[3].trim();
    let alt = parts[4].trim();

    if rf.is_empty() || rf == "." || alt.is_empty() || alt == "." {
        return Ok(0);
    }

    let pos = pos_str.parse::<u32>().unwrap_or(0);
    let chr_id = match chr_name_to_id(chr) {
        Some(v) => v,
        None => return Ok(0),
    };

    let alt_alleles: Vec<&str> = alt.split(',').collect();

    if alt_alleles.len() > 1 {
        multiallelic_count.fetch_add(1, Ordering::Relaxed);
    }

    let key = make_position_key(chr_id, pos, rf);

    if entries_map.contains_key(&key) {
        duplicates_skipped.fetch_add(1, Ordering::Relaxed);

        if debug {
            eprintln!("[ani-build] Overwriting duplicate: {}:{} {}", chr, pos, rf);
        }
    }

    let id = parts[2].trim();
    let qual = parts[5].trim();
    let filter = parts[6].trim();

    let ref_ofs = append_cstr(pool, rf);
    let alt_ofs = append_cstr(pool, alt);

    let id_ofs = if id != "." && !id.is_empty() {
        append_cstr(pool, id)
    } else {
        append_cstr(pool, ".")
    };

    let qual_ofs = if qual != "." && !qual.is_empty() {
        append_cstr(pool, qual)
    } else {
        append_cstr(pool, ".")
    };

    let filter_ofs = if filter != "." && !filter.is_empty() {
        append_cstr(pool, filter)
    } else {
        append_cstr(pool, ".")
    };

    let info_str = parts[7].trim();
    let info_ofs = if info_str != "." && !info_str.is_empty() {
        let encoded = url_encode_info_value(info_str);
        append_cstr(pool, &encoded)
    } else {
        append_cstr(pool, ".")
    };

    let entry = AniEntry {
        chr_id,
        pos,
        ref_ofs: ref_ofs as u32,
        alt_ofs: alt_ofs as u32,
        id_ofs: id_ofs as u32,
        qual_ofs: qual_ofs as u32,
        filter_ofs: filter_ofs as u32,
        info_ofs: info_ofs as u32,
        info_len: 0,
    };

    entries_map.insert(key, (entry, *insertion_order));
    *insertion_order += 1;

    Ok(alt_alleles.len())
}

fn process_tab_line_multiallelic(
    line: &str,
    schema: &TabSchema,
    entries_map: &mut FxHashMap<u64, (AniEntry, usize)>,
    pool: &mut Vec<u8>,
    insertion_order: &mut usize,
    duplicates_skipped: &AtomicUsize,
    multiallelic_count: &AtomicUsize,
    debug: bool,
) -> Result<usize> {
    let parts: Vec<&str> = line.split('\t').collect();

    let chr = parts
        .get(schema.chrom_idx)
        .ok_or_else(|| anyhow::anyhow!("Missing CHROM"))?;
    let pos_str = parts
        .get(schema.pos_idx)
        .ok_or_else(|| anyhow::anyhow!("Missing POS"))?;
    let pos = pos_str.parse::<u32>().unwrap_or(0);

    let rf = schema
        .ref_idx
        .and_then(|i| parts.get(i))
        .unwrap_or(&".")
        .trim();
    let alt = schema
        .alt_idx
        .and_then(|i| parts.get(i))
        .unwrap_or(&".")
        .trim();

    if rf.is_empty() || rf == "." || alt.is_empty() || alt == "." {
        return Ok(0);
    }

    let chr_id = match chr_name_to_id(chr) {
        Some(v) => v,
        None => return Ok(0),
    };

    let alt_alleles: Vec<&str> = alt.split(',').collect();

    if alt_alleles.len() > 1 {
        multiallelic_count.fetch_add(1, Ordering::Relaxed);
    }

    let key = make_position_key(chr_id, pos, rf);

    if entries_map.contains_key(&key) {
        duplicates_skipped.fetch_add(1, Ordering::Relaxed);

        if debug {
            eprintln!(
                "[ani-build-tab] Overwriting duplicate: {}:{} {}",
                chr, pos, rf
            );
        }
    }

    let id = schema
        .id_idx
        .and_then(|i| parts.get(i))
        .map(|s| s.trim())
        .filter(|s| !s.is_empty() && *s != ".")
        .unwrap_or(".");

    let qual = schema
        .qual_idx
        .and_then(|i| parts.get(i))
        .map(|s| s.trim())
        .filter(|s| !s.is_empty() && *s != ".")
        .unwrap_or(".");

    let filt = schema
        .filter_idx
        .and_then(|i| parts.get(i))
        .map(|s| s.trim())
        .filter(|s| !s.is_empty() && *s != ".")
        .unwrap_or(".");

    let ref_ofs = append_cstr(pool, rf);
    let alt_ofs = append_cstr(pool, alt);
    let id_ofs = append_cstr(pool, id);
    let qual_ofs = append_cstr(pool, qual);
    let filter_ofs = append_cstr(pool, filt);

    let info_start_ofs = pool.len();

    if let Some(info_idx) = schema.info_start {
        if let Some(info_str) = parts.get(info_idx) {
            let info_str = info_str.trim();
            if !info_str.is_empty() && info_str != "." {
                let encoded = url_encode_info_value(info_str);
                pool.extend_from_slice(encoded.as_bytes());
            }
        }
    }

    for col in &schema.info_cols {
        if let Some(val_str) = parts.get(col.index) {
            let val_str = val_str.trim();
            if !val_str.is_empty() && val_str != "." {
                if pool.len() > info_start_ofs {
                    pool.push(b';');
                }

                pool.extend_from_slice(col.key.as_bytes());
                pool.push(b'=');

                let encoded = url_encode_info_value(val_str);
                pool.extend_from_slice(encoded.as_bytes());
            }
        }
    }

    pool.push(0);

    let entry = AniEntry {
        chr_id,
        pos,
        ref_ofs: ref_ofs as u32,
        alt_ofs: alt_ofs as u32,
        id_ofs: id_ofs as u32,
        qual_ofs: qual_ofs as u32,
        filter_ofs: filter_ofs as u32,
        info_ofs: info_start_ofs as u32,
        info_len: 0,
    };

    entries_map.insert(key, (entry, *insertion_order));
    *insertion_order += 1;

    Ok(alt_alleles.len())
}

fn make_position_key(chr: u8, pos: u32, rf: &str) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h: u64 = ((chr as u64) << 56) | ((pos as u64) << 24);

    let mut hasher = fxhash::FxHasher::default();
    rf.hash(&mut hasher);
    let rf_hash = hasher.finish();

    h ^= (rf_hash & 0x00FFFFFF);
    h
}

fn finalize_ani_index(
    rows: Vec<(u64, AniEntry)>,
    pool: Vec<u8>,
    output: &Path,
    timing: bool,
) -> Result<()> {
    use kira_kv_engine::{BuildConfig, Builder};

    if rows.is_empty() {
        anyhow::bail!("No valid entries to index");
    }

    let n = rows.len();

    let mph_start = std::time::Instant::now();

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
            "[ani-build] MPH construction: {:.3}s",
            mph_start.elapsed().as_secs_f64()
        );
    }

    let mut entries = vec![
        AniEntry {
            chr_id: 0,
            pos: 0,
            ref_ofs: 0,
            alt_ofs: 0,
            id_ofs: 0,
            qual_ofs: 0,
            filter_ofs: 0,
            info_ofs: 0,
            info_len: 0,
        };
        n
    ];

    for (k, entry) in rows {
        let key_bytes = k.to_le_bytes();
        let idx = mph.index(&key_bytes) as usize;
        if idx < n {
            entries[idx] = entry;
        }
    }

    let header = AniHeader {
        magic: ANI_MAGIC,
        version: ANI_VERSION,
        n_entries: n as u64,
        mph_m: mph.m as u64,
        mph_salt: mph.salt,
        off_mph_g: 0,
        off_entries: 0,
        off_strings: 0,
    };

    let hdr_size = mem::size_of::<AniHeader>();
    let g_size = mph.g.len() * 4;
    let ent_size = n * mem::size_of::<AniEntry>();
    let str_size = pool.len();

    let mut header = header;
    header.off_mph_g = hdr_size as u64;
    header.off_entries = (hdr_size + g_size) as u64;
    header.off_strings = (hdr_size + g_size + ent_size) as u64;

    let write_start = std::time::Instant::now();
    let mut file = File::create(output)?;

    let hdr_bytes: &[u8] =
        unsafe { std::slice::from_raw_parts(&header as *const _ as *const u8, hdr_size) };
    file.write_all(hdr_bytes)?;

    for g in &mph.g {
        file.write_all(&g.to_le_bytes())?;
    }

    for e in &entries {
        let e_bytes: &[u8] =
            unsafe { std::slice::from_raw_parts(e as *const _ as *const u8, mem::size_of_val(e)) };
        file.write_all(e_bytes)?;
    }

    file.write_all(&pool)?;

    if timing {
        eprintln!(
            "[ani-build] Write to disk: {:.3}s",
            write_start.elapsed().as_secs_f64()
        );
        eprintln!(
            "[ani-build] Total ANI size: {} bytes",
            hdr_size + g_size + ent_size + str_size
        );
    }

    Ok(())
}
