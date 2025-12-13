use anyhow::Result;
use fxhash::FxHashMap;
use kira_kv_engine::{BuildConfig, Builder};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::mem;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::annotate::structs::{AniEntry, AniHeader, TabSchema, ANI_MAGIC, ANI_VERSION};
use crate::util::{append_cstr, chr_name_to_id, url_encode_info_value};
use crate::vcf::simd::SimdVcfParser;
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

        let processed = process_vcf_line_multiallelic_simd(
            line.as_bytes(),
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
        eprintln!("  Parse time: {:.3}s", start.elapsed().as_secs_f64());
    }

    let mut rows: Vec<(u64, AniEntry, usize)> = entries_map
        .into_iter()
        .map(|(key, (entry, order))| (key, entry, order))
        .collect();

    rows.sort_unstable_by_key(|(_, _, order)| *order);

    let rows: Vec<(u64, AniEntry)> = rows.into_iter().map(|(k, e, _)| (k, e)).collect();

    finalize_ani_index(rows, pool, output, timing)?;

    Ok(())
}

pub fn build_ani_index_from_tab(input: &Path, output: &Path, columns: Option<&str>) -> Result<()> {
    let timing = std::env::var("KIRA_BT_TIMING").is_ok();
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok();
    let start = std::time::Instant::now();

    let schema = TabSchema::parse(input, columns)?;

    if debug {
        eprintln!("[ani-build-tab] Schema detected:");
        eprintln!("  CHROM={}, POS={}", schema.chrom_idx, schema.pos_idx);
        if let Some(i) = schema.ref_idx {
            eprintln!("  REF={}", i);
        }
        if let Some(i) = schema.alt_idx {
            eprintln!("  ALT={}", i);
        }
        eprintln!("  INFO columns: {}", schema.info_cols.len());
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
        eprintln!("  Parse time: {:.3}s", start.elapsed().as_secs_f64());
    }

    let mut rows: Vec<(u64, AniEntry, usize)> = entries_map
        .into_iter()
        .map(|(key, (entry, order))| (key, entry, order))
        .collect();

    rows.sort_unstable_by_key(|(_, _, order)| *order);

    let rows: Vec<(u64, AniEntry)> = rows.into_iter().map(|(k, e, _)| (k, e)).collect();

    finalize_ani_index(rows, pool, output, timing)?;

    Ok(())
}

fn process_vcf_line_multiallelic_simd(
    line: &[u8],
    entries_map: &mut FxHashMap<u64, (AniEntry, usize)>,
    pool: &mut Vec<u8>,
    insertion_order: &mut usize,
    duplicates_skipped: &AtomicUsize,
    multiallelic_count: &AtomicUsize,
    debug: bool,
) -> Result<usize> {
    let fields = match SimdVcfParser::parse_fields(line) {
        Some(f) => f,
        None => return Ok(0),
    };

    let chr = fields.chrom;
    let chr_id = match chr_name_to_id(chr) {
        Some(v) => v,
        None => return Ok(0),
    };

    let pos = match fields.pos.parse::<u32>() {
        Ok(p) => p,
        Err(_) => return Ok(0),
    };

    let rf = fields.ref_allele.trim();
    let alt = fields.alt.trim();

    if rf.is_empty() || rf == "." || alt.is_empty() || alt == "." {
        return Ok(0);
    }

    let alt_alleles: Vec<&str> = alt.split(',').collect();

    if alt_alleles.len() > 1 {
        multiallelic_count.fetch_add(1, Ordering::Relaxed);
    }

    let id = if fields.id != "." && !fields.id.is_empty() {
        fields.id
    } else {
        "."
    };

    let qual = if fields.qual != "." && !fields.qual.is_empty() {
        fields.qual
    } else {
        "."
    };

    let filter = if fields.filter != "." && !fields.filter.is_empty() {
        fields.filter
    } else {
        "."
    };

    let ref_ofs = append_cstr(pool, rf);
    let id_ofs = append_cstr(pool, id);
    let qual_ofs = append_cstr(pool, qual);
    let filter_ofs = append_cstr(pool, filter);

    let info_str = fields.info.trim();

    for (alt_idx, alt_single) in alt_alleles.iter().enumerate() {
        let key = make_position_key(chr_id, pos, rf, alt_single);

        if entries_map.contains_key(&key) {
            duplicates_skipped.fetch_add(1, Ordering::Relaxed);

            if debug {
                eprintln!(
                    "[ani-build] Overwriting duplicate: {}:{} {} {}",
                    chr, pos, rf, alt_single
                );
            }
        } else if debug && *insertion_order < 10 {
            eprintln!(
                "[ani-build] Creating entry #{}: {}:{} {} {} (key={:016x})",
                insertion_order, chr, pos, rf, alt_single, key
            );
        }

        let alt_ofs = append_cstr(pool, alt_single);

        let info_ofs = if info_str != "." && !info_str.is_empty() && alt_alleles.len() > 1 {
            let split_info = split_info_for_allele(info_str, alt_idx, alt_alleles.len());
            let encoded = url_encode_info_value(&split_info);
            append_cstr(pool, &encoded) as u32
        } else if info_str != "." && !info_str.is_empty() {
            let encoded = url_encode_info_value(info_str);
            append_cstr(pool, &encoded) as u32
        } else {
            append_cstr(pool, ".") as u32
        };

        let entry = AniEntry {
            chr_id,
            pos,
            ref_ofs: ref_ofs as u32,
            alt_ofs: alt_ofs as u32,
            id_ofs: id_ofs as u32,
            qual_ofs: qual_ofs as u32,
            filter_ofs: filter_ofs as u32,
            info_ofs,
            info_len: 0,
        };

        entries_map.insert(key, (entry, *insertion_order));
        *insertion_order += 1;
    }

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
    let id_ofs = append_cstr(pool, id);
    let qual_ofs = append_cstr(pool, qual);
    let filter_ofs = append_cstr(pool, filt);

    let mut base_info_str = String::new();

    if let Some(info_idx) = schema.info_start {
        if let Some(info_str) = parts.get(info_idx) {
            let info_str = info_str.trim();
            if !info_str.is_empty() && info_str != "." {
                base_info_str.push_str(info_str);
            }
        }
    }

    for col in &schema.info_cols {
        if let Some(val_str) = parts.get(col.index) {
            let val_str = val_str.trim();
            if !val_str.is_empty() && val_str != "." {
                if !base_info_str.is_empty() {
                    base_info_str.push(';');
                }

                base_info_str.push_str(&col.key);
                base_info_str.push('=');
                base_info_str.push_str(val_str);
            }
        }
    }

    for (alt_idx, alt_single) in alt_alleles.iter().enumerate() {
        let key = make_position_key(chr_id, pos, rf, alt_single);

        if entries_map.contains_key(&key) {
            duplicates_skipped.fetch_add(1, Ordering::Relaxed);

            if debug {
                eprintln!(
                    "[ani-build-tab] Overwriting duplicate: {}:{} {} {}",
                    chr, pos, rf, alt_single
                );
            }
        }

        let alt_ofs = append_cstr(pool, alt_single);

        if debug {
            eprintln!(
                "[ani-build-tab] Processing alt_idx={} alt_single={} num_alts={} base_info='{}'",
                alt_idx,
                alt_single,
                alt_alleles.len(),
                if base_info_str.len() > 50 {
                    &base_info_str[..50]
                } else {
                    &base_info_str
                }
            );
        }

        let info_ofs = if !base_info_str.is_empty() && base_info_str != "." {
            let final_info = if alt_alleles.len() > 1 {
                if debug {
                    eprintln!("[ani-build-tab] Calling split_info_for_allele");
                }
                split_info_for_allele(&base_info_str, alt_idx, alt_alleles.len())
            } else {
                base_info_str.clone()
            };
            if debug {
                eprintln!(
                    "[ani-build-tab] final_info='{}'",
                    if final_info.len() > 50 {
                        &final_info[..50]
                    } else {
                        &final_info
                    }
                );
            }
            let encoded = url_encode_info_value(&final_info);
            append_cstr(pool, &encoded) as u32
        } else {
            append_cstr(pool, ".") as u32
        };

        let entry = AniEntry {
            chr_id,
            pos,
            ref_ofs: ref_ofs as u32,
            alt_ofs: alt_ofs as u32,
            id_ofs: id_ofs as u32,
            qual_ofs: qual_ofs as u32,
            filter_ofs: filter_ofs as u32,
            info_ofs,
            info_len: 0,
        };

        entries_map.insert(key, (entry, *insertion_order));
        *insertion_order += 1;
    }

    Ok(alt_alleles.len())
}

fn split_info_for_allele(info: &str, alt_idx: usize, num_alts: usize) -> String {
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok();
    let mut result = Vec::new();

    for field in info.split(';') {
        if field.is_empty() {
            continue;
        }

        if let Some(eq_pos) = field.find('=') {
            let key = &field[..eq_pos];
            let value = &field[eq_pos + 1..];

            if value.contains(',') {
                let values: Vec<&str> = value.split(',').collect();

                if debug {
                    eprintln!(
                        "[split] key={} values.len={} num_alts={} alt_idx={} value={}",
                        key,
                        values.len(),
                        num_alts,
                        alt_idx,
                        value
                    );
                }

                if values.len() == num_alts {
                    let selected = values[alt_idx];
                    result.push(format!("{}={}", key, selected));
                    if debug {
                        eprintln!("[split] Number=A: selected={}", selected);
                    }
                } else if values.len() == num_alts + 1 {
                    let ref_val = values[0];
                    let alt_val = values[alt_idx + 1];
                    result.push(format!("{}={},{}", key, ref_val, alt_val));
                    if debug {
                        eprintln!("[split] Number=R: ref={} alt={}", ref_val, alt_val);
                    }
                } else {
                    result.push(field.to_string());
                    if debug {
                        eprintln!("[split] Unknown number, keeping as-is");
                    }
                }
            } else {
                result.push(field.to_string());
            }
        } else {
            result.push(field.to_string());
        }
    }

    if result.is_empty() {
        ".".to_string()
    } else {
        result.join(";")
    }
}

fn make_position_key(chr_id: u8, pos: u32, rf: &str, alt: &str) -> u64 {
    use fxhash::hash64;
    let base_h = (chr_id as u64) << 32 | (pos as u64);
    let ref_hash = hash64(rf.as_bytes());
    let alt_hash = hash64(alt.as_bytes());
    let mut h = base_h;
    h ^= ref_hash;
    h ^= alt_hash;

    let debug = std::env::var("KIRA_BT_DEBUG_HASH").is_ok();
    if debug {
        eprintln!(
            "[HASH-BUILD] {}:{} {:?} {:?} -> base={:016x} ref={:016x} alt={:016x} key={:016x}",
            chr_id, pos, rf, alt, base_h, ref_hash, alt_hash, h
        );
    }

    h
}

fn finalize_ani_index(
    rows: Vec<(u64, AniEntry)>,
    pool: Vec<u8>,
    output: &Path,
    timing: bool,
) -> Result<()> {
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
        eprintln!("[ani-build] First 5 keys for verification:");
        for (i, (k, e)) in rows.iter().enumerate().take(5) {
            eprintln!("  [{}] key={:016x} chr={} pos={}", i, k, e.chr_id, e.pos);
        }
    }

    let keys_bytes: Vec<[u8; 8]> = rows.iter().map(|(k, _)| k.to_le_bytes()).collect();

    if timing {
        let mut key_set = std::collections::HashSet::new();
        let mut dup_count = 0;
        for (i, key) in keys_bytes.iter().enumerate() {
            let k = u64::from_le_bytes(*key);
            if !key_set.insert(k) {
                dup_count += 1;
                if dup_count <= 5 {
                    let (_, e) = &rows[i];
                    eprintln!(
                        "[ani-build] WARNING: Duplicate key {:016x} for chr={} pos={}",
                        k, e.chr_id, e.pos
                    );
                }
            }
        }
        if dup_count > 0 {
            eprintln!(
                "[ani-build] ERROR: Found {} duplicate keys! MPH will not work correctly!",
                dup_count
            );
        }
    }

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

    let entries: Vec<AniEntry> = rows.into_iter().map(|(_, entry)| entry).collect();

    if timing {
        eprintln!("[ani-build] Verifying MPH lookups (checking 100 random entries):");
        let mut errors = 0;
        let check_count = 100.min(entries.len());
        let step = entries.len() / check_count;

        for i in (0..entries.len()).step_by(step.max(1)).take(check_count) {
            let key = u64::from_le_bytes(keys_bytes[i]);
            let idx = mph.index(&keys_bytes[i]) as usize;
            let retrieved = &entries[idx];
            let expected = &entries[i];

            if retrieved.chr_id != expected.chr_id || retrieved.pos != expected.pos {
                errors += 1;
                if errors <= 5 {
                    eprintln!(
                        "  [{}] ERROR: key={:016x} -> mph_idx={} -> chr={} pos={} (expected chr={} pos={})",
                        i, key, idx, retrieved.chr_id, retrieved.pos, expected.chr_id, expected.pos
                    );
                }
            } else if i < 3 {
                eprintln!(
                    "  [{}] OK: key={:016x} -> mph_idx={} -> chr={} pos={}",
                    i, key, idx, retrieved.chr_id, retrieved.pos
                );
            }
        }

        if errors > 0 {
            eprintln!(
                "[ani-build] ERROR: MPH verification FAILED with {} errors out of {} checks!",
                errors, check_count
            );
            eprintln!("[ani-build] This means the index will NOT work correctly!");
        } else {
            eprintln!(
                "[ani-build] MPH verification: All {} checks passed ✓",
                check_count
            );
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
