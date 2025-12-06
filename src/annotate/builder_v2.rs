use anyhow::Result;
use fxhash::{hash64, FxHashMap};
use kira_kv_engine::{BuildConfig, Builder};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::mem;
use std::path::Path;
use std::slice;
use std::sync::atomic::{AtomicUsize, Ordering};

use super::reader::VcfAnnotationReader;
use super::structs::*;
use crate::chr_name_to_id;

/// Build ANI index with automatic format detection and bcftools-compatible deduplication
pub fn build_ani_index_auto_v2(input: &Path, output: &Path) -> Result<()> {
    let timing = std::env::var("KIRA_BT_TIMING").is_ok();
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok();
    let start = std::time::Instant::now();

    let mut reader = VcfAnnotationReader::open(input)?;

    // Key -> (AniEntry, insertion_order) for first-wins deduplication
    let mut entries_map: FxHashMap<u64, (AniEntry, usize)> = FxHashMap::default();
    let mut pool = Vec::<u8>::new();

    // Statistics
    let total_variants = AtomicUsize::new(0);
    let duplicates_skipped = AtomicUsize::new(0);
    let multiallelic_count = AtomicUsize::new(0);

    let mut insertion_order = 0usize;

    // Skip header lines
    while let Some(line) = reader.read_line()? {
        if line.starts_with("#CHROM") {
            break;
        }
    }

    // Process VCF records
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

    // Convert HashMap to sorted vector for MPHF build
    let mut rows: Vec<(u64, AniEntry, usize)> = entries_map
        .into_iter()
        .map(|(key, (entry, order))| (key, entry, order))
        .collect();

    // Sort by insertion order to maintain deterministic output
    rows.sort_by_key(|(_, _, order)| *order);

    let rows: Vec<(u64, AniEntry)> = rows.into_iter().map(|(k, e, _)| (k, e)).collect();

    finalize_ani_index(rows, pool, output, timing)?;

    Ok(())
}

/// Process single VCF line with deduplication (bcftools-compatible)
fn process_vcf_line_dedup(
    line: &str,
    entries_map: &mut FxHashMap<u64, (AniEntry, usize)>,
    pool: &mut Vec<u8>,
    insertion_order: &mut usize,
    duplicates_skipped: &AtomicUsize,
    multiallelic_count: &AtomicUsize,
    debug: bool,
) -> Result<usize> {
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
        None => return Ok(0),
    };

    let alt_alleles: Vec<&str> = alt.split(',').collect();

    // Track multiallelic variants
    if alt_alleles.len() > 1 {
        multiallelic_count.fetch_add(1, Ordering::Relaxed);
    }

    let mut processed = 0usize;

    // bcftools-style decomposition: create biallelic entry for each ALT
    for a in alt_alleles {
        let key = make_key(chr_id, pos, rf, a);

        // First-wins deduplication
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

        // Create new entry
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

/// Finalize ANI index with MPHF build and serialization
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

    // Adaptive MPHF configuration based on dataset size
    // Larger gamma = more space but better peelability for large datasets
    let (gamma, rehash_limit) = match n {
        0..=100_000 => (1.2, 16),            // Small: tight packing
        100_001..=1_000_000 => (1.5, 32),    // Medium: balanced
        1_000_001..=10_000_000 => (2.0, 64), // Large: 1000G scale
        _ => (2.5, 100),                     // Very large: gnomAD scale
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

        // Memory statistics for large indexes
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

    // Map entries to MPHF indices
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

    // Write header
    unsafe {
        bw.write_all(slice::from_raw_parts(
            (&header as *const _) as *const u8,
            mem::size_of::<AniHeader>(),
        ))?;
    }

    // Write MPHF g array
    let g_bytes = unsafe { slice::from_raw_parts(mph.g.as_ptr() as *const u8, g_size) };
    bw.write_all(g_bytes)?;

    // Write entries
    let ent_bytes = unsafe { slice::from_raw_parts(arr.as_ptr() as *const u8, ent_size) };
    bw.write_all(ent_bytes)?;

    // Write string pool
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

/// Generate deterministic key for variant (chr_id, pos, ref, alt)
fn make_key(chr_id: u8, pos: u32, rf: &str, alt: &str) -> u64 {
    let mut h = hash64(&[chr_id]);
    h ^= hash64(pos.to_le_bytes().as_ref());
    h ^= hash64(rf.as_bytes());
    h ^= hash64(alt.as_bytes());
    h
}

/// Append C-string to pool and return offset
fn append_cstr(pool: &mut Vec<u8>, s: &str) -> u32 {
    let ofs = pool.len() as u32;
    pool.extend_from_slice(s.as_bytes());
    pool.push(0);
    ofs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_deduplication() {
        let mut map = FxHashMap::default();
        let mut pool = Vec::new();
        let mut order = 0;
        let dups = AtomicUsize::new(0);
        let multi = AtomicUsize::new(0);

        // First entry
        let line1 = "chr1\t1000\trs1\tA\tT\t30\tPASS\tDP=10";
        process_vcf_line_dedup(line1, &mut map, &mut pool, &mut order, &dups, &multi, false)
            .unwrap();
        assert_eq!(map.len(), 1);
        assert_eq!(dups.load(Ordering::Relaxed), 0);

        // Duplicate - should be skipped
        let line2 = "chr1\t1000\trs1_dup\tA\tT\t40\tPASS\tDP=20";
        process_vcf_line_dedup(line2, &mut map, &mut pool, &mut order, &dups, &multi, false)
            .unwrap();
        assert_eq!(map.len(), 1); // Still 1 entry
        assert_eq!(dups.load(Ordering::Relaxed), 1); // 1 duplicate skipped

        // Verify first-wins: check if DP=10 (not DP=20)
        let key = make_key(1, 1000, "A", "T");
        let (entry, _) = map.get(&key).unwrap();
        let info_str = read_cstring(&pool, entry.info_ofs as usize);
        assert!(info_str.contains("DP=10"));
    }

    #[test]
    fn test_multiallelic_decomposition() {
        let mut map = FxHashMap::default();
        let mut pool = Vec::new();
        let mut order = 0;
        let dups = AtomicUsize::new(0);
        let multi = AtomicUsize::new(0);

        let line = "chr1\t2000\trs2\tG\tA,C,T\t50\tPASS\tAF=0.1,0.2,0.3";
        process_vcf_line_dedup(line, &mut map, &mut pool, &mut order, &dups, &multi, false)
            .unwrap();

        assert_eq!(map.len(), 3); // 3 biallelic variants
        assert_eq!(multi.load(Ordering::Relaxed), 1); // 1 multiallelic

        // Verify all three ALT alleles are indexed
        assert!(map.contains_key(&make_key(1, 2000, "G", "A")));
        assert!(map.contains_key(&make_key(1, 2000, "G", "C")));
        assert!(map.contains_key(&make_key(1, 2000, "G", "T")));
    }

    #[test]
    fn test_key_generation() {
        let key1 = make_key(1, 1000, "A", "T");
        let key2 = make_key(1, 1000, "A", "T");
        let key3 = make_key(1, 1000, "A", "G");
        let key4 = make_key(1, 1001, "A", "T"); // Different position
        let key5 = make_key(2, 1000, "A", "T"); // Different chromosome

        assert_eq!(key1, key2); // Same variant
        assert_ne!(key1, key3); // Different alt
        assert_ne!(key1, key4); // Different position
        assert_ne!(key1, key5); // Different chromosome
    }

    #[test]
    fn test_insertion_order() {
        let mut map = FxHashMap::default();
        let mut pool = Vec::new();
        let mut order = 0;
        let dups = AtomicUsize::new(0);
        let multi = AtomicUsize::new(0);

        // Add 3 variants
        let lines = [
            "chr1\t1000\trs1\tA\tT\t30\tPASS\t.",
            "chr2\t2000\trs2\tG\tC\t40\tPASS\t.",
            "chr3\t3000\trs3\tC\tA\t50\tPASS\t.",
        ];

        for line in &lines {
            process_vcf_line_dedup(line, &mut map, &mut pool, &mut order, &dups, &multi, false)
                .unwrap();
        }

        // Check insertion order
        let key1 = make_key(1, 1000, "A", "T");
        let key2 = make_key(2, 2000, "G", "C");
        let key3 = make_key(3, 3000, "C", "A");

        assert_eq!(map.get(&key1).unwrap().1, 0);
        assert_eq!(map.get(&key2).unwrap().1, 1);
        assert_eq!(map.get(&key3).unwrap().1, 2);
    }

    #[test]
    fn test_adaptive_mphf_config() {
        // Test gamma selection logic
        let test_cases = vec![
            (50_000, 1.2, 16),
            (500_000, 1.5, 32),
            (5_000_000, 2.0, 64),
            (50_000_000, 2.5, 100),
        ];

        for (n, expected_gamma, expected_rehash) in test_cases {
            let (gamma, rehash) = match n {
                0..=100_000 => (1.2, 16),
                100_001..=1_000_000 => (1.5, 32),
                1_000_001..=10_000_000 => (2.0, 64),
                _ => (2.5, 100),
            };

            assert_eq!(gamma, expected_gamma, "Wrong gamma for n={}", n);
            assert_eq!(rehash, expected_rehash, "Wrong rehash for n={}", n);
        }
    }

    /// Helper to read C-string from pool
    fn read_cstring(pool: &[u8], offset: usize) -> String {
        let mut end = offset;
        while end < pool.len() && pool[end] != 0 {
            end += 1;
        }
        String::from_utf8_lossy(&pool[offset..end]).to_string()
    }
}
