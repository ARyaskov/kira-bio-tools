pub mod entry_processing;
pub mod finalization;
pub mod multiallelic;
pub mod pool;
pub mod tab_processing;
pub mod vcf_processing;

pub use finalization::finalize_ani_index;
pub use pool::StringPool;
pub use tab_processing::process_tab_line_multiallelic;
pub use vcf_processing::process_vcf_line_multiallelic_simd;

use anyhow::Result;
use fxhash::FxHashMap;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::annotate::builder_v2::entry_processing::EntryEntry;
use crate::annotate::structs::ani::{AniEntry, ContigDict};
use crate::annotate::structs::bundle::FieldNumber;
use crate::annotate::structs::tab::TabSchema;
use crate::util::{extract_info_key, extract_info_number};
use crate::vcf::UnifiedVcfReader;

pub fn build_ani_index_auto_v2(input: &Path, output: &Path) -> Result<()> {
    let timing = std::env::var("KIRA_BT_TIMING").is_ok();
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok();
    let start = std::time::Instant::now();

    let mut reader = UnifiedVcfReader::open(input)?;
    let headers = reader.header()?;

    // Pre-seed the contig dict from `##contig=` lines in the source VCF — gives
    // a canonical id-order matching the source genome's declared contigs.
    // Variants on contigs not declared in the header still get an id via
    // first-seen insert inside vcf_processing.
    let mut contigs = ContigDict::from_header_lines(headers.iter().map(String::as_str));
    let mut entries_map: FxHashMap<u64, EntryEntry> = FxHashMap::default();
    let mut pool = StringPool::new();

    let field_meta = extract_info_metadata(&headers);
    let format_meta = extract_format_metadata(&headers);
    let expected_sample_count = extract_sample_count_from_headers(&headers);

    save_vcf_headers_to_pool(&headers, &mut pool)?;
    append_header_end_marker(&mut pool);

    let total_variants = AtomicUsize::new(0);
    let duplicates_skipped = AtomicUsize::new(0);
    let collisions_detected = AtomicUsize::new(0);
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
            &mut contigs,
            &mut entries_map,
            &mut pool,
            &mut insertion_order,
            &duplicates_skipped,
            &collisions_detected,
            &multiallelic_count,
            &field_meta,
            &format_meta,
            expected_sample_count,
            debug,
        )?;

        count += processed;
        if count > 0 && count % report_interval == 0 {
            report_progress(&total_variants, &entries_map, &duplicates_skipped);
        }
    }

    if timing {
        print_final_stats(
            &total_variants,
            &duplicates_skipped,
            &collisions_detected,
            &multiallelic_count,
            &entries_map,
            &contigs,
            &pool,
            &start,
        );
    }

    let mut rows: Vec<(u64, AniEntry, usize)> = entries_map
        .into_iter()
        .map(|(key, ee)| (key, ee.entry, ee.order))
        .collect();

    rows.sort_unstable_by_key(|(_, _, order)| *order);
    let rows: Vec<(u64, AniEntry)> = rows.into_iter().map(|(k, e, _)| (k, e)).collect();

    finalize_ani_index(rows, pool, &contigs, output, timing)?;

    Ok(())
}

pub fn build_ani_index_from_tab(input: &Path, output: &Path, columns: Option<&str>) -> Result<()> {
    let timing = std::env::var("KIRA_BT_TIMING").is_ok();
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok();
    let start = std::time::Instant::now();

    let mut schema = TabSchema::parse(input, columns)?;

    if schema.ref_idx.is_some() && schema.alt_idx.is_some() {
        let alt_idx = schema.alt_idx;

        for col in schema.info_cols.iter_mut() {
            if col.number.is_none() {
                let inferred = infer_number_from_data(input, alt_idx, col.index);
                if let Some(number) = inferred {
                    if debug {
                        eprintln!(
                            "[ani-build-tab] Inferred Number={:?} for column {}",
                            number, col.key
                        );
                    }
                    col.number = Some(number);
                }
            }

            if let Some(number) = col.number {
                schema
                    .field_metadata
                    .entry(col.key.clone())
                    .or_insert(number);
            }
        }
    }

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
        for col in &schema.info_cols {
            eprintln!("    {} -> Number={:?}", col.key, col.number);
        }
    }

    let file = File::open(input)?;
    let reader = BufReader::new(file);

    let mut contigs = ContigDict::default();
    let mut entries_map: FxHashMap<u64, EntryEntry> = FxHashMap::default();
    let mut pool = StringPool::new();

    save_tab_headers_to_pool(&schema, &mut pool)?;
    append_header_end_marker(&mut pool);

    let total_variants = AtomicUsize::new(0);
    let duplicates_skipped = AtomicUsize::new(0);
    let collisions_detected = AtomicUsize::new(0);
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
            &mut contigs,
            &mut entries_map,
            &mut pool,
            &mut insertion_order,
            &duplicates_skipped,
            &collisions_detected,
            &multiallelic_count,
            debug,
        )?;

        count += processed;
        if count > 0 && count % report_interval == 0 {
            report_progress(&total_variants, &entries_map, &duplicates_skipped);
        }
    }

    if timing {
        print_final_stats(
            &total_variants,
            &duplicates_skipped,
            &collisions_detected,
            &multiallelic_count,
            &entries_map,
            &contigs,
            &pool,
            &start,
        );
    }

    let mut rows: Vec<(u64, AniEntry, usize)> = entries_map
        .into_iter()
        .map(|(key, ee)| (key, ee.entry, ee.order))
        .collect();

    rows.sort_unstable_by_key(|(_, _, order)| *order);
    let rows: Vec<(u64, AniEntry)> = rows.into_iter().map(|(k, e, _)| (k, e)).collect();

    finalize_ani_index(rows, pool, &contigs, output, timing)?;

    Ok(())
}

fn save_tab_headers_to_pool(schema: &TabSchema, pool: &mut StringPool) -> Result<()> {
    if schema.ref_idx.is_none() && schema.alt_idx.is_none() && schema.to_idx.is_some() {
        pool.append_cstr("##KIRA_BT_ANI_INTERVALS");
    }

    for col in &schema.info_cols {
        let number_str = match col.number {
            Some(FieldNumber::Zero) => "0",
            Some(FieldNumber::A) => "A",
            Some(FieldNumber::R) => "R",
            Some(FieldNumber::G) => "G",
            Some(FieldNumber::One) => "1",
            Some(FieldNumber::Many) => ".",
            None => ".",
        };

        // bcftools never guesses Type from the column name; the .tab format
        // doesn't carry per-column types. We always declare `Type=String`
        // here — the canonical place to override is a sibling `.hdr` file
        // (parsed by `TabSchema::parse_header_file`) where the user provides
        // explicit `##INFO=<…,Type=…>` lines that take precedence at lookup
        // time (`load_field_metadata` in cpu_v2::field_metadata).
        let header = format!(
            "##INFO=<ID={key},Number={number_str},Type=String,Description=\"Annotation field\">",
            key = col.key
        );
        pool.append_cstr(&header);
    }

    Ok(())
}

fn report_progress(
    total_variants: &AtomicUsize,
    entries_map: &FxHashMap<u64, EntryEntry>,
    duplicates_skipped: &AtomicUsize,
) {
    let total = total_variants.load(Ordering::Relaxed);
    let unique = entries_map.len();
    let dups = duplicates_skipped.load(Ordering::Relaxed);
    let dup_rate = (dups as f64 / total as f64) * 100.0;

    eprintln!(
        "[ani-build] Progress: {total} variants → {unique} unique entries ({dups} dups, {dup_rate:.1}%)"
    );
}

fn print_final_stats(
    total_variants: &AtomicUsize,
    duplicates_skipped: &AtomicUsize,
    collisions_detected: &AtomicUsize,
    multiallelic_count: &AtomicUsize,
    entries_map: &FxHashMap<u64, EntryEntry>,
    contigs: &ContigDict,
    pool: &StringPool,
    start: &std::time::Instant,
) {
    let total = total_variants.load(Ordering::Relaxed);
    let dups = duplicates_skipped.load(Ordering::Relaxed);
    let collisions = collisions_detected.load(Ordering::Relaxed);
    let multi = multiallelic_count.load(Ordering::Relaxed);

    eprintln!("[ani-build] Parse complete:");
    eprintln!("  Total variants:     {total}");
    eprintln!("  Unique entries:     {}", entries_map.len());
    eprintln!("  Contigs:            {}", contigs.len());
    eprintln!(
        "  Duplicates skipped: {dups} ({rate:.1}%)",
        rate = (dups as f64 / total.max(1) as f64) * 100.0
    );
    if collisions > 0 {
        eprintln!(
            "  HASH COLLISIONS:    {collisions}  ← review previous warnings; \
             previous entries were dropped"
        );
    }
    eprintln!(
        "  Multiallelic:       {multi} ({rate:.1}%)",
        rate = (multi as f64 / total.max(1) as f64) * 100.0
    );
    let hits = pool.intern_hits();
    let saved_mb = pool.intern_bytes_saved() as f64 / (1024.0 * 1024.0);
    eprintln!(
        "  Intern hits:        {hits} ({saved_mb:.1} MB saved before deflate)"
    );
    eprintln!("  Parse time: {:.3}s", start.elapsed().as_secs_f64());
}

fn infer_number_from_data(
    file_path: &Path,
    alt_idx: Option<usize>,
    col_idx: usize,
) -> Option<FieldNumber> {
    let file = File::open(file_path).ok()?;
    let reader = BufReader::new(file);

    let mut value_counts = Vec::new();
    let mut alt_counts = Vec::new();

    for line in reader.lines().take(100) {
        let line = line.ok()?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        let val = parts.get(col_idx)?;
        if val.trim().is_empty() || val.trim() == "." {
            continue;
        }

        let num_values = val.split(',').count();
        value_counts.push(num_values);

        if let Some(alt_idx) = alt_idx {
            if let Some(alt) = parts.get(alt_idx) {
                let num_alts = alt.split(',').count();
                alt_counts.push(num_alts);
            }
        }
    }

    if value_counts.is_empty() {
        return None;
    }

    if alt_counts.is_empty() {
        return Some(if value_counts[0] == 1 {
            FieldNumber::One
        } else {
            FieldNumber::Many
        });
    }

    let mut count_a = 0;
    let mut count_r = 0;

    for (val_count, alt_count) in value_counts.iter().zip(alt_counts.iter()) {
        if *val_count == *alt_count {
            count_a += 1;
        } else if *val_count == *alt_count + 1 {
            count_r += 1;
        }
    }

    if count_a > count_r && count_a > value_counts.len() / 2 {
        Some(FieldNumber::A)
    } else if count_r > count_a && count_r > value_counts.len() / 2 {
        Some(FieldNumber::R)
    } else if value_counts[0] == 1 {
        Some(FieldNumber::One)
    } else {
        Some(FieldNumber::Many)
    }
}

fn save_vcf_headers_to_pool(headers: &[String], pool: &mut StringPool) -> Result<()> {
    for h in headers {
        if h.starts_with("##INFO=")
            || h.starts_with("##FORMAT=")
            || h.starts_with("##FILTER=")
            || h.starts_with("#CHROM")
        {
            pool.append_cstr(h);
        }
    }

    Ok(())
}

fn append_header_end_marker(pool: &mut StringPool) {
    use crate::annotate::structs::ani::ANI_HEADER_END;
    pool.append_cstr(ANI_HEADER_END);
}

fn extract_info_metadata(headers: &[String]) -> HashMap<String, FieldNumber> {
    let mut meta = HashMap::new();
    for h in headers {
        if !h.starts_with("##INFO=") {
            continue;
        }
        if let Some(key) = extract_info_key(h) {
            if let Some(number) = extract_info_number(h) {
                meta.insert(key, number);
            }
        }
    }
    meta
}

fn extract_format_metadata(headers: &[String]) -> HashMap<String, FieldNumber> {
    let mut meta = HashMap::new();
    for h in headers {
        if !h.starts_with("##FORMAT=") {
            continue;
        }
        if let Some(key) = extract_info_key(h) {
            if let Some(number) = extract_info_number(h) {
                meta.insert(key, number);
            }
        }
    }
    meta
}

fn extract_sample_count_from_headers(headers: &[String]) -> usize {
    for h in headers {
        if h.starts_with("#CHROM") {
            let parts: Vec<&str> = h.split('\t').collect();
            if parts.len() > 9 {
                return parts.len() - 9;
            }
            break;
        }
    }
    0
}

#[cfg(test)]
#[path = "../../../tests/unit/annotate_builder_v2.rs"]
mod tests;
