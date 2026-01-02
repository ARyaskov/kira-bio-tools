pub mod entry_processing;
pub mod finalization;
pub mod multiallelic;
pub mod tab_processing;
pub mod vcf_processing;

pub use finalization::finalize_ani_index;
pub use tab_processing::process_tab_line_multiallelic;
pub use vcf_processing::process_vcf_line_multiallelic_simd;

use anyhow::Result;
use fxhash::FxHashMap;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::annotate::structs::ani::AniEntry;
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

    let mut entries_map: FxHashMap<u64, (AniEntry, usize)> = FxHashMap::default();
    let mut pool = Vec::<u8>::new();

    let field_meta = extract_info_metadata(&headers);
    let expected_sample_count = extract_sample_count_from_headers(&headers);

    save_vcf_headers_to_pool(&headers, &mut pool)?;
    append_header_end_marker(&mut pool);

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
            &field_meta,
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
            &multiallelic_count,
            &entries_map,
            &start,
        );
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

    let mut schema = TabSchema::parse(input, columns)?;

    if schema.ref_idx.is_some() && schema.alt_idx.is_some() {
        let chrom_idx = schema.chrom_idx;
        let pos_idx = schema.pos_idx;
        let ref_idx = schema.ref_idx;
        let alt_idx = schema.alt_idx;

        for col in schema.info_cols.iter_mut() {
            if col.number.is_none() {
                let inferred =
                    infer_number_from_data(input, chrom_idx, pos_idx, ref_idx, alt_idx, col.index);
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

    let mut entries_map: FxHashMap<u64, (AniEntry, usize)> = FxHashMap::default();
    let mut pool = Vec::<u8>::new();

    save_tab_headers_to_pool(&schema, &mut pool)?;
    append_header_end_marker(&mut pool);

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
            report_progress(&total_variants, &entries_map, &duplicates_skipped);
        }
    }

    if timing {
        print_final_stats(
            &total_variants,
            &duplicates_skipped,
            &multiallelic_count,
            &entries_map,
            &start,
        );
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

fn save_tab_headers_to_pool(schema: &TabSchema, pool: &mut Vec<u8>) -> Result<()> {
    use crate::util::append_cstr;

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

        let type_str = if col.key.starts_with('I') || col.key.ends_with("INT") {
            "Integer"
        } else if col.key.starts_with('F') || col.key.ends_with("FLT") || col.key.ends_with("FLOAT")
        {
            "Float"
        } else if col.key.starts_with('S')
            || col.key.ends_with("STR")
            || col.key.ends_with("STRING")
        {
            "String"
        } else {
            "String"
        };

        let header = format!(
            "##INFO=<ID={},Number={},Type={},Description=\"Annotation field\">",
            col.key, number_str, type_str
        );
        append_cstr(pool, &header);
    }

    Ok(())
}

fn report_progress(
    total_variants: &AtomicUsize,
    entries_map: &FxHashMap<u64, (AniEntry, usize)>,
    duplicates_skipped: &AtomicUsize,
) {
    let total = total_variants.load(Ordering::Relaxed);
    let unique = entries_map.len();
    let dups = duplicates_skipped.load(Ordering::Relaxed);
    let dup_rate = (dups as f64 / total as f64) * 100.0;

    eprintln!(
        "[ani-build] Progress: {} variants → {} unique entries ({} dups, {:.1}%)",
        total, unique, dups, dup_rate
    );
}

fn print_final_stats(
    total_variants: &AtomicUsize,
    duplicates_skipped: &AtomicUsize,
    multiallelic_count: &AtomicUsize,
    entries_map: &FxHashMap<u64, (AniEntry, usize)>,
    start: &std::time::Instant,
) {
    let total = total_variants.load(Ordering::Relaxed);
    let dups = duplicates_skipped.load(Ordering::Relaxed);
    let multi = multiallelic_count.load(Ordering::Relaxed);

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

fn infer_number_from_data(
    file_path: &Path,
    chrom_idx: usize,
    pos_idx: usize,
    ref_idx: Option<usize>,
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

fn save_vcf_headers_to_pool(headers: &[String], pool: &mut Vec<u8>) -> Result<()> {
    use crate::util::append_cstr;

    for h in headers {
        if h.starts_with("##INFO=")
            || h.starts_with("##FORMAT=")
            || h.starts_with("##FILTER=")
            || h.starts_with("#CHROM")
        {
            append_cstr(pool, h);
        }
    }

    Ok(())
}

fn append_header_end_marker(pool: &mut Vec<u8>) {
    use crate::annotate::structs::ani::ANI_HEADER_END;
    use crate::util::append_cstr;

    append_cstr(pool, ANI_HEADER_END);
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
mod tests {
    use super::build_ani_index_from_tab;
    use crate::annotate::structs::ani::AniIndex;
    use std::fs;

    #[test]
    fn test_tab_infers_numbers_for_split_info() {
        let dir = std::env::temp_dir();
        let tab_path = dir.join("kira_tab_infer_test.tab");
        let ani_path = dir.join("kira_tab_infer_test.ani");

        let tab = "1\t1\tC\tA,T\t0,1.1\t1.1,0,2.2\n";
        fs::write(&tab_path, tab).unwrap();

        build_ani_index_from_tab(&tab_path, &ani_path, Some("CHROM,POS,REF,ALT,FA,FR")).unwrap();

        let ani = AniIndex::open(&ani_path).unwrap();
        let bundle = ani.lookup_exact("1", 1, "C", "T").unwrap();

        let fa = bundle
            .info
            .iter()
            .find(|f| f.key == "FA")
            .map(|f| f.values.clone())
            .unwrap();
        let fr = bundle
            .info
            .iter()
            .find(|f| f.key == "FR")
            .map(|f| f.values.clone())
            .unwrap();

        assert_eq!(fa, vec!["1.1"]);
        assert_eq!(fr, vec!["1.1", "2.2"]);

        let _ = fs::remove_file(&tab_path);
        let _ = fs::remove_file(&ani_path);
    }
}
