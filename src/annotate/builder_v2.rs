use anyhow::Result;
use fxhash::FxHashMap;
use kira_kv_engine::Mphf;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::mem;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::annotate::structs::{AniEntry, AniHeader, ANI_MAGIC, ANI_VERSION};
use crate::chr_name_to_id;
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

    let mut line_no = 0usize;
    for line in reader.lines() {
        let line = line?;
        line_no += 1;

        if line.trim().is_empty() || line.starts_with('#') {
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

#[derive(Debug, Clone)]
struct TabColumn {
    index: usize,
    key: String,
    is_append: bool,
}

#[derive(Debug, Clone)]
struct TabSchema {
    chrom_idx: usize,
    pos_idx: usize,
    ref_idx: Option<usize>,
    alt_idx: Option<usize>,
    id_idx: Option<usize>,
    qual_idx: Option<usize>,
    filter_idx: Option<usize>,
    info_start: Option<usize>,
    info_cols: Vec<TabColumn>,
}

impl TabSchema {
    fn parse(path: &Path, columns: Option<&str>) -> Result<Self> {
        if let Some(cols) = columns {
            return Self::from_column_spec(cols);
        }

        let file = File::open(path)?;
        let mut reader = BufReader::new(file);
        let mut first_line = String::new();
        reader.read_line(&mut first_line)?;

        if first_line.starts_with('#') {
            let header = first_line.trim_start_matches('#').trim();
            Self::from_header(header)
        } else {
            let ncols = first_line.split('\t').count();
            if ncols >= 9 {
                Self::from_column_spec("CHROM,POS,REF,ALT,ID,QUAL,FILTER,INFO")
            } else if ncols >= 5 {
                Self::from_column_spec("CHROM,POS,REF,ALT,ID")
            } else if ncols >= 4 {
                Self::from_column_spec("CHROM,POS,REF,ALT")
            } else {
                anyhow::bail!("Cannot detect TAB schema from {} columns", ncols)
            }
        }
    }

    fn from_header(header: &str) -> Result<Self> {
        let parts: Vec<&str> = header.split('\t').map(|s| s.trim()).collect();
        let spec = parts.join(",");
        Self::from_column_spec(&spec)
    }

    fn from_column_spec(spec: &str) -> Result<Self> {
        let parts: Vec<&str> = spec.split(',').map(|s| s.trim()).collect();

        let mut chrom_idx = None;
        let mut pos_idx = None;
        let mut ref_idx = None;
        let mut alt_idx = None;
        let mut id_idx = None;
        let mut qual_idx = None;
        let mut filter_idx = None;
        let mut info_start = None;
        let mut info_cols = Vec::new();

        for (i, part) in parts.iter().enumerate() {
            let (is_append, clean_part) = if part.starts_with('+') {
                (true, &part[1..])
            } else if part.starts_with('-') {
                (false, &part[1..])
            } else {
                (false, *part)
            };

            match clean_part {
                "CHROM" => chrom_idx = Some(i),
                "POS" => pos_idx = Some(i),
                "REF" => ref_idx = Some(i),
                "ALT" => alt_idx = Some(i),
                "ID" => id_idx = Some(i),
                "QUAL" => qual_idx = Some(i),
                "FILTER" => filter_idx = Some(i),
                "INFO" => info_start = Some(i),
                _ if clean_part.starts_with("INFO/") => {
                    let key = clean_part.strip_prefix("INFO/").unwrap();
                    info_cols.push(TabColumn {
                        index: i,
                        key: key.to_string(),
                        is_append,
                    });
                }
                _ if clean_part.starts_with("FMT/") || clean_part.starts_with("FORMAT/") => {
                    let key = if let Some(k) = clean_part.strip_prefix("FMT/") {
                        k
                    } else {
                        clean_part.strip_prefix("FORMAT/").unwrap()
                    };
                    info_cols.push(TabColumn {
                        index: i,
                        key: key.to_string(),
                        is_append,
                    });
                }
                _ => {
                    info_cols.push(TabColumn {
                        index: i,
                        key: clean_part.to_string(),
                        is_append,
                    });
                }
            }
        }

        let chrom_idx = chrom_idx.ok_or_else(|| anyhow::anyhow!("CHROM column required"))?;
        let pos_idx = pos_idx.ok_or_else(|| anyhow::anyhow!("POS column required"))?;

        Ok(Self {
            chrom_idx,
            pos_idx,
            ref_idx,
            alt_idx,
            id_idx,
            qual_idx,
            filter_idx,
            info_start,
            info_cols,
        })
    }
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
                pool.extend_from_slice(info_str.as_bytes());
            }
        }
    } else {
        let mut info_parts = Vec::new();

        for col in &schema.info_cols {
            if let Some(val) = parts.get(col.index) {
                let val = val.trim();
                if val.is_empty() || val == "." {
                    continue;
                }

                let encoded_val = url_encode_info_value(val);

                let field_str = if col.is_append {
                    format!("+{}={}", col.key, encoded_val)
                } else {
                    format!("{}={}", col.key, encoded_val)
                };
                info_parts.push(field_str);
            }
        }

        if !info_parts.is_empty() {
            pool.extend_from_slice(info_parts.join(";").as_bytes());
        }
    }

    pool.push(0);
    let info_ofs = info_start_ofs as u32;
    let info_len = (pool.len() - info_start_ofs - 1) as u32;

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

    Ok(1)
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
    let pos = parts[1].parse::<u32>().unwrap_or(0);
    let id = if parts[2] == "." { "." } else { parts[2] };
    let rf = parts[3];
    let alt = parts[4];
    let qual = if parts[5] == "." { "." } else { parts[5] };
    let filter = if parts[6] == "." { "." } else { parts[6] };
    let info = parts[7];

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
                "[ani-build-vcf] Overwriting duplicate: {}:{} {}",
                chr, pos, rf
            );
        }
    }

    let ref_ofs = append_cstr(pool, rf);
    let alt_ofs = append_cstr(pool, alt);
    let id_ofs = append_cstr(pool, id);
    let qual_ofs = append_cstr(pool, qual);
    let filter_ofs = append_cstr(pool, filter);
    let info_ofs = append_cstr(pool, info);

    let entry = AniEntry {
        chr_id,
        pos,
        ref_ofs,
        alt_ofs,
        id_ofs,
        qual_ofs,
        filter_ofs,
        info_ofs,
        info_len: 0,
    };

    entries_map.insert(key, (entry, *insertion_order));
    *insertion_order += 1;

    Ok(1)
}

fn url_encode_info_value(val: &str) -> String {
    let mut result = String::with_capacity(val.len() * 2);

    for c in val.chars() {
        match c {
            '=' => result.push_str("%3D"),
            ';' => result.push_str("%3B"),
            ',' => result.push_str("%2C"),
            '%' => result.push_str("%25"),
            ' ' => result.push_str("%20"),
            '\t' => result.push_str("%09"),
            '\n' => result.push_str("%0A"),
            '\r' => result.push_str("%0D"),
            _ => result.push(c),
        }
    }

    result
}

fn make_position_key(chr_id: u8, pos: u32, rf: &str) -> u64 {
    use fxhash::hash64;

    let mut h = hash64(&[chr_id]);
    h ^= hash64(pos.to_le_bytes().as_ref());
    h ^= hash64(rf.as_bytes());
    h
}

fn append_cstr(pool: &mut Vec<u8>, s: &str) -> u32 {
    let ofs = pool.len() as u32;
    pool.extend_from_slice(s.as_bytes());
    pool.push(0);
    ofs
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

    let entries: Vec<AniEntry> = rows.into_iter().map(|(_, e)| e).collect();

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
