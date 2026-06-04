use indexmap::IndexMap;
use std::collections::HashMap;

use crate::annotate::cpu_v2::field_metadata::is_missing_value;
use crate::annotate::cpu_v2::vcf_parsing::{
    ParsedVcfRecord, parse_vcf_record_simd, patch_samples_from_line,
};
use crate::annotate::structs::ani::AniIndex;
use crate::annotate::structs::annotate_mode::AnnotateMode;
use crate::annotate::structs::bundle::{AnnotationBundle, FieldNumber};
use crate::util::{chr_name_to_id, fast_hash64};
use std::sync::atomic::{AtomicU64, Ordering};

pub struct BundleTimingAccum {
    read_ns: AtomicU64,
    info_ns: AtomicU64,
    optional_ns: AtomicU64,
    samples_ns: AtomicU64,
}

impl BundleTimingAccum {
    pub fn new() -> Self {
        Self {
            read_ns: AtomicU64::new(0),
            info_ns: AtomicU64::new(0),
            optional_ns: AtomicU64::new(0),
            samples_ns: AtomicU64::new(0),
        }
    }

    pub fn add(&self, t: &crate::annotate::structs::ani::lookup::BundleTiming) {
        self.read_ns
            .fetch_add((t.read_s * 1e9) as u64, Ordering::Relaxed);
        self.info_ns
            .fetch_add((t.info_s * 1e9) as u64, Ordering::Relaxed);
        self.optional_ns
            .fetch_add((t.optional_s * 1e9) as u64, Ordering::Relaxed);
        self.samples_ns
            .fetch_add((t.samples_s * 1e9) as u64, Ordering::Relaxed);
    }

    pub fn snapshot_seconds(&self) -> (f64, f64, f64, f64) {
        (
            self.read_ns.load(Ordering::Relaxed) as f64 / 1e9,
            self.info_ns.load(Ordering::Relaxed) as f64 / 1e9,
            self.optional_ns.load(Ordering::Relaxed) as f64 / 1e9,
            self.samples_ns.load(Ordering::Relaxed) as f64 / 1e9,
        )
    }
}

fn split_mapped_ref(raw: &str) -> (&str, &str) {
    if let Some((src, dst)) = raw.split_once("=>") {
        (src, dst)
    } else {
        (raw, raw)
    }
}

fn is_fixed_column(raw: &str) -> bool {
    raw.eq_ignore_ascii_case("ID")
        || raw.eq_ignore_ascii_case("QUAL")
        || raw.eq_ignore_ascii_case("FILTER")
        || raw.eq_ignore_ascii_case("ALT")
        || raw.eq_ignore_ascii_case("FMT")
        || raw.eq_ignore_ascii_case("FORMAT")
}

fn is_format_ref(raw: &str) -> bool {
    let upper = raw.to_ascii_uppercase();
    upper == "FMT" || upper == "FORMAT" || upper.starts_with("FMT/") || upper.starts_with("FORMAT/")
}

pub fn column_targets_format(raw: &str) -> bool {
    let (_, dst) = split_mapped_ref(raw);
    is_format_ref(dst)
}

fn column_targets_info(raw: &str) -> bool {
    let (_, dst) = split_mapped_ref(raw);
    !is_fixed_column(dst) && !is_format_ref(dst)
}

fn columns_are_info_only(column_modes: &[(String, AnnotateMode)], format_overwrite_all: bool) -> bool {
    !format_overwrite_all && column_modes.iter().all(|(k, _)| column_targets_info(k))
}

fn can_use_entry_info_fast_path(
    info_overwrite_all: bool,
    format_overwrite_all: bool,
    column_modes: &[(String, AnnotateMode)],
) -> bool {
    info_overwrite_all
        && !format_overwrite_all
        && !column_modes.is_empty()
        && column_modes.iter().all(|(k, mode)| {
            let (src, dst) = split_mapped_ref(k);
            *mode == AnnotateMode::default_mode()
                && column_targets_info(k)
                && info_key_from_ref(src) == info_key_from_ref(dst)
        })
}

/// Build the set of INFO tag names the merge will read from each bundle.
/// Returns `None` when the full bundle is required (overwrite-all modes or
/// non-INFO targets).
fn compute_info_filter<'a>(
    info_overwrite_all: bool,
    format_overwrite_all: bool,
    column_modes: &'a [(String, AnnotateMode)],
) -> Option<std::collections::HashSet<&'a str>> {
    if info_overwrite_all || format_overwrite_all {
        return None;
    }
    let mut set: std::collections::HashSet<&'a str> =
        std::collections::HashSet::with_capacity(column_modes.len());
    for (k, _mode) in column_modes {
        let (src_raw, dst_raw) = split_mapped_ref(k);
        if !column_targets_info(k) {
            return None;
        }
        let dst_upper = dst_raw.to_ascii_uppercase();
        if dst_upper == "ID"
            || dst_upper == "QUAL"
            || dst_upper == "FILTER"
            || dst_upper == "ALT"
            || dst_upper == "FMT"
            || dst_upper == "FORMAT"
            || dst_upper.starts_with("FMT/")
            || dst_upper.starts_with("FORMAT/")
        {
            return None;
        }
        set.insert(info_key_from_ref(src_raw));
    }
    Some(set)
}

fn info_key_from_ref(raw: &str) -> &str {
    raw.strip_prefix("INFO/").unwrap_or(raw)
}

fn format_key_from_ref(raw: &str) -> &str {
    raw.strip_prefix("FMT/")
        .or_else(|| raw.strip_prefix("FORMAT/"))
        .unwrap_or(raw)
}

/// Returns the byte offset of `child` within `parent` if `child` is a
/// proper substring slice of `parent`; `None` otherwise.
#[inline]
fn offset_within(parent: &str, child: &str) -> Option<usize> {
    let p = parent.as_ptr() as usize;
    let c = child.as_ptr() as usize;
    let p_end = p.checked_add(parent.len())?;
    let c_end = c.checked_add(child.len())?;
    if c >= p && c_end <= p_end {
        Some(c - p)
    } else {
        None
    }
}

fn tail_after_info<'a>(line: &'a str, parsed: &ParsedVcfRecord<'a>) -> &'a str {
    let Some(start) = offset_within(line, parsed.info) else {
        return "";
    };
    let end = start + parsed.info.len();
    line.get(end..).unwrap_or("")
}

fn build_info_replaced_line(line: &str, parsed: &ParsedVcfRecord, new_info: &str) -> String {
    let Some(start) = offset_within(line, parsed.info) else {
        return line.to_string();
    };
    let end = start + parsed.info.len();
    let Some(prefix) = line.get(..start) else {
        return line.to_string();
    };
    let Some(tail) = line.get(end..) else {
        return line.to_string();
    };
    let info = if new_info.is_empty() { "." } else { new_info };
    let mut out = String::with_capacity(prefix.len() + info.len() + tail.len());
    out.push_str(prefix);
    out.push_str(info);
    out.push_str(tail);
    out
}

fn serialize_info_map(info_map: &IndexMap<String, String>) -> String {
    if info_map.is_empty() {
        return ".".to_string();
    }
    let mut len = 0usize;
    for (k, v) in info_map.iter() {
        len += k.len();
        if !v.is_empty() {
            len += 1 + v.len();
        }
    }
    len += info_map.len().saturating_sub(1);
    let mut out = String::with_capacity(len);
    for (i, (k, v)) in info_map.iter().enumerate() {
        if i > 0 {
            out.push(';');
        }
        if v.is_empty() {
            out.push_str(k);
        } else {
            out.push_str(k);
            out.push('=');
            out.push_str(v);
        }
    }
    if out.is_empty() { ".".to_string() } else { out }
}

fn exact_mph_entry_idx(
    ani: &AniIndex,
    chr_id: u32,
    pos: u32,
    ref_allele: &str,
    alt: &str,
    mph_hits: &[Option<usize>],
) -> Option<usize> {
    let idx = mph_hits.first().copied().flatten()?;
    let e = ani.entries.get(idx)?;
    if e.chr_id != chr_id || e.pos != pos {
        return None;
    }
    let rf_str = ani.read_cstring(e.ref_ofs as usize);
    if rf_str.as_ref() != ref_allele {
        return None;
    }
    let alt_str = ani.read_cstring(e.alt_ofs as usize);
    if alt_str.as_ref() != alt {
        return None;
    }
    Some(idx)
}

fn read_scalar_from_bundle(bundle: &AnnotationBundle, src_ref: &str) -> Option<String> {
    let src_upper = src_ref.to_ascii_uppercase();
    if src_upper == "ID" {
        return bundle.id.clone();
    }
    if src_upper == "QUAL" {
        return bundle.qual.clone();
    }
    if src_upper == "FILTER" {
        return bundle.filter.clone();
    }
    if src_upper == "ALT" {
        return Some(bundle.alt.clone());
    }
    if src_upper.starts_with("INFO/") || !src_ref.contains('/') {
        let key = info_key_from_ref(src_ref);
        return bundle
            .info
            .iter()
            .find(|f| f.key == key)
            .map(|f| join_values(&f.values));
    }
    None
}

fn apply_format_field_mappings(
    format: &mut Option<String>,
    samples: &mut Vec<String>,
    bundle: &AnnotationBundle,
    column_modes: &[(String, AnnotateMode)],
    sample_map: &[Option<usize>],
) {
    let Some(db_format) = bundle.format_str.as_ref() else {
        return;
    };
    if samples.is_empty() {
        return;
    }

    let db_keys: Vec<&str> = db_format.split(':').collect();
    let mut out_keys: Vec<String> = format
        .as_ref()
        .map(|s| {
            if s.is_empty() {
                Vec::new()
            } else {
                s.split(':').map(|v| v.to_string()).collect()
            }
        })
        .unwrap_or_default();
    let original_keys = out_keys.clone();
    let mut out_sample_values: Vec<Vec<String>> = samples
        .iter()
        .map(|s| {
            if s.is_empty() {
                Vec::new()
            } else {
                s.split(':').map(|v| v.to_string()).collect()
            }
        })
        .collect();

    for (raw_key, mode) in column_modes {
        let (src_ref, dst_ref) = split_mapped_ref(raw_key);
        if !is_format_ref(dst_ref) {
            continue;
        }
        if dst_ref.eq_ignore_ascii_case("FMT") || dst_ref.eq_ignore_ascii_case("FORMAT") {
            continue;
        }
        if !is_format_ref(src_ref) {
            continue;
        }

        let src_key = format_key_from_ref(src_ref);
        let dst_key = format_key_from_ref(dst_ref);
        let src_db_idx = db_keys.iter().position(|k| *k == src_key);
        let Some(src_db_idx) = src_db_idx else {
            continue;
        };
        let dst_out_idx = match out_keys.iter().position(|k| k == dst_key) {
            Some(i) => i,
            None => {
                out_keys.push(dst_key.to_string());
                for vals in &mut out_sample_values {
                    vals.push(".".to_string());
                }
                out_keys.len() - 1
            }
        };

        for (sample_idx, vals) in out_sample_values.iter_mut().enumerate() {
            if dst_out_idx >= vals.len() {
                vals.resize(dst_out_idx + 1, ".".to_string());
            }
            let existing = vals.get(dst_out_idx).map(|s| s.as_str()).unwrap_or(".");
            let db_idx = sample_map.get(sample_idx).and_then(|v| *v);
            if db_idx.is_none() {
                continue;
            }
            let db_val = db_idx
                .and_then(|idx| bundle.format_samples.get(idx))
                .and_then(|s| s.split(':').nth(src_db_idx))
                .unwrap_or(".");
            if let Some(v) = merge_scalar_field(existing, db_val, *mode, ",") {
                vals[dst_out_idx] = v;
            }
        }
    }

    let keep: Vec<usize> = out_keys
        .iter()
        .enumerate()
        .filter_map(|(idx, key)| {
            if original_keys.iter().any(|k| k == key)
                || out_sample_values.iter().any(|vals| {
                    vals.get(idx)
                        .map(|v| !is_missing_value(v))
                        .unwrap_or(false)
                })
            {
                Some(idx)
            } else {
                None
            }
        })
        .collect();
    if keep.len() != out_keys.len() {
        out_keys = keep.iter().map(|&idx| out_keys[idx].clone()).collect();
        for vals in &mut out_sample_values {
            let next: Vec<String> = keep
                .iter()
                .map(|&idx| vals.get(idx).cloned().unwrap_or_else(|| ".".to_string()))
                .collect();
            *vals = next;
        }
    }

    if out_keys.is_empty() {
        *format = None;
    } else {
        *format = Some(join_keys(&out_keys));
    }
    *samples = out_sample_values
        .into_iter()
        .map(|vals| normalize_sample_values(&vals))
        .collect();
}

pub fn annotate_line(
    line: &str,
    ani: &AniIndex,
    field_meta: &HashMap<String, FieldNumber>,
    column_modes: &[(String, AnnotateMode)],
    sample_map: &[Option<usize>],
    info_overwrite_all: bool,
    format_overwrite_all: bool,
) -> String {
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok();

    let want_format =
        format_overwrite_all || column_modes.iter().any(|(k, _)| column_targets_format(k));
    let info_only = columns_are_info_only(column_modes, format_overwrite_all);

    let Some(mut parsed) = parse_vcf_record_simd(line, want_format) else {
        return line.to_string();
    };
    if want_format {
        patch_samples_from_line(&mut parsed, line);
    }
    let raw_tail = if !want_format {
        Some(tail_after_info(line, &parsed))
    } else {
        None
    };

    let raw_samples = if want_format && sample_map.iter().any(|v| v.is_none()) {
        let mut fields = line.split('\t');
        for _ in 0..9 {
            let _ = fields.next();
        }
        Some(fields.collect::<Vec<&str>>())
    } else {
        None
    };

    let chrom = &parsed.chrom;
    let pos = parsed.pos;
    let ref_allele = &parsed.ref_allele;
    let alt_alleles: Vec<&str> = parsed.alt.split(',').collect();

    if debug {
        eprintln!(
            "[ANNOTATE] {}:{} {} -> {:?}",
            chrom, pos, ref_allele, alt_alleles
        );
    }

    let Some(chr_id) = ani
        .contig_id(chrom)
        .or_else(|| chr_name_to_id(chrom).map(u32::from))
    else {
        return line.to_string();
    };
    let ref_hash = fast_hash64(ref_allele.as_bytes());

    let need_info = info_overwrite_all || column_modes.iter().any(|(k, _)| column_targets_info(k));
    let need_format = want_format;

    let mut bundles = Vec::with_capacity(alt_alleles.len());
    for (vcf_idx, alt) in alt_alleles.iter().enumerate() {
        if let Some(bundle) = ani.lookup_exact_by_chr_id_pos_index_opts(
            chr_id,
            pos,
            ref_allele,
            ref_hash,
            alt,
            field_meta,
            need_info,
            need_format,
        ) {
            if debug {
                eprintln!(
                    "[ANNOTATE] Found bundle for ALT {}: INFO fields={}",
                    alt,
                    bundle.info.len()
                );
            }
            bundles.push((vcf_idx, bundle));
        }
    }

    if bundles.is_empty() {
        return line.to_string();
    }

    if info_only {
        return annotate_info_only_with_bundles(
            line,
            &parsed,
            &alt_alleles,
            &bundles,
            field_meta,
            column_modes,
            info_overwrite_all,
        )
        .unwrap_or_else(|| line.to_string());
    }

    annotate_record_with_alts(
        &parsed,
        &alt_alleles,
        &bundles,
        field_meta,
        column_modes,
        sample_map,
        raw_tail,
        raw_samples.as_deref(),
        info_overwrite_all,
        format_overwrite_all,
        debug,
    )
}

/// Per-line key tuple emitted by the worker thread's batch-pre-lookup pass.
pub struct LineKeys<'a> {
    pub parsed: ParsedVcfRecord<'a>,
    pub chr_id: Option<u32>,
    pub alt_alleles: smallvec::SmallVec<[&'a str; 2]>,
    pub keys: smallvec::SmallVec<[u64; 2]>,
}

/// Pure parse + key-extraction pass. Used by the batched worker to compute
/// all MPH keys for a batch up-front so a single SIMD lookup can be used.
pub fn extract_line_keys<'a>(
    line: &'a str,
    ani: &AniIndex,
    want_format: bool,
) -> Option<LineKeys<'a>> {
    use crate::annotate::structs::ani::make_variant_key;

    // Chrom pre-filter (single memchr, no full parse).
    let bytes = line.as_bytes();
    let tab1 = memchr::memchr(b'\t', bytes)?;
    // SAFETY: input was &str, prefix slice of UTF-8 bytes is valid UTF-8.
    let chrom = unsafe { std::str::from_utf8_unchecked(&bytes[..tab1]) };
    let chr_id_early = ani
        .contig_id(chrom)
        .or_else(|| chr_name_to_id(chrom).map(u32::from));
    match chr_id_early {
        Some(id) if !ani.chrom_has_entries(id) => {
            return None;
        }
        None => {
            return None;
        }
        Some(_) => {}
    }

    let mut parsed = parse_vcf_record_simd(line, want_format)?;
    if want_format {
        patch_samples_from_line(&mut parsed, line);
    }
    let chr_id = ani
        .contig_id(parsed.chrom)
        .or_else(|| chr_name_to_id(parsed.chrom).map(u32::from));
    let alt_alleles: smallvec::SmallVec<[&str; 2]> = parsed.alt.split(',').collect();
    let mut keys: smallvec::SmallVec<[u64; 2]> = smallvec::SmallVec::with_capacity(alt_alleles.len());
    if let Some(chr_id) = chr_id {
        for alt in &alt_alleles {
            keys.push(make_variant_key(
                chr_id,
                parsed.pos,
                parsed.ref_allele.as_bytes(),
                alt.as_bytes(),
            ));
        }
    }
    Some(LineKeys {
        parsed,
        chr_id,
        alt_alleles,
        keys,
    })
}

/// Pre-batched annotation: identical semantics to [`annotate_line`] but
/// consumes precomputed MPH lookup results.
pub fn annotate_line_with_mph_hints<'a>(
    line: &'a str,
    line_keys: &LineKeys<'a>,
    mph_hits: &[Option<usize>],
    ani: &AniIndex,
    field_meta: &HashMap<String, FieldNumber>,
    column_modes: &[(String, AnnotateMode)],
    sample_map: &[Option<usize>],
    info_overwrite_all: bool,
    format_overwrite_all: bool,
) -> Option<String> {
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok();
    let parsed = &line_keys.parsed;
    let alt_alleles = &line_keys.alt_alleles;

    let want_format =
        format_overwrite_all || column_modes.iter().any(|(k, _)| column_targets_format(k));
    let info_only = columns_are_info_only(column_modes, format_overwrite_all);
    let entry_info_fast =
        can_use_entry_info_fast_path(info_overwrite_all, format_overwrite_all, column_modes);
    let need_info = info_overwrite_all || column_modes.iter().any(|(k, _)| column_targets_info(k));
    let need_format = want_format;
    let info_filter: Option<std::collections::HashSet<&str>> =
        compute_info_filter(info_overwrite_all, format_overwrite_all, column_modes);

    let raw_samples = if want_format && sample_map.iter().any(|v| v.is_none()) {
        let mut fields = line.split('\t');
        for _ in 0..9 {
            let _ = fields.next();
        }
        Some(fields.collect::<Vec<&str>>())
    } else {
        None
    };

    let Some(chr_id) = line_keys.chr_id else {
        return None;
    };
    let ref_allele = parsed.ref_allele;
    let pos = parsed.pos;
    let ref_hash = fast_hash64(ref_allele.as_bytes());

    if entry_info_fast && alt_alleles.len() == 1 {
        if let Some(entry_idx) =
            exact_mph_entry_idx(ani, chr_id, pos, ref_allele, alt_alleles[0], mph_hits)
        {
            let info = ani
                .entry_info_string(entry_idx)
                .unwrap_or_else(|| ".".to_string());
            if info == parsed.info {
                return None;
            }
            return Some(build_info_replaced_line(line, parsed, &info));
        }
    }

    let mut bundles = Vec::with_capacity(alt_alleles.len());
    for (vcf_idx, alt) in alt_alleles.iter().enumerate() {
        let bundle = mph_hits
            .get(vcf_idx)
            .copied()
            .flatten()
            .and_then(|mph_idx| {
                let entries = &ani.entries;
                if mph_idx >= entries.len() {
                    return None;
                }
                let e = &entries[mph_idx];
                if e.chr_id != chr_id || e.pos != pos {
                    return None;
                }
                let rf_str = ani.read_cstring(e.ref_ofs as usize);
                if rf_str.as_ref() != ref_allele {
                    return None;
                }
                let alt_str = ani.read_cstring(e.alt_ofs as usize);
                if alt_str.as_ref() != *alt {
                    return None;
                }
                Some(ani.build_bundle_from_entry_idx_opts_with_meta_filtered(
                    mph_idx,
                    field_meta,
                    need_info,
                    need_format,
                    info_filter.as_ref(),
                ))
            })
            .or_else(|| {
                ani.lookup_exact_by_chr_id_pos_index_opts_filtered(
                    chr_id,
                    pos,
                    ref_allele,
                    ref_hash,
                    alt,
                    field_meta,
                    need_info,
                    need_format,
                    info_filter.as_ref(),
                )
            });

        if let Some(bundle) = bundle {
            bundles.push((vcf_idx, bundle));
        }
    }

    if bundles.is_empty() {
        return None;
    }

    if info_only {
        return annotate_info_only_with_bundles(
            line,
            parsed,
            alt_alleles,
            &bundles,
            field_meta,
            column_modes,
            info_overwrite_all,
        );
    }

    Some(annotate_record_with_alts(
        parsed,
        alt_alleles,
        &bundles,
        field_meta,
        column_modes,
        sample_map,
        if !want_format {
            Some(tail_after_info(line, parsed))
        } else {
            None
        },
        raw_samples.as_deref(),
        info_overwrite_all,
        format_overwrite_all,
        debug,
    ))
}

/// Timing-instrumented annotate. Returns `None` when the line is unchanged.
pub fn annotate_line_with_timing(
    line: &str,
    ani: &AniIndex,
    field_meta: &HashMap<String, FieldNumber>,
    column_modes: &[(String, AnnotateMode)],
    sample_map: &[Option<usize>],
    info_overwrite_all: bool,
    format_overwrite_all: bool,
    acc: &BundleTimingAccum,
) -> Option<String> {
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok();

    let want_format =
        format_overwrite_all || column_modes.iter().any(|(k, _)| column_targets_format(k));
    let info_only = columns_are_info_only(column_modes, format_overwrite_all);

    let mut parsed = parse_vcf_record_simd(line, want_format)?;
    if want_format {
        patch_samples_from_line(&mut parsed, line);
    }
    let raw_tail = if !want_format {
        Some(tail_after_info(line, &parsed))
    } else {
        None
    };

    let raw_samples = if want_format && sample_map.iter().any(|v| v.is_none()) {
        let mut fields = line.split('\t');
        for _ in 0..9 {
            let _ = fields.next();
        }
        Some(fields.collect::<Vec<&str>>())
    } else {
        None
    };

    let chrom = &parsed.chrom;
    let pos = parsed.pos;
    let ref_allele = &parsed.ref_allele;
    let alt_alleles: Vec<&str> = parsed.alt.split(',').collect();

    if debug {
        eprintln!(
            "[ANNOTATE] {}:{} {} -> {:?}",
            chrom, pos, ref_allele, alt_alleles
        );
    }

    let chr_id = ani
        .contig_id(chrom)
        .or_else(|| chr_name_to_id(chrom).map(u32::from))?;
    let ref_hash = fast_hash64(ref_allele.as_bytes());

    let need_info = info_overwrite_all || column_modes.iter().any(|(k, _)| column_targets_info(k));
    let need_format = want_format;

    let mut bundles = Vec::with_capacity(alt_alleles.len());
    for (vcf_idx, alt) in alt_alleles.iter().enumerate() {
        if let Some(bundle) = ani.lookup_exact_by_chr_id_pos_index_opts(
            chr_id,
            pos,
            ref_allele,
            ref_hash,
            alt,
            field_meta,
            need_info,
            need_format,
        ) {
            if debug {
                eprintln!(
                    "[ANNOTATE] Found bundle for ALT {}: INFO fields={}",
                    alt,
                    bundle.info.len()
                );
            }
            bundles.push((vcf_idx, bundle));
        } else if let Some((bundle, t)) = ani.lookup_exact_by_chr_id_timed_opts(
            chr_id,
            pos,
            ref_allele,
            ref_hash,
            alt,
            need_info,
            need_format,
        ) {
            if debug {
                eprintln!(
                    "[ANNOTATE] Found bundle for ALT {}: INFO fields={}",
                    alt,
                    bundle.info.len()
                );
            }
            acc.add(&t);
            bundles.push((vcf_idx, bundle));
        }
    }

    if bundles.is_empty() {
        return None;
    }

    if info_only {
        return annotate_info_only_with_bundles(
            line,
            &parsed,
            &alt_alleles,
            &bundles,
            field_meta,
            column_modes,
            info_overwrite_all,
        );
    }

    Some(annotate_record_with_alts(
        &parsed,
        &alt_alleles,
        &bundles,
        field_meta,
        column_modes,
        sample_map,
        raw_tail,
        raw_samples.as_deref(),
        info_overwrite_all,
        format_overwrite_all,
        debug,
    ))
}

fn annotate_record(
    parsed: &ParsedVcfRecord,
    bundles: &[(usize, AnnotationBundle)],
    field_meta: &HashMap<String, FieldNumber>,
    column_modes: &[(String, AnnotateMode)],
    sample_map: &[Option<usize>],
    info_overwrite_all: bool,
    format_overwrite_all: bool,
    debug: bool,
) -> String {
    let alt_alleles: Vec<&str> = parsed.alt.split(',').collect();
    annotate_record_with_alts(
        parsed,
        &alt_alleles,
        bundles,
        field_meta,
        column_modes,
        sample_map,
        None,
        None,
        info_overwrite_all,
        format_overwrite_all,
        debug,
    )
}

/// INFO-only annotation. Returns `Some(new_line)` when the record changed,
/// `None` on a no-op merge so the writer can reuse the original line.
fn annotate_info_only_with_bundles(
    line: &str,
    parsed: &ParsedVcfRecord,
    alt_alleles: &[&str],
    bundles: &[(usize, AnnotationBundle)],
    field_meta: &HashMap<String, FieldNumber>,
    column_modes: &[(String, AnnotateMode)],
    info_overwrite_all: bool,
) -> Option<String> {
    use crate::annotate::cpu_v2::merge_info::merge_info_fields;

    let base_info = if info_overwrite_all && !bundles.is_empty() {
        "."
    } else {
        parsed.info
    };
    let merged = merge_info_fields(
        base_info,
        bundles,
        &None,
        parsed.ref_allele,
        alt_alleles,
        field_meta,
        column_modes,
    );
    let info_was_reset = info_overwrite_all && !bundles.is_empty() && parsed.info != ".";
    if !merged.dirty && !info_was_reset {
        return None;
    }
    let new_info = serialize_info_map(&merged.map);
    Some(build_info_replaced_line(line, parsed, &new_info))
}

fn annotate_record_with_alts(
    parsed: &ParsedVcfRecord,
    alt_alleles: &[&str],
    bundles: &[(usize, AnnotationBundle)],
    field_meta: &HashMap<String, FieldNumber>,
    column_modes: &[(String, AnnotateMode)],
    sample_map: &[Option<usize>],
    raw_tail: Option<&str>,
    raw_samples: Option<&[&str]>,
    info_overwrite_all: bool,
    format_overwrite_all: bool,
    debug: bool,
) -> String {
    use crate::annotate::cpu_v2::merge_info::merge_info_fields;

    let mut new_id = parsed.id.to_string();
    let mut new_qual = parsed.qual.to_string();
    let mut new_filter = parsed.filter.to_string();
    let mut new_alt = parsed.alt.to_string();
    let format_modified =
        format_overwrite_all || column_modes.iter().any(|(k, _)| column_targets_format(k));
    let mut new_format: Option<String> = if format_modified {
        parsed.format.as_ref().map(|f| f.raw.to_string())
    } else {
        None
    };
    let mut new_samples: Vec<String> = if format_modified {
        parsed
            .samples
            .iter()
            .map(|s| normalize_sample_raw(s.raw).into_owned())
            .collect()
    } else {
        Vec::new()
    };

    let base_info = if info_overwrite_all && !bundles.is_empty() {
        "."
    } else {
        parsed.info
    };

    let merged_info = merge_info_fields(
        base_info,
        bundles,
        &None,
        &parsed.ref_allele,
        &alt_alleles,
        field_meta,
        column_modes,
    );

    let new_info = serialize_info_map(&merged_info.map);

    if !bundles.is_empty() {
        let bundle = &bundles[0].1;

        for (key, mode) in column_modes {
            let (src_ref, dst_ref) = split_mapped_ref(key);
            match dst_ref.to_ascii_uppercase().as_str() {
                "ID" => {
                    if let Some(db_id) = read_scalar_from_bundle(bundle, src_ref) {
                        if let Some(v) = merge_scalar_field(&new_id, &db_id, *mode, ";") {
                            new_id = v;
                        }
                    }
                }
                "QUAL" => {
                    if let Some(db_qual) = read_scalar_from_bundle(bundle, src_ref) {
                        if let Some(v) = merge_scalar_field(&new_qual, &db_qual, *mode, ";") {
                            new_qual = v;
                        }
                    }
                }
                "FILTER" => {
                    if let Some(db_filter) = read_scalar_from_bundle(bundle, src_ref) {
                        if let Some(v) = merge_scalar_field(&new_filter, &db_filter, *mode, ";") {
                            new_filter = v;
                        }
                    }
                }
                "ALT" => {
                    if let Some(db_alt) = read_scalar_from_bundle(bundle, src_ref) {
                        if let Some(v) = merge_scalar_field(&new_alt, &db_alt, *mode, ",") {
                            new_alt = v;
                        }
                    }
                }
                "FMT" | "FORMAT" => {
                    let (fmt, samples) = merge_all_format(
                        parsed,
                        bundle,
                        *mode,
                        sample_map,
                        raw_samples,
                        format_overwrite_all,
                        debug,
                    );
                    new_format = fmt;
                    new_samples = samples;
                }
                _ => {}
            }
        }
        if format_modified {
            apply_format_field_mappings(
                &mut new_format,
                &mut new_samples,
                bundle,
                column_modes,
                sample_map,
            );
        }
    }

    let mut out = String::new();
    out.push_str(&parsed.chrom);
    out.push('\t');
    out.push_str(&parsed.pos.to_string());
    out.push('\t');
    out.push_str(&new_id);
    out.push('\t');
    out.push_str(&parsed.ref_allele);
    out.push('\t');
    out.push_str(&new_alt);
    out.push('\t');
    out.push_str(&new_qual);
    out.push('\t');
    out.push_str(&new_filter);
    out.push('\t');
    out.push_str(&new_info);

    if format_modified {
        if let Some(fmt) = new_format {
            out.push('\t');
            out.push_str(&fmt);
            for s in &new_samples {
                out.push('\t');
                out.push_str(&s);
            }
        }
    } else if let Some(tail) = raw_tail {
        out.push_str(tail);
    } else if let Some(fmt) = parsed.format {
        out.push('\t');
        out.push_str(fmt.raw);
        for s in &parsed.samples {
            out.push('\t');
            out.push_str(s.raw);
        }
    }

    out
}

pub fn annotate_record_with_bundles(
    parsed: &ParsedVcfRecord,
    bundles: &[(usize, AnnotationBundle)],
    field_meta: &HashMap<String, FieldNumber>,
    column_modes: &[(String, AnnotateMode)],
    sample_map: &[Option<usize>],
    info_overwrite_all: bool,
    format_overwrite_all: bool,
    debug: bool,
) -> String {
    annotate_record(
        parsed,
        bundles,
        field_meta,
        column_modes,
        sample_map,
        info_overwrite_all,
        format_overwrite_all,
        debug,
    )
}

pub fn annotate_record_with_bundles_and_info(
    parsed: &ParsedVcfRecord,
    bundles: &[(usize, AnnotationBundle)],
    field_meta: &HashMap<String, FieldNumber>,
    column_modes: &[(String, AnnotateMode)],
    sample_map: &[Option<usize>],
    info_overwrite_all: bool,
    format_overwrite_all: bool,
    debug: bool,
    merged_info: Option<&str>,
) -> String {
    annotate_record_with_alts_and_info(
        parsed,
        &parsed.alt.split(',').collect::<Vec<&str>>(),
        bundles,
        field_meta,
        column_modes,
        sample_map,
        None,
        None,
        info_overwrite_all,
        format_overwrite_all,
        debug,
        merged_info,
    )
}

fn annotate_record_with_alts_and_info(
    parsed: &ParsedVcfRecord,
    alt_alleles: &[&str],
    bundles: &[(usize, AnnotationBundle)],
    field_meta: &HashMap<String, FieldNumber>,
    column_modes: &[(String, AnnotateMode)],
    sample_map: &[Option<usize>],
    raw_tail: Option<&str>,
    raw_samples: Option<&[&str]>,
    info_overwrite_all: bool,
    format_overwrite_all: bool,
    debug: bool,
    merged_info: Option<&str>,
) -> String {
    use crate::annotate::cpu_v2::merge_info::merge_info_fields;

    let mut new_id = parsed.id.to_string();
    let mut new_qual = parsed.qual.to_string();
    let mut new_filter = parsed.filter.to_string();
    let mut new_alt = parsed.alt.to_string();
    let format_modified =
        format_overwrite_all || column_modes.iter().any(|(k, _)| column_targets_format(k));
    let mut new_format: Option<String> = if format_modified {
        parsed.format.as_ref().map(|f| f.raw.to_string())
    } else {
        None
    };
    let mut new_samples: Vec<String> = if format_modified {
        parsed
            .samples
            .iter()
            .map(|s| normalize_sample_raw(s.raw).into_owned())
            .collect()
    } else {
        Vec::new()
    };

    let new_info = if let Some(info) = merged_info {
        if info.is_empty() {
            ".".to_string()
        } else {
            info.to_string()
        }
    } else {
        let base_info = if info_overwrite_all && !bundles.is_empty() {
            "."
        } else {
            parsed.info
        };
        let merged = merge_info_fields(
            base_info,
            bundles,
            &None,
            parsed.ref_allele,
            alt_alleles,
            field_meta,
            column_modes,
        );
        serialize_info_map(&merged.map)
    };

    if !bundles.is_empty() {
        let bundle = &bundles[0].1;

        for (key, mode) in column_modes {
            let (src_ref, dst_ref) = split_mapped_ref(key);
            match dst_ref.to_ascii_uppercase().as_str() {
                "ID" => {
                    if let Some(db_id) = read_scalar_from_bundle(bundle, src_ref) {
                        if let Some(v) = merge_scalar_field(&new_id, &db_id, *mode, ";") {
                            new_id = v;
                        }
                    }
                }
                "QUAL" => {
                    if let Some(db_qual) = read_scalar_from_bundle(bundle, src_ref) {
                        if let Some(v) = merge_scalar_field(&new_qual, &db_qual, *mode, ";") {
                            new_qual = v;
                        }
                    }
                }
                "FILTER" => {
                    if let Some(db_filter) = read_scalar_from_bundle(bundle, src_ref) {
                        if let Some(v) = merge_scalar_field(&new_filter, &db_filter, *mode, ";") {
                            new_filter = v;
                        }
                    }
                }
                "ALT" => {
                    if let Some(db_alt) = read_scalar_from_bundle(bundle, src_ref) {
                        if let Some(v) = merge_scalar_field(&new_alt, &db_alt, *mode, ",") {
                            new_alt = v;
                        }
                    }
                }
                "FMT" | "FORMAT" => {
                    let (fmt, samples) = merge_all_format(
                        parsed,
                        bundle,
                        *mode,
                        sample_map,
                        raw_samples,
                        format_overwrite_all,
                        debug,
                    );
                    new_format = fmt;
                    new_samples = samples;
                }
                _ => {}
            }
        }
        if format_modified {
            apply_format_field_mappings(
                &mut new_format,
                &mut new_samples,
                bundle,
                column_modes,
                sample_map,
            );
        }
    }

    let mut out = String::new();
    out.push_str(&parsed.chrom);
    out.push('\t');
    out.push_str(&parsed.pos.to_string());
    out.push('\t');
    out.push_str(&new_id);
    out.push('\t');
    out.push_str(&parsed.ref_allele);
    out.push('\t');
    out.push_str(&new_alt);
    out.push('\t');
    out.push_str(&new_qual);
    out.push('\t');
    out.push_str(&new_filter);
    out.push('\t');
    out.push_str(&new_info);

    if format_modified {
        if let Some(fmt) = new_format {
            out.push('\t');
            out.push_str(&fmt);
            for s in &new_samples {
                out.push('\t');
                out.push_str(&s);
            }
        }
    } else if let Some(tail) = raw_tail {
        out.push_str(tail);
    } else if let Some(fmt) = parsed.format {
        out.push('\t');
        out.push_str(fmt.raw);
        for s in &parsed.samples {
            out.push('\t');
            out.push_str(s.raw);
        }
    }

    out
}

fn should_replace(existing: &str, new_val: &str, mode: AnnotateMode) -> bool {
    let existing_missing = is_missing_value(existing);
    let new_missing = is_missing_value(new_val);

    if new_missing && !mode.carry_over_missing {
        return false;
    }

    if mode.replace_all {
        return true;
    }

    if mode.replace_missing && existing_missing {
        return true;
    }

    if mode.replace_non_missing && !existing_missing {
        return true;
    }

    false
}

fn merge_scalar_field(
    existing: &str,
    new_val: &str,
    mode: AnnotateMode,
    sep: &str,
) -> Option<String> {
    let existing_missing = is_missing_value(existing);
    let new_missing = is_missing_value(new_val);
    if mode.set_or_append {
        if new_missing && !mode.carry_over_missing {
            return None;
        }
        if existing_missing {
            return Some(new_val.to_string());
        }
        if new_missing {
            return None;
        }
        return Some(format!("{existing}{sep}{new_val}"));
    }
    if should_replace(existing, new_val, mode) {
        return Some(new_val.to_string());
    }
    None
}

fn merge_all_format(
    parsed: &ParsedVcfRecord,
    bundle: &AnnotationBundle,
    mode: AnnotateMode,
    sample_map: &[Option<usize>],
    raw_samples: Option<&[&str]>,
    overwrite_all: bool,
    debug: bool,
) -> (Option<String>, Vec<String>) {
    if parsed.samples.is_empty() {
        return (parsed.format.as_ref().map(|f| f.raw.to_string()), Vec::new());
    }

    let db_format = match &bundle.format_str {
        Some(f) => f,
        None => {
            return (
                parsed.format.as_ref().map(|f| f.raw.to_string()),
                parsed.samples.iter().map(|s| s.raw.to_string()).collect(),
            );
        }
    };

    let db_keys: Vec<&str> = db_format.split(':').collect();
    let input_keys: Vec<String> = parsed
        .format
        .as_ref()
        .map(|f| f.keys().map(str::to_string).collect())
        .unwrap_or_default();

    const MISSING_IDX: usize = usize::MAX;
    let mut final_keys: Vec<String> = Vec::new();
    let mut input_idx_of_final: Vec<usize> = Vec::new();
    let mut db_idx_of_final: Vec<usize> = Vec::new();

    if overwrite_all {
        final_keys.reserve(db_keys.len());
        input_idx_of_final.reserve(db_keys.len());
        db_idx_of_final.reserve(db_keys.len());
        for (i, k) in db_keys.iter().enumerate() {
            final_keys.push((*k).to_string());
            let input_idx = input_keys
                .iter()
                .position(|ik| ik == k)
                .unwrap_or(MISSING_IDX);
            input_idx_of_final.push(input_idx);
            db_idx_of_final.push(i);
        }
    } else {
        final_keys.reserve(input_keys.len() + db_keys.len());
        input_idx_of_final.reserve(input_keys.len() + db_keys.len());
        db_idx_of_final.reserve(input_keys.len() + db_keys.len());
        for (i, k) in input_keys.iter().enumerate() {
            final_keys.push(k.clone());
            input_idx_of_final.push(i);
            let db_idx = db_keys
                .iter()
                .position(|dk| *dk == k)
                .unwrap_or(MISSING_IDX);
            db_idx_of_final.push(db_idx);
        }
        for (db_i, k) in db_keys.iter().enumerate() {
            if input_keys.iter().any(|ik| ik == k) {
                continue;
            }
            final_keys.push((*k).to_string());
            input_idx_of_final.push(MISSING_IDX);
            db_idx_of_final.push(db_i);
        }
    }

    if debug {
        eprintln!(
            "[FORMAT] input_keys={:?} db_keys={:?} final_keys={:?}",
            input_keys, db_keys, final_keys
        );
        eprintln!("[FORMAT] db samples: {:?}", bundle.format_samples);
    }

    let mut sample_order_identity = false;
    if overwrite_all && parsed.samples.len() == bundle.format_samples.len() {
        sample_order_identity = sample_map
            .iter()
            .enumerate()
            .all(|(i, v)| v.map(|idx| idx == i).unwrap_or(false));
    }

    let final_keys_is_input = !overwrite_all
        && final_keys.len() == input_keys.len()
        && final_keys
            .iter()
            .zip(input_keys.iter())
            .all(|(a, b)| a == b);
    let final_keys_is_db = overwrite_all
        && final_keys.len() == db_keys.len()
        && final_keys.iter().zip(db_keys.iter()).all(|(a, b)| *a == *b);

    if final_keys_is_input {
        for i in 0..final_keys.len() {
            input_idx_of_final[i] = i;
        }
    }
    if final_keys_is_db {
        for i in 0..final_keys.len() {
            db_idx_of_final[i] = i;
        }
    }

    let mut result_samples: Vec<String> = Vec::with_capacity(parsed.samples.len());

    for (input_idx, input_sample) in parsed.samples.iter().enumerate() {
        let mut db_idx = if sample_order_identity {
            Some(input_idx)
        } else {
            sample_map.get(input_idx).and_then(|v| *v)
        };
        if let Some(idx) = db_idx {
            if idx >= bundle.format_samples.len() {
                db_idx = None;
            }
        }

        if debug {
            eprintln!("[FORMAT] Sample {} -> DB idx {:?}", input_idx, db_idx);
        }

        let input_subfields: smallvec::SmallVec<[&str; 16]> =
            input_sample.raw.split(':').collect();

        let raw_split = if db_idx.is_none() {
            raw_samples
                .and_then(|rs| rs.get(input_idx))
                .map(|raw| raw.split(':').collect::<Vec<&str>>())
        } else {
            None
        };

        if !overwrite_all && final_keys_is_input && db_idx.is_none() {
            result_samples.push(normalize_sample_raw(input_sample.raw).into_owned());
            continue;
        }

        let mut sample_values: Vec<String> = Vec::with_capacity(final_keys.len());

        if overwrite_all && final_keys_is_db {
            if let Some(idx) = db_idx {
                let mut it = bundle
                    .format_samples
                    .get(idx)
                    .map(|s| s.split(':'))
                    .into_iter()
                    .flatten();
                for _ in 0..final_keys.len() {
                    let v = it.next().unwrap_or(".");
                    sample_values.push(v.to_string());
                }
            } else {
                for k_i in 0..final_keys.len() {
                    let key_idx = input_idx_of_final[k_i];
                    let mut v = match key_idx {
                        MISSING_IDX => None,
                        idx => {
                            if let Some(raw_vals) = raw_split.as_ref() {
                                raw_vals.get(idx).copied()
                            } else {
                                input_subfields.get(idx).copied()
                            }
                        }
                    };
                    if v.map_or(true, |s| s.is_empty()) && key_idx != MISSING_IDX {
                        v = if let Some(raw_vals) = raw_split.as_ref() {
                            raw_vals.get(key_idx).copied()
                        } else {
                            input_value_from_raw(raw_samples, input_idx, key_idx)
                        };
                    }
                    sample_values.push(v.unwrap_or(".").to_string());
                }
            }
        } else {
            for (k_i, _) in final_keys.iter().enumerate() {
                let key_idx = input_idx_of_final[k_i];
                let mut input_val = match key_idx {
                    MISSING_IDX => None,
                    idx => {
                        if let Some(raw_vals) = raw_split.as_ref() {
                            raw_vals.get(idx).copied()
                        } else {
                            input_subfields.get(idx).copied()
                        }
                    }
                };
                if input_val.map_or(true, |s| s.is_empty()) && key_idx != MISSING_IDX {
                    input_val = if let Some(raw_vals) = raw_split.as_ref() {
                        raw_vals.get(key_idx).copied()
                    } else {
                        input_value_from_raw(raw_samples, input_idx, key_idx)
                    };
                }
                let db_val = db_idx.and_then(|idx| {
                    let vals = bundle.format_samples.get(idx)?;
                    let key_idx = db_idx_of_final[k_i];
                    if key_idx == MISSING_IDX {
                        return None;
                    }
                    vals.split(':').nth(key_idx)
                });

                if overwrite_all {
                    let final_val = if db_idx.is_some() {
                        db_val.unwrap_or(".")
                    } else {
                        input_val.unwrap_or(".")
                    };
                    sample_values.push(final_val.to_string());
                } else {
                    let final_val = merge_format_value_str(input_val, db_val, mode);
                    sample_values.push(final_val);
                }
            }
        }

        result_samples.push(normalize_sample_values(&sample_values));
    }

    let format_out = if final_keys_is_db {
        db_format.clone()
    } else {
        join_keys(&final_keys)
    };
    (Some(format_out), result_samples)
}

fn merge_format_value_str(
    input: Option<&str>,
    db: Option<&str>,
    mode: AnnotateMode,
) -> String {
    let input_val = match input {
        Some(v) if !v.is_empty() => v,
        _ => ".",
    };
    let db_val = match db {
        Some(v) if !v.is_empty() => v,
        _ => ".",
    };

    let input_missing = is_missing_value(input_val);
    let db_missing = is_missing_value(db_val);

    if mode.replace_all {
        if !db_missing || mode.carry_over_missing {
            return db_val.to_string();
        }
    }

    if mode.replace_missing && input_missing {
        if !db_missing || mode.carry_over_missing {
            return db_val.to_string();
        }
    }

    if mode.replace_non_missing && !input_missing {
        if !db_missing || mode.carry_over_missing {
            return db_val.to_string();
        }
    }

    if mode.set_or_append {
        return merge_scalar_field(input_val, db_val, mode, ",")
            .unwrap_or_else(|| input_val.to_string());
    }

    input_val.to_string()
}

fn normalize_sample_values(values: &[String]) -> String {
    if values.is_empty() {
        return ".".to_string();
    }

    let mut out = String::new();
    for (i, v) in values.iter().enumerate() {
        if i > 0 {
            out.push(':');
        }
        if v.is_empty() {
            out.push('.');
        } else {
            out.push_str(v);
        }
    }
    out
}

fn normalize_sample_values_strs(values: &[&str]) -> String {
    if values.is_empty() {
        return ".".to_string();
    }

    let mut out = String::new();
    for (i, v) in values.iter().enumerate() {
        if i > 0 {
            out.push(':');
        }
        if v.is_empty() {
            out.push('.');
        } else {
            out.push_str(v);
        }
    }
    out
}

/// Normalises an already-joined sample string (e.g. `"0/1::99"` → `"0/1:.:99"`).
fn normalize_sample_raw(raw: &str) -> std::borrow::Cow<'_, str> {
    if raw.is_empty() {
        return std::borrow::Cow::Borrowed(".");
    }
    let mut prev_is_colon = true;
    let mut needs_rewrite = false;
    for b in raw.as_bytes() {
        if *b == b':' {
            if prev_is_colon {
                needs_rewrite = true;
                break;
            }
            prev_is_colon = true;
        } else {
            prev_is_colon = false;
        }
    }
    if prev_is_colon {
        needs_rewrite = true;
    }
    if !needs_rewrite {
        return std::borrow::Cow::Borrowed(raw);
    }
    let parts: Vec<&str> = raw.split(':').collect();
    std::borrow::Cow::Owned(normalize_sample_values_strs(&parts))
}

fn input_value_from_raw<'a>(
    raw_samples: Option<&'a [&'a str]>,
    input_idx: usize,
    key_idx: usize,
) -> Option<&'a str> {
    let raw = raw_samples?.get(input_idx)?;
    raw.split(':').nth(key_idx)
}

fn join_keys(values: &[String]) -> String {
    join_values(values)
}

fn join_values(values: &[String]) -> String {
    if values.is_empty() {
        return String::new();
    }
    let mut len = values.iter().map(|v| v.len()).sum::<usize>();
    len += values.len().saturating_sub(1);
    let mut out = String::with_capacity(len);
    for (i, v) in values.iter().enumerate() {
        if i > 0 {
            out.push(':');
        }
        out.push_str(v);
    }
    out
}

#[cfg(test)]
#[path = "../../../tests/unit/annotate_cpu_v2_annotation.rs"]
mod tests;
