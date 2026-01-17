use std::collections::HashMap;

use crate::annotate::cpu_v2::field_metadata::is_missing_value;
use crate::annotate::cpu_v2::vcf_parsing::{
    parse_vcf_record_simd, patch_samples_from_line, ParsedVcfRecord,
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

    let mut want_format = format_overwrite_all
        || column_modes
            .iter()
            .any(|(k, _)| k.eq_ignore_ascii_case("FMT") || k.eq_ignore_ascii_case("FORMAT"));

    if !want_format {
        let mut tabs = 0;
        for &b in line.as_bytes() {
            if b == b'\t' {
                tabs += 1;
                if tabs >= 8 {
                    want_format = true;
                    break;
                }
            }
        }
    }

    let Some(mut parsed) = parse_vcf_record_simd(line, want_format) else {
        return line.to_string();
    };
    patch_samples_from_line(&mut parsed, line);

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

    let chr_id = match chr_name_to_id(chrom) {
        Some(id) => id,
        None => return line.to_string(),
    };
    let ref_hash = fast_hash64(ref_allele.as_bytes());

    let need_info = info_overwrite_all
        || column_modes.iter().any(|(k, _)| {
            !(k.eq_ignore_ascii_case("ID")
                || k.eq_ignore_ascii_case("QUAL")
                || k.eq_ignore_ascii_case("FILTER")
                || k.eq_ignore_ascii_case("FMT")
                || k.eq_ignore_ascii_case("FORMAT"))
        });
    let need_format = want_format;

    let mut bundles = Vec::with_capacity(alt_alleles.len());
    for (vcf_idx, alt) in alt_alleles.iter().enumerate() {
        if let Some(bundle) = ani.lookup_exact_by_chr_id_opts(
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
            bundles.push((vcf_idx, bundle));
        }
    }

    if bundles.is_empty() {
        return line.to_string();
    }

    annotate_record_with_alts(
        &parsed,
        &alt_alleles,
        &bundles,
        field_meta,
        column_modes,
        sample_map,
        raw_samples.as_deref(),
        info_overwrite_all,
        format_overwrite_all,
        debug,
    )
}

pub fn annotate_line_with_timing(
    line: &str,
    ani: &AniIndex,
    field_meta: &HashMap<String, FieldNumber>,
    column_modes: &[(String, AnnotateMode)],
    sample_map: &[Option<usize>],
    info_overwrite_all: bool,
    format_overwrite_all: bool,
    acc: &BundleTimingAccum,
) -> String {
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok();

    let mut want_format = format_overwrite_all
        || column_modes
            .iter()
            .any(|(k, _)| k.eq_ignore_ascii_case("FMT") || k.eq_ignore_ascii_case("FORMAT"));

    if !want_format {
        let mut tabs = 0;
        for &b in line.as_bytes() {
            if b == b'\t' {
                tabs += 1;
                if tabs >= 8 {
                    want_format = true;
                    break;
                }
            }
        }
    }

    let Some(mut parsed) = parse_vcf_record_simd(line, want_format) else {
        return line.to_string();
    };
    patch_samples_from_line(&mut parsed, line);

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

    let chr_id = match chr_name_to_id(chrom) {
        Some(id) => id,
        None => return line.to_string(),
    };
    let ref_hash = fast_hash64(ref_allele.as_bytes());

    let need_info = info_overwrite_all
        || column_modes.iter().any(|(k, _)| {
            !(k.eq_ignore_ascii_case("ID")
                || k.eq_ignore_ascii_case("QUAL")
                || k.eq_ignore_ascii_case("FILTER")
                || k.eq_ignore_ascii_case("FMT")
                || k.eq_ignore_ascii_case("FORMAT"))
        });
    let need_format = want_format;

    let mut bundles = Vec::with_capacity(alt_alleles.len());
    for (vcf_idx, alt) in alt_alleles.iter().enumerate() {
        if let Some((bundle, t)) = ani.lookup_exact_by_chr_id_timed_opts(
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
        return line.to_string();
    }

    annotate_record_with_alts(
        &parsed,
        &alt_alleles,
        &bundles,
        field_meta,
        column_modes,
        sample_map,
        raw_samples.as_deref(),
        info_overwrite_all,
        format_overwrite_all,
        debug,
    )
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
        info_overwrite_all,
        format_overwrite_all,
        debug,
    )
}

fn annotate_record_with_alts(
    parsed: &ParsedVcfRecord,
    alt_alleles: &[&str],
    bundles: &[(usize, AnnotationBundle)],
    field_meta: &HashMap<String, FieldNumber>,
    column_modes: &[(String, AnnotateMode)],
    sample_map: &[Option<usize>],
    raw_samples: Option<&[&str]>,
    info_overwrite_all: bool,
    format_overwrite_all: bool,
    debug: bool,
) -> String {
    use crate::annotate::cpu_v2::merge_info::merge_info_fields;

    let mut new_id = parsed.id.clone();
    let mut new_qual = parsed.qual.clone();
    let mut new_filter = parsed.filter.clone();
    let mut new_format: Option<String> = parsed.format.as_ref().map(|f| join_keys(&f.keys));
    let mut new_samples: Vec<String> = parsed
        .samples
        .iter()
        .map(|s| normalize_sample_values(&s.raw))
        .collect();
    let format_modified = format_overwrite_all
        || column_modes
            .iter()
            .any(|(k, _)| k.eq_ignore_ascii_case("FMT") || k.eq_ignore_ascii_case("FORMAT"));
    if !format_modified {
        if let Some(raw) = raw_samples {
            if new_samples.len() != raw.len() {
                new_samples = raw
                    .iter()
                    .map(|s| {
                        if s.is_empty() {
                            ".".to_string()
                        } else {
                            (*s).to_string()
                        }
                    })
                    .collect();
            }
        }
    }

    let base_info = if info_overwrite_all && !bundles.is_empty() {
        "."
    } else {
        parsed.info.as_str()
    };

    let info_map = merge_info_fields(
        base_info,
        bundles,
        &None,
        &alt_alleles,
        field_meta,
        column_modes,
    );

    let new_info = if info_map.is_empty() {
        ".".to_string()
    } else {
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
        out
    };

    if !bundles.is_empty() {
        let bundle = &bundles[0].1;

        for (key, mode) in column_modes {
            match key.as_str() {
                "ID" => {
                    if let Some(ref db_id) = bundle.id {
                        if should_replace(&parsed.id, db_id, *mode) {
                            new_id = db_id.clone();
                        }
                    }
                }
                "QUAL" => {
                    if let Some(ref db_qual) = bundle.qual {
                        if should_replace(&parsed.qual, db_qual, *mode) {
                            new_qual = db_qual.clone();
                        }
                    }
                }
                "FILTER" => {
                    if let Some(ref db_filter) = bundle.filter {
                        if should_replace(&parsed.filter, db_filter, *mode) {
                            new_filter = db_filter.clone();
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
    out.push_str(&parsed.alt);
    out.push('\t');
    out.push_str(&new_qual);
    out.push('\t');
    out.push_str(&new_filter);
    out.push('\t');
    out.push_str(&new_info);

    if let Some(fmt) = new_format {
        out.push('\t');
        out.push_str(&fmt);
        for s in &new_samples {
            out.push('\t');
            out.push_str(&s);
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
        return (parsed.format.as_ref().map(|f| f.keys.join(":")), Vec::new());
    }

    let db_format = match &bundle.format_str {
        Some(f) => f,
        None => {
            if overwrite_all {
                return (
                    parsed.format.as_ref().map(|f| join_keys(&f.keys)),
                    parsed.samples.iter().map(|s| join_values(&s.raw)).collect(),
                );
            }
            return (
                parsed.format.as_ref().map(|f| join_keys(&f.keys)),
                parsed.samples.iter().map(|s| join_values(&s.raw)).collect(),
            );
        }
    };

    let db_keys: Vec<&str> = db_format.split(':').collect();
    let input_keys: Vec<String> = parsed
        .format
        .as_ref()
        .map(|f| f.keys.clone())
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

        let raw_split = if db_idx.is_none() {
            raw_samples
                .and_then(|rs| rs.get(input_idx))
                .map(|raw| raw.split(':').collect::<Vec<&str>>())
        } else {
            None
        };

        if !overwrite_all && final_keys_is_input && db_idx.is_none() {
            result_samples.push(normalize_sample_values(&input_sample.raw));
            continue;
        }

        let mut sample_values: Vec<&str> = Vec::with_capacity(final_keys.len());

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
                    sample_values.push(v);
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
                                input_sample.raw.get(idx).map(|s| s.as_str())
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
                    sample_values.push(v.unwrap_or("."));
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
                            input_sample.raw.get(idx).map(|s| s.as_str())
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
                    sample_values.push(final_val);
                } else {
                    let final_val = merge_format_value_str(input_val, db_val, mode);
                    sample_values.push(final_val);
                }
            }
        }

        result_samples.push(normalize_sample_values_strs(&sample_values));
    }

    let format_out = if final_keys_is_db {
        db_format.clone()
    } else {
        join_keys(&final_keys)
    };
    (Some(format_out), result_samples)
}

fn merge_format_value_str<'a>(
    input: Option<&'a str>,
    db: Option<&'a str>,
    mode: AnnotateMode,
) -> &'a str {
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
            return db_val;
        }
    }

    if mode.replace_missing && input_missing {
        if !db_missing || mode.carry_over_missing {
            return db_val;
        }
    }

    input_val
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
mod tests {
    use super::*;
    use crate::annotate::cpu_v2::{ParsedFormat, ParsedSample};

    #[test]
    fn test_empty_sample_column_is_dot() {
        let vals = vec![String::new()];
        assert_eq!(normalize_sample_values(&vals), ".");
    }

    #[test]
    fn test_missing_gt_with_format_fields() {
        let vals = vec![String::new(), String::new(), String::new(), String::new()];
        assert_eq!(normalize_sample_values(&vals), ".:.:.:.");
    }

    #[test]
    fn test_missing_subfields_are_dots() {
        let vals = vec![
            "0/0".to_string(),
            "".to_string(),
            "1.1".to_string(),
            "".to_string(),
        ];
        assert_eq!(normalize_sample_values(&vals), "0/0:.:1.1:.");
    }

    #[test]
    fn test_format_samples_mapped_by_name() {
        let parsed = ParsedVcfRecord {
            chrom: "1".to_string(),
            pos: 1,
            id: ".".to_string(),
            ref_allele: "A".to_string(),
            alt: "C".to_string(),
            qual: ".".to_string(),
            filter: ".".to_string(),
            info: ".".to_string(),
            format: Some(ParsedFormat {
                keys: vec![
                    "GT".to_string(),
                    "FINT".to_string(),
                    "FFLT".to_string(),
                    "FSTR".to_string(),
                ],
            }),
            samples: vec![
                ParsedSample {
                    raw: vec![
                        "0/0".to_string(),
                        "11".to_string(),
                        "1.1".to_string(),
                        "AAA".to_string(),
                    ],
                },
                ParsedSample {
                    raw: vec![
                        "0/1".to_string(),
                        "22".to_string(),
                        "2.2".to_string(),
                        "BBB".to_string(),
                    ],
                },
                ParsedSample {
                    raw: vec![
                        "0/0".to_string(),
                        "33".to_string(),
                        "3.3".to_string(),
                        "CCC".to_string(),
                    ],
                },
            ],
        };

        let bundle = AnnotationBundle {
            alt: "C".to_string(),
            id: None,
            qual: None,
            filter: None,
            info: Vec::new(),
            format_str: Some("GT:FINT:FFLT:FSTR".to_string()),
            format_samples: vec![
                "1/1:88:8.8:BBB_DB".to_string(), // db sample B
                "0/1:77:7.7:AAA_DB".to_string(), // db sample A
            ],
        };

        let sample_map = vec![Some(1), Some(0), None]; // input A,B,C -> db A,B,missing

        let (fmt, samples) = merge_all_format(
            &parsed,
            &bundle,
            AnnotateMode::default_mode(),
            &sample_map,
            None,
            true,
            false,
        );

        assert_eq!(fmt, Some("GT:FINT:FFLT:FSTR".to_string()));
        assert_eq!(samples[0], "0/1:77:7.7:AAA_DB");
        assert_eq!(samples[1], "1/1:88:8.8:BBB_DB");
        assert_eq!(samples[2], "0/0:33:3.3:CCC");
    }
}
