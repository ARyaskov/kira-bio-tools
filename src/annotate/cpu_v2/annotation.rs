use std::collections::HashMap;

use crate::annotate::cpu_v2::field_metadata::is_missing_value;
use crate::annotate::cpu_v2::vcf_parsing::{parse_vcf_record, ParsedVcfRecord};
use crate::annotate::structs::ani::AniIndex;
use crate::annotate::structs::annotate_mode::AnnotateMode;
use crate::annotate::structs::bundle::{AnnotationBundle, FieldNumber};

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

    let Some(parsed) = parse_vcf_record(line) else {
        return line.to_string();
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

    let mut bundles = Vec::new();
    for (vcf_idx, alt) in alt_alleles.iter().enumerate() {
        if let Some(bundle) = ani.lookup_exact(chrom, pos, ref_allele, alt) {
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

    annotate_record(
        &parsed,
        &bundles,
        field_meta,
        column_modes,
        sample_map,
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
    use crate::annotate::cpu_v2::merge_info::merge_info_fields;

    let mut new_id = parsed.id.clone();
    let mut new_qual = parsed.qual.clone();
    let mut new_filter = parsed.filter.clone();
    let mut new_format: Option<String> = parsed.format.as_ref().map(|f| f.keys.join(":"));
    let mut new_samples: Vec<String> = parsed
        .samples
        .iter()
        .map(|s| normalize_sample_values(&s.raw))
        .collect();

    let alt_alleles: Vec<&str> = parsed.alt.split(',').collect();

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
        info_map
            .iter()
            .map(|(k, v)| {
                if v.is_empty() {
                    k.clone()
                } else {
                    format!("{}={}", k, v)
                }
            })
            .collect::<Vec<_>>()
            .join(";")
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
        for s in new_samples {
            out.push('\t');
            out.push_str(&s);
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

fn merge_all_format(
    parsed: &ParsedVcfRecord,
    bundle: &AnnotationBundle,
    mode: AnnotateMode,
    sample_map: &[Option<usize>],
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
                    parsed.format.as_ref().map(|f| f.keys.join(":")),
                    parsed.samples.iter().map(|s| s.raw.join(":")).collect(),
                );
            }
            return (
                parsed.format.as_ref().map(|f| f.keys.join(":")),
                parsed.samples.iter().map(|s| s.raw.join(":")).collect(),
            );
        }
    };

    let db_keys: Vec<&str> = db_format.split(':').collect();
    let input_keys: Vec<String> = parsed
        .format
        .as_ref()
        .map(|f| f.keys.clone())
        .unwrap_or_default();

    let mut final_keys: Vec<String> = Vec::new();
    let mut key_set = std::collections::HashSet::new();

    if overwrite_all {
        for k in &db_keys {
            if key_set.insert(k.to_string()) {
                final_keys.push(k.to_string());
            }
        }
    } else {
        for k in &input_keys {
            if key_set.insert(k.clone()) {
                final_keys.push(k.clone());
            }
        }
        for k in &db_keys {
            if key_set.insert(k.to_string()) {
                final_keys.push(k.to_string());
            }
        }
    }

    if debug {
        eprintln!(
            "[FORMAT] input_keys={:?} db_keys={:?} final_keys={:?}",
            input_keys, db_keys, final_keys
        );
        eprintln!("[FORMAT] db samples: {:?}", bundle.format_samples);
    }

    let db_sample_values: Vec<Vec<&str>> = bundle
        .format_samples
        .iter()
        .map(|s| s.split(':').collect())
        .collect();

    let input_index: std::collections::HashMap<&str, usize> = input_keys
        .iter()
        .enumerate()
        .map(|(i, k)| (k.as_str(), i))
        .collect();
    let db_index: std::collections::HashMap<&str, usize> =
        db_keys.iter().enumerate().map(|(i, k)| (*k, i)).collect();

    let mut result_samples: Vec<String> = Vec::new();

    for (input_idx, input_sample) in parsed.samples.iter().enumerate() {
        let db_idx = sample_map.get(input_idx).and_then(|v| *v);

        if debug {
            eprintln!("[FORMAT] Sample {} -> DB idx {:?}", input_idx, db_idx);
        }

        let mut sample_values: Vec<String> = Vec::new();

        for key in &final_keys {
            let input_val = input_index
                .get(key.as_str())
                .and_then(|idx| input_sample.raw.get(*idx).map(|s| s.as_str()));
            let db_val = db_idx.and_then(|idx| {
                db_sample_values.get(idx).and_then(|vals| {
                    let key_idx = db_index.get(key.as_str())?;
                    vals.get(*key_idx).map(|&s| s)
                })
            });

            if overwrite_all {
                let final_val = if db_idx.is_some() {
                    db_val.unwrap_or(".").to_string()
                } else {
                    input_val.unwrap_or(".").to_string()
                };
                sample_values.push(final_val);
            } else {
                let final_val = merge_format_value(input_val, db_val, mode);
                sample_values.push(final_val);
            }
        }

        result_samples.push(normalize_sample_values(&sample_values));
    }

    (Some(final_keys.join(":")), result_samples)
}

fn merge_format_value(input: Option<&str>, db: Option<&str>, mode: AnnotateMode) -> String {
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
            true,
            false,
        );

        assert_eq!(fmt, Some("GT:FINT:FFLT:FSTR".to_string()));
        assert_eq!(samples[0], "0/1:77:7.7:AAA_DB");
        assert_eq!(samples[1], "1/1:88:8.8:BBB_DB");
        assert_eq!(samples[2], "0/0:33:3.3:CCC");
    }
}
