use indexmap::IndexMap;
use std::collections::HashMap;

use super::field_metadata::{infer_field_type, is_missing_value};
use super::merge_info_helpers::*;
use crate::annotate::structs::annotate_mode::AnnotateMode;
use crate::annotate::structs::bundle::{AnnotationBundle, FieldNumber, StructuredInfoField};

pub fn merge_info_fields(
    existing_info: &str,
    bundles: &[(usize, AnnotationBundle)],
    multiallelic_bundle: &Option<AnnotationBundle>,
    vcf_alt_alleles: &[&str],
    field_meta: &HashMap<String, FieldNumber>,
    column_specs: &[(String, AnnotateMode)],
) -> IndexMap<String, String> {
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok();

    if debug {
        eprintln!("[MERGE] Starting merge_info_fields");
        eprintln!("[MERGE] Existing info: '{}'", existing_info);
        eprintln!("[MERGE] VCF alleles: {:?}", vcf_alt_alleles);
        eprintln!("[MERGE] Column specs: {:?}", column_specs);
    }

    let mut info_map = parse_existing_info(existing_info);

    let mut bundle_ref_values: HashMap<String, String> = HashMap::new();
    for (_vcf_idx, bundle) in bundles {
        for field in &bundle.info {
            if let Some(field_number) = field_meta.get(&field.key) {
                if *field_number == FieldNumber::R {
                    if let Some(ref_val) = field.values.first() {
                        bundle_ref_values.insert(field.key.clone(), ref_val.clone());
                    }
                }
            }
        }
    }

    if let Some(ref bundle) = multiallelic_bundle {
        for field in &bundle.info {
            if let Some(field_number) = field_meta.get(&field.key) {
                if *field_number == FieldNumber::R && !bundle_ref_values.contains_key(&field.key) {
                    if let Some(ref_val) = field.values.first() {
                        bundle_ref_values.insert(field.key.clone(), ref_val.clone());
                    }
                }
            }
        }
    }

    if debug {
        eprintln!("[MERGE] Bundle REF values: {:?}", bundle_ref_values);
        eprintln!("[MERGE] Parsed existing info: {:?}", info_map);
    }

    for (spec_key, mode) in column_specs {
        let key = spec_key.strip_prefix("INFO/").unwrap_or(spec_key);

        if debug {
            eprintln!("[MERGE] Processing field: {} (mode: {})", key, mode);
        }

        let field_number = field_meta.get(key).copied().unwrap_or(FieldNumber::One);
        let vcf_has_field = info_map.contains_key(key);
        let existing_val = info_map.get(key).cloned().unwrap_or_default();

        if field_number == FieldNumber::Zero {
            let has_flag = bundles
                .iter()
                .any(|(_, bundle)| bundle.info.iter().any(|f| f.key == key));
            let has_annotation_data = has_flag;

            if debug {
                eprintln!(
                    "[MERGE] Field {} (Flag): vcf_has_field={}, existing_val='{}', has_annotation_data={}",
                    key, vcf_has_field, existing_val, has_annotation_data
                );
            }

            let should_transfer =
                if mode.replace_missing && !mode.replace_non_missing && !mode.replace_all {
                    has_annotation_data
                } else {
                    mode.should_transfer(
                        !has_annotation_data,
                        vcf_has_field,
                        is_missing_value(&existing_val),
                    )
                };

            if should_transfer && has_flag {
                info_map.insert(key.to_string(), String::new());
            }

            continue;
        }

        if field_number != FieldNumber::A && field_number != FieldNumber::R {
            let mut annotated_val: Option<String> = None;
            for (_vcf_idx, bundle) in bundles {
                if let Some(field) = bundle.info.iter().find(|f| f.key == key) {
                    let joined = field.values.join(",");
                    if !is_missing_value(&joined) {
                        annotated_val = Some(joined);
                        break;
                    }
                }
            }

            let has_annotation_data = annotated_val
                .as_ref()
                .map(|v| !is_missing_value(v))
                .unwrap_or(false);

            if debug {
                eprintln!(
                    "[MERGE] Field {} (Record): vcf_has_field={}, existing_val='{}', annotated_val={:?}",
                    key, vcf_has_field, existing_val, annotated_val
                );
            }

            let should_transfer =
                if mode.replace_missing && !mode.replace_non_missing && !mode.replace_all {
                    has_annotation_data
                } else {
                    mode.should_transfer(
                        !has_annotation_data,
                        vcf_has_field,
                        is_missing_value(&existing_val),
                    )
                };

            if should_transfer {
                if let Some(val) = annotated_val {
                    info_map.insert(key.to_string(), val);
                }
            }

            continue;
        }

        let vcf_values = collect_annotations_for_field(
            key,
            bundles,
            multiallelic_bundle,
            vcf_alt_alleles,
            field_meta,
        );

        let has_annotation_data = vcf_values.iter().any(|v| {
            if let Some(val) = v {
                !is_missing_value(val)
            } else {
                false
            }
        });

        if debug {
            eprintln!("[MERGE] Field {}: vcf_values={:?}, vcf_has_field={}, existing_val='{}', has_annotation_data={}",
                      key, vcf_values, vcf_has_field, existing_val, has_annotation_data);
        }

        let should_transfer =
            if mode.replace_missing && !mode.replace_non_missing && !mode.replace_all {
                has_annotation_data
            } else {
                mode.should_transfer(
                    !has_annotation_data,
                    vcf_has_field,
                    is_missing_value(&existing_val),
                )
            };

        if debug {
            eprintln!(
                "[MERGE] Field {}: should_transfer={} (special +TAG logic: {})",
                key,
                should_transfer,
                mode.replace_missing && !mode.replace_non_missing && !mode.replace_all
            );
        }

        if !should_transfer {
            continue;
        }

        let existing_parts: Vec<&str> = if !existing_val.is_empty() && existing_val != "." {
            existing_val.split(',').collect()
        } else {
            vec![]
        };

        let field_type = infer_field_type(key);
        let is_integer = field_type == "Integer";

        if debug {
            eprintln!(
                "[MERGE] Field {}: field_number={:?}, existing_parts={:?}",
                key, field_number, existing_parts
            );
        }

        let final_values: Vec<String> = match field_number {
            FieldNumber::A => merge_field_number_a(
                &vcf_values,
                &existing_parts,
                *mode,
                is_integer,
                vcf_has_field,
            ),
            FieldNumber::R => merge_field_number_r(
                &vcf_values,
                &bundle_ref_values,
                key,
                &existing_parts,
                *mode,
                is_integer,
                vcf_has_field,
            ),
            _ => vcf_values.iter().filter_map(|v| v.clone()).collect(),
        };

        if debug {
            eprintln!("[MERGE] Field {}: final_values={:?}", key, final_values);
        }

        if !final_values.is_empty() && !final_values.iter().all(|v| is_missing_value(v)) {
            info_map.insert(key.to_string(), final_values.join(","));
            if debug {
                eprintln!("[MERGE] Field {}: inserted into info_map", key);
            }
        } else if debug {
            eprintln!("[MERGE] Field {}: not inserted (empty or all missing)", key);
        }
    }

    if debug {
        eprintln!("[MERGE] Final info_map: {:?}", info_map);
    }

    info_map
}

fn collect_annotations_for_field(
    key: &str,
    bundles: &[(usize, AnnotationBundle)],
    multiallelic_bundle: &Option<AnnotationBundle>,
    vcf_alt_alleles: &[&str],
    field_meta: &HashMap<String, FieldNumber>,
) -> Vec<Option<String>> {
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok();
    let mut values = vec![None; vcf_alt_alleles.len()];
    let field_number = field_meta.get(key).copied().unwrap_or(FieldNumber::One);

    if debug {
        eprintln!(
            "[COLLECT] Field {} (Number={:?}) for alleles {:?}",
            key, field_number, vcf_alt_alleles
        );
    }

    for (vcf_idx, bundle) in bundles {
        if let Some(field) = bundle.info.iter().find(|f| f.key == key) {
            if debug {
                eprintln!(
                    "[COLLECT] Found exact match for allele {}: {} -> {:?}",
                    vcf_idx, vcf_alt_alleles[*vcf_idx], field.values
                );
            }

            match field_number {
                FieldNumber::A => {
                    if let Some(val) = field.values.first() {
                        if *vcf_idx < values.len() {
                            values[*vcf_idx] = Some(val.clone());
                        }
                    }
                }
                FieldNumber::R => {
                    if let Some(alt_val) = field.values.get(1) {
                        if *vcf_idx < values.len() {
                            values[*vcf_idx] = Some(alt_val.clone());

                            if debug {
                                eprintln!(
                                    "[COLLECT] R field exact match: REF={:?}, ALT[{}]={}",
                                    field.values.first(),
                                    vcf_idx,
                                    alt_val
                                );
                            }
                        }
                    }
                }
                _ => {
                    if let Some(val) = field.values.first() {
                        if *vcf_idx < values.len() {
                            values[*vcf_idx] = Some(val.clone());
                        }
                    }
                }
            }
        }
    }

    if let Some(ref bundle) = multiallelic_bundle {
        if let Some(field) = bundle.info.iter().find(|f| f.key == key) {
            if debug {
                eprintln!(
                    "[COLLECT] Processing multiallelic for field {}: alt={}, values={:?}",
                    key, bundle.alt, field.values
                );
            }

            let db_alts: Vec<&str> = bundle.alt.split(',').collect();

            match field_number {
                FieldNumber::A => {
                    for (vcf_idx, vcf_alt) in vcf_alt_alleles.iter().enumerate() {
                        if values[vcf_idx].is_none() {
                            if let Some(db_idx) =
                                db_alts.iter().position(|&db_alt| db_alt == *vcf_alt)
                            {
                                if let Some(val) = field.values.get(db_idx) {
                                    values[vcf_idx] = Some(val.clone());
                                    if debug {
                                        eprintln!("[COLLECT] Multiallelic A match: VCF[{}]={} -> DB[{}]={}", vcf_idx, vcf_alt, db_idx, val);
                                    }
                                }
                            }
                        }
                    }
                }
                FieldNumber::R => {
                    for (vcf_idx, vcf_alt) in vcf_alt_alleles.iter().enumerate() {
                        if values[vcf_idx].is_none() {
                            if let Some(db_idx) =
                                db_alts.iter().position(|&db_alt| db_alt == *vcf_alt)
                            {
                                if let Some(val) = field.values.get(db_idx + 1) {
                                    values[vcf_idx] = Some(val.clone());
                                    if debug {
                                        eprintln!("[COLLECT] Multiallelic R match: VCF[{}]={} -> DB[{}]={}", vcf_idx, vcf_alt, db_idx + 1, val);
                                    }
                                }
                            }
                        }
                    }
                }
                _ => {
                    if let Some(val) = field.values.first() {
                        for i in 0..values.len() {
                            if values[i].is_none() {
                                values[i] = Some(val.clone());
                            }
                        }
                    }
                }
            }
        }
    }

    if debug {
        eprintln!("[COLLECT] Final values for {}: {:?}", key, values);
    }

    values
}

fn merge_field_number_a(
    vcf_values: &[Option<String>],
    existing_parts: &[&str],
    mode: AnnotateMode,
    is_integer: bool,
    _vcf_has_field: bool,
) -> Vec<String> {
    let mut result = Vec::with_capacity(vcf_values.len());

    for (i, vcf_val) in vcf_values.iter().enumerate() {
        let existing = existing_parts.get(i).copied().unwrap_or(".");

        if mode.should_append() {
            match vcf_val {
                Some(val) if !is_missing_value(val) => {
                    result.push(merge_append_values(existing, val));
                }
                _ => result.push(existing.to_string()),
            }
            continue;
        }

        let annotated = vcf_val.as_ref().filter(|v| !is_missing_value(v));

        let overwrite = annotated.is_some() && (is_missing_value(existing) || is_integer);

        if overwrite {
            result.push(annotated.unwrap().clone());
        } else {
            result.push(existing.to_string());
        }
    }

    result
}

fn merge_field_number_r(
    vcf_values: &[Option<String>],
    bundle_ref_values: &HashMap<String, String>,
    key: &str,
    existing_parts: &[&str],
    mode: AnnotateMode,
    is_integer: bool,
    _vcf_has_field: bool,
) -> Vec<String> {
    let mut result = Vec::with_capacity(vcf_values.len() + 1);

    let ref_value = bundle_ref_values
        .get(key)
        .cloned()
        .unwrap_or_else(|| ".".to_string());

    result.push(ref_value);

    for (i, vcf_val) in vcf_values.iter().enumerate() {
        let existing = existing_parts.get(i + 1).copied().unwrap_or(".");

        if mode.should_append() {
            match vcf_val {
                Some(val) if !is_missing_value(val) => {
                    result.push(merge_append_values(existing, val));
                }
                _ => result.push(existing.to_string()),
            }
            continue;
        }

        let annotated = vcf_val.as_ref().filter(|v| !is_missing_value(v));

        let overwrite = annotated.is_some() && (is_missing_value(existing) || is_integer);

        if overwrite {
            result.push(annotated.unwrap().clone());
        } else {
            result.push(existing.to_string());
        }
    }

    result
}

fn parse_existing_info(info: &str) -> IndexMap<String, String> {
    let mut map = IndexMap::new();
    if info == "." || info.is_empty() {
        return map;
    }

    for pair in info.split(';') {
        if let Some(eq_pos) = pair.find('=') {
            let key = pair[..eq_pos].to_string();
            let value = pair[eq_pos + 1..].to_string();
            map.insert(key, value);
        } else {
            map.insert(pair.to_string(), String::new());
        }
    }

    map
}

pub fn format_info_string(info_map: &IndexMap<String, String>, field_order: &[String]) -> String {
    if info_map.is_empty() {
        return ".".to_string();
    }

    let mut ordered_keys: Vec<String> = Vec::new();
    let mut unordered_keys: Vec<String> = Vec::new();

    for key in info_map.keys() {
        if field_order.contains(key) {
            ordered_keys.push(key.clone());
        } else {
            unordered_keys.push(key.clone());
        }
    }

    ordered_keys.sort_by_key(|k| {
        field_order
            .iter()
            .position(|f| f == k)
            .unwrap_or(usize::MAX)
    });
    unordered_keys.sort();

    let mut all_keys = ordered_keys;
    all_keys.extend(unordered_keys);

    let parts: Vec<String> = all_keys
        .into_iter()
        .filter_map(|k| {
            let v = info_map.get(&k)?;
            if v.is_empty() {
                Some(k)
            } else {
                Some(format!("{}={}", k, v))
            }
        })
        .collect();

    if parts.is_empty() {
        ".".to_string()
    } else {
        parts.join(";")
    }
}
