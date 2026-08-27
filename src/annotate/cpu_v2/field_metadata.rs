use anyhow::Result;
use std::collections::{HashMap, HashSet};

use crate::annotate::structs::ani::{ANI_HEADER_END, AniIndex};
use crate::annotate::structs::bundle::FieldNumber;
use crate::util::{
    choose_best_number, extract_info_key, extract_info_number, url_decode_info_value,
};

pub fn load_and_infer_metadata(
    ani: &AniIndex,
    debug: bool,
) -> Result<HashMap<String, FieldNumber>> {
    let mut field_meta = load_field_metadata(ani, debug)?;

    if field_meta.is_empty() {
        if debug {
            eprintln!("[annotate] No metadata in ANI headers, inferring from data...");
        }

        let field_names: Vec<String> = ani
            .entries
            .iter()
            .take(10)
            .flat_map(|e| {
                let info_str = ani.read_cstring(e.info_ofs as usize);
                let decoded = url_decode_info_value(info_str.as_ref());
                decoded
                    .split(';')
                    .filter_map(|kv| kv.split('=').next().map(|k| k.to_string()))
                    .collect::<Vec<_>>()
            })
            .collect::<HashSet<_>>()
            .into_iter()
            .collect();

        field_meta = infer_field_metadata_from_data(ani, &field_names, debug);
    }

    Ok(field_meta)
}

fn load_field_metadata(ani: &AniIndex, debug: bool) -> Result<HashMap<String, FieldNumber>> {
    let mut metadata = HashMap::new();
    let headers = iter_ani_header_lines(ani);

    if debug {
        eprintln!("[DEBUG] ANI strings total length: {}", ani.strings_len());
        eprintln!(
            "[DEBUG] Found {} header lines in ANI strings",
            headers.len()
        );
    }

    for line in headers {
        if !line.starts_with("##INFO=") {
            continue;
        }

        if let Some(key) = extract_info_key(&line) {
            if let Some(number) = extract_info_number(&line) {
                metadata.insert(key.clone(), number);

                if debug {
                    eprintln!(
                        "[DEBUG] Loaded metadata from header: {} -> {:?}",
                        key, number
                    );
                }
            }
        }
    }

    if debug {
        eprintln!(
            "[DEBUG] Total metadata entries loaded from headers: {}",
            metadata.len()
        );
    }

    Ok(metadata)
}

pub fn iter_ani_header_lines(ani: &AniIndex) -> Vec<String> {
    let mut headers = Vec::new();
    let mut saw_header = false;
    let mut idx = 0usize;
    while idx < ani.strings_len() {
        let line_ref = ani.read_cstring(idx);
        let line = line_ref.as_ref();

        if line.is_empty() {
            idx += 1;
            continue;
        }

        if line == ANI_HEADER_END {
            break;
        }

        let is_header = line.starts_with("##INFO=")
            || line.starts_with("##FORMAT=")
            || line.starts_with("##FILTER=")
            || line.starts_with("#CHROM");

        if is_header {
            headers.push(line.to_string());
            saw_header = true;
        } else if saw_header {
            break;
        }

        idx += line.len() + 1;
    }

    headers
}

fn infer_field_metadata_from_data(
    ani: &AniIndex,
    field_names: &[String],
    debug: bool,
) -> HashMap<String, FieldNumber> {
    let field_set: HashSet<&str> = field_names.iter().map(|s| s.as_str()).collect();
    let mut candidates: HashMap<String, Vec<FieldNumber>> = HashMap::new();

    for e in ani.entries.iter().take(1000) {
        let info = ani.read_cstring(e.info_ofs as usize);
        let decoded = url_decode_info_value(info.as_ref());

        for kv in decoded.split(';') {
            if kv.is_empty() || !kv.contains('=') {
                continue;
            }

            let mut parts = kv.splitn(2, '=');
            let key = parts.next().unwrap();
            let value = match parts.next() {
                Some(v) => v,
                None => continue,
            };

            if !field_set.contains(key) {
                continue;
            }

            let value_parts: Vec<&str> = value.split(',').collect();
            let num_values = value_parts.len();

            let inferred = if num_values == 1 {
                FieldNumber::One
            } else {
                FieldNumber::Many
            };

            candidates
                .entry(key.to_string())
                .or_default()
                .push(inferred);
        }
    }

    let mut metadata = HashMap::new();
    for (key, numbers) in candidates {
        let best = choose_best_number(&numbers);
        metadata.insert(key.clone(), best);

        if debug {
            eprintln!(
                "[DEBUG] Inferred metadata: {} -> {:?} (from {} samples)",
                key,
                best,
                numbers.len()
            );
        }
    }

    metadata
}

pub fn is_missing_value(val: &str) -> bool {
    val.is_empty() || val == "."
}

/// Default-safe field type. Used only when there's no explicit `Type=` in
/// the embedded `.ani` header. bcftools requires an explicit type — guessing
/// from the key name mis-handles real VEP/SnpEff fields (IMPACT, FATHMM, …),
/// so we default to the only universally-safe type: String.
pub fn infer_field_type(_key: &str) -> &'static str {
    "String"
}

/// Parses `##FORMAT=<ID=...,Number=...,Type=...,Description=...>` header lines
/// from the embedded `.ani` headers, returning a `key → FieldNumber` map for
/// FORMAT-side annotations. Previously only `##INFO=` lines were parsed,
/// which forced FORMAT annotations through the (now-removed) letter
/// heuristic.
pub fn load_format_metadata(
    ani: &AniIndex,
    debug: bool,
) -> Result<HashMap<String, FieldNumber>> {
    let mut metadata = HashMap::new();
    for line in iter_ani_header_lines(ani) {
        if !line.starts_with("##FORMAT=") {
            continue;
        }
        if let (Some(key), Some(number)) =
            (extract_info_key(&line), extract_info_number(&line))
        {
            if debug {
                eprintln!("[DEBUG] FORMAT meta: {key} -> {number:?}");
            }
            metadata.insert(key, number);
        }
    }
    Ok(metadata)
}
