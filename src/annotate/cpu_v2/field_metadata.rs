use anyhow::Result;
use std::collections::{HashMap, HashSet};

use crate::annotate::structs::ani::AniIndex;
use crate::annotate::structs::bundle::FieldNumber;
use crate::util::{
    choose_best_number, extract_info_key, extract_info_number, read_cstring, url_decode_info_value,
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
                let info_str = read_cstring(&ani.strings, e.info_ofs as usize);
                let decoded = url_decode_info_value(info_str);
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
    let strings_str = std::str::from_utf8(&ani.strings).unwrap_or("");

    if debug {
        eprintln!("[DEBUG] ANI strings total length: {}", ani.strings.len());
        let preview_len = strings_str.len().min(500);
        eprintln!(
            "[DEBUG] First {} chars of strings: {:?}",
            preview_len,
            &strings_str[..preview_len]
        );
        let header_count = strings_str
            .split('\0')
            .filter(|s| s.starts_with("##INFO="))
            .count();
        eprintln!(
            "[DEBUG] Found {} ##INFO headers in ANI strings",
            header_count
        );
    }

    for line in strings_str.split('\0') {
        if !line.starts_with("##INFO=") {
            continue;
        }

        if let Some(key) = extract_info_key(line) {
            if let Some(number) = extract_info_number(line) {
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

fn infer_field_metadata_from_data(
    ani: &AniIndex,
    field_names: &[String],
    debug: bool,
) -> HashMap<String, FieldNumber> {
    let mut candidates: HashMap<String, Vec<FieldNumber>> = HashMap::new();

    for entry in ani.entries.iter().take(1000) {
        let info_str = read_cstring(&ani.strings, entry.info_ofs as usize);
        let decoded = url_decode_info_value(info_str);

        for kv in decoded.split(';') {
            if kv.is_empty() || kv == "." {
                continue;
            }

            let mut parts = kv.splitn(2, '=');
            let key = parts.next().unwrap();
            let value = match parts.next() {
                Some(v) => v,
                None => continue,
            };

            if !field_names.contains(&key.to_string()) {
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
                .or_insert_with(Vec::new)
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

pub fn infer_field_type(key: &str) -> &'static str {
    if key.starts_with('I') || key.ends_with("INT") {
        return "Integer";
    }
    if key.starts_with('F') || key.ends_with("FLT") || key.ends_with("FLOAT") {
        return "Float";
    }
    if key.starts_with('S') || key.ends_with("STR") || key.ends_with("STRING") {
        return "String";
    }
    "String"
}
