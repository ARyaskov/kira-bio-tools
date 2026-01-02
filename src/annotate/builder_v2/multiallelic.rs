use std::collections::HashMap;

use crate::annotate::structs::bundle::FieldNumber;

pub fn split_info_for_allele(
    info: &str,
    alt_idx: usize,
    num_alts: usize,
    field_meta: &HashMap<String, FieldNumber>,
) -> String {
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok();
    let mut result = Vec::new();

    for field in info.split(';') {
        if field.is_empty() {
            continue;
        }

        if let Some(eq_pos) = field.find('=') {
            let key = &field[..eq_pos];
            let value = &field[eq_pos + 1..];

            let parts: Vec<&str> = value.split(',').collect();
            let number = field_meta.get(key).copied();

            match number {
                Some(FieldNumber::A) => {
                    if parts.len() == num_alts {
                        if let Some(v) = parts.get(alt_idx) {
                            result.push(format!("{}={}", key, v));
                        }
                    } else {
                        result.push(field.to_string());
                    }
                }
                Some(FieldNumber::R) => {
                    if parts.len() == num_alts + 1 {
                        if let Some(ref_val) = parts.first() {
                            if let Some(alt_val) = parts.get(alt_idx + 1) {
                                result.push(format!("{}={},{}", key, ref_val, alt_val));
                            }
                        }
                    } else {
                        result.push(field.to_string());
                    }
                }
                _ => {
                    result.push(field.to_string());
                }
            }
        } else {
            result.push(field.to_string());
        }
    }

    let final_result = result.join(";");

    if debug && result.len() > 0 {
        eprintln!(
            "[split_info] alt_idx={} num_alts={} -> '{}'",
            alt_idx,
            num_alts,
            if final_result.len() > 80 {
                &final_result[..80]
            } else {
                &final_result
            }
        );
    }

    final_result
}
