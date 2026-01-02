use anyhow::Result;
use fxhash::FxHashMap;
use std::sync::atomic::{AtomicUsize, Ordering};

use super::entry_processing::{insert_or_update_entry, make_position_key};
use super::multiallelic::split_info_for_allele;
use crate::annotate::builder_v2::StringPool;
use crate::annotate::structs::ani::{AniEntry, ANI_STR_NONE};
use crate::annotate::structs::tab::TabSchema;
use crate::util::{chr_name_to_id, url_encode_info_value};

pub fn process_tab_line_multiallelic(
    line: &str,
    schema: &TabSchema,
    entries_map: &mut FxHashMap<u64, (AniEntry, usize)>,
    pool: &mut StringPool,
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

    let ref_ofs = pool.append_cstr(rf);
    let id_ofs = pool.append_cstr(id);
    let qual_ofs = pool.append_cstr(qual);
    let filter_ofs = pool.append_cstr(filt);

    let mut base_info_str = String::new();

    if let Some(info_idx) = schema.info_start {
        if let Some(info_str) = parts.get(info_idx) {
            let info_str = info_str.trim();
            if !info_str.is_empty() && info_str != "." {
                base_info_str.push_str(info_str);
            }
        }
    }

    for col in &schema.info_cols {
        if let Some(val_str) = parts.get(col.index) {
            let val_str = val_str.trim();
            if !val_str.is_empty() && val_str != "." {
                if !base_info_str.is_empty() {
                    base_info_str.push(';');
                }

                base_info_str.push_str(&col.key);
                base_info_str.push('=');
                base_info_str.push_str(val_str);
            }
        }
    }

    for (alt_idx, alt_single) in alt_alleles.iter().enumerate() {
        let key = make_position_key(chr_id, pos, rf, alt_single);

        if debug {
            eprintln!(
                "[TAB-INDEX] {}:{} {}>{} alt_idx={} key={:016x}",
                chr, pos, rf, alt_single, alt_idx, key
            );
        }

        let alt_ofs = pool.append_cstr(alt_single);

        let info_ofs = if !base_info_str.is_empty() && base_info_str != "." {
            let final_info = if alt_alleles.len() > 1 {
                split_info_for_allele(
                    &base_info_str,
                    alt_idx,
                    alt_alleles.len(),
                    &schema.field_metadata,
                )
            } else {
                base_info_str.clone()
            };
            if debug {
                eprintln!(
                    "[TAB-INDEX] final_info for alt_idx={}: '{}'",
                    alt_idx,
                    if final_info.len() > 100 {
                        &final_info[..100]
                    } else {
                        &final_info
                    }
                );
            }
            let encoded = url_encode_info_value(&final_info);
            pool.append_cstr(&encoded) as u32
        } else {
            pool.append_cstr(".") as u32
        };

        let entry = AniEntry {
            chr_id,
            pos,
            ref_ofs: ref_ofs as u32,
            alt_ofs: alt_ofs as u32,
            id_ofs: id_ofs as u32,
            qual_ofs: qual_ofs as u32,
            filter_ofs: filter_ofs as u32,
            info_ofs,
            info_len: 0,
            format_ofs: ANI_STR_NONE,
            samples_ofs: ANI_STR_NONE,
        };

        insert_or_update_entry(
            key,
            entry,
            entries_map,
            insertion_order,
            duplicates_skipped,
            debug,
            chr,
            pos,
            rf,
            alt_single,
        );
    }

    Ok(alt_alleles.len())
}
