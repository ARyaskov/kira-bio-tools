use anyhow::Result;
use fxhash::FxHashMap;
use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};

use super::entry_processing::{EntryEntry, insert_or_update_entry, make_position_key};
use super::multiallelic::split_info_for_allele;
use crate::annotate::builder_v2::StringPool;
use crate::annotate::structs::ani::ANI_STR_NONE;
use crate::annotate::structs::ani::AniEntry;
use crate::annotate::structs::ani::ContigDict;
use crate::annotate::structs::bundle::FieldNumber;
use crate::util::url_encode_info_value;
use crate::vcf::simd::SimdVcfParser;

pub fn process_vcf_line_multiallelic_simd(
    line: &[u8],
    contigs: &mut ContigDict,
    entries_map: &mut FxHashMap<u64, EntryEntry>,
    pool: &mut StringPool,
    insertion_order: &mut usize,
    duplicates_skipped: &AtomicUsize,
    collisions_detected: &AtomicUsize,
    multiallelic_count: &AtomicUsize,
    field_meta: &HashMap<String, FieldNumber>,
    format_meta: &HashMap<String, FieldNumber>,
    expected_sample_count: usize,
    debug: bool,
) -> Result<usize> {
    let Some(parsed) = SimdVcfParser::parse_line(line) else {
        return Ok(0);
    };

    let Ok(pos) = parsed.pos.parse::<u32>() else {
        return Ok(0);
    };
    // Insert-on-first-seen: source VCF body wins the contig id assignment.
    // (When a `##contig=` block was parsed up-front into `contigs`, this is a
    // pure dict lookup with no insertion.)
    let chr_id = contigs.insert(parsed.chrom);

    let alt_alleles: Vec<&str> = parsed.alt.split(',').collect();

    if alt_alleles.len() > 1 {
        multiallelic_count.fetch_add(1, Ordering::Relaxed);
    }

    let ref_ofs = pool.append_cstr(parsed.ref_allele.trim());
    let id_ofs = pool.append_cstr(parsed.id);
    let qual_ofs = pool.append_cstr(parsed.qual);
    let filter_ofs = pool.append_cstr(parsed.filter);

    let mut samples = parsed.samples.clone();
    if expected_sample_count > 0
        && (samples.len() != expected_sample_count || samples.iter().any(|s| s.is_empty()))
    {
        if let Ok(line_str) = std::str::from_utf8(line) {
            let cols: Vec<&str> = line_str.split('\t').collect();
            if cols.len() > 9 {
                samples = cols[9..].to_vec();
            } else {
                samples.clear();
            }
        }
    }

    let format_ofs = if samples.is_empty() {
        ANI_STR_NONE
    } else {
        match parsed.format {
            Some(fmt) if fmt != "." => pool.append_cstr(fmt) as u32,
            _ => ANI_STR_NONE,
        }
    };

    for (alt_idx, alt_single) in alt_alleles.iter().enumerate() {
        let alt_single = alt_single.trim();
        if debug && pos <= 5 {
            eprintln!(
                "[VCF-INDEX] {}:{} {}>{} info='{}' alts={}",
                parsed.chrom, pos, parsed.ref_allele, alt_single, parsed.info, parsed.alt
            );
        }
        let key = make_position_key(chr_id, pos, parsed.ref_allele.trim(), alt_single);

        let alt_ofs = pool.append_cstr(alt_single);
        let samples_ofs = if samples.is_empty() {
            ANI_STR_NONE
        } else {
            let entry_samples = if alt_alleles.len() > 1 {
                split_format_samples_for_allele(
                    &samples,
                    parsed.format,
                    alt_idx,
                    alt_alleles.len(),
                    format_meta,
                )
            } else {
                samples.iter().map(|s| (*s).to_string()).collect()
            };
            let joined = entry_samples.join("\t");
            pool.append_cstr(&joined) as u32
        };

        let (info_ofs, info_len) = if !parsed.info.is_empty() && parsed.info != "." {
            let final_info = if alt_alleles.len() > 1 {
                split_info_for_allele(parsed.info, alt_idx, alt_alleles.len(), field_meta)
            } else {
                parsed.info.to_string()
            };
            let encoded = url_encode_info_value(&final_info);
            (pool.append_cstr(&encoded) as u32, encoded.len() as u32)
        } else {
            (pool.append_cstr(".") as u32, 1)
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
            info_len,
            format_ofs,
            samples_ofs,
        };

        insert_or_update_entry(
            key,
            entry,
            entries_map,
            insertion_order,
            duplicates_skipped,
            collisions_detected,
            debug,
            parsed.chrom,
            pos,
            parsed.ref_allele,
            alt_single,
        );
    }

    Ok(alt_alleles.len())
}

fn split_format_samples_for_allele(
    samples: &[&str],
    format: Option<&str>,
    alt_idx: usize,
    alt_count: usize,
    format_meta: &HashMap<String, FieldNumber>,
) -> Vec<String> {
    let Some(format) = format else {
        return samples.iter().map(|s| (*s).to_string()).collect();
    };
    let keys: Vec<&str> = format.split(':').collect();
    if keys.is_empty() || format_meta.is_empty() {
        return samples.iter().map(|s| (*s).to_string()).collect();
    }
    let numbers: Vec<Option<FieldNumber>> = keys
        .iter()
        .map(|k| format_meta.get(*k).copied())
        .collect();
    if !numbers
        .iter()
        .any(|n| matches!(n, Some(FieldNumber::A | FieldNumber::R | FieldNumber::G)))
    {
        return samples.iter().map(|s| (*s).to_string()).collect();
    }
    samples
        .iter()
        .map(|sample| {
            let vals: Vec<&str> = sample.split(':').collect();
            let mut out = Vec::with_capacity(keys.len());
            for (idx, number) in numbers.iter().enumerate() {
                let raw = vals.get(idx).copied().unwrap_or(".");
                out.push(split_format_value_for_allele(
                    raw,
                    *number,
                    alt_idx,
                    alt_count,
                ));
            }
            out.join(":")
        })
        .collect()
}

fn split_format_value_for_allele(
    raw: &str,
    number: Option<FieldNumber>,
    alt_idx: usize,
    alt_count: usize,
) -> String {
    if raw.is_empty() || raw == "." {
        return raw.to_string();
    }
    let Some(number) = number else {
        return raw.to_string();
    };
    match number {
        FieldNumber::A => raw
            .split(',')
            .nth(alt_idx)
            .unwrap_or(".")
            .to_string(),
        FieldNumber::R => {
            let values: Vec<&str> = raw.split(',').collect();
            let ref_val = values.first().copied().unwrap_or(".");
            let alt_val = values.get(alt_idx + 1).copied().unwrap_or(".");
            format!("{ref_val},{alt_val}")
        }
        FieldNumber::G => {
            let values: Vec<&str> = raw.split(',').collect();
            let allele = alt_idx + 1;
            let idx00 = genotype_index(0, 0);
            let idx01 = genotype_index(0, allele);
            let idx11 = genotype_index(allele, allele);
            if values.len() != (alt_count + 1) * (alt_count + 2) / 2 {
                return raw.to_string();
            }
            format!(
                "{},{},{}",
                values.get(idx00).copied().unwrap_or("."),
                values.get(idx01).copied().unwrap_or("."),
                values.get(idx11).copied().unwrap_or(".")
            )
        }
        _ => raw.to_string(),
    }
}

fn genotype_index(a: usize, b: usize) -> usize {
    let lo = a.min(b);
    let hi = a.max(b);
    hi * (hi + 1) / 2 + lo
}

#[cfg(test)]
#[path = "../../../tests/unit/annotate_builder_v2_vcf_processing.rs"]
mod tests;
