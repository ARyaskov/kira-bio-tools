use anyhow::Result;
use fxhash::FxHashMap;
use std::sync::atomic::{AtomicUsize, Ordering};

use super::entry_processing::{insert_or_update_entry, make_position_key, parse_chrom_and_pos};
use super::multiallelic::split_info_for_allele;
use crate::annotate::structs::ani::AniEntry;
use crate::util::{append_cstr, url_encode_info_value};
use crate::vcf::simd::SimdVcfParser;

pub fn process_vcf_line_multiallelic_simd(
    line: &[u8],
    entries_map: &mut FxHashMap<u64, (AniEntry, usize)>,
    pool: &mut Vec<u8>,
    insertion_order: &mut usize,
    duplicates_skipped: &AtomicUsize,
    multiallelic_count: &AtomicUsize,
    debug: bool,
) -> Result<usize> {
    let Some(parsed) = SimdVcfParser::parse_line(line) else {
        return Ok(0);
    };

    let (chr_id, pos) = match parse_chrom_and_pos(parsed.chrom, parsed.pos) {
        Some(v) => v,
        None => return Ok(0),
    };

    let alt_alleles: Vec<&str> = parsed.alt.split(',').collect();

    if alt_alleles.len() > 1 {
        multiallelic_count.fetch_add(1, Ordering::Relaxed);
    }

    let ref_ofs = append_cstr(pool, parsed.ref_allele.trim());
    let id_ofs = append_cstr(pool, parsed.id);
    let qual_ofs = append_cstr(pool, parsed.qual);
    let filter_ofs = append_cstr(pool, parsed.filter);

    for (alt_idx, alt_single) in alt_alleles.iter().enumerate() {
        let alt_single = alt_single.trim();
        let key = make_position_key(chr_id, pos, parsed.ref_allele.trim(), alt_single);

        let alt_ofs = append_cstr(pool, alt_single);

        let info_ofs = if !parsed.info.is_empty() && parsed.info != "." {
            let final_info = if alt_alleles.len() > 1 {
                split_info_for_allele(parsed.info, alt_idx, alt_alleles.len())
            } else {
                parsed.info.to_string()
            };
            let encoded = url_encode_info_value(&final_info);
            append_cstr(pool, &encoded) as u32
        } else {
            append_cstr(pool, ".") as u32
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
        };

        insert_or_update_entry(
            key,
            entry,
            entries_map,
            insertion_order,
            duplicates_skipped,
            debug,
            parsed.chrom,
            pos,
            parsed.ref_allele,
            alt_single,
        );
    }

    Ok(alt_alleles.len())
}
