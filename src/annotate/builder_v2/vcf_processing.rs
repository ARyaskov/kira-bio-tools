use anyhow::Result;
use fxhash::FxHashMap;
use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};

use super::entry_processing::{insert_or_update_entry, make_position_key, parse_chrom_and_pos};
use super::multiallelic::split_info_for_allele;
use crate::annotate::builder_v2::StringPool;
use crate::annotate::structs::ani::ANI_STR_NONE;
use crate::annotate::structs::ani::AniEntry;
use crate::annotate::structs::bundle::FieldNumber;
use crate::util::url_encode_info_value;
use crate::vcf::simd::SimdVcfParser;

pub fn process_vcf_line_multiallelic_simd(
    line: &[u8],
    entries_map: &mut FxHashMap<u64, (AniEntry, usize)>,
    pool: &mut StringPool,
    insertion_order: &mut usize,
    duplicates_skipped: &AtomicUsize,
    multiallelic_count: &AtomicUsize,
    field_meta: &HashMap<String, FieldNumber>,
    expected_sample_count: usize,
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

    let (format_ofs, samples_ofs) = if samples.is_empty() {
        (ANI_STR_NONE, ANI_STR_NONE)
    } else {
        let fmt_ofs = match parsed.format {
            Some(fmt) if fmt != "." => pool.append_cstr(fmt) as u32,
            _ => ANI_STR_NONE,
        };
        let joined = samples.join("\t");
        let samp_ofs = pool.append_cstr(&joined) as u32;
        (fmt_ofs, samp_ofs)
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
            debug,
            parsed.chrom,
            pos,
            parsed.ref_allele,
            alt_single,
        );
    }

    Ok(alt_alleles.len())
}

#[cfg(test)]
mod tests {
    use super::process_vcf_line_multiallelic_simd;
    use crate::annotate::builder_v2::StringPool;
    use crate::annotate::structs::ani::ANI_STR_NONE;
    use crate::annotate::structs::bundle::FieldNumber;
    use crate::util::read_cstring;
    use fxhash::FxHashMap;
    use std::collections::HashMap;
    use std::sync::atomic::AtomicUsize;

    #[test]
    fn test_samples_are_not_truncated() {
        let line = b"1\t3000002\tid\tC\tT\t99\tq99\tFLAG;IINT=88,99;IFLT=8.8,9.9;ISTR=888,999\tGT:FINT:FFLT:FSTR\t1|1:88,99:8.8,9.9:888,999\t0|1:77:7.7:77";
        let mut entries_map: FxHashMap<u64, (crate::annotate::structs::ani::AniEntry, usize)> =
            FxHashMap::default();
        let mut pool = StringPool::new();
        let mut insertion_order = 0usize;
        let duplicates_skipped = AtomicUsize::new(0);
        let multiallelic_count = AtomicUsize::new(0);
        let field_meta: HashMap<String, FieldNumber> = HashMap::new();

        let processed = process_vcf_line_multiallelic_simd(
            line,
            &mut entries_map,
            &mut pool,
            &mut insertion_order,
            &duplicates_skipped,
            &multiallelic_count,
            &field_meta,
            2,
            false,
        )
        .unwrap();

        assert_eq!(processed, 1);
        let entry = entries_map.values().next().unwrap().0;
        assert_ne!(entry.samples_ofs, ANI_STR_NONE);
        let pool_bytes = pool.materialize().unwrap();
        let samples = read_cstring(&pool_bytes, entry.samples_ofs as usize);
        assert_eq!(samples, "1|1:88,99:8.8,9.9:888,999\t0|1:77:7.7:77");
    }
}
