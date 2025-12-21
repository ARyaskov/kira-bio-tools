use std::collections::HashMap;

use super::column_spec::ColumnSpec;
use super::lookup::{lookup_annotations, update_id_and_filter};
use super::merge_info::{format_info_string, merge_info_fields};
use super::vcf_parsing::parse_vcf_record;
use crate::annotate::structs::ani::AniIndex;
use crate::annotate::structs::bundle::FieldNumber;

pub fn annotate_line(
    line: &str,
    ani: &AniIndex,
    field_meta: &HashMap<String, FieldNumber>,
    field_order: &[String],
    column_specs: &[ColumnSpec],
) -> String {
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok();

    if debug {
        eprintln!("[ANNOTATE] Processing line: {}", line);
        eprintln!(
            "[ANNOTATE] Column specs: {:?}",
            column_specs
                .iter()
                .map(|c| format!("{}{}", c.mode, c.key))
                .collect::<Vec<_>>()
        );
    }

    let Some(parsed) = parse_vcf_record(line) else {
        if debug {
            eprintln!("[ANNOTATE] Failed to parse VCF line");
        }
        return line.to_string();
    };

    if debug {
        eprintln!(
            "[ANNOTATE] Parsed: {}:{} {} -> {:?}",
            parsed.chrom, parsed.pos, parsed.ref_allele, parsed.vcf_alt_alleles
        );
        eprintln!("[ANNOTATE] Existing INFO: {:?}", parsed.info_map);
    }

    let lookup_result = lookup_annotations(
        ani,
        parsed.chrom,
        parsed.pos,
        parsed.ref_allele,
        &parsed.vcf_alt_alleles,
    );

    if debug {
        eprintln!(
            "[ANNOTATE] Lookup result: Found {} exact matches, multiallelic: {}",
            lookup_result.found_bundles.len(),
            lookup_result.multiallelic_bundle.is_some()
        );

        for (idx, bundle) in &lookup_result.found_bundles {
            eprintln!(
                "[ANNOTATE] Exact match [{}]: alt={}, info fields={}",
                idx,
                bundle.alt,
                bundle.info.len()
            );
        }

        if let Some(ref bundle) = lookup_result.multiallelic_bundle {
            eprintln!(
                "[ANNOTATE] Multiallelic: alt={}, info fields={}",
                bundle.alt,
                bundle.info.len()
            );
        }
    }

    let (updated_id, updated_filter) = update_id_and_filter(
        &lookup_result.found_bundles,
        &lookup_result.multiallelic_bundle,
        parsed.updated_id.clone(),
        parsed.updated_filter.clone(),
        &parsed.updated_id,
    );

    let column_specs_tuples: Vec<(String, super::super::structs::annotate_mode::AnnotateMode)> =
        column_specs
            .iter()
            .map(|cs| (cs.key.clone(), cs.mode))
            .collect();

    let existing_info = format_info_string(&parsed.info_map, &[]);
    if debug {
        eprintln!("[ANNOTATE] Existing info string: '{}'", existing_info);
    }

    let info_map = merge_info_fields(
        &existing_info,
        &lookup_result.found_bundles,
        &lookup_result.multiallelic_bundle,
        &parsed.vcf_alt_alleles,
        field_meta,
        &column_specs_tuples,
    );

    let final_info = format_info_string(&info_map, field_order);

    if debug {
        eprintln!("[ANNOTATE] Final INFO: '{}'", final_info);
    }

    format!(
        "{}\t{}\t{}\t{}\t{}\t.\t{}\t{}",
        parsed.chrom,
        parsed.pos,
        updated_id,
        parsed.ref_allele,
        parsed.vcf_alt_alleles.join(","),
        updated_filter,
        final_info
    )
}
