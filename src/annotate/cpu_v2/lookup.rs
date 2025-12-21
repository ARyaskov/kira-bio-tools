use std::collections::HashMap;

use crate::annotate::structs::ani::AniIndex;
use crate::annotate::AnnotationBundle;

pub struct LookupResult {
    pub found_bundles: Vec<(usize, AnnotationBundle)>,
    pub multiallelic_bundle: Option<AnnotationBundle>,
}

pub fn lookup_annotations<'a>(
    ani: &AniIndex,
    chrom: &str,
    pos: u32,
    ref_allele: &str,
    vcf_alt_alleles: &[&'a str],
) -> LookupResult {
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok();
    let mut found_bundles: Vec<(usize, AnnotationBundle)> = Vec::new();

    if debug {
        eprintln!("[LOOKUP] ===== Starting lookup =====");
        eprintln!("[LOOKUP] Position: {}:{}", chrom, pos);
        eprintln!("[LOOKUP] REF: '{}'", ref_allele);
        eprintln!("[LOOKUP] VCF ALTs: {:?}", vcf_alt_alleles);
        eprintln!("[LOOKUP] Total entries in ANI: {}", ani.entries.len());

        // Debug first few entries
        for (i, entry) in ani.entries.iter().take(10).enumerate() {
            eprintln!(
                "[LOOKUP] Entry[{}]: chr_id={}, pos={}, ref_ofs={}, alt_ofs={}, info_ofs={}",
                i, entry.chr_id, entry.pos, entry.ref_ofs, entry.alt_ofs, entry.info_ofs
            );
        }
    }

    // Try exact lookup for each VCF allele
    for (vcf_idx, vcf_alt) in vcf_alt_alleles.iter().enumerate() {
        let vcf_alt_trimmed = vcf_alt.trim();

        if debug {
            eprintln!(
                "[LOOKUP] Trying exact lookup for allele[{}]: '{}'",
                vcf_idx, vcf_alt_trimmed
            );
        }

        if let Some(bundle) = ani.lookup(chrom, pos, ref_allele.trim(), vcf_alt_trimmed) {
            if debug {
                eprintln!(
                    "[LOOKUP] ✓ Found exact match for allele[{}]: {}",
                    vcf_idx, vcf_alt_trimmed
                );
                eprintln!("[LOOKUP]   Bundle INFO fields: {}", bundle.info.len());
                for field in &bundle.info {
                    eprintln!("[LOOKUP]     Field {}: {:?}", field.key, field.values);
                }
            }
            found_bundles.push((vcf_idx, bundle));
        } else if debug {
            eprintln!(
                "[LOOKUP] ✗ No exact match for allele[{}]: '{}'",
                vcf_idx, vcf_alt_trimmed
            );
        }
    }

    // No need for multiallelic_bundle since we're doing exact lookups
    let multiallelic_bundle = None;

    if debug {
        eprintln!("[LOOKUP] ===== Lookup summary =====");
        eprintln!("[LOOKUP] Exact matches: {}", found_bundles.len());
        eprintln!("[LOOKUP] Multiallelic bundle: false");
        eprintln!("[LOOKUP] =============================");
    }

    LookupResult {
        found_bundles,
        multiallelic_bundle,
    }
}

pub fn update_id_and_filter(
    found_bundles: &[(usize, AnnotationBundle)],
    multiallelic_bundle: &Option<AnnotationBundle>,
    mut updated_id: String,
    mut updated_filter: String,
    original_id: &str,
) -> (String, String) {
    for (_vcf_idx, bundle) in found_bundles {
        if let Some(id) = &bundle.id {
            if updated_id == "." || updated_id.is_empty() {
                updated_id = id.clone();
            } else if !original_id.contains(id.as_str()) {
                updated_id = format!("{};{}", updated_id, id);
            }
        }

        if let Some(filt) = &bundle.filter {
            if updated_filter == "." || updated_filter.is_empty() {
                updated_filter = filt.clone();
            }
        }
    }

    if let Some(ref bundle) = multiallelic_bundle {
        if let Some(id) = &bundle.id {
            if updated_id == "." || updated_id.is_empty() {
                updated_id = id.clone();
            }
        }
        if let Some(filt) = &bundle.filter {
            if updated_filter == "." || updated_filter.is_empty() {
                updated_filter = filt.clone();
            }
        }
    }

    (updated_id, updated_filter)
}

pub fn build_alt_mapping(alt_string: &str) -> HashMap<&str, usize> {
    let db_alts: Vec<&str> = alt_string.split(',').collect();
    let mut alt_to_db_idx: HashMap<&str, usize> = HashMap::new();
    for (i, alt) in db_alts.iter().enumerate() {
        alt_to_db_idx.insert(*alt, i);
    }
    alt_to_db_idx
}
