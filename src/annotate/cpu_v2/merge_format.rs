use std::collections::HashMap;
use crate::annotate::cpu_v2::vcf_parsing::{ParsedFormat, ParsedSample};
use crate::annotate::structs::bundle::{AnnotationBundle, FieldNumber};

pub fn merge_format_and_samples(
    input_format: Option<&ParsedFormat>,
    input_samples: &[ParsedSample],
    _bundles: &[(usize, AnnotationBundle)],
    _allele_map: &[Option<usize>],
    _field_meta: &HashMap<String, FieldNumber>,
) -> (Option<String>, Vec<String>) {
    if input_samples.is_empty() {
        return (None, Vec::new());
    }

    let format_str = input_format.map(|f| f.keys.join(":"));
    let sample_strs: Vec<String> = input_samples.iter()
        .map(|s| s.raw.join(":"))
        .collect();

    (format_str, sample_strs)
}