use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FieldNumber {
    Zero,
    One,
    Many,
    A,
    R,
    G,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FieldType {
    Int,
    Float,
    Str,
    Flag,
}

#[derive(Debug, Clone)]
pub struct StructuredInfoField {
    pub key: String,
    pub number: FieldNumber,
    pub ty: FieldType,
    pub values: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct AnnotationBundle {
    pub alt: String,
    pub id: Option<String>,
    pub qual: Option<String>,
    pub filter: Option<String>,
    pub info: Vec<StructuredInfoField>,
}

impl AnnotationBundle {
    pub fn match_alleles(&self, vcf_alt_alleles: &[&str]) -> Vec<Option<usize>> {
        let db_alts: Vec<&str> = self.alt.split(',').collect();

        let mut db_map: HashMap<&str, Vec<usize>> = HashMap::new();
        for (idx, alt) in db_alts.iter().enumerate() {
            db_map.entry(*alt).or_default().push(idx);
        }

        vcf_alt_alleles
            .iter()
            .map(|vcf_alt| {
                db_map
                    .get(vcf_alt)
                    .and_then(|indices| indices.first().copied())
            })
            .collect()
    }

    pub fn remap_field_values(
        &self,
        field: &StructuredInfoField,
        allele_map: &[Option<usize>],
    ) -> Vec<String> {
        match field.number {
            FieldNumber::A => allele_map
                .iter()
                .map(|opt_idx| {
                    opt_idx
                        .and_then(|idx| field.values.get(idx))
                        .cloned()
                        .unwrap_or_else(|| ".".to_string())
                })
                .collect(),
            FieldNumber::R => {
                let mut result = Vec::new();

                if let Some(ref_val) = field.values.first() {
                    result.push(ref_val.clone());
                } else {
                    result.push(".".to_string());
                }

                for opt_idx in allele_map {
                    let val = opt_idx
                        .and_then(|idx| field.values.get(idx + 1))
                        .cloned()
                        .unwrap_or_else(|| ".".to_string());
                    result.push(val);
                }

                result
            }
            _ => field.values.clone(),
        }
    }
}

pub fn parse_info_field(info_str: &str) -> Vec<StructuredInfoField> {
    let mut fields = Vec::new();

    if info_str == "." || info_str.is_empty() {
        return fields;
    }

    // URL-decode the info string first
    let decoded_info = url_decode_info_value(info_str);

    for pair in decoded_info.split(';') {
        if let Some(eq_pos) = pair.find('=') {
            let key = &pair[..eq_pos];
            let value = &pair[eq_pos + 1..];
            let values: Vec<String> = value.split(',').map(|s| s.to_string()).collect();

            let number = if values.len() == 1 {
                FieldNumber::One
            } else {
                FieldNumber::Many
            };

            let ty = infer_field_type(key);

            fields.push(StructuredInfoField {
                key: key.to_string(),
                number,
                ty,
                values,
            });
        } else {
            fields.push(StructuredInfoField {
                key: pair.to_string(),
                number: FieldNumber::Zero,
                ty: FieldType::Flag,
                values: vec![],
            });
        }
    }

    fields
}

pub fn infer_structured_info_fields(
    alt_alleles: &[&str],
    info_str: &str,
) -> Vec<StructuredInfoField> {
    let mut fields = Vec::new();

    if info_str == "." || info_str.is_empty() {
        return fields;
    }

    // URL-decode the info string first
    let decoded_info = url_decode_info_value(info_str);

    for pair in decoded_info.split(';') {
        if let Some(eq_pos) = pair.find('=') {
            let key = &pair[..eq_pos];
            let value = &pair[eq_pos + 1..];
            let values: Vec<String> = value.split(',').map(|s| s.to_string()).collect();

            let number = if values.len() == alt_alleles.len() {
                FieldNumber::A
            } else if values.len() == alt_alleles.len() + 1 {
                FieldNumber::R
            } else if values.len() == 1 {
                FieldNumber::One
            } else {
                FieldNumber::Many
            };

            let ty = infer_field_type(key);

            fields.push(StructuredInfoField {
                key: key.to_string(),
                number,
                ty,
                values,
            });
        } else {
            fields.push(StructuredInfoField {
                key: pair.to_string(),
                number: FieldNumber::Zero,
                ty: FieldType::Flag,
                values: vec![],
            });
        }
    }

    fields
}

fn url_decode_info_value(encoded: &str) -> String {
    encoded
        .replace("%3D", "=")
        .replace("%2C", ",")
        .replace("%3B", ";")
        .replace("%2F", "/")
        .replace("%20", " ")
        .replace("%25", "%")
}

fn infer_field_type(key: &str) -> FieldType {
    if key.ends_with("_AF") || key.ends_with("_FREQ") || key == "AF" || key == "FREQ" {
        return FieldType::Float;
    }

    if key.ends_with("_AC")
        || key.ends_with("_AN")
        || key.ends_with("_DP")
        || key == "AC"
        || key == "AN"
        || key == "DP"
        || key == "NS"
    {
        return FieldType::Int;
    }

    FieldType::Str
}
