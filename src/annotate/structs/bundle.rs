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
                    if let Some(idx) = opt_idx {
                        if let Some(val) = field.values.get(idx + 1) {
                            result.push(val.clone());
                        } else {
                            result.push(".".to_string());
                        }
                    } else {
                        result.push(".".to_string());
                    }
                }

                result
            }
            _ => field.values.clone(),
        }
    }
}

pub fn parse_info_field(info: &str) -> Vec<StructuredInfoField> {
    use crate::util::url_decode_info_value;

    let decoded_info = url_decode_info_value(info);
    let mut fields = Vec::new();

    for kv in decoded_info.split(';') {
        if kv.is_empty() || kv == "." {
            continue;
        }

        if let Some(eq_pos) = kv.find('=') {
            let key = &kv[..eq_pos];
            let value_part = &kv[eq_pos + 1..];

            let values: Vec<String> = value_part.split(',').map(|s| s.to_string()).collect();

            let ty = infer_field_type(key);

            let number = match values.len() {
                0 => FieldNumber::Zero,
                1 => FieldNumber::One,
                _ => FieldNumber::Many,
            };

            fields.push(StructuredInfoField {
                key: key.to_string(),
                number,
                ty,
                values,
            });
        } else {
            fields.push(StructuredInfoField {
                key: kv.to_string(),
                number: FieldNumber::Zero,
                ty: FieldType::Flag,
                values: vec![],
            });
        }
    }

    fields
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
