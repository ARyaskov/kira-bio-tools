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
pub struct StructuredInfoField<'a> {
    pub key: &'a str,
    pub number: FieldNumber,
    pub ty: FieldType,
    pub values: Vec<&'a str>,
}

#[derive(Debug)]
pub struct AnnotationBundle<'a> {
    pub alt: &'a str,
    pub id: Option<&'a str>,
    pub qual: Option<&'a str>,
    pub filter: Option<&'a str>,
    pub info: Vec<StructuredInfoField<'a>>,
}

impl<'a> AnnotationBundle<'a> {
    pub fn match_alleles(&self, vcf_alt_alleles: &[&str]) -> Vec<Option<usize>> {
        let db_alts: Vec<&str> = self.alt.split(',').collect();

        let mut db_map: HashMap<&str, Vec<usize>> = HashMap::new();
        for (idx, alt) in db_alts.iter().enumerate() {
            db_map.entry(*alt).or_insert_with(Vec::new).push(idx);
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
        field: &StructuredInfoField<'a>,
        allele_map: &[Option<usize>],
    ) -> Vec<String> {
        match field.number {
            FieldNumber::A => allele_map
                .iter()
                .map(|opt_idx| {
                    opt_idx
                        .and_then(|idx| field.values.get(idx))
                        .map(|v| v.to_string())
                        .unwrap_or_else(|| ".".to_string())
                })
                .collect(),
            FieldNumber::R => {
                let mut result = Vec::new();

                if let Some(ref_val) = field.values.first() {
                    result.push(ref_val.to_string());
                } else {
                    result.push(".".to_string());
                }

                for opt_idx in allele_map {
                    if let Some(idx) = opt_idx {
                        let alt_idx = idx + 1;
                        if let Some(val) = field.values.get(alt_idx) {
                            result.push(val.to_string());
                        } else {
                            result.push(".".to_string());
                        }
                    } else {
                        result.push(".".to_string());
                    }
                }

                result
            }
            _ => field.values.iter().map(|v| v.to_string()).collect(),
        }
    }
}

pub fn parse_info_field<'a>(raw: &'a str) -> Vec<StructuredInfoField<'a>> {
    if raw == "." || raw.is_empty() {
        return vec![];
    }

    let mut out = Vec::new();

    for part in raw.split(';') {
        if part.is_empty() {
            continue;
        }

        let mut kv = part.splitn(2, '=');
        let key = kv.next().unwrap();
        let vals = kv.next().unwrap_or("");

        let number = if vals.is_empty() {
            FieldNumber::Zero
        } else if !vals.contains(',') {
            FieldNumber::One
        } else {
            FieldNumber::Many
        };

        let ty = if number == FieldNumber::Zero {
            FieldType::Flag
        } else if vals.parse::<i64>().is_ok() {
            FieldType::Int
        } else if vals.parse::<f64>().is_ok() {
            FieldType::Float
        } else {
            FieldType::Str
        };

        let values: Vec<&str> = if number == FieldNumber::Zero {
            vec![]
        } else {
            vals.split(',').collect()
        };

        out.push(StructuredInfoField {
            key,
            number,
            ty,
            values,
        });
    }

    out
}

pub fn infer_field_number(header_key: &str, value_count: usize, n_alt: usize) -> FieldNumber {
    if value_count == 0 {
        FieldNumber::Zero
    } else if value_count == 1 {
        FieldNumber::One
    } else if value_count == n_alt {
        FieldNumber::A
    } else if value_count == n_alt + 1 {
        FieldNumber::R
    } else {
        FieldNumber::Many
    }
}
