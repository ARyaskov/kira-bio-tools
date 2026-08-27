use std::collections::HashMap;
use url::form_urlencoded;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FieldNumber {
    Zero,
    One,
    A,
    R,
    G,
    Many,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FieldType {
    Integer,
    Float,
    String,
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
    pub format_str: Option<String>,
    pub format_samples: Vec<String>,
    /// REF allele as recorded in the source `.ani` for this entry. Needed by
    /// the bcftools-style vcmp allele matcher in `cpu_v2::vcmp` so that
    /// differently-padded indels between source and target VCFs still
    /// resolve to the same biological event. `String::new()` for the rare
    /// legacy code paths that don't carry it.
    pub db_ref: String,
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
            .map(|vcf_alt| db_map.get(vcf_alt).and_then(|v| v.first().copied()))
            .collect()
    }

    pub fn remap_field_values(
        &self,
        field: &StructuredInfoField,
        allele_map: &[Option<usize>],
    ) -> Vec<String> {
        match field.number {
            FieldNumber::A => remap_a(&field.values, allele_map),
            FieldNumber::R => remap_r(&field.values, allele_map),
            FieldNumber::G => remap_g_diploid(&field.values, allele_map),
            _ => field.values.clone(),
        }
    }
}

pub fn parse_info_field(info_str: &str) -> Vec<StructuredInfoField> {
    infer_structured_info_fields(info_str, &HashMap::new())
}

pub fn infer_structured_info_fields(
    info_str: &str,
    field_meta: &HashMap<String, FieldNumber>,
) -> Vec<StructuredInfoField> {
    let mut fields = Vec::new();

    if info_str == "." || info_str.is_empty() {
        return fields;
    }

    let decoded_info = url_decode_info_value(info_str);

    for pair in decoded_info.split(';') {
        if pair.is_empty() {
            continue;
        }

        if let Some((key, value)) = pair.split_once('=') {
            let values: Vec<String> = value.split(',').map(|s| s.to_string()).collect();

            let number = field_meta.get(key).copied().unwrap_or_else(|| {
                if values.len() == 1 {
                    FieldNumber::One
                } else {
                    FieldNumber::Many
                }
            });

            let ty = infer_field_type(key);

            fields.push(StructuredInfoField {
                key: key.to_string(),
                number,
                ty,
                values,
            });
        } else {
            let key = pair;
            let number = field_meta.get(key).copied().unwrap_or(FieldNumber::Zero);

            fields.push(StructuredInfoField {
                key: key.to_string(),
                number,
                ty: FieldType::Flag,
                values: vec![],
            });
        }
    }

    fields
}

fn url_decode_info_value(s: &str) -> String {
    if s.contains('%') {
        let decoded: String = form_urlencoded::parse(s.as_bytes())
            .map(|(k, v)| {
                if v.is_empty() {
                    k.into_owned()
                } else {
                    format!("{}={}", k, v)
                }
            })
            .collect::<Vec<_>>()
            .join("&");
        decoded.replace('&', ";")
    } else {
        s.to_string()
    }
}

fn infer_field_type(_key: &str) -> FieldType {
    // bcftools never guesses INFO/FORMAT type from the key name — it requires
    // an explicit `##INFO=<…,Type=…>` header. The legacy heuristic
    // ("starts-with I → Integer, F → Float") miscategorised real SnpEff/VEP
    // keys like IMPACT (String), FATHMM (String), FUNSEQ (String) and
    // silently corrupted output. Default to String; proper types are loaded
    // from the embedded `.ani` header by `load_field_metadata` in
    // cpu_v2::field_metadata.
    FieldType::String
}

fn remap_a(values: &[String], allele_map: &[Option<usize>]) -> Vec<String> {
    allele_map
        .iter()
        .map(|opt_idx| match opt_idx {
            Some(db_idx) => values
                .get(*db_idx)
                .cloned()
                .unwrap_or_else(|| ".".to_string()),
            None => ".".to_string(),
        })
        .collect()
}

fn remap_r(values: &[String], allele_map: &[Option<usize>]) -> Vec<String> {
    let mut result = Vec::new();

    result.push(values.first().cloned().unwrap_or_else(|| ".".to_string()));

    for opt_idx in allele_map {
        match opt_idx {
            Some(db_idx) => {
                let alt_val = values
                    .get(*db_idx + 1)
                    .cloned()
                    .unwrap_or_else(|| ".".to_string());
                result.push(alt_val);
            }
            None => result.push(".".to_string()),
        }
    }

    result
}

fn remap_g_diploid(values: &[String], allele_map: &[Option<usize>]) -> Vec<String> {
    let mut result = Vec::new();
    let max_alt = allele_map.len();

    for i in 0..=max_alt {
        for j in i..=max_alt {
            let db_i = if i == 0 {
                Some(0)
            } else {
                allele_map.get(i - 1).copied().flatten().map(|x| x + 1)
            };
            let db_j = if j == 0 {
                Some(0)
            } else {
                allele_map.get(j - 1).copied().flatten().map(|x| x + 1)
            };

            match (db_i, db_j) {
                (Some(a), Some(b)) => {
                    let g_idx = if a <= b {
                        b * (b + 1) / 2 + a
                    } else {
                        a * (a + 1) / 2 + b
                    };
                    result.push(
                        values
                            .get(g_idx)
                            .cloned()
                            .unwrap_or_else(|| ".".to_string()),
                    );
                }
                _ => result.push(".".to_string()),
            }
        }
    }

    result
}
