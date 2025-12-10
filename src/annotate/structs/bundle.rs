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

pub fn encode_structured_info(f: &[StructuredInfoField]) -> String {
    let mut out = String::new();

    for (i, fld) in f.iter().enumerate() {
        if i > 0 {
            out.push(';');
        }

        match fld.number {
            FieldNumber::Zero => {
                out.push_str(fld.key);
            }

            FieldNumber::One => {
                out.push_str(fld.key);
                out.push('=');
                if let Some(v) = fld.values.first() {
                    out.push_str(v);
                }
            }

            FieldNumber::Many | FieldNumber::A | FieldNumber::R | FieldNumber::G => {
                out.push_str(fld.key);
                out.push('=');
                out.push_str(&fld.values.join(","));
            }
        }
    }

    out
}
