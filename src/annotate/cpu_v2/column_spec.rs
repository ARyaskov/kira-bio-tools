use super::super::structs::annotate_mode::AnnotateMode;

#[derive(Debug, Clone)]
pub struct ColumnSpec {
    pub key: String,
    pub dst_key: String,
    pub mode: AnnotateMode,
}

impl ColumnSpec {
    fn canonical_ref(raw: &str) -> String {
        if let Some(rest) = raw.strip_prefix("FORMAT/") {
            return format!("FMT/{rest}");
        }
        raw.to_string()
    }

    fn is_fixed_column(raw: &str) -> bool {
        raw.eq_ignore_ascii_case("ID")
            || raw.eq_ignore_ascii_case("QUAL")
            || raw.eq_ignore_ascii_case("FILTER")
            || raw.eq_ignore_ascii_case("INFO")
            || raw.eq_ignore_ascii_case("FMT")
            || raw.eq_ignore_ascii_case("FORMAT")
            || raw.eq_ignore_ascii_case("ALT")
    }

    pub fn parse(spec: &str) -> Self {
        let (mode, rest) = AnnotateMode::parse(spec);

        let (src_key, dst_key) = if rest.contains(":=") {
            let parts: Vec<&str> = rest.splitn(2, ":=").collect();
            if parts.len() == 2 {
                let dst = Self::canonical_ref(parts[0]);
                let mut src = Self::canonical_ref(parts[1]);
                if !src.contains('/')
                    && !Self::is_fixed_column(&src)
                    && dst.to_ascii_uppercase().starts_with("FMT/")
                {
                    src = format!("FMT/{src}");
                } else if !src.contains('/')
                    && !Self::is_fixed_column(&src)
                    && dst.to_ascii_uppercase().starts_with("INFO/")
                {
                    src = format!("INFO/{src}");
                }
                (src, dst)
            } else {
                let key = Self::canonical_ref(rest);
                (key.clone(), key)
            }
        } else {
            let key = Self::canonical_ref(rest);
            (key.clone(), key)
        };

        let runtime_key = if src_key == dst_key {
            src_key.clone()
        } else {
            format!("{src_key}=>{dst_key}")
        };

        Self {
            key: runtime_key,
            dst_key,
            mode,
        }
    }

    pub fn parse_all(columns: &[String]) -> Vec<Self> {
        columns
            .iter()
            .filter(|c| {
                let upper = c.to_uppercase();
                !upper.starts_with("CHROM")
                    && !upper.starts_with("POS")
                    && !upper.starts_with("REF")
                    && !upper.starts_with("FROM")
                    && !upper.starts_with("TO")
                    && !upper.starts_with("BEG")
                    && !upper.starts_with("END")
                    && *c != "-"
            })
            .map(|c| Self::parse(c))
            .collect()
    }
}
