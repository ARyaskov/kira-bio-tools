use super::super::structs::annotate_mode::AnnotateMode;

#[derive(Debug, Clone)]
pub struct ColumnSpec {
    pub key: String,
    pub dst_key: String,
    pub mode: AnnotateMode,
}

impl ColumnSpec {
    pub fn parse(spec: &str) -> Self {
        let (mode, rest) = AnnotateMode::parse(spec);

        let (src_key, dst_key) = if rest.contains(":=") {
            let parts: Vec<&str> = rest.splitn(2, ":=").collect();
            if parts.len() == 2 {
                let dst = parts[0].strip_prefix("INFO/").unwrap_or(parts[0]);
                let src = parts[1].strip_prefix("INFO/").unwrap_or(parts[1]);
                (src.to_string(), dst.to_string())
            } else {
                let key = rest.strip_prefix("INFO/").unwrap_or(rest).to_string();
                (key.clone(), key)
            }
        } else {
            let key = rest.strip_prefix("INFO/").unwrap_or(rest).to_string();
            (key.clone(), key)
        };

        Self {
            key: src_key,
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
                    && !upper.starts_with("ALT")
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
