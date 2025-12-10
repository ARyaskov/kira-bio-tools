use anyhow::Result;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct TabColumn {
    pub index: usize,
    pub key: String,
    pub is_append: bool,
}

#[derive(Debug, Clone)]
pub struct TabSchema {
    pub chrom_idx: usize,
    pub pos_idx: usize,
    pub ref_idx: Option<usize>,
    pub alt_idx: Option<usize>,
    pub id_idx: Option<usize>,
    pub qual_idx: Option<usize>,
    pub filter_idx: Option<usize>,
    pub info_start: Option<usize>,
    pub info_cols: Vec<TabColumn>,
}

impl TabSchema {
    pub fn parse(path: &Path, columns: Option<&str>) -> Result<Self> {
        if let Some(cols) = columns {
            return Self::from_column_spec(cols);
        }

        let file = File::open(path)?;
        let mut reader = BufReader::new(file);
        let mut first_line = String::new();
        reader.read_line(&mut first_line)?;

        if first_line.starts_with('#') {
            let header = first_line.trim_start_matches('#').trim();
            Self::from_header(header)
        } else {
            let ncols = first_line.split('\t').count();
            if ncols >= 9 {
                Self::from_column_spec("CHROM,POS,REF,ALT,ID,QUAL,FILTER,INFO")
            } else if ncols >= 5 {
                Self::from_column_spec("CHROM,POS,REF,ALT,ID")
            } else if ncols >= 4 {
                Self::from_column_spec("CHROM,POS,REF,ALT")
            } else {
                anyhow::bail!("Cannot detect TAB schema from {} columns", ncols)
            }
        }
    }

    fn from_header(header: &str) -> Result<Self> {
        let parts: Vec<&str> = header.split('\t').map(|s| s.trim()).collect();
        let spec = parts.join(",");
        Self::from_column_spec(&spec)
    }

    fn from_column_spec(spec: &str) -> Result<Self> {
        let parts: Vec<&str> = spec.split(',').map(|s| s.trim()).collect();

        let mut chrom_idx = None;
        let mut pos_idx = None;
        let mut ref_idx = None;
        let mut alt_idx = None;
        let mut id_idx = None;
        let mut qual_idx = None;
        let mut filter_idx = None;
        let mut info_start = None;
        let mut info_cols = Vec::new();

        for (i, part) in parts.iter().enumerate() {
            let (is_append, clean_part) = if part.starts_with('+') {
                (true, &part[1..])
            } else if part.starts_with('-') {
                (false, &part[1..])
            } else {
                (false, *part)
            };

            match clean_part {
                "CHROM" => chrom_idx = Some(i),
                "POS" => pos_idx = Some(i),
                "REF" => ref_idx = Some(i),
                "ALT" => alt_idx = Some(i),
                "ID" => id_idx = Some(i),
                "QUAL" => qual_idx = Some(i),
                "FILTER" => filter_idx = Some(i),
                "INFO" => info_start = Some(i),
                _ if clean_part.starts_with("INFO/") => {
                    let key = clean_part.strip_prefix("INFO/").unwrap();
                    info_cols.push(TabColumn {
                        index: i,
                        key: key.to_string(),
                        is_append,
                    });
                }
                _ if clean_part.starts_with("FMT/") || clean_part.starts_with("FORMAT/") => {
                    let key = if let Some(k) = clean_part.strip_prefix("FMT/") {
                        k
                    } else {
                        clean_part.strip_prefix("FORMAT/").unwrap()
                    };
                    info_cols.push(TabColumn {
                        index: i,
                        key: key.to_string(),
                        is_append,
                    });
                }
                _ => {
                    info_cols.push(TabColumn {
                        index: i,
                        key: clean_part.to_string(),
                        is_append,
                    });
                }
            }
        }

        let chrom_idx = chrom_idx.ok_or_else(|| anyhow::anyhow!("CHROM column required"))?;
        let pos_idx = pos_idx.ok_or_else(|| anyhow::anyhow!("POS column required"))?;

        Ok(Self {
            chrom_idx,
            pos_idx,
            ref_idx,
            alt_idx,
            id_idx,
            qual_idx,
            filter_idx,
            info_start,
            info_cols,
        })
    }
}
