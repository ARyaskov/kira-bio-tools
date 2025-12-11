use anyhow::Result;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use super::bundle::FieldNumber;

#[derive(Debug, Clone)]
pub struct TabColumn {
    pub index: usize,
    pub key: String,
    pub is_append: bool,
    pub number: Option<FieldNumber>,
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
    pub field_metadata: HashMap<String, FieldNumber>,
}

impl TabSchema {
    pub fn parse(path: &Path, columns: Option<&str>) -> Result<Self> {
        let field_metadata = Self::parse_header_file(path)?;

        if let Some(cols) = columns {
            return Self::from_column_spec(cols, field_metadata);
        }

        let file = File::open(path)?;
        let mut reader = BufReader::new(file);
        let mut first_line = String::new();
        reader.read_line(&mut first_line)?;

        if first_line.starts_with('#') && !first_line.starts_with("##") {
            let header = first_line.trim_start_matches('#').trim();
            Self::from_header(header, field_metadata)
        } else {
            let ncols = first_line.split('\t').count();
            if ncols >= 9 {
                Self::from_column_spec("CHROM,POS,REF,ALT,ID,QUAL,FILTER,INFO", field_metadata)
            } else if ncols >= 5 {
                Self::from_column_spec("CHROM,POS,REF,ALT,ID", field_metadata)
            } else if ncols >= 4 {
                Self::from_column_spec("CHROM,POS,REF,ALT", field_metadata)
            } else {
                anyhow::bail!("Cannot detect TAB schema from {} columns", ncols)
            }
        }
    }

    fn parse_header_file(tab_path: &Path) -> Result<HashMap<String, FieldNumber>> {
        let mut metadata = HashMap::new();

        let hdr_path = tab_path.with_extension("hdr");
        if !hdr_path.exists() {
            return Ok(metadata);
        }

        let file = File::open(&hdr_path)?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line?;
            if !line.starts_with("##INFO=") {
                continue;
            }

            if let Some(key) = Self::extract_info_key(&line) {
                if let Some(number) = Self::extract_info_number(&line) {
                    metadata.insert(key, number);
                }
            }
        }

        Ok(metadata)
    }

    fn extract_info_key(line: &str) -> Option<String> {
        if let Some(start) = line.find("ID=") {
            let rest = &line[start + 3..];
            if let Some(end) = rest.find(',') {
                return Some(rest[..end].to_string());
            }
        }
        None
    }

    fn extract_info_number(line: &str) -> Option<FieldNumber> {
        if let Some(start) = line.find("Number=") {
            let rest = &line[start + 7..];
            if let Some(end) = rest.find(',') {
                let number_str = &rest[..end];
                return match number_str {
                    "0" => Some(FieldNumber::Zero),
                    "1" => Some(FieldNumber::One),
                    "A" => Some(FieldNumber::A),
                    "R" => Some(FieldNumber::R),
                    "G" => Some(FieldNumber::G),
                    "." => Some(FieldNumber::Many),
                    _ => {
                        if number_str.parse::<i32>().is_ok() {
                            Some(FieldNumber::Many)
                        } else {
                            None
                        }
                    }
                };
            }
        }
        None
    }

    fn from_header(header: &str, field_metadata: HashMap<String, FieldNumber>) -> Result<Self> {
        let parts: Vec<&str> = header.split('\t').map(|s| s.trim()).collect();
        let spec = parts.join(",");
        Self::from_column_spec(&spec, field_metadata)
    }

    fn from_column_spec(spec: &str, field_metadata: HashMap<String, FieldNumber>) -> Result<Self> {
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
                    let number = field_metadata.get(key).copied();
                    info_cols.push(TabColumn {
                        index: i,
                        key: key.to_string(),
                        is_append,
                        number,
                    });
                }
                _ if clean_part.starts_with("FMT/") || clean_part.starts_with("FORMAT/") => {
                    let key = if let Some(k) = clean_part.strip_prefix("FMT/") {
                        k
                    } else {
                        clean_part.strip_prefix("FORMAT/").unwrap()
                    };
                    let number = field_metadata.get(key).copied();
                    info_cols.push(TabColumn {
                        index: i,
                        key: key.to_string(),
                        is_append,
                        number,
                    });
                }
                _ => {
                    let number = field_metadata.get(clean_part).copied();
                    info_cols.push(TabColumn {
                        index: i,
                        key: clean_part.to_string(),
                        is_append,
                        number,
                    });
                }
            }
        }

        let chrom_idx = chrom_idx.ok_or_else(|| anyhow::anyhow!("CHROM column missing"))?;
        let pos_idx = pos_idx.ok_or_else(|| anyhow::anyhow!("POS column missing"))?;

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
            field_metadata,
        })
    }
}
