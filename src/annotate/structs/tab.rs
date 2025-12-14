use anyhow::Result;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use super::annotate_mode::AnnotateMode;
use super::bundle::FieldNumber;
use crate::util::{choose_best_number, extract_info_key, extract_info_number};

#[derive(Debug, Clone)]
pub struct TabColumn {
    pub index: usize,
    pub key: String,
    pub dst_key: String,
    pub mode: AnnotateMode,
    pub number: Option<FieldNumber>,
}

impl TabColumn {
    pub fn new(index: usize, key: String, mode: AnnotateMode, number: Option<FieldNumber>) -> Self {
        Self {
            index,
            dst_key: key.clone(),
            key,
            mode,
            number,
        }
    }

    pub fn with_rename(
        index: usize,
        src_key: String,
        dst_key: String,
        mode: AnnotateMode,
        number: Option<FieldNumber>,
    ) -> Self {
        Self {
            index,
            key: src_key,
            dst_key,
            mode,
            number,
        }
    }

    pub fn is_add_if_missing(&self) -> bool {
        self.mode.replace_missing
    }

    pub fn is_append(&self) -> bool {
        self.mode.set_or_append
    }

    pub fn should_transfer(
        &self,
        src_is_missing: bool,
        dst_exists: bool,
        dst_is_missing: bool,
    ) -> bool {
        self.mode
            .should_transfer(src_is_missing, dst_exists, dst_is_missing)
    }
}

#[derive(Debug, Clone)]
pub struct TabSchema {
    pub chrom_idx: usize,
    pub pos_idx: usize,
    pub ref_idx: Option<usize>,
    pub alt_idx: Option<usize>,
    pub id_idx: Option<usize>,
    pub id_mode: AnnotateMode,
    pub qual_idx: Option<usize>,
    pub qual_mode: AnnotateMode,
    pub filter_idx: Option<usize>,
    pub filter_mode: AnnotateMode,
    pub info_start: Option<usize>,
    pub info_all: bool,
    pub info_all_except: Vec<String>,
    pub info_cols: Vec<TabColumn>,
    pub field_metadata: HashMap<String, FieldNumber>,
    pub match_id: bool,
    pub match_end: bool,
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

            if let Some(key) = extract_info_key(&line) {
                if let Some(number) = extract_info_number(&line) {
                    metadata.insert(key, number);
                }
            }
        }

        Ok(metadata)
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
        let mut id_mode = AnnotateMode::default_mode();
        let mut qual_idx = None;
        let mut qual_mode = AnnotateMode::default_mode();
        let mut filter_idx = None;
        let mut filter_mode = AnnotateMode::default_mode();
        let mut info_start = None;
        let mut info_all = false;
        let mut info_all_except = Vec::new();
        let mut info_cols = Vec::new();
        let mut match_id = false;
        let mut match_end = false;

        for (i, part) in parts.iter().enumerate() {
            if part.is_empty() || *part == "-" {
                continue;
            }

            let (mode, clean_part) = AnnotateMode::parse(part);

            if clean_part.contains(":=") {
                let rename_parts: Vec<&str> = clean_part.splitn(2, ":=").collect();
                if rename_parts.len() == 2 {
                    let dst = rename_parts[0];
                    let src = rename_parts[1];
                    let src_clean = src
                        .strip_prefix("INFO/")
                        .or_else(|| src.strip_prefix("FMT/"))
                        .or_else(|| src.strip_prefix("FORMAT/"))
                        .unwrap_or(src);
                    let dst_clean = dst
                        .strip_prefix("INFO/")
                        .or_else(|| dst.strip_prefix("FMT/"))
                        .or_else(|| dst.strip_prefix("FORMAT/"))
                        .unwrap_or(dst);
                    let number = field_metadata.get(src_clean).copied();
                    info_cols.push(TabColumn::with_rename(
                        i,
                        src_clean.to_string(),
                        dst_clean.to_string(),
                        mode,
                        number,
                    ));
                    continue;
                }
            }

            let base_part = clean_part
                .strip_prefix("INFO/")
                .or_else(|| clean_part.strip_prefix("FMT/"))
                .or_else(|| clean_part.strip_prefix("FORMAT/"))
                .unwrap_or(clean_part);

            if base_part.starts_with('^') {
                let except_tag = &base_part[1..];
                let except_clean = except_tag.strip_prefix("INFO/").unwrap_or(except_tag);
                info_all_except.push(except_clean.to_string());
                continue;
            }

            match base_part.to_uppercase().as_str() {
                "CHROM" => chrom_idx = Some(i),
                "POS" => pos_idx = Some(i),
                "FROM" | "BEG" => {
                    if pos_idx.is_none() {
                        pos_idx = Some(i);
                    }
                }
                "TO" | "END" => {}
                "REF" => ref_idx = Some(i),
                "ALT" => alt_idx = Some(i),
                "ID" => {
                    if mode.match_value {
                        match_id = true;
                    } else {
                        id_idx = Some(i);
                        id_mode = mode;
                    }
                }
                "QUAL" => {
                    qual_idx = Some(i);
                    qual_mode = mode;
                }
                "FILTER" => {
                    filter_idx = Some(i);
                    filter_mode = mode;
                }
                "INFO" => {
                    if mode.match_value {
                        if base_part.to_uppercase().contains("END") {
                            match_end = true;
                        }
                    } else {
                        info_start = Some(i);
                        info_all = true;
                    }
                }
                _ => {
                    if mode.match_value {
                        if base_part.to_uppercase() == "INFO/END"
                            || base_part.to_uppercase() == "END"
                        {
                            match_end = true;
                        }
                        continue;
                    }
                    let number = field_metadata.get(base_part).copied();
                    info_cols.push(TabColumn::new(i, base_part.to_string(), mode, number));
                }
            }
        }

        let chrom_idx = chrom_idx.ok_or_else(|| anyhow::anyhow!("Missing CHROM column"))?;
        let pos_idx = pos_idx.ok_or_else(|| anyhow::anyhow!("Missing POS column"))?;

        Ok(Self {
            chrom_idx,
            pos_idx,
            ref_idx,
            alt_idx,
            id_idx,
            id_mode,
            qual_idx,
            qual_mode,
            filter_idx,
            filter_mode,
            info_start,
            info_all,
            info_all_except,
            info_cols,
            field_metadata,
            match_id,
            match_end,
        })
    }
}
