use std::collections::HashMap;

use crate::vcf::structs::VcfParsedRecord;

pub fn parse_vcf_full_line(line: &str) -> Option<VcfParsedRecord> {
    if line.starts_with('#') {
        return None;
    }
    let cols: Vec<&str> = line.trim_end().split('\t').collect();
    if cols.len() < 8 {
        return None;
    }

    let filter = cols[6].to_string();
    let mut info = HashMap::new();

    for item in cols[7].split(';') {
        if let Some((k, v)) = item.split_once('=') {
            info.insert(k.to_string(), v.to_string());
        }
    }

    Some(VcfParsedRecord {
        chrom: cols[0].to_string(),
        pos: cols[1].parse().ok()?,
        filter,
        info,
        raw_line: line.to_string(),
    })
}

pub fn extract_contig_id(line: &str) -> Option<String> {
    let start = line.find("ID=")? + 3;
    let rest = &line[start..];
    let end = rest.find(|c| c == ',' || c == '>')?;
    Some(rest[..end].to_string())
}
