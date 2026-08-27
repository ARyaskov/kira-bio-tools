use fxhash::FxHashMap;

#[derive(Default, Clone)]
pub struct BcfHeaderDict {
    pub raw_lines: Vec<String>,
    pub contigs: Vec<String>,
    pub samples: Vec<String>,
    pub info: Vec<HdrField>,
    pub format: Vec<HdrField>,
    pub filter: Vec<HdrField>,
    pub contig_idx: FxHashMap<String, u32>,
    pub info_idx: FxHashMap<String, u32>,
    pub format_idx: FxHashMap<String, u32>,
    pub filter_idx: FxHashMap<String, u32>,
}

#[derive(Clone)]
pub struct HdrField {
    pub id: String,
    pub number: String,
    pub typ: String,
    pub description: String,
    pub idx: u32,
}

pub fn parse_header_to_dict(headers: &[String]) -> BcfHeaderDict {
    let mut d = BcfHeaderDict::default();
    let mut contig_next: u32 = 0;
    let mut info_format_filter_next: u32 = 1;
    let mut have_pass = false;

    let extract_idx = |s: &str| -> Option<u32> {
        s.split(',').find_map(|kv| kv.strip_prefix("IDX="))
            .and_then(|v| v.trim_end_matches('>').parse().ok())
    };

    for h in headers {
        d.raw_lines.push(h.clone());
        if let Some(rest) = h.strip_prefix("##contig=") {
            if let Some(id) = extract_struct_id(rest) {
                let idx = extract_idx(rest).unwrap_or(contig_next);
                contig_next = contig_next.max(idx + 1);
                d.contig_idx.insert(id.clone(), idx);
                d.contigs.push(id);
            }
        } else if let Some(rest) = h.strip_prefix("##INFO=<") {
            if let Some(mut f) = parse_struct_field(rest) {
                let idx = extract_idx(rest).unwrap_or(info_format_filter_next);
                info_format_filter_next = info_format_filter_next.max(idx + 1);
                f.idx = idx;
                d.info_idx.insert(f.id.clone(), idx);
                d.info.push(f);
            }
        } else if let Some(rest) = h.strip_prefix("##FORMAT=<") {
            if let Some(mut f) = parse_struct_field(rest) {
                let idx = extract_idx(rest).unwrap_or(info_format_filter_next);
                info_format_filter_next = info_format_filter_next.max(idx + 1);
                f.idx = idx;
                d.format_idx.insert(f.id.clone(), idx);
                d.format.push(f);
            }
        } else if let Some(rest) = h.strip_prefix("##FILTER=<") {
            if let Some(mut f) = parse_struct_field(rest) {
                let idx = extract_idx(rest).unwrap_or(if f.id == "PASS" { 0 } else { info_format_filter_next });
                if f.id != "PASS" {
                    info_format_filter_next = info_format_filter_next.max(idx + 1);
                }
                f.idx = idx;
                if f.id == "PASS" { have_pass = true; }
                d.filter_idx.insert(f.id.clone(), idx);
                d.filter.push(f);
            }
        } else if h.starts_with("#CHROM") {
            let cols: Vec<&str> = h.split('\t').collect();
            if cols.len() > 9 {
                d.samples = cols[9..].iter().map(|s| s.to_string()).collect();
            }
        }
    }
    if !have_pass {
        d.filter_idx.insert("PASS".into(), 0);
        d.filter.insert(0, HdrField {
            id: "PASS".into(), number: ".".into(), typ: ".".into(),
            description: "All filters passed".into(), idx: 0,
        });
    }
    d
}

fn extract_struct_id(rest: &str) -> Option<String> {
    let body = rest.strip_prefix('<')?.strip_suffix('>')?;
    for kv in body.split(',') {
        if let Some(v) = kv.strip_prefix("ID=") { return Some(v.to_string()); }
    }
    None
}

fn parse_struct_field(rest: &str) -> Option<HdrField> {
    let body = rest.strip_suffix('>')?;
    let mut id = String::new();
    let mut number = ".".to_string();
    let mut typ = "String".to_string();
    let mut desc = String::new();
    for kv in body.split(',') {
        if let Some(v) = kv.strip_prefix("ID=") { id = v.to_string(); }
        else if let Some(v) = kv.strip_prefix("Number=") { number = v.to_string(); }
        else if let Some(v) = kv.strip_prefix("Type=") { typ = v.to_string(); }
        else if let Some(v) = kv.strip_prefix("Description=") { desc = v.trim_matches('"').to_string(); }
    }
    if id.is_empty() { return None; }
    Some(HdrField { id, number, typ, description: desc, idx: 0 })
}

/// Serialize header dict back to a VCF header text block (with IDX= tags inserted).
pub fn serialize_header(d: &BcfHeaderDict) -> String {
    let mut out = String::new();
    let mut had_fileformat = false;
    for line in &d.raw_lines {
        if line.starts_with("##fileformat=") { had_fileformat = true; }
        if line.starts_with("#CHROM") { continue; }
        if line.starts_with("##contig=") {
            if let Some(id) = line.strip_prefix("##contig=").and_then(extract_struct_id) {
                if let Some(idx) = d.contig_idx.get(&id) {
                    out.push_str(&inject_idx(line, *idx));
                    out.push('\n');
                    continue;
                }
            }
        }
        if line.starts_with("##INFO=<") || line.starts_with("##FORMAT=<") || line.starts_with("##FILTER=<") {
            if let Some(id) = extract_struct_id(line.split_once('<').map(|p| p.1).unwrap_or("")) {
                let idx_opt = if line.starts_with("##INFO=") { d.info_idx.get(&id) }
                    else if line.starts_with("##FORMAT=") { d.format_idx.get(&id) }
                    else { d.filter_idx.get(&id) };
                if let Some(idx) = idx_opt {
                    out.push_str(&inject_idx(line, *idx));
                    out.push('\n');
                    continue;
                }
            }
        }
        out.push_str(line);
        out.push('\n');
    }
    if !had_fileformat { out = format!("##fileformat=VCFv4.2\n{}", out); }
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    if !d.samples.is_empty() {
        out.push_str("\tFORMAT");
        for s in &d.samples { out.push('\t'); out.push_str(s); }
    }
    out.push('\n');
    out
}

fn inject_idx(line: &str, idx: u32) -> String {
    if line.contains("IDX=") { return line.to_string(); }
    if let Some(p) = line.rfind('>') {
        let mut s = String::with_capacity(line.len() + 12);
        s.push_str(&line[..p]);
        s.push_str(&format!(",IDX={}>", idx));
        s
    } else { line.to_string() }
}

#[cfg(test)]
#[path = "../../tests/unit/bcf_header.rs"]
mod tests;
