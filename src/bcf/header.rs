use fxhash::FxHashMap;

use crate::vcf::header::{parse_struct_line, samples_from_chrom_line};

/// BCF header dictionaries. INFO, FORMAT and FILTER share one id space keyed
/// by tag name (htslib `BCF_DT_ID`, with PASS fixed at 0); contigs have their
/// own (`BCF_DT_CTG`).
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
    info_by_idx: Vec<u32>,
    format_by_idx: Vec<u32>,
    filter_by_idx: Vec<u32>,
    contig_by_idx: Vec<u32>,
}

#[derive(Clone)]
pub struct HdrField {
    pub id: String,
    pub number: String,
    pub typ: String,
    pub description: String,
    pub idx: u32,
}

const NONE: u32 = u32::MAX;

fn set_by_idx(table: &mut Vec<u32>, idx: u32, pos: usize) {
    let i = idx as usize;
    if table.len() <= i {
        table.resize(i + 1, NONE);
    }
    if table[i] == NONE {
        table[i] = pos as u32;
    }
}

impl BcfHeaderDict {
    pub fn info_field(&self, idx: u32) -> Option<&HdrField> {
        self.info_by_idx.get(idx as usize).and_then(|&p| if p == NONE { None } else { self.info.get(p as usize) })
    }

    pub fn format_field(&self, idx: u32) -> Option<&HdrField> {
        self.format_by_idx.get(idx as usize).and_then(|&p| if p == NONE { None } else { self.format.get(p as usize) })
    }

    pub fn filter_field(&self, idx: u32) -> Option<&HdrField> {
        self.filter_by_idx.get(idx as usize).and_then(|&p| if p == NONE { None } else { self.filter.get(p as usize) })
    }

    pub fn contig_name(&self, rid: u32) -> Option<&str> {
        self.contig_by_idx
            .get(rid as usize)
            .and_then(|&p| if p == NONE { None } else { self.contigs.get(p as usize).map(String::as_str) })
    }

    fn rebuild_tables(&mut self) {
        self.info_by_idx.clear();
        self.format_by_idx.clear();
        self.filter_by_idx.clear();
        self.contig_by_idx.clear();
        for (p, f) in self.info.iter().enumerate() {
            set_by_idx(&mut self.info_by_idx, f.idx, p);
        }
        for (p, f) in self.format.iter().enumerate() {
            set_by_idx(&mut self.format_by_idx, f.idx, p);
        }
        for (p, f) in self.filter.iter().enumerate() {
            set_by_idx(&mut self.filter_by_idx, f.idx, p);
        }
        for (p, name) in self.contigs.iter().enumerate() {
            if let Some(&idx) = self.contig_idx.get(name) {
                set_by_idx(&mut self.contig_by_idx, idx, p);
            }
        }
    }
}

struct IdSpace {
    by_name: FxHashMap<String, u32>,
    used: std::collections::HashSet<u32>,
    next: u32,
}

impl IdSpace {
    fn new(start: u32) -> Self {
        Self {
            by_name: FxHashMap::default(),
            used: std::collections::HashSet::new(),
            next: start,
        }
    }

    fn assign(&mut self, name: &str, explicit: Option<u32>) -> u32 {
        if let Some(&idx) = self.by_name.get(name) {
            return idx;
        }
        let idx = match explicit {
            Some(i) => i,
            None => {
                while self.used.contains(&self.next) {
                    self.next += 1;
                }
                self.next
            }
        };
        self.used.insert(idx);
        self.by_name.insert(name.to_string(), idx);
        idx
    }
}

pub fn parse_header_to_dict(headers: &[String]) -> BcfHeaderDict {
    let mut d = BcfHeaderDict::default();
    let mut ids = IdSpace::new(0);
    let mut contig_ids = IdSpace::new(0);
    // PASS is always the first entry of the shared dictionary.
    let pass_idx = ids.assign("PASS", Some(0));
    let mut have_pass_line = false;

    for h in headers {
        d.raw_lines.push(h.clone());
        if let Some((kind, kvs)) = parse_struct_line(h) {
            let id = kvs.iter().find(|(k, _)| *k == "ID").map(|(_, v)| v.to_string());
            let explicit = kvs.iter().find(|(k, _)| *k == "IDX").and_then(|(_, v)| v.parse::<u32>().ok());
            let Some(id) = id else { continue };
            match kind {
                "contig" => {
                    let idx = contig_ids.assign(&id, explicit);
                    if !d.contig_idx.contains_key(&id) {
                        d.contig_idx.insert(id.clone(), idx);
                        d.contigs.push(id);
                    }
                }
                "INFO" | "FORMAT" | "FILTER" => {
                    let mut f = field_from_kvs(&id, &kvs);
                    f.idx = ids.assign(&id, explicit);
                    match kind {
                        "INFO" => {
                            if !d.info_idx.contains_key(&id) {
                                d.info_idx.insert(id, f.idx);
                                d.info.push(f);
                            }
                        }
                        "FORMAT" => {
                            if !d.format_idx.contains_key(&id) {
                                d.format_idx.insert(id, f.idx);
                                d.format.push(f);
                            }
                        }
                        _ => {
                            if id == "PASS" {
                                have_pass_line = true;
                            }
                            if !d.filter_idx.contains_key(&id) {
                                d.filter_idx.insert(id, f.idx);
                                d.filter.push(f);
                            }
                        }
                    }
                }
                _ => {}
            }
        } else if h.starts_with("#CHROM") {
            d.samples = samples_from_chrom_line(h);
        }
    }
    if !have_pass_line {
        d.filter_idx.insert("PASS".into(), pass_idx);
        d.filter.insert(
            0,
            HdrField {
                id: "PASS".into(),
                number: ".".into(),
                typ: ".".into(),
                description: "All filters passed".into(),
                idx: pass_idx,
            },
        );
    }
    d.rebuild_tables();
    d
}

fn field_from_kvs(id: &str, kvs: &[(&str, &str)]) -> HdrField {
    let mut f = HdrField {
        id: id.to_string(),
        number: ".".into(),
        typ: "String".into(),
        description: String::new(),
        idx: 0,
    };
    for (k, v) in kvs {
        match *k {
            "Number" => f.number = v.to_string(),
            "Type" => f.typ = v.to_string(),
            "Description" => f.description = v.to_string(),
            _ => {}
        }
    }
    f
}

/// Serialize the header back to text with `IDX=` tags on every dictionary
/// line, a PASS filter line if the source lacked one, and the `#CHROM` line.
pub fn serialize_header(d: &BcfHeaderDict) -> String {
    let mut out = String::new();
    let mut had_fileformat = false;
    let has_pass_line = d.raw_lines.iter().any(|l| {
        parse_struct_line(l).is_some_and(|(k, kvs)| k == "FILTER" && kvs.iter().any(|(a, b)| *a == "ID" && *b == "PASS"))
    });
    let pass_line = format!(
        "##FILTER=<ID=PASS,Description=\"All filters passed\",IDX={}>",
        d.filter_idx.get("PASS").copied().unwrap_or(0)
    );
    for line in &d.raw_lines {
        if line.starts_with("##fileformat=") {
            had_fileformat = true;
            out.push_str(line);
            out.push('\n');
            if !has_pass_line {
                out.push_str(&pass_line);
                out.push('\n');
            }
            continue;
        }
        if line.starts_with("#CHROM") {
            continue;
        }
        if let Some((kind, kvs)) = parse_struct_line(line) {
            let id = kvs.iter().find(|(k, _)| *k == "ID").map(|(_, v)| *v);
            let idx = match (kind, id) {
                ("contig", Some(id)) => d.contig_idx.get(id),
                ("INFO", Some(id)) => d.info_idx.get(id),
                ("FORMAT", Some(id)) => d.format_idx.get(id),
                ("FILTER", Some(id)) => d.filter_idx.get(id),
                _ => None,
            };
            if let Some(idx) = idx {
                out.push_str(&inject_idx(line, *idx));
                out.push('\n');
                continue;
            }
        }
        out.push_str(line);
        out.push('\n');
    }
    if !had_fileformat {
        let mut pre = String::from("##fileformat=VCFv4.2\n");
        if !has_pass_line {
            pre.push_str(&pass_line);
            pre.push('\n');
        }
        out = format!("{pre}{out}");
    }
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    if !d.samples.is_empty() {
        out.push_str("\tFORMAT");
        for s in &d.samples {
            out.push('\t');
            out.push_str(s);
        }
    }
    out.push('\n');
    out
}

fn inject_idx(line: &str, idx: u32) -> String {
    if let Some((_, kvs)) = parse_struct_line(line) {
        if kvs.iter().any(|(k, _)| *k == "IDX") {
            return line.to_string();
        }
    }
    if let Some(p) = line.rfind('>') {
        let mut s = String::with_capacity(line.len() + 12);
        s.push_str(&line[..p]);
        s.push_str(&format!(",IDX={}>", idx));
        s
    } else {
        line.to_string()
    }
}

#[cfg(test)]
#[path = "../../tests/unit/bcf_header.rs"]
mod tests;
