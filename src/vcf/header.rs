//! VCF header metadata: contig dictionary, INFO/FORMAT field definitions, samples.

use fxhash::FxHashMap;

/// Contig name <-> id table in header order. Ids are dense `u32`, assigned in
/// order of declaration; contigs seen in data but absent from the header are
/// appended in order of first appearance (htslib does the same).
#[derive(Debug, Clone, Default)]
pub struct ContigDict {
    names: Vec<String>,
    lengths: Vec<Option<u64>>,
    ids: FxHashMap<String, u32>,
}

impl ContigDict {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_header_lines<'a, I>(lines: I) -> Self
    where
        I: IntoIterator<Item = &'a str>,
    {
        let mut d = Self::new();
        for line in lines {
            if let Some((id, len)) = parse_contig_line(line) {
                d.insert_with_length(&id, len);
            }
        }
        d
    }

    #[inline]
    pub fn id(&self, name: &str) -> Option<u32> {
        self.ids.get(name).copied()
    }

    pub fn insert(&mut self, name: &str) -> u32 {
        self.insert_with_length(name, None)
    }

    pub fn insert_with_length(&mut self, name: &str, len: Option<u64>) -> u32 {
        if let Some(&id) = self.ids.get(name) {
            if len.is_some() && self.lengths[id as usize].is_none() {
                self.lengths[id as usize] = len;
            }
            return id;
        }
        let id = self.names.len() as u32;
        self.names.push(name.to_string());
        self.lengths.push(len);
        self.ids.insert(name.to_string(), id);
        id
    }

    #[inline]
    pub fn name(&self, id: u32) -> Option<&str> {
        self.names.get(id as usize).map(String::as_str)
    }

    #[inline]
    pub fn length(&self, id: u32) -> Option<u64> {
        self.lengths.get(id as usize).copied().flatten()
    }

    pub fn len(&self) -> usize {
        self.names.len()
    }

    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    pub fn names(&self) -> &[String] {
        &self.names
    }

    pub fn iter(&self) -> impl Iterator<Item = (u32, &str, Option<u64>)> {
        self.names
            .iter()
            .enumerate()
            .map(|(i, n)| (i as u32, n.as_str(), self.lengths[i]))
    }
}

/// Declared cardinality of an INFO/FORMAT field.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FieldNumber {
    Fixed(u32),
    A,
    R,
    G,
    Dot,
}

impl FieldNumber {
    pub fn parse(s: &str) -> Self {
        match s.trim() {
            "A" => Self::A,
            "R" => Self::R,
            "G" => Self::G,
            "." => Self::Dot,
            n => n.parse::<u32>().map(Self::Fixed).unwrap_or(Self::Dot),
        }
    }

    pub fn is_per_allele(self) -> bool {
        matches!(self, Self::A | Self::R | Self::G)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FieldType {
    Integer,
    Float,
    String,
    Character,
    Flag,
}

impl FieldType {
    pub fn parse(s: &str) -> Self {
        match s.trim() {
            "Integer" => Self::Integer,
            "Float" => Self::Float,
            "Character" => Self::Character,
            "Flag" => Self::Flag,
            _ => Self::String,
        }
    }
}

#[derive(Clone, Debug)]
pub struct FieldDef {
    pub id: String,
    pub number: FieldNumber,
    pub ty: FieldType,
    pub description: String,
}

/// Parsed view of a VCF header.
#[derive(Clone, Debug, Default)]
pub struct HeaderInfo {
    pub info: FxHashMap<String, FieldDef>,
    pub format: FxHashMap<String, FieldDef>,
    pub filters: Vec<String>,
    pub samples: Vec<String>,
    pub contigs: ContigDict,
}

impl HeaderInfo {
    pub fn parse<S: AsRef<str>>(lines: &[S]) -> Self {
        let mut h = Self::default();
        for l in lines {
            let line = l.as_ref();
            if let Some((kind, kvs)) = parse_struct_line(line) {
                match kind {
                    "INFO" | "FORMAT" => {
                        let mut def = FieldDef {
                            id: String::new(),
                            number: FieldNumber::Dot,
                            ty: FieldType::String,
                            description: String::new(),
                        };
                        for (k, v) in kvs {
                            match k {
                                "ID" => def.id = v.to_string(),
                                "Number" => def.number = FieldNumber::parse(v),
                                "Type" => def.ty = FieldType::parse(v),
                                "Description" => def.description = v.to_string(),
                                _ => {}
                            }
                        }
                        if def.id.is_empty() {
                            continue;
                        }
                        let map = if kind == "INFO" { &mut h.info } else { &mut h.format };
                        map.insert(def.id.clone(), def);
                    }
                    "FILTER" => {
                        if let Some((_, v)) = kvs.iter().find(|(k, _)| *k == "ID") {
                            h.filters.push(v.to_string());
                        }
                    }
                    "contig" => {
                        if let Some((id, len)) = parse_contig_line(line) {
                            h.contigs.insert_with_length(&id, len);
                        }
                    }
                    _ => {}
                }
            } else if line.starts_with("#CHROM") {
                h.samples = samples_from_chrom_line(line);
            }
        }
        h
    }

    pub fn info_number(&self, key: &str) -> FieldNumber {
        self.info.get(key).map(|d| d.number).unwrap_or(FieldNumber::Dot)
    }

    pub fn format_number(&self, key: &str) -> FieldNumber {
        self.format.get(key).map(|d| d.number).unwrap_or(FieldNumber::Dot)
    }

    pub fn info_type(&self, key: &str) -> Option<FieldType> {
        self.info.get(key).map(|d| d.ty)
    }

    pub fn format_type(&self, key: &str) -> Option<FieldType> {
        self.format.get(key).map(|d| d.ty)
    }
}

/// `##KIND=<k=v,...>` -> ("KIND", [(k, v)]) with quotes stripped from values.
pub fn parse_struct_line(line: &str) -> Option<(&str, Vec<(&str, &str)>)> {
    let rest = line.strip_prefix("##")?;
    let eq = rest.find('=')?;
    let kind = &rest[..eq];
    let body = rest[eq + 1..].strip_prefix('<')?.strip_suffix('>')?;
    let mut out = Vec::new();
    for kv in split_top_commas(body) {
        let Some((k, v)) = kv.split_once('=') else { continue };
        let v = v.trim_matches('"');
        out.push((k.trim(), v));
    }
    Some((kind, out))
}

/// Split on commas that are outside quotes and angle brackets.
pub fn split_top_commas(s: &str) -> Vec<&str> {
    let mut out = Vec::new();
    let bytes = s.as_bytes();
    let mut start = 0;
    let mut depth = 0i32;
    let mut in_quote = false;
    for (i, &b) in bytes.iter().enumerate() {
        match b {
            b'"' => in_quote = !in_quote,
            b'<' if !in_quote => depth += 1,
            b'>' if !in_quote => depth -= 1,
            b',' if !in_quote && depth == 0 => {
                out.push(&s[start..i]);
                start = i + 1;
            }
            _ => {}
        }
    }
    out.push(&s[start..]);
    out
}

/// `##contig=<ID=...,length=...>` -> (ID, length).
pub fn parse_contig_line(line: &str) -> Option<(String, Option<u64>)> {
    let (kind, kvs) = parse_struct_line(line)?;
    if kind != "contig" {
        return None;
    }
    let mut id = None;
    let mut len = None;
    for (k, v) in kvs {
        match k {
            "ID" => id = Some(v.to_string()),
            "length" => len = v.parse::<u64>().ok(),
            _ => {}
        }
    }
    id.map(|i| (i, len))
}

pub fn samples_from_chrom_line(line: &str) -> Vec<String> {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() > 9 {
        cols[9..].iter().map(|s| s.to_string()).collect()
    } else {
        Vec::new()
    }
}

pub fn extract_samples<S: AsRef<str>>(headers: &[S]) -> Vec<String> {
    headers
        .iter()
        .map(|h| h.as_ref())
        .find(|h| h.starts_with("#CHROM"))
        .map(samples_from_chrom_line)
        .unwrap_or_default()
}

#[cfg(test)]
#[path = "../../tests/unit/vcf_header.rs"]
mod tests;
