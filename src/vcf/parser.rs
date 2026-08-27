use std::collections::HashMap;

use crate::vcf::structs::{VcfFields, VcfParsedRecord};

pub struct VcfParser<'a> {
    line: &'a str,
    pos: usize,
}

impl<'a> VcfParser<'a> {
    #[inline]
    pub fn new(line: &'a str) -> Self {
        Self { line, pos: 0 }
    }

    #[inline]
    pub fn parse_standard_fields(&mut self) -> Option<VcfFields<'a>> {
        let fields = VcfFields {
            chrom: self.next_field()?,
            pos: self.next_field()?,
            id: self.next_field()?,
            ref_allele: self.next_field()?,
            alt: self.next_field()?,
            qual: self.next_field()?,
            filter: self.next_field()?,
            info: self.next_field()?,
        };

        Some(fields)
    }

    #[inline]
    fn next_field(&mut self) -> Option<&'a str> {
        let start = self.pos;
        let bytes = self.line.as_bytes();

        while self.pos < bytes.len() && bytes[self.pos] != b'\t' {
            self.pos += 1;
        }

        let field = &self.line[start..self.pos];

        if self.pos < bytes.len() {
            self.pos += 1;
        }

        Some(field)
    }

    #[inline]
    pub fn rest(&self) -> &'a str {
        &self.line[self.pos..]
    }
}

pub struct InfoParser<'a> {
    data: &'a str,
}

impl<'a> InfoParser<'a> {
    pub fn new(data: &'a str) -> Self {
        Self { data }
    }

    pub fn iter(&self) -> InfoIterator<'a> {
        InfoIterator {
            data: self.data,
            pos: 0,
        }
    }
}

pub struct InfoIterator<'a> {
    data: &'a str,
    pos: usize,
}

impl<'a> Iterator for InfoIterator<'a> {
    type Item = (&'a str, Option<&'a str>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.data.len() {
            return None;
        }

        let start = self.pos;
        let bytes = self.data.as_bytes();

        while self.pos < bytes.len() && bytes[self.pos] != b';' {
            self.pos += 1;
        }

        let segment = &self.data[start..self.pos];

        if self.pos < bytes.len() {
            self.pos += 1;
        }

        if let Some(eq_pos) = segment.find('=') {
            Some((&segment[..eq_pos], Some(&segment[eq_pos + 1..])))
        } else {
            Some((segment, None))
        }
    }
}

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
