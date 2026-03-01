use crate::vcf::VcfParser;
use anyhow::Result;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

#[inline]
pub fn normalize(ref_allele: &str, alt_allele: &str) -> (usize, usize) {
    let rb = ref_allele.as_bytes();
    let ab = alt_allele.as_bytes();

    let mut prefix = 0;
    while prefix < rb.len() && prefix < ab.len() && rb[prefix] == ab[prefix] {
        prefix += 1;
    }

    let mut suffix = 0;
    let mut ri = rb.len();
    let mut ai = ab.len();
    while ri > prefix && ai > prefix && rb[ri - 1] == ab[ai - 1] {
        ri -= 1;
        ai -= 1;
        suffix += 1;
    }

    (prefix, suffix)
}

pub fn turbo_norm_vcf(input: &Path, output: &Path) -> Result<()> {
    let in_file = File::open(input)?;
    let mut reader = BufReader::new(in_file);
    let out_file = File::create(output)?;
    let mut out = BufWriter::new(out_file);

    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        let raw = line.trim_end_matches(['\n', '\r']);
        if raw.is_empty() || raw.starts_with('#') {
            writeln!(out, "{raw}")?;
            continue;
        }

        let mut parser = VcfParser::new(raw);
        if let Some(fields) = parser.parse_standard_fields() {
            let ref_allele = fields.ref_allele.as_bytes();
            let alt_allele = fields.alt.as_bytes();
            let (prefix, suffix) = normalize(fields.ref_allele, fields.alt);
            let nr = &ref_allele[prefix..ref_allele.len() - suffix];
            let na = &alt_allele[prefix..alt_allele.len() - suffix];

            write!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                fields.chrom,
                fields.pos,
                fields.id,
                String::from_utf8_lossy(nr),
                String::from_utf8_lossy(na),
                fields.qual,
                fields.filter,
                fields.info
            )?;
            let rest = parser.rest();
            if !rest.is_empty() {
                write!(out, "\t{rest}")?;
            }
            writeln!(out)?;
        } else {
            writeln!(out, "{raw}")?;
        }
    }
    out.flush()?;
    Ok(())
}
