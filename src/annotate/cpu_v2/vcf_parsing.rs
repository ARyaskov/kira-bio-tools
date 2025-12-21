use crate::vcf::VcfParser;
use indexmap::IndexMap;

pub struct ParsedVcfRecord<'a> {
    pub chrom: &'a str,
    pub pos: u32,
    pub ref_allele: &'a str,
    pub vcf_alt_alleles: Vec<&'a str>,
    pub updated_id: String,
    pub updated_filter: String,
    pub info_map: IndexMap<String, String>,
    pub rest: &'a str,
}

pub fn parse_vcf_record(line: &str) -> Option<ParsedVcfRecord> {
    let mut parser = VcfParser::new(line);
    let rec = parser.parse_standard_fields()?;

    let pos = rec.pos.parse::<u32>().unwrap_or(0);
    let vcf_alt_alleles: Vec<&str> = rec.alt.split(',').collect();

    let mut info_map: IndexMap<String, String> = IndexMap::new();
    for kv in rec.info.split(';') {
        if kv.is_empty() || kv == "." {
            continue;
        }
        let mut parts = kv.splitn(2, '=');
        let k = parts.next().unwrap();
        let v = parts.next().unwrap_or("");
        info_map.insert(k.to_string(), v.to_string());
    }

    Some(ParsedVcfRecord {
        chrom: rec.chrom,
        pos,
        ref_allele: rec.ref_allele,
        vcf_alt_alleles,
        updated_id: rec.id.to_string(),
        updated_filter: rec.filter.to_string(),
        info_map,
        rest: parser.rest(),
    })
}

pub fn format_vcf_output(
    chrom: &str,
    pos: u32,
    ref_allele: &str,
    vcf_alt_alleles: &[&str],
    rest: &str,
    updated_id: String,
    updated_filter: String,
    info_str: String,
) -> String {
    let mut fields = vec![
        chrom.to_string(),
        pos.to_string(),
        updated_id,
        ref_allele.to_string(),
        vcf_alt_alleles.join(","),
        ".".to_string(),
        updated_filter,
        info_str,
    ];

    if !rest.is_empty() {
        fields.push(rest.to_string());
    }

    fields.join("\t")
}
