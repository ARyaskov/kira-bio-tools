#[derive(Debug, Clone)]
pub struct AnnotatedVcfRecord {
    pub chrom: String,
    pub pos: u32,
    pub id: String,
    pub ref_allele: String,
    pub alt: String,
    pub qual: String,
    pub filter: String,
    pub info: String,
    pub format: Option<String>,
    pub samples: Vec<String>,
}

pub fn format_vcf_output(rec: &AnnotatedVcfRecord) -> String {
    let mut out = Vec::with_capacity(9 + rec.samples.len());

    out.push(rec.chrom.clone());
    out.push(rec.pos.to_string());
    out.push(rec.id.clone());
    out.push(rec.ref_allele.clone());
    out.push(rec.alt.clone());
    out.push(rec.qual.clone());
    out.push(rec.filter.clone());
    out.push(rec.info.clone());

    if let Some(fmt) = &rec.format {
        out.push(fmt.clone());
        out.extend(rec.samples.clone());
    }

    out.join("\t")
}
