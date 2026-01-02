#[derive(Debug, Clone)]
pub struct ParsedFormat {
    pub keys: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct ParsedSample {
    pub raw: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct ParsedVcfRecord {
    pub chrom: String,
    pub pos: u32,
    pub id: String,
    pub ref_allele: String,
    pub alt: String,
    pub qual: String,
    pub filter: String,
    pub info: String,
    pub format: Option<ParsedFormat>,
    pub samples: Vec<ParsedSample>,
}

use crate::vcf::simd::SimdVcfParser;

pub fn parse_vcf_record(line: &str) -> Option<ParsedVcfRecord> {
    let mut fields = line.split('\t');

    let chrom = fields.next()?.to_string();
    let pos = fields.next()?.parse::<u32>().ok()?;
    let id = fields.next()?.to_string();
    let ref_allele = fields.next()?.to_string();
    let alt = fields.next()?.to_string();
    let qual = fields.next()?.to_string();
    let filter = fields.next()?.to_string();
    let info = fields.next()?.to_string();

    let format_str = fields.next();
    let sample_fields: Vec<&str> = fields.collect();

    let (format, samples) = if let Some(fmt_str) = format_str {
        let format_keys: Vec<String> = fmt_str.split(':').map(|s| s.to_string()).collect();
        let format_obj = ParsedFormat { keys: format_keys };

        let parsed_samples: Vec<ParsedSample> = sample_fields
            .iter()
            .map(|sample_data| {
                let sample_values: Vec<String> =
                    sample_data.split(':').map(|s| s.to_string()).collect();
                ParsedSample { raw: sample_values }
            })
            .collect();

        (Some(format_obj), parsed_samples)
    } else {
        (None, Vec::new())
    };

    Some(ParsedVcfRecord {
        chrom,
        pos,
        id,
        ref_allele,
        alt,
        qual,
        filter,
        info,
        format,
        samples,
    })
}

pub fn parse_vcf_record_simd(line: &str, want_format: bool) -> Option<ParsedVcfRecord> {
    let bytes = line.as_bytes();

    if want_format {
        let parsed = SimdVcfParser::parse_line(bytes)?;

        let format = parsed.format.map(|fmt_str| ParsedFormat {
            keys: fmt_str.split(':').map(|s| s.to_string()).collect(),
        });

        let samples = parsed
            .samples
            .iter()
            .map(|sample_data| ParsedSample {
                raw: sample_data.split(':').map(|s| s.to_string()).collect(),
            })
            .collect();

        return Some(ParsedVcfRecord {
            chrom: parsed.chrom.to_string(),
            pos: parsed.position()?,
            id: parsed.id.to_string(),
            ref_allele: parsed.ref_allele.to_string(),
            alt: parsed.alt.to_string(),
            qual: parsed.qual.to_string(),
            filter: parsed.filter.to_string(),
            info: parsed.info.to_string(),
            format,
            samples,
        });
    }

    let parsed = SimdVcfParser::parse_fields(bytes)?;
    Some(ParsedVcfRecord {
        chrom: parsed.chrom.to_string(),
        pos: parsed.position()?,
        id: parsed.id.to_string(),
        ref_allele: parsed.ref_allele.to_string(),
        alt: parsed.alt.to_string(),
        qual: parsed.qual.to_string(),
        filter: parsed.filter.to_string(),
        info: parsed.info.to_string(),
        format: None,
        samples: Vec::new(),
    })
}

impl ParsedVcfRecord {
    pub fn to_line(&self) -> String {
        let mut result = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.chrom,
            self.pos,
            self.id,
            self.ref_allele,
            self.alt,
            self.qual,
            self.filter,
            self.info
        );

        if let Some(format) = &self.format {
            result.push('\t');
            result.push_str(&format.keys.join(":"));

            for sample in &self.samples {
                result.push('\t');
                result.push_str(&sample.raw.join(":"));
            }
        }

        result
    }
}
