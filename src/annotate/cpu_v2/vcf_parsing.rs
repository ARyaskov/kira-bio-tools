//! Zero-copy VCF record parsing for the annotate hot loop.
//!
//! Lifetime-borrowed view over the input line; samples are stored as raw
//! `&'a str` slices and only split-and-reassembled by the FORMAT-annotation
//! path when actually modified.

use crate::vcf::simd::SimdVcfParser;

/// Colon-separated `FORMAT` key list, e.g. `"GT:DP:GQ"`.
#[derive(Debug, Clone, Copy)]
pub struct ParsedFormat<'a> {
    pub raw: &'a str,
}

impl<'a> ParsedFormat<'a> {
    /// Iterator over individual key names.
    #[inline]
    pub fn keys(&self) -> impl Iterator<Item = &'a str> {
        self.raw.split(':')
    }

    /// 0-based index of `key` in the FORMAT, or `None`.
    #[inline]
    pub fn position(&self, key: &str) -> Option<usize> {
        self.keys().position(|k| k == key)
    }

    /// Owned `Vec<String>` materialisation.
    #[inline]
    pub fn keys_owned(&self) -> Vec<String> {
        self.keys().map(str::to_string).collect()
    }
}

/// Colon-separated single-sample value string, e.g. `"0/1:30:99"`.
#[derive(Debug, Clone, Copy)]
pub struct ParsedSample<'a> {
    pub raw: &'a str,
}

impl<'a> ParsedSample<'a> {
    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = &'a str> {
        self.raw.split(':')
    }

    /// `idx`-th subfield, or `None` if the sample is shorter than expected.
    #[inline]
    pub fn value(&self, idx: usize) -> Option<&'a str> {
        self.iter().nth(idx)
    }

    /// Whether the sample is empty (no subfields, or all subfields empty).
    #[inline]
    pub fn is_empty_subfields(&self) -> bool {
        self.raw.is_empty() || self.iter().all(str::is_empty)
    }
}

#[derive(Debug, Clone)]
pub struct ParsedVcfRecord<'a> {
    pub chrom: &'a str,
    pub pos: u32,
    pub id: &'a str,
    pub ref_allele: &'a str,
    pub alt: &'a str,
    pub qual: &'a str,
    pub filter: &'a str,
    pub info: &'a str,
    pub format: Option<ParsedFormat<'a>>,
    pub samples: Vec<ParsedSample<'a>>,
}

pub fn parse_vcf_record<'a>(line: &'a str) -> Option<ParsedVcfRecord<'a>> {
    let mut fields = line.split('\t');

    let chrom = fields.next()?;
    let pos = fields.next()?.parse::<u32>().ok()?;
    let id = fields.next()?;
    let ref_allele = fields.next()?;
    let alt = fields.next()?;
    let qual = fields.next()?;
    let filter = fields.next()?;
    let info = fields.next()?;

    let format_str = fields.next();
    let sample_fields: Vec<&str> = fields.collect();

    let (format, samples) = match format_str {
        Some(fmt) => (
            Some(ParsedFormat { raw: fmt }),
            sample_fields
                .into_iter()
                .map(|raw| ParsedSample { raw })
                .collect(),
        ),
        None => (None, Vec::new()),
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

pub fn parse_vcf_record_simd<'a>(line: &'a str, want_format: bool) -> Option<ParsedVcfRecord<'a>> {
    let bytes = line.as_bytes();

    if want_format {
        let parsed = SimdVcfParser::parse_line(bytes)?;
        let mut format = parsed.format.map(|fmt| ParsedFormat { raw: fmt });
        let mut samples: Vec<ParsedSample<'a>> = parsed
            .samples
            .iter()
            .map(|raw| ParsedSample { raw })
            .collect();

        // SIMD parser can produce shorter sample lists on unusual whitespace;
        // fall back to a scalar re-split when so.
        let needs_fallback = samples.is_empty()
            || samples.iter().any(|s| s.is_empty_subfields())
            || sample_count_from_line(line) != samples.len();
        if needs_fallback {
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() > 9 {
                let fmt_str = cols[8];
                format = if fmt_str.is_empty() || fmt_str == "." {
                    None
                } else {
                    Some(ParsedFormat { raw: fmt_str })
                };
                samples = cols[9..]
                    .iter()
                    .map(|raw| ParsedSample { raw })
                    .collect();
            }
        }

        return Some(ParsedVcfRecord {
            chrom: parsed.chrom,
            pos: parsed.position()?,
            id: parsed.id,
            ref_allele: parsed.ref_allele,
            alt: parsed.alt,
            qual: parsed.qual,
            filter: parsed.filter,
            info: parsed.info,
            format,
            samples,
        });
    }

    let parsed = SimdVcfParser::parse_fields(bytes)?;
    Some(ParsedVcfRecord {
        chrom: parsed.chrom,
        pos: parsed.position()?,
        id: parsed.id,
        ref_allele: parsed.ref_allele,
        alt: parsed.alt,
        qual: parsed.qual,
        filter: parsed.filter,
        info: parsed.info,
        format: None,
        samples: Vec::new(),
    })
}

/// Rewires the sample list if the SIMD parser fell short. Resulting borrows
/// still point into the original `line: &'a str`.
pub fn patch_samples_from_line<'a>(parsed: &mut ParsedVcfRecord<'a>, line: &'a str) {
    let sample_count = sample_count_from_line(line);
    if !parsed.samples.is_empty()
        && sample_count == parsed.samples.len()
        && !parsed.samples.iter().any(ParsedSample::is_empty_subfields)
    {
        return;
    }

    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() > 9 {
        let fmt_str = cols[8];
        parsed.format = if fmt_str.is_empty() || fmt_str == "." {
            None
        } else {
            Some(ParsedFormat { raw: fmt_str })
        };
        parsed.samples = cols[9..]
            .iter()
            .map(|raw| ParsedSample { raw: *raw })
            .collect();
    } else {
        parsed.format = None;
        parsed.samples.clear();
    }
}

#[inline]
fn sample_count_from_line(line: &str) -> usize {
    let mut tabs = 0usize;
    for &b in line.as_bytes() {
        if b == b'\t' {
            tabs += 1;
        }
    }
    tabs.saturating_sub(8)
}

impl<'a> ParsedVcfRecord<'a> {
    /// Serialises the record back to a tab-separated VCF line.
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
            result.push_str(format.raw);
            for sample in &self.samples {
                result.push('\t');
                result.push_str(sample.raw);
            }
        }

        result
    }
}
