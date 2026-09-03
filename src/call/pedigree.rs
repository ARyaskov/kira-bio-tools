//! Pedigree / sample-groups / ploidy file parsers for `bcftools call`.

use anyhow::{Context, Result, bail};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use super::mcall::{SampleGroup, TrioFamily};

#[derive(Clone, Debug)]
pub struct PloidyRegion {
    pub chrom: String,
    pub beg: u32,
    pub end: u32,
    pub sex: String,
    pub ploidy: u8,
}

/// Parse `--ploidy-file FILE` table (`chrom beg end sex ploidy`); `*` is a
/// wildcard in the first three columns, and lines match in file order.
pub fn parse_ploidy_file<P: AsRef<Path>>(p: P) -> Result<Vec<PloidyRegion>> {
    let mut out = Vec::new();
    for line in BufReader::new(File::open(p.as_ref())?).lines() {
        let l = line?;
        let t = l.trim();
        if t.is_empty() || t.starts_with('#') { continue; }
        let parts: Vec<&str> = t.split_whitespace().collect();
        if parts.len() < 5 { bail!("--ploidy-file: expected 5 fields, got {}", parts.len()); }
        out.push(PloidyRegion {
            chrom: parts[0].to_string(),
            beg: parse_bound(parts[1], 0).context("ploidy-file beg")?,
            end: parse_bound(parts[2], u32::MAX).context("ploidy-file end")?,
            sex: parts[3].to_string(),
            ploidy: parts[4].parse().context("ploidy-file ploidy")?,
        });
    }
    Ok(out)
}

fn parse_bound(field: &str, wildcard: u32) -> Result<u32> {
    if field == "*" { return Ok(wildcard); }
    Ok(field.parse()?)
}

/// Compute per-sample ploidy at a given site (chrom:pos) given:
///  - ploidy regions (`chrom beg end sex ploidy`)
///  - sex map (sample → "M"/"F"/...)
///  - sample list with default ploidy fallback
pub fn ploidy_at_site(
    chrom: &str, pos: u32,
    samples: &[String],
    sex_map: &HashMap<String, String>,
    regions: &[PloidyRegion],
    default_ploidy: u8,
) -> Vec<u8> {
    samples.iter().map(|name| {
        let sex = sex_map.get(name).map(|s| s.as_str()).unwrap_or("");
        // A bare 0/1/2 in the sample file is the ploidy itself, not a sex.
        match sex {
            "0" => return 0,
            "1" => return 1,
            "2" => return 2,
            _ => {}
        }
        for r in regions {
            // `*` applies to every chromosome / every sex, as in bcftools.
            let chrom_matches = r.chrom == "*" || r.chrom == chrom;
            let sex_matches = r.sex == "*" || r.sex == sex;
            if chrom_matches && pos >= r.beg && pos <= r.end && sex_matches {
                return r.ploidy;
            }
        }
        default_ploidy
    }).collect()
}

/// Parse PLINK-style PED file: family father mother sex affection (6 cols).
/// Returns one TrioFamily per (father, mother, child) trio found.
pub fn parse_ped<P: AsRef<Path>>(p: P, samples: &[String]) -> Result<Vec<TrioFamily>> {
    let idx: HashMap<&str, usize> = samples.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
    let mut out = Vec::new();
    for line in BufReader::new(File::open(p.as_ref())?).lines() {
        let l = line?;
        let t = l.trim();
        if t.is_empty() || t.starts_with('#') { continue; }
        let parts: Vec<&str> = t.split_whitespace().collect();
        if parts.len() < 6 { continue; }
        let (_fam, iid, f, m, sex, _aff) = (parts[0], parts[1], parts[2], parts[3], parts[4], parts[5]);
        if f == "0" || m == "0" { continue; }
        let child = idx.get(iid).copied();
        let father = idx.get(f).copied();
        let mother = idx.get(m).copied();
        if child.is_none() || father.is_none() || mother.is_none() { continue; }
        let is_son = sex == "1";
        out.push(TrioFamily { father, mother, child, is_son });
    }
    Ok(out)
}

/// Parse `--group-samples FILE`: each line "group_name sample_name".
pub fn parse_groups<P: AsRef<Path>>(p: P, samples: &[String]) -> Result<Vec<SampleGroup>> {
    let idx: HashMap<&str, usize> = samples.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
    let mut by_name: HashMap<String, Vec<usize>> = HashMap::new();
    for line in BufReader::new(File::open(p.as_ref())?).lines() {
        let l = line?;
        let t = l.trim();
        if t.is_empty() || t.starts_with('#') { continue; }
        // `sample<TAB>group`, as `bcftools call -G` reads it.
        let parts: Vec<&str> = t.split_whitespace().collect();
        if parts.len() < 2 { continue; }
        let sample = parts[0];
        let group = parts[1].to_string();
        if let Some(&si) = idx.get(sample) {
            by_name.entry(group).or_default().push(si);
        }
    }
    Ok(by_name.into_iter().map(|(name, sample_idxs)| SampleGroup { name, sample_idxs }).collect())
}

/// Parse `--samples-file FILE` sex map from either a `sample sex` table or a PED.
/// The sex label is kept verbatim, since ploidy files use `M`/`F` or `1`/`2`.
pub fn parse_sex_file<P: AsRef<Path>>(p: P) -> Result<HashMap<String, String>> {
    let mut m = HashMap::new();
    for line in BufReader::new(File::open(p.as_ref())?).lines() {
        let l = line?;
        let t = l.trim();
        if t.is_empty() || t.starts_with('#') { continue; }
        let parts: Vec<&str> = t.split_whitespace().collect();
        match parts.len() {
            0 | 1 => {}
            2..=4 => { m.insert(parts[0].to_string(), parts[1].to_string()); }
            // PED: sex 1/2 becomes the M/F the ploidy file is keyed on.
            _ => {
                let sex = match parts[4] { "1" => "M", "2" => "F", other => other };
                m.insert(parts[1].to_string(), sex.to_string());
            }
        }
    }
    Ok(m)
}

/// Parse `--prior-freqs AN_AC` spec: e.g. `AN,AC` (use INFO/AN, INFO/AC) or
/// just `AF` (use INFO/AF directly).
pub fn parse_prior_freqs(spec: &str) -> Result<PriorFreqsSpec> {
    if let Some((an, ac)) = spec.split_once(',') {
        return Ok(PriorFreqsSpec::AnAc { an: an.to_string(), ac: ac.to_string() });
    }
    Ok(PriorFreqsSpec::Af(spec.to_string()))
}

#[derive(Clone, Debug)]
pub enum PriorFreqsSpec {
    Af(String),
    AnAc { an: String, ac: String },
}

impl PriorFreqsSpec {
    pub fn extract_af(&self, info: &str, n_alleles: usize) -> Option<Vec<f64>> {
        match self {
            PriorFreqsSpec::Af(tag) => {
                for kv in info.split(';') {
                    if let Some((k, v)) = kv.split_once('=') {
                        if k == tag {
                            let alts: Vec<f64> = v.split(',').filter_map(|s| s.parse().ok()).collect();
                            let mut af = vec![0.0f64; n_alleles];
                            let alt_sum: f64 = alts.iter().sum();
                            af[0] = (1.0 - alt_sum).max(0.0);
                            for (i, &a) in alts.iter().enumerate() {
                                if i + 1 < n_alleles { af[i + 1] = a; }
                            }
                            return Some(af);
                        }
                    }
                }
                None
            }
            PriorFreqsSpec::AnAc { an, ac } => {
                let mut an_val: u32 = 0;
                let mut ac_vec: Vec<u32> = Vec::new();
                for kv in info.split(';') {
                    if let Some((k, v)) = kv.split_once('=') {
                        if k == an { an_val = v.parse().unwrap_or(0); }
                        else if k == ac { ac_vec = v.split(',').filter_map(|s| s.parse().ok()).collect(); }
                    }
                }
                if an_val == 0 { return None; }
                let mut af = vec![0.0f64; n_alleles];
                let total_alt: u32 = ac_vec.iter().sum();
                af[0] = ((an_val - total_alt) as f64) / (an_val as f64);
                for (i, &c) in ac_vec.iter().enumerate() {
                    if i + 1 < n_alleles { af[i + 1] = (c as f64) / (an_val as f64); }
                }
                Some(af)
            }
        }
    }
}

#[cfg(test)]
#[path = "../../tests/unit/call_pedigree.rs"]
mod tests;
