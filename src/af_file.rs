//! `--AF-file` of `roh`/`cnv`: a `CHR\tPOS\tREF,ALT\tAF` table (plain or
//! bgzipped) that acts as a targets file (sites absent from it are skipped)
//! and supplies the allele frequency when the alleles match.

use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use anyhow::{Result, bail};

pub struct AfFile {
    fname: String,
    /// First line seen per `(chrom, 1-based pos)`: alleles, the raw AF
    /// column (parsed only once the alleles match, as bcftools does) and the
    /// raw line for error messages.
    sites: HashMap<(String, u32), (Vec<String>, String, String)>,
}

impl AfFile {
    pub fn load(path: &Path) -> Result<Self> {
        let mut rd = crate::util::open_text_maybe_gz(path)?;
        let mut sites = HashMap::new();
        let mut line = String::new();
        loop {
            line.clear();
            if rd.read_line(&mut line)? == 0 {
                break;
            }
            let l = line.trim_end_matches(['\n', '\r']);
            if l.is_empty() || l.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = l.split('\t').collect();
            if cols.len() < 3 {
                bail!("Could not parse the targets file {}: {}", path.display(), l);
            }
            let Ok(pos) = cols[1].trim().parse::<u32>() else {
                bail!("Could not parse the targets file {}: {}", path.display(), l);
            };
            let alleles: Vec<String> = cols[2].split(',').map(str::to_string).collect();
            let af = cols.get(3).map(|s| s.to_string()).unwrap_or_default();
            sites
                .entry((cols[0].to_string(), pos))
                .or_insert((alleles, af, l.to_string()));
        }
        Ok(Self {
            fname: path.display().to_string(),
            sites,
        })
    }

    /// Targets semantics: is the site listed at all?
    pub fn has_site(&self, chrom: &str, pos1: u32) -> bool {
        self.sites.contains_key(&(chrom.to_string(), pos1))
    }

    /// `read_AF`: `Ok(None)` when the alleles differ or the AF is missing.
    pub fn read_af(&self, chrom: &str, pos1: u32, alleles: &[String]) -> Result<Option<f64>> {
        let Some((als, af, line)) = self.sites.get(&(chrom.to_string(), pos1)) else {
            return Ok(None);
        };
        if als.len() < 2 {
            bail!(
                "Expected two comma-separated alleles (REF,ALT) in the third column of {}, found:\n\t{}",
                self.fname,
                line
            );
        }
        if als.len() != alleles.len() {
            return Ok(None);
        }
        if als.iter().zip(alleles).any(|(a, b)| a != b) {
            return Ok(None);
        }
        // strtod: leading number, then the rest must be blank.
        let b = af.as_bytes();
        let mut end = 0;
        while end < b.len()
            && (b[end].is_ascii_digit() || matches!(b[end], b'.' | b'-' | b'+' | b'e' | b'E'))
        {
            end += 1;
        }
        let head = &af[..end];
        let rest = af[end..].chars().next();
        let value = head.parse::<f64>().ok();
        match (value, rest) {
            (Some(v), None) | (Some(v), Some(' ')) | (Some(v), Some('\t')) => {
                if !(0.0..=1.0).contains(&v) {
                    bail!("Could not parse AF: [{line}]");
                }
                Ok(Some(v))
            }
            _ => {
                if b.first() == Some(&b'.') && (b.len() == 1 || b[1].is_ascii_whitespace()) {
                    return Ok(None);
                }
                bail!("Could not parse: [{line}]")
            }
        }
    }
}
