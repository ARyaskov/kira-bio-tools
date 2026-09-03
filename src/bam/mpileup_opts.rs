//! Mpileup config presets, annotate spec, flag filters, sample filtering.

use anyhow::{Result, bail};

#[derive(Clone, Debug, Default)]
pub struct AnnotateSpec {
    pub fmt_ad: bool,
    pub fmt_dp: bool,
    pub fmt_qs: bool,
    pub fmt_sp: bool,
    pub fmt_adf: bool,
    pub fmt_adr: bool,
    pub fmt_scr: bool,
    pub info_ad: bool,
    pub info_adf: bool,
    pub info_adr: bool,
    pub info_scr: bool,
    pub fmt_pl: bool,
    pub fmt_gq: bool,
    /// FORMAT/DV: non-reference depth (deprecated bcftools tag).
    pub fmt_dv: bool,
    /// FORMAT/DP4: ref-fwd, ref-rev, alt-fwd, alt-rev.
    pub fmt_dp4: bool,
    pub info_dp4: bool,
}

impl AnnotateSpec {
    /// bcftools syntax: comma list of tags, `FORMAT/`/`FMT/`/`INFO/` prefixes,
    /// a leading `-` removes a tag; unknown tags are reported and ignored.
    pub fn parse(spec: Option<&str>) -> Result<Self> {
        let mut s = Self { fmt_ad: true, fmt_dp: true, fmt_pl: true, ..Default::default() };
        let Some(spec) = spec else { return Ok(s); };
        s = Self::default();
        s.fmt_pl = true;
        for tok in spec.split(',') {
            let t = tok.trim();
            let (on, name) = match t.strip_prefix('-') {
                Some(n) => (false, n),
                None => (true, t),
            };
            let upper = name.to_uppercase();
            let flag: Option<&mut bool> = match upper.as_str() {
                "FORMAT/AD" | "FMT/AD" | "AD" | "FORMAT/DPR" | "FMT/DPR" | "DPR" => Some(&mut s.fmt_ad),
                "FORMAT/DP" | "FMT/DP" | "DP" => Some(&mut s.fmt_dp),
                "FORMAT/QS" | "FMT/QS" | "QS" => Some(&mut s.fmt_qs),
                "FORMAT/SP" | "FMT/SP" | "SP" => Some(&mut s.fmt_sp),
                "FORMAT/ADF" | "FMT/ADF" | "ADF" => Some(&mut s.fmt_adf),
                "FORMAT/ADR" | "FMT/ADR" | "ADR" => Some(&mut s.fmt_adr),
                "FORMAT/SCR" | "FMT/SCR" | "SCR" => Some(&mut s.fmt_scr),
                "FORMAT/PL" | "FMT/PL" | "PL" => Some(&mut s.fmt_pl),
                "FORMAT/GQ" | "FMT/GQ" | "GQ" => Some(&mut s.fmt_gq),
                "FORMAT/DV" | "FMT/DV" | "DV" => Some(&mut s.fmt_dv),
                "FORMAT/DP4" | "FMT/DP4" | "DP4" => Some(&mut s.fmt_dp4),
                "INFO/AD" | "INFO/DPR" => Some(&mut s.info_ad),
                "INFO/ADF" => Some(&mut s.info_adf),
                "INFO/ADR" => Some(&mut s.info_adr),
                "INFO/SCR" => Some(&mut s.info_scr),
                "INFO/DP4" => Some(&mut s.info_dp4),
                "" => None,
                other => {
                    if other.starts_with("INFO/") || other.starts_with("FORMAT/") || other.starts_with("FMT/") {
                        eprintln!("[mpileup] warning: annotation {t:?} is not supported and is ignored");
                        None
                    } else {
                        bail!("--annotate: unknown tag {t:?}");
                    }
                }
            };
            if let Some(f) = flag { *f = on; }
        }
        Ok(s)
    }
}

/// Preset configurations matching bcftools `-X` flag.
#[derive(Clone, Debug)]
pub struct PresetConfig {
    pub min_mq: Option<u32>,
    pub min_bq: Option<u32>,
    pub no_baq: Option<bool>,
    pub indel_size: Option<u32>,
    pub max_depth: Option<u32>,
    pub gap_frac: Option<f64>,
    pub tandem_qual: Option<u32>,
}

impl PresetConfig {
    pub fn parse(name: &str) -> Result<Self> {
        Ok(match name {
            "illumina-1.18" | "illumina" => Self {
                min_mq: Some(1), min_bq: Some(13), no_baq: Some(false),
                indel_size: Some(110), max_depth: Some(250),
                gap_frac: Some(0.002), tandem_qual: Some(500),
            },
            "ont" | "nanopore" => Self {
                min_mq: Some(7), min_bq: Some(5), no_baq: Some(true),
                indel_size: Some(110), max_depth: Some(250),
                gap_frac: Some(0.05), tandem_qual: Some(100),
            },
            "pacbio-ccs" | "pacbio" | "hifi" => Self {
                min_mq: Some(10), min_bq: Some(13), no_baq: Some(true),
                indel_size: Some(110), max_depth: Some(250),
                gap_frac: Some(0.005), tandem_qual: Some(500),
            },
            other => bail!("-X: unknown preset {other:?}; expected illumina-1.18|ont|pacbio-ccs"),
        })
    }
}

/// SAM flag filters (per-read). Maps to bcftools `--ef/--df/--if/--nf`.
#[derive(Clone, Debug, Default)]
pub struct FlagFilters {
    pub require_flags: u16,
    pub exclude_flags: u16,
    /// `--skip-all-set`: drop reads with every listed bit set.
    pub skip_all_set: u16,
    /// `--skip-all-unset`: drop reads with every listed bit unset.
    pub skip_all_unset: u16,
    pub include_proper_pair_only: bool,
    pub exclude_proper_pair: bool,
}

impl FlagFilters {
    pub fn from_args(ef: Option<u32>, df: Option<u32>, if_: Option<u32>, nf: Option<u32>) -> Self {
        Self {
            exclude_flags: ef.unwrap_or(0) as u16 | df.unwrap_or(0) as u16,
            require_flags: if_.unwrap_or(0) as u16 | nf.unwrap_or(0) as u16,
            skip_all_set: 0,
            skip_all_unset: 0,
            include_proper_pair_only: false,
            exclude_proper_pair: false,
        }
    }

    pub fn from_full(ef: Option<u32>, df: Option<u32>, if_: Option<u32>, nf: Option<u32>,
                     rf: Option<&str>, ff: Option<&str>) -> Result<Self> {
        let mut s = Self::from_args(ef, df, if_, nf);
        if let Some(t) = rf { s.require_flags |= parse_sam_flags(t)?; }
        if let Some(t) = ff { s.exclude_flags |= parse_sam_flags(t)?; }
        Ok(s)
    }

    pub fn with_skip_all(mut self, all_set: Option<&str>, all_unset: Option<&str>) -> Result<Self> {
        if let Some(t) = all_set { self.skip_all_set |= parse_sam_flags(t)?; }
        if let Some(t) = all_unset { self.skip_all_unset |= parse_sam_flags(t)?; }
        Ok(self)
    }

    pub fn is_active(&self) -> bool {
        self.require_flags != 0 || self.exclude_flags != 0 || self.skip_all_set != 0 || self.skip_all_unset != 0
    }

    pub fn passes(&self, flags: u16) -> bool {
        if self.require_flags != 0 && (flags & self.require_flags) != self.require_flags { return false; }
        if self.exclude_flags != 0 && (flags & self.exclude_flags) != 0 { return false; }
        if self.skip_all_set != 0 && (flags & self.skip_all_set) == self.skip_all_set { return false; }
        if self.skip_all_unset != 0 && (flags & self.skip_all_unset) == 0 { return false; }
        true
    }
}

pub fn parse_sam_flags(s: &str) -> Result<u16> {
    let t = s.trim();
    if let Some(rest) = t.strip_prefix("0x").or_else(|| t.strip_prefix("0X")) {
        return Ok(u16::from_str_radix(rest, 16)?);
    }
    if t.chars().all(|c| c.is_ascii_digit()) {
        return Ok(t.parse::<u16>()?);
    }
    let mut acc: u16 = 0;
    for tok in t.split(',') {
        acc |= match tok.trim().to_ascii_uppercase().as_str() {
            "PAIRED" => 0x1,
            "PROPER_PAIR" | "PROPERPAIR" => 0x2,
            "UNMAP" | "UNMAPPED" => 0x4,
            "MUNMAP" | "MATE_UNMAPPED" => 0x8,
            "REVERSE" | "REV" => 0x10,
            "MREVERSE" | "MATE_REVERSE" => 0x20,
            "READ1" | "FIRST" => 0x40,
            "READ2" | "SECOND" => 0x80,
            "SECONDARY" => 0x100,
            "QCFAIL" | "FAIL" => 0x200,
            "DUP" | "DUPLICATE" => 0x400,
            "SUPPLEMENTARY" | "SUPP" => 0x800,
            "" => 0,
            other => bail!("unknown SAM flag {other:?}"),
        };
    }
    Ok(acc)
}

/// Parse `-s LIST` (comma-list) or `-S FILE` (one per line) to filter samples by name.
pub fn parse_samples_filter(spec: Option<&str>, file: Option<&std::path::Path>) -> Result<Option<Vec<String>>> {
    let mut names: Vec<String> = Vec::new();
    if let Some(s) = spec {
        for tok in s.split(',') {
            let t = tok.trim();
            if !t.is_empty() { names.push(t.to_string()); }
        }
    }
    if let Some(p) = file {
        use std::io::BufRead;
        let f = std::fs::File::open(p)?;
        for line in std::io::BufReader::new(f).lines() {
            let l = line?;
            let t = l.trim();
            if !t.is_empty() && !t.starts_with('#') { names.push(t.to_string()); }
        }
    }
    if names.is_empty() { Ok(None) } else { Ok(Some(names)) }
}

#[cfg(test)]
#[path = "../../tests/unit/bam_mpileup_opts.rs"]
mod tests;
