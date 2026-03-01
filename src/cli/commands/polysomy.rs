use anyhow::Result;
use std::collections::BTreeMap;

use crate::VcfReader;
use crate::cli::args::PolysomyArgs;

pub fn cmd_polysomy(args: PolysomyArgs) -> Result<()> {
    let cfg = parse_args(&args.bcftools_args);
    let mut r = VcfReader::open(&args.input)?;
    let headers = r.header()?;
    let sample_names = parse_sample_names(&headers);
    let sample_idx = cfg
        .sample_name
        .as_deref()
        .and_then(|n| sample_names.iter().position(|s| s == n))
        .unwrap_or(0);

    let mut by_chrom = BTreeMap::<String, ChromStats>::new();

    while let Some(rec) = r.next_record()? {
        if rec.alt == "." {
            continue;
        }
        if !cfg.region.contains(&rec.chrom, rec.pos) {
            continue;
        }
        let af = site_alt_fraction(&rec, sample_idx).unwrap_or(0.0);
        let maf = af.min(1.0 - af);
        if maf < cfg.min_minor_af {
            continue;
        }

        let st = by_chrom.entry(rec.chrom.clone()).or_default();
        st.n_sites += 1;
        st.alt_sum += af;
        if af < cfg.baf_low || af > cfg.baf_high {
            st.abnormal_baf += 1;
        }

        match gt_class(&rec, sample_idx) {
            GtClass::HomRef => st.hom_ref += 1,
            GtClass::Het => st.het += 1,
            GtClass::HomAlt => st.hom_alt += 1,
            GtClass::Missing => st.missing += 1,
        }
    }

    println!("#CHROM\tN_SITES\tHET\tHOM_REF\tHOM_ALT\tMEAN_AF\tABN_FRAC\tHET_RATE\tSCORE\tSTATUS");
    for (chrom, st) in by_chrom {
        let mean = if st.n_sites > 0 {
            st.alt_sum / st.n_sites as f64
        } else {
            0.0
        };
        let abn_frac = if st.n_sites > 0 {
            st.abnormal_baf as f64 / st.n_sites as f64
        } else {
            0.0
        };
        let called = st.hom_ref + st.het + st.hom_alt;
        let het_rate = if called > 0 {
            st.het as f64 / called as f64
        } else {
            0.0
        };
        let score = if st.n_sites > 0 {
            ((abn_frac * 1.4) + (mean - 0.5).abs() + (0.5 - het_rate).abs())
                * (st.n_sites as f64).sqrt()
        } else {
            0.0
        };

        let status = if st.n_sites < cfg.min_sites {
            "LOW_DATA"
        } else if abn_frac >= cfg.polysomy_frac
            || mean < cfg.baf_low
            || mean > cfg.baf_high
            || score > 3.0
        {
            "POLYSOMY_CANDIDATE"
        } else {
            "DIPLOID_LIKE"
        };

        if cfg.verbosity == 0 && status == "LOW_DATA" {
            continue;
        }

        if cfg.include_intermediate {
            println!(
                "{chrom}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.3}\t{}",
                st.n_sites, st.het, st.hom_ref, st.hom_alt, mean, abn_frac, het_rate, score, status
            );
        } else {
            println!(
                "{chrom}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.3}\t{}",
                st.n_sites, st.het, st.hom_ref, st.hom_alt, mean, abn_frac, het_rate, score, status
            );
        }
    }

    Ok(())
}

#[derive(Default)]
struct ChromStats {
    n_sites: usize,
    het: usize,
    hom_ref: usize,
    hom_alt: usize,
    missing: usize,
    alt_sum: f64,
    abnormal_baf: usize,
}

struct PolyCfg {
    min_minor_af: f64,
    min_sites: usize,
    verbosity: u8,
    sample_name: Option<String>,
    region: Region,
    baf_low: f64,
    baf_high: f64,
    polysomy_frac: f64,
    include_intermediate: bool,
}

#[derive(Clone)]
struct Region {
    chrom: Option<String>,
    start: u32,
    end: u32,
}

impl Region {
    fn all() -> Self {
        Self {
            chrom: None,
            start: 1,
            end: u32::MAX,
        }
    }

    fn contains(&self, chrom: &str, pos: u32) -> bool {
        if let Some(c) = &self.chrom {
            if c != chrom {
                return false;
            }
        }
        pos >= self.start && pos <= self.end
    }
}

fn parse_args(args: &[String]) -> PolyCfg {
    let mut cfg = PolyCfg {
        min_minor_af: 0.0,
        min_sites: 10,
        verbosity: 1,
        sample_name: None,
        region: Region::all(),
        baf_low: 0.35,
        baf_high: 0.65,
        polysomy_frac: 0.6,
        include_intermediate: false,
    };
    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "-s" => {
                i += 1;
                if let Some(v) = args.get(i) {
                    cfg.sample_name = Some(v.clone());
                }
            }
            "-r" => {
                i += 1;
                if let Some(v) = args.get(i) {
                    cfg.region = parse_region(v);
                }
            }
            "-m" => {
                i += 1;
                if let Some(v) = args.get(i).and_then(|x| x.parse::<f64>().ok()) {
                    cfg.min_minor_af = v.clamp(0.0, 0.5);
                }
            }
            "-f" => {
                i += 1;
                if let Some(v) = args.get(i).and_then(|x| x.parse::<f64>().ok()) {
                    cfg.min_minor_af = cfg.min_minor_af.max(v.clamp(0.0, 0.5));
                }
            }
            "-b" => {
                i += 1;
                if let Some(v) = args.get(i).and_then(|x| x.parse::<f64>().ok()) {
                    cfg.baf_low = v.clamp(0.0, 1.0);
                }
            }
            "-c" => {
                i += 1;
                if let Some(v) = args.get(i).and_then(|x| x.parse::<f64>().ok()) {
                    cfg.baf_high = v.clamp(0.0, 1.0);
                }
            }
            "-p" => {
                i += 1;
                if let Some(v) = args.get(i).and_then(|x| x.parse::<f64>().ok()) {
                    cfg.polysomy_frac = v.clamp(0.0, 1.0);
                }
            }
            "-i" => {
                cfg.include_intermediate = true;
            }
            "-v" => {
                i += 1;
                if let Some(v) = args.get(i).and_then(|x| x.parse::<u8>().ok()) {
                    cfg.verbosity = v;
                }
            }
            _ => {}
        }
        i += 1;
    }
    cfg
}

enum GtClass {
    HomRef,
    Het,
    HomAlt,
    Missing,
}

fn site_alt_fraction(rec: &crate::vcf::structs::VcfRecord, sample_idx: usize) -> Option<f64> {
    let fmt = rec.format.as_deref()?;
    let keys = fmt.split(':').collect::<Vec<_>>();
    let sample = rec.samples.get(sample_idx)?.as_str();
    let vals = sample.split(':').collect::<Vec<_>>();

    if let Some(ad_i) = keys.iter().position(|k| *k == "AD") {
        if let Some(ad) = vals.get(ad_i) {
            let parts = ad
                .split(',')
                .filter_map(|x| x.parse::<f64>().ok())
                .collect::<Vec<_>>();
            if parts.len() >= 2 {
                let ref_dp = parts[0];
                let alt_dp = parts[1..].iter().sum::<f64>();
                let d = ref_dp + alt_dp;
                if d > 0.0 {
                    return Some(alt_dp / d);
                }
            }
        }
    }

    if let Some(gt_i) = keys.iter().position(|k| *k == "GT") {
        if let Some(gt) = vals.get(gt_i) {
            return gt_to_fraction(gt);
        }
    }

    None
}

fn gt_class(rec: &crate::vcf::structs::VcfRecord, sample_idx: usize) -> GtClass {
    let Some(fmt) = rec.format.as_deref() else {
        return GtClass::Missing;
    };
    let keys = fmt.split(':').collect::<Vec<_>>();
    let Some(gt_i) = keys.iter().position(|k| *k == "GT") else {
        return GtClass::Missing;
    };
    let Some(sample) = rec.samples.get(sample_idx) else {
        return GtClass::Missing;
    };
    let vals = sample.split(':').collect::<Vec<_>>();
    let Some(gt) = vals.get(gt_i) else {
        return GtClass::Missing;
    };

    match gt_to_fraction(gt) {
        Some(x) if x <= 0.01 => GtClass::HomRef,
        Some(x) if x >= 0.99 => GtClass::HomAlt,
        Some(_) => GtClass::Het,
        None => GtClass::Missing,
    }
}

fn parse_sample_names(headers: &[String]) -> Vec<String> {
    headers
        .iter()
        .find(|h| h.starts_with("#CHROM\t"))
        .map(|h| {
            h.split('\t')
                .skip(9)
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
        })
        .unwrap_or_default()
}

fn parse_region(s: &str) -> Region {
    if let Some((chrom, rest)) = s.split_once(':') {
        if let Some((a, b)) = rest.split_once('-') {
            let start = a.parse::<u32>().unwrap_or(1);
            let end = b.parse::<u32>().unwrap_or(u32::MAX);
            return Region {
                chrom: Some(chrom.to_string()),
                start,
                end,
            };
        }
        return Region {
            chrom: Some(chrom.to_string()),
            start: 1,
            end: u32::MAX,
        };
    }
    Region {
        chrom: Some(s.to_string()),
        start: 1,
        end: u32::MAX,
    }
}

fn gt_to_fraction(gt: &str) -> Option<f64> {
    let alleles = gt.split(['/', '|']).collect::<Vec<_>>();
    if alleles.is_empty() {
        return None;
    }
    let mut alt = 0f64;
    let mut n = 0f64;
    for a in alleles {
        if a == "." {
            continue;
        }
        if let Ok(v) = a.parse::<i32>() {
            n += 1.0;
            if v > 0 {
                alt += 1.0;
            }
        }
    }
    if n == 0.0 { None } else { Some(alt / n) }
}
