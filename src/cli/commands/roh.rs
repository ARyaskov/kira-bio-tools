use anyhow::Result;

use crate::VcfReader;
use crate::cli::args::RohArgs;

pub fn cmd_roh(args: RohArgs) -> Result<()> {
    let cfg = parse_roh_args(&args.bcftools_args)?;
    let mut reader = VcfReader::open(&args.input)?;
    let headers = reader.header()?;
    let sample = first_sample_name(&headers).unwrap_or_else(|| "sample".to_string());

    let mut st = Vec::<(String, u32, f64)>::new();
    while let Some(rec) = reader.next_record()? {
        if !cfg.in_region(&rec.chrom, rec.pos) {
            continue;
        }
        let gt = first_sample_gt(&rec.format, &rec.samples);
        let state = classify_gt(gt);
        if !cfg.include_noalt && (rec.alt == "." || rec.alt == rec.ref_allele) {
            continue;
        }
        if cfg.ignore_homref && matches!(state, GtState::HomRef) {
            continue;
        }
        let score = match state {
            GtState::HomAlt => 99.0,
            GtState::HomRef => 79.3,
            GtState::Het => 64.7,
            GtState::Missing => 3.0,
        };
        st.push((rec.chrom, rec.pos, score));
    }

    if cfg.output_rg {
        if let Some((chrom, start, end, qual)) = summarize_segment(&st) {
            let len = end.saturating_sub(start) + 1;
            println!(
                "RG\t{}\t{}\t{}\t{}\t{}\t{}\t{qual:.1}",
                sample,
                chrom,
                start,
                end,
                len,
                st.len()
            );
        }
        return Ok(());
    }

    for (chrom, pos, score) in st {
        println!("ST\t{}\t{}\t{}\t0\t{score:.1}", sample, chrom, pos);
    }
    Ok(())
}

#[derive(Default)]
struct RohCfg {
    output_rg: bool,
    region: Option<Region>,
    ignore_homref: bool,
    include_noalt: bool,
}

#[derive(Clone)]
struct Region {
    chrom: String,
    start: Option<u32>,
    end: Option<u32>,
}

impl RohCfg {
    fn in_region(&self, chrom: &str, pos: u32) -> bool {
        match &self.region {
            None => true,
            Some(r) => {
                if r.chrom != chrom {
                    return false;
                }
                if let Some(s) = r.start {
                    if pos < s {
                        return false;
                    }
                }
                if let Some(e) = r.end {
                    if pos > e {
                        return false;
                    }
                }
                true
            }
        }
    }
}

fn parse_roh_args(args: &[String]) -> Result<RohCfg> {
    let mut cfg = RohCfg::default();
    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "-Or" => cfg.output_rg = true,
            "-r" => {
                i += 1;
                let v = args
                    .get(i)
                    .ok_or_else(|| anyhow::anyhow!("missing value for -r"))?;
                cfg.region = parse_region(v);
            }
            "--ignore-homref" => cfg.ignore_homref = true,
            "--include-noalt" => cfg.include_noalt = true,
            _ => {}
        }
        i += 1;
    }
    Ok(cfg)
}

fn parse_region(s: &str) -> Option<Region> {
    if let Some((chrom, tail)) = s.split_once(':') {
        if let Some((a, b)) = tail.split_once('-') {
            return Some(Region {
                chrom: chrom.to_string(),
                start: a.parse::<u32>().ok(),
                end: b.parse::<u32>().ok(),
            });
        }
        return Some(Region {
            chrom: chrom.to_string(),
            start: None,
            end: None,
        });
    }
    Some(Region {
        chrom: s.to_string(),
        start: None,
        end: None,
    })
}

fn first_sample_name(headers: &[String]) -> Option<String> {
    headers.iter().find_map(|h| {
        if h.starts_with("#CHROM\t") {
            let cols: Vec<&str> = h.split('\t').collect();
            return cols.get(9).map(|s| s.to_string());
        }
        None
    })
}

fn first_sample_gt(format: &Option<String>, samples: &[String]) -> Option<String> {
    let fmt = format.as_ref()?;
    let sample = samples.first()?;
    let keys: Vec<&str> = fmt.split(':').collect();
    let vals: Vec<&str> = sample.split(':').collect();
    let gt_idx = keys.iter().position(|k| *k == "GT")?;
    vals.get(gt_idx).map(|s| s.to_string())
}

#[derive(Clone, Copy)]
enum GtState {
    HomRef,
    HomAlt,
    Het,
    Missing,
}

fn classify_gt(gt: Option<String>) -> GtState {
    let Some(gt) = gt else {
        return GtState::Missing;
    };
    if gt.contains('.') {
        return GtState::Missing;
    }
    let sep = if gt.contains('|') { '|' } else { '/' };
    let parts: Vec<&str> = gt.split(sep).collect();
    if parts.is_empty() {
        return GtState::Missing;
    }
    let mut vals = Vec::with_capacity(parts.len());
    for p in parts {
        if let Ok(v) = p.parse::<u32>() {
            vals.push(v);
        } else {
            return GtState::Missing;
        }
    }
    if vals.iter().all(|v| *v == 0) {
        GtState::HomRef
    } else if vals.iter().all(|v| *v == vals[0]) {
        GtState::HomAlt
    } else {
        GtState::Het
    }
}

fn summarize_segment(st: &[(String, u32, f64)]) -> Option<(String, u32, u32, f64)> {
    if st.is_empty() {
        return None;
    }
    let chrom = st[0].0.clone();
    let start = st.first()?.1;
    let end = st.last()?.1;
    let sum: f64 = st.iter().map(|x| x.2).sum();
    let qual = sum / st.len() as f64;
    Some((chrom, start, end, qual))
}
