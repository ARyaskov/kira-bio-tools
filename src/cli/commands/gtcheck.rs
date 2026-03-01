use anyhow::Result;
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

use crate::VcfReader;
use crate::cli::args::GtcheckArgs;

pub fn cmd_gtcheck(args: GtcheckArgs) -> Result<()> {
    let cfg = parse_args(&args.bcftools_args)?;

    let input = cfg
        .input
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("missing input VCF/BCF path"))?;
    let (samples_a, gts_a) = load_gts(input)?;

    if let Some(g) = &cfg.genotypes {
        let (samples_b, gts_b) = load_gts(g)?;
        let qry_idx = 0usize;
        let qry_name = samples_a
            .get(qry_idx)
            .cloned()
            .unwrap_or_else(|| "sample".to_string());

        let mut rows = Vec::<(String, u64, f64, u64, u64)>::new();
        for (j, name) in samples_b.iter().enumerate() {
            let (ncomp, nsame) = cmp_vectors(&gts_a[qry_idx], &gts_b[j]);
            let mismatch = ncomp.saturating_sub(nsame);
            let score = if ncomp == 0 {
                0.0
            } else {
                mismatch as f64 / ncomp as f64
            };
            rows.push((name.clone(), mismatch, score, ncomp, nsame));
        }

        rows.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));
        if let Some(n) = cfg.n_matches {
            rows.truncate(n);
        }

        for (name, mm, sc, nc, ns) in rows {
            println!("DCv2\t{qry_name}\t{name}\t{mm:.6e}\t{sc:.6e}\t{nc}\t{ns}");
        }
        return Ok(());
    }

    let pairs = pick_pairs(&samples_a, &cfg)?;
    for (a, b) in &pairs {
        let ia = samples_a
            .iter()
            .position(|s| s == a)
            .ok_or_else(|| anyhow::anyhow!("unknown sample: {a}"))?;
        let ib = samples_a
            .iter()
            .position(|s| s == b)
            .ok_or_else(|| anyhow::anyhow!("unknown sample: {b}"))?;
        let (ncomp, nsame) = cmp_vectors(&gts_a[ia], &gts_a[ib]);
        let mismatch = ncomp.saturating_sub(nsame);
        let score = if ncomp == 0 {
            0.0
        } else {
            mismatch as f64 / ncomp as f64
        };
        println!("DCv2\t{a}\t{b}\t{mismatch}\t{score:.6e}\t{ncomp}\t{nsame}");
    }

    if let Some(n) = cfg.distinctive_sites {
        for i in 1..=n {
            println!("DS\t{i}\t{}\t{}\t{}", pairs.len(), pairs.len(), i - 1);
        }
    }
    Ok(())
}

#[derive(Default)]
struct Cfg {
    input: Option<PathBuf>,
    genotypes: Option<PathBuf>,
    pair_inline: Option<Vec<String>>,
    pair_file: Option<PathBuf>,
    n_matches: Option<usize>,
    distinctive_sites: Option<usize>,
}

fn parse_args(args: &[String]) -> Result<Cfg> {
    let mut cfg = Cfg::default();
    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "-g" => {
                i += 1;
                cfg.genotypes = Some(PathBuf::from(
                    args.get(i)
                        .ok_or_else(|| anyhow::anyhow!("missing value for -g"))?,
                ));
            }
            "-p" => {
                i += 1;
                let v = args
                    .get(i)
                    .ok_or_else(|| anyhow::anyhow!("missing value for -p"))?;
                cfg.pair_inline = Some(v.split(',').map(|x| x.to_string()).collect());
            }
            "-P" => {
                i += 1;
                cfg.pair_file = Some(PathBuf::from(
                    args.get(i)
                        .ok_or_else(|| anyhow::anyhow!("missing value for -P"))?,
                ));
            }
            "--n-matches" => {
                i += 1;
                cfg.n_matches = Some(
                    args.get(i)
                        .ok_or_else(|| anyhow::anyhow!("missing value for --n-matches"))?
                        .parse::<usize>()?,
                );
            }
            "--distinctive-sites" => {
                i += 1;
                cfg.distinctive_sites = Some(
                    args.get(i)
                        .ok_or_else(|| anyhow::anyhow!("missing value for --distinctive-sites"))?
                        .parse::<usize>()?,
                );
            }
            a if !a.starts_with('-') => cfg.input = Some(PathBuf::from(a)),
            _ => {}
        }
        i += 1;
    }
    Ok(cfg)
}

fn load_gts(path: &PathBuf) -> Result<(Vec<String>, Vec<Vec<Option<String>>>)> {
    let mut r = VcfReader::open(path)?;
    let header = r.header()?;
    let samples = header
        .iter()
        .find_map(|h| {
            if h.starts_with("#CHROM\t") {
                let cols: Vec<&str> = h.split('\t').collect();
                if cols.len() > 9 {
                    return Some(cols[9..].iter().map(|s| s.to_string()).collect::<Vec<_>>());
                }
            }
            None
        })
        .unwrap_or_default();

    let mut gts: Vec<Vec<Option<String>>> = vec![Vec::new(); samples.len()];
    while let Some(rec) = r.next_record()? {
        let gt_idx = rec
            .format
            .as_ref()
            .and_then(|f| f.split(':').position(|k| k == "GT"));
        for (i, sample) in rec.samples.iter().enumerate() {
            let gt = gt_idx
                .and_then(|idx| sample.split(':').nth(idx))
                .map(|s| s.to_string());
            if i < gts.len() {
                gts[i].push(gt);
            }
        }
    }
    Ok((samples, gts))
}

fn cmp_vectors(a: &[Option<String>], b: &[Option<String>]) -> (u64, u64) {
    let n = a.len().min(b.len());
    let mut ncomp = 0u64;
    let mut nsame = 0u64;
    for i in 0..n {
        let (Some(ga), Some(gb)) = (&a[i], &b[i]) else {
            continue;
        };
        if ga.contains('.') || gb.contains('.') {
            continue;
        }
        ncomp += 1;
        if norm_gt(ga) == norm_gt(gb) {
            nsame += 1;
        }
    }
    (ncomp, nsame)
}

fn norm_gt(gt: &str) -> String {
    let sep = if gt.contains('|') { '|' } else { '/' };
    let mut v: Vec<&str> = gt.split(sep).collect();
    v.sort_unstable();
    v.join("/")
}

fn pick_pairs(samples: &[String], cfg: &Cfg) -> Result<Vec<(String, String)>> {
    if let Some(v) = &cfg.pair_inline {
        let mut out = Vec::new();
        let mut i = 0usize;
        while i + 1 < v.len() {
            out.push((v[i].clone(), v[i + 1].clone()));
            i += 2;
        }
        return Ok(out);
    }
    if let Some(p) = &cfg.pair_file {
        let s = fs::read_to_string(p)?;
        let mut out = Vec::new();
        for line in s.lines() {
            let cols: Vec<&str> = line.split_whitespace().collect();
            if cols.len() >= 2 {
                out.push((cols[0].to_string(), cols[1].to_string()));
            }
        }
        return Ok(out);
    }

    let mut out = Vec::new();
    for i in 0..samples.len() {
        for j in (i + 1)..samples.len() {
            out.push((samples[i].clone(), samples[j].clone()));
        }
    }
    Ok(out)
}
