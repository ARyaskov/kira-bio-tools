use anyhow::{Context, Result, anyhow};
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};

use crate::VcfReader;

#[derive(Clone, Debug, Default)]
pub struct CnvConfig {
    pub input: PathBuf,
    pub output_dir: PathBuf,
    pub query_sample: String,
    pub control_sample: Option<String>,
    pub af_file: Option<PathBuf>,
    pub regions: Vec<RegionSpec>,
    pub targets: Vec<RegionSpec>,
    pub regions_overlap: u8,
    pub targets_overlap: u8,
    pub baf_weight: f64,
    pub lrr_weight: f64,
    pub baf_dev_query: f64,
    pub baf_dev_control: f64,
    pub lrr_dev_query: f64,
    pub lrr_dev_control: f64,
    pub optimize: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct RegionSpec {
    chr: String,
    start: Option<u32>,
    end: Option<u32>,
}

/// CNV states for the Viterbi HMM.
/// 0 = full loss (CN0), 1 = hemizygous loss (CN1), 2 = diploid (CN2),
/// 3 = duplication (CN3), 4 = amplification (CN4+).
const N_STATES: usize = 5;

#[derive(Clone, Debug)]
struct SamplePoint {
    chr: String,
    pos: u32,
    baf: f64,
    lrr: f64,
    cn: u8,
    probs: [f64; N_STATES],
    is_het: bool,
}

#[derive(Clone, Debug)]
struct SampleSummaryRow {
    chr: String,
    start: u32,
    end: u32,
    cn: u8,
    quality: f64,
    n_sites: usize,
    n_hets: usize,
}

pub fn run_from_args(args: &[String]) -> Result<()> {
    let cfg = parse_args(args)?;
    run(cfg)
}

pub fn run(cfg: CnvConfig) -> Result<()> {
    fs::create_dir_all(&cfg.output_dir)
        .with_context(|| format!("failed to create output dir {:?}", cfg.output_dir))?;

    let af_map = if let Some(path) = &cfg.af_file {
        load_af_file(path)?
    } else {
        HashMap::new()
    };

    let mut reader = VcfReader::open(&cfg.input)?;
    let header = reader.header()?;
    let sample_names = extract_sample_names(&header)?;
    let query_idx = sample_names
        .iter()
        .position(|s| s == &cfg.query_sample)
        .ok_or_else(|| anyhow!("query sample not found: {}", cfg.query_sample))?;
    let control_idx = cfg
        .control_sample
        .as_ref()
        .map(|name| {
            sample_names
                .iter()
                .position(|s| s == name)
                .ok_or_else(|| anyhow!("control sample not found: {}", name))
        })
        .transpose()?;

    let mut query_points = Vec::new();
    let mut control_points = Vec::new();

    while let Some(rec) = reader.next_record()? {
        if !record_matches_regions(&rec.chrom, rec.pos, &cfg.regions)
            || !record_matches_regions(&rec.chrom, rec.pos, &cfg.targets)
        {
            continue;
        }

        let format = match &rec.format {
            Some(f) => f,
            None => continue,
        };
        let fmt_index = build_format_index(format);
        let baf_i = match fmt_index.get("BAF") {
            Some(i) => *i,
            None => continue,
        };
        let lrr_i = match fmt_index.get("LRR") {
            Some(i) => *i,
            None => continue,
        };
        let gt_i = fmt_index.get("GT").copied();

        let af = read_af_value(&rec.info).or_else(|| {
            af_map
                .get(&(
                    rec.chrom.clone(),
                    rec.pos,
                    rec.ref_allele.clone(),
                    rec.alt.clone(),
                ))
                .copied()
        });

        if let Some(point) = parse_sample_point(
            &rec.chrom,
            rec.pos,
            rec.samples.get(query_idx),
            gt_i,
            baf_i,
            lrr_i,
            af,
            cfg.baf_weight,
            cfg.lrr_weight,
            cfg.baf_dev_query,
            cfg.lrr_dev_query,
            cfg.optimize,
        ) {
            query_points.push(point);
        }

        if let Some(ci) = control_idx {
            if let Some(point) = parse_sample_point(
                &rec.chrom,
                rec.pos,
                rec.samples.get(ci),
                gt_i,
                baf_i,
                lrr_i,
                af,
                cfg.baf_weight,
                cfg.lrr_weight,
                cfg.baf_dev_control,
                cfg.lrr_dev_control,
                cfg.optimize,
            ) {
                control_points.push(point);
            }
        }
    }

    // Real Viterbi HMM per chromosome (replaces the previous wrapper that
    // forced CN=2 everywhere). t_same=0.9999 → mean segment ~10^4 sites.
    run_viterbi_per_chrom(&mut query_points, 0.9999, cfg.baf_weight, cfg.lrr_weight, cfg.baf_dev_query, cfg.lrr_dev_query);
    if cfg.control_sample.is_some() {
        run_viterbi_per_chrom(&mut control_points, 0.9999, cfg.baf_weight, cfg.lrr_weight, cfg.baf_dev_control, cfg.lrr_dev_control);
    }

    write_sample_outputs(
        &cfg.output_dir,
        &cfg.query_sample,
        &query_points,
        cfg.baf_dev_query,
        &cfg,
    )?;
    if let Some(control) = &cfg.control_sample {
        write_sample_outputs(
            &cfg.output_dir,
            control,
            &control_points,
            cfg.baf_dev_control,
            &cfg,
        )?;
        write_pair_summary(
            &cfg.output_dir,
            &cfg.query_sample,
            control,
            &query_points,
            &control_points,
            cfg.baf_dev_query,
            cfg.baf_dev_control,
            &cfg,
        )?;
    }

    Ok(())
}

fn parse_args(args: &[String]) -> Result<CnvConfig> {
    let mut cfg = CnvConfig {
        output_dir: PathBuf::from("cnv"),
        regions_overlap: 1,
        targets_overlap: 0,
        baf_weight: 1.0,
        lrr_weight: 0.2,
        baf_dev_query: 0.04,
        baf_dev_control: 0.04,
        lrr_dev_query: 0.2,
        lrr_dev_control: 0.2,
        ..Default::default()
    };

    let mut i = 0usize;
    while i < args.len() {
        let arg = &args[i];
        let next = |i: &mut usize| -> Result<&String> {
            *i += 1;
            args.get(*i)
                .ok_or_else(|| anyhow!("missing value for {}", arg))
        };

        match arg.as_str() {
            "-s" | "--query-sample" => cfg.query_sample = next(&mut i)?.clone(),
            "-c" | "--control-sample" => cfg.control_sample = Some(next(&mut i)?.clone()),
            "-o" | "--output-dir" => cfg.output_dir = PathBuf::from(next(&mut i)?),
            "-f" | "--AF-file" => cfg.af_file = Some(PathBuf::from(next(&mut i)?)),
            "-r" | "--regions" => cfg.regions.extend(parse_region_list(next(&mut i)?)),
            "-R" | "--regions-file" => cfg.regions.extend(parse_region_file(next(&mut i)?)?),
            "-t" | "--targets" => cfg.targets.extend(parse_region_list(next(&mut i)?)),
            "-T" | "--targets-file" => cfg.targets.extend(parse_region_file(next(&mut i)?)?),
            "--regions-overlap" => cfg.regions_overlap = next(&mut i)?.parse().unwrap_or(1),
            "--targets-overlap" => cfg.targets_overlap = next(&mut i)?.parse().unwrap_or(0),
            "-b" | "--BAF-weight" => cfg.baf_weight = next(&mut i)?.parse().unwrap_or(1.0),
            "-l" | "--LRR-weight" => cfg.lrr_weight = next(&mut i)?.parse().unwrap_or(0.2),
            "-d" | "--BAF-dev" => {
                let (a, b) = parse_pair(next(&mut i)?, 0.04, 0.04);
                cfg.baf_dev_query = a;
                cfg.baf_dev_control = b;
            }
            "-k" | "--LRR-dev" => {
                let (a, b) = parse_pair(next(&mut i)?, 0.2, 0.2);
                cfg.lrr_dev_query = a;
                cfg.lrr_dev_control = b;
            }
            "-O" | "--optimize" => cfg.optimize = Some(next(&mut i)?.parse().unwrap_or(1.0)),
            "-a" | "--aberrant" | "-e" | "--err-probability" | "--err-prob" | "-L"
            | "--LRR-smooth-win" | "-P" | "--same-prob" | "-x" | "--xy-prob" | "-p"
            | "--plot-threshold" | "-v" | "--verbosity" => {
                let _ = next(&mut i)?;
            }
            _ if arg.starts_with('-') => {}
            _ => cfg.input = PathBuf::from(arg),
        }
        i += 1;
    }

    if cfg.input.as_os_str().is_empty() {
        return Err(anyhow!("missing input VCF/BCF file for cnv"));
    }
    if cfg.query_sample.is_empty() {
        return Err(anyhow!("missing query sample, use -s/--query-sample"));
    }
    Ok(cfg)
}

fn parse_pair(s: &str, d1: f64, d2: f64) -> (f64, f64) {
    let mut it = s.split(',');
    let a = it.next().and_then(|x| x.parse::<f64>().ok()).unwrap_or(d1);
    let b = it.next().and_then(|x| x.parse::<f64>().ok()).unwrap_or(d2);
    (a, b)
}

fn parse_region_list(s: &str) -> Vec<RegionSpec> {
    s.split(',').filter_map(parse_region).collect()
}

fn parse_region_file(path: &str) -> Result<Vec<RegionSpec>> {
    let file = File::open(path).with_context(|| format!("failed to open region file {}", path))?;
    let mut out = Vec::new();
    for line in BufReader::new(file).lines() {
        let line = line?;
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        if t.contains('\t') {
            let cols: Vec<&str> = t.split('\t').collect();
            if cols.len() >= 3 {
                let chr = cols[0].to_string();
                let start = cols[1].parse::<u32>().ok().map(|x| x.saturating_add(1));
                let end = cols[2].parse::<u32>().ok();
                out.push(RegionSpec { chr, start, end });
                continue;
            }
        }
        if let Some(r) = parse_region(t) {
            out.push(r);
        }
    }
    Ok(out)
}

fn parse_region(s: &str) -> Option<RegionSpec> {
    if s.contains(':') {
        let (chr, beg, end) = crate::regions::parse_region_spec(s).ok()?;
        return Some(RegionSpec { chr, start: Some(beg), end: (end != u32::MAX).then_some(end) });
    }
    Some(RegionSpec {
        chr: s.to_string(),
        start: None,
        end: None,
    })
}

fn record_matches_regions(chrom: &str, pos: u32, regions: &[RegionSpec]) -> bool {
    if regions.is_empty() {
        return true;
    }
    regions.iter().any(|r| {
        if r.chr != chrom {
            return false;
        }
        let s = r.start.unwrap_or(0);
        let e = r.end.unwrap_or(u32::MAX);
        pos >= s && pos <= e
    })
}

fn extract_sample_names(header: &[String]) -> Result<Vec<String>> {
    if !header.iter().any(|l| l.starts_with("#CHROM")) {
        return Err(anyhow!("VCF header does not contain #CHROM line"));
    }
    let samples = crate::vcf::header::extract_samples(header);
    if samples.is_empty() {
        return Err(anyhow!("VCF does not contain sample columns"));
    }
    Ok(samples)
}

fn build_format_index(format: &str) -> HashMap<String, usize> {
    format
        .split(':')
        .enumerate()
        .map(|(i, k)| (k.to_string(), i))
        .collect()
}

fn parse_sample_point(
    chrom: &str,
    pos: u32,
    sample: Option<&String>,
    gt_i: Option<usize>,
    baf_i: usize,
    lrr_i: usize,
    af: Option<f64>,
    baf_weight: f64,
    lrr_weight: f64,
    baf_dev: f64,
    lrr_dev: f64,
    optimize: Option<f64>,
) -> Option<SamplePoint> {
    let sample = sample?;
    let fields: Vec<&str> = sample.split(':').collect();
    let baf = fields.get(baf_i)?.parse::<f64>().ok()?;
    let lrr = fields.get(lrr_i)?.parse::<f64>().ok()?;
    let gt = gt_i.and_then(|i| fields.get(i)).copied().unwrap_or("./.");
    let is_het = is_heterozygous(gt);
    let probs = state_probabilities(
        baf, lrr, af, baf_weight, lrr_weight, baf_dev, lrr_dev, optimize,
    );
    let cn = probs
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i as u8)
        .unwrap_or(2);

    Some(SamplePoint {
        chr: chrom.to_string(),
        pos,
        baf,
        lrr,
        cn,
        probs,
        is_het,
    })
}

fn is_heterozygous(gt: &str) -> bool {
    let sep = if gt.contains('|') { '|' } else { '/' };
    let mut it = gt.split(sep);
    let a = it.next().unwrap_or(".");
    let b = it.next().unwrap_or(".");
    a != "." && b != "." && a != b
}

/// Log-emission of (BAF, LRR) for each of the 5 CN states.
fn log_emission(
    baf: f64,
    lrr: f64,
    baf_weight: f64,
    lrr_weight: f64,
    baf_dev: f64,
    lrr_dev: f64,
) -> [f64; N_STATES] {
    // BAF target pairs (two symmetric peaks per state).
    let baf_targets: [[f64; 2]; N_STATES] = [
        [0.0, 1.0],           // CN0 — no signal but extremes
        [0.0, 1.0],           // CN1 — hemizygous: only AA or BB possible
        [0.5, 0.5],           // CN2 — diploid AB
        [1.0 / 3.0, 2.0 / 3.0], // CN3 — duplication AAB/ABB
        [0.25, 0.75],         // CN4 — amp AAAB/ABBB
    ];
    let lrr_targets = [-2.0, -0.5, 0.0, 0.35, 0.65];

    let mut s = [0.0f64; N_STATES];
    let bd = baf_dev.max(1e-6);
    let ld = lrr_dev.max(1e-6);
    for state in 0..N_STATES {
        let db = (baf - baf_targets[state][0])
            .abs()
            .min((baf - baf_targets[state][1]).abs());
        let dl = (lrr - lrr_targets[state]).abs();
        s[state] = -0.5
            * (baf_weight * (db / bd).powi(2)
                + lrr_weight * (dl / ld).powi(2));
    }
    s
}

fn state_probabilities(
    baf: f64,
    lrr: f64,
    _af: Option<f64>,
    baf_weight: f64,
    lrr_weight: f64,
    baf_dev: f64,
    lrr_dev: f64,
    optimize: Option<f64>,
) -> [f64; N_STATES] {
    let log_e = log_emission(baf, lrr, baf_weight, lrr_weight, baf_dev, lrr_dev);
    let mut s = [0.0f64; N_STATES];
    for i in 0..N_STATES { s[i] = log_e[i].exp(); }

    if let Some(o) = optimize {
        let ab = (1.0 - o).clamp(0.0, 1.0);
        let priors = [0.05 * ab + 1e-6, 0.15 * ab + 1e-6, 1.0, 0.25 * ab + 1e-6, 0.05 * ab + 1e-6];
        for state in 0..N_STATES { s[state] *= priors[state]; }
    }
    let z: f64 = s.iter().sum();
    if z > 0.0 { for v in &mut s { *v /= z; } }
    else { s = [1.0 / N_STATES as f64; N_STATES]; }
    s
}

/// Viterbi HMM over a single chromosome's points.
/// `t_same` ∈ (0,1): probability of staying in same state at adjacent sites
/// (default ≈ 0.9999 → mean segment length ~10⁴ sites).
fn viterbi_chromosome(points: &mut [SamplePoint], t_same: f64, baf_weight: f64, lrr_weight: f64, baf_dev: f64, lrr_dev: f64) {
    let n = points.len();
    if n == 0 { return; }
    let t_self = t_same.clamp(1e-6, 1.0 - 1e-6).ln();
    let t_switch = ((1.0 - t_same).max(1e-9) / (N_STATES as f64 - 1.0)).ln();

    // log-prior toward diploid
    let log_prior: [f64; N_STATES] = [
        (0.02_f64).ln(),
        (0.08_f64).ln(),
        (0.80_f64).ln(),
        (0.08_f64).ln(),
        (0.02_f64).ln(),
    ];

    let mut dp = vec![[f64::NEG_INFINITY; N_STATES]; n];
    let mut bt = vec![[0u8; N_STATES]; n];

    let e0 = log_emission(points[0].baf, points[0].lrr, baf_weight, lrr_weight, baf_dev, lrr_dev);
    for s in 0..N_STATES { dp[0][s] = log_prior[s] + e0[s]; }

    for t in 1..n {
        let e = log_emission(points[t].baf, points[t].lrr, baf_weight, lrr_weight, baf_dev, lrr_dev);
        for j in 0..N_STATES {
            let mut best = f64::NEG_INFINITY;
            let mut arg = 0u8;
            for i in 0..N_STATES {
                let trans = if i == j { t_self } else { t_switch };
                let v = dp[t - 1][i] + trans + e[j];
                if v > best { best = v; arg = i as u8; }
            }
            dp[t][j] = best;
            bt[t][j] = arg;
        }
    }

    // Backtrack
    let mut path = vec![0u8; n];
    let (mut best_last, mut best_v) = (0u8, f64::NEG_INFINITY);
    for s in 0..N_STATES {
        if dp[n - 1][s] > best_v { best_v = dp[n - 1][s]; best_last = s as u8; }
    }
    path[n - 1] = best_last;
    for t in (1..n).rev() {
        path[t - 1] = bt[t][path[t] as usize];
    }

    // Write back to points: cn = path[t] (0..=4); probs = posterior of state via softmax of dp
    for t in 0..n {
        points[t].cn = path[t];
        let max = dp[t].iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let mut probs = [0.0f64; N_STATES];
        let mut z = 0.0;
        for s in 0..N_STATES { probs[s] = (dp[t][s] - max).exp(); z += probs[s]; }
        if z > 0.0 { for v in &mut probs { *v /= z; } }
        points[t].probs = probs;
    }
}

fn read_af_value(info: &str) -> Option<f64> {
    for item in info.split(';') {
        if let Some(rest) = item.strip_prefix("AF=") {
            let val = rest.split(',').next().unwrap_or(rest);
            if let Ok(v) = val.parse::<f64>() {
                return Some(v);
            }
        }
    }
    None
}

fn load_af_file(path: &Path) -> Result<HashMap<(String, u32, String, String), f64>> {
    let file = File::open(path).with_context(|| format!("failed to open AF file {:?}", path))?;
    let reader: Box<dyn Read> = if matches!(
        path.extension().and_then(|x| x.to_str()),
        Some("gz" | "bgz" | "bgzf")
    ) {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let mut out = HashMap::new();
    for line in BufReader::new(reader).lines() {
        let line = line?;
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = t.split('\t').collect();
        if cols.len() < 4 {
            continue;
        }
        let pos = match cols[1].parse::<u32>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let (r, a) = match cols[2].split_once(',') {
            Some(v) => v,
            None => continue,
        };
        let af = match cols[3].parse::<f64>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        out.insert((cols[0].to_string(), pos, r.to_string(), a.to_string()), af);
    }
    Ok(out)
}

fn write_sample_outputs(
    out_dir: &Path,
    sample: &str,
    points: &[SamplePoint],
    _baf_dev: f64,
    cfg: &CnvConfig,
) -> Result<()> {
    let mut dat = File::create(out_dir.join(format!("dat.{}.tab", sample)))?;
    writeln!(dat, "# [1]Chromosome\t[2]Position\t[3]BAF\t[4]LRR")?;
    for p in points {
        writeln!(dat, "{}\t{}\t{:.3}\t{:.3}", p.chr, p.pos, p.baf, p.lrr)?;
    }

    let mut cn = File::create(out_dir.join(format!("cn.{}.tab", sample)))?;
    writeln!(
        cn,
        "# [1]Chromosome\t[2]Position\t[3]CN\t[4]P(CN0)\t[5]P(CN1)\t[6]P(CN2)\t[7]P(CN3)\t[8]P(CN4)"
    )?;
    for p in points {
        writeln!(
            cn,
            "{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            p.chr, p.pos, p.cn, p.probs[0], p.probs[1], p.probs[2], p.probs[3], p.probs[4]
        )?;
    }

    let rows = build_summary_rows(points);
    let mut sum = File::create(out_dir.join(format!("summary.{}.tab", sample)))?;
    writeln!(
        sum,
        "# RG, Regions [2]Chromosome\t[3]Start\t[4]End\t[5]Copy Number state\t[6]Quality\t[7]nSites\t[8]nHETs"
    )?;
    writeln!(sum, "# This file was produced by: kira-bt cnv(native)")?;
    writeln!(
        sum,
        "# The command line was:\tkira-bt cnv -s {}{} -o {} {}",
        cfg.query_sample,
        cfg.control_sample
            .as_ref()
            .map(|c| format!(" -c {}", c))
            .unwrap_or_default(),
        cfg.output_dir.display(),
        cfg.input.display()
    )?;
    writeln!(sum, "#")?;
    writeln!(
        sum,
        "# RG, Regions\t[2]Chromosome\t[3]Start\t[4]End\t[5]Copy number:{}\t[6]Quality\t[7]nSites\t[8]nHETs",
        sample
    )?;
    for r in rows {
        writeln!(
            sum,
            "RG\t{}\t{}\t{}\t{}\t{:.1}\t{}\t{}",
            r.chr, r.start, r.end, r.cn, r.quality, r.n_sites, r.n_hets
        )?;
    }

    Ok(())
}

fn write_pair_summary(
    out_dir: &Path,
    query: &str,
    control: &str,
    q_points: &[SamplePoint],
    c_points: &[SamplePoint],
    q_baf_dev: f64,
    c_baf_dev: f64,
    _cfg: &CnvConfig,
) -> Result<()> {
    let q_rows = build_summary_rows(q_points);
    let mut f = File::create(out_dir.join("summary.tab"))?;
    writeln!(f, "# This file was produced by: kira-bt cnv(native)")?;
    writeln!(
        f,
        "# RG, Regions\t[2]Chromosome\t[3]Start\t[4]End\t[5]Copy number:{}\t[6]Copy number:{}\t[7]Quality\t[8]nSites in (5)\t[9]nHETs in (5)\t[10]nSites in (6)\t[11]nHETs in(6)",
        query, control
    )?;
    let _ = (q_baf_dev, c_baf_dev);
    for q in &q_rows {
        let c = region_summary_for(c_points, &q.chr, q.start, q.end).unwrap_or(SampleSummaryRow {
            chr: q.chr.clone(),
            start: q.start,
            end: q.end,
            cn: 2,
            quality: 0.0,
            n_sites: 0,
            n_hets: 0,
        });
        let qual = (q.quality + c.quality) * 0.5;
        writeln!(
            f,
            "RG\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{}\t{}\t{}\t{}",
            q.chr, q.start, q.end, q.cn, c.cn, qual, q.n_sites, q.n_hets, c.n_sites, c.n_hets
        )?;
    }
    Ok(())
}

fn build_summary_rows(points: &[SamplePoint]) -> Vec<SampleSummaryRow> {
    if points.is_empty() {
        return Vec::new();
    }
    let mut rows = Vec::new();
    let mut cur_chr = points[0].chr.as_str();
    let mut cur_cn = points[0].cn;
    let mut start = points[0].pos;
    let mut end = points[0].pos;
    let mut n_sites = 1usize;
    let mut n_hets = usize::from(points[0].is_het);
    let mut prob_sum = points[0].probs[cur_cn as usize];

    for p in points.iter().skip(1) {
        let split = p.chr != cur_chr || p.cn != cur_cn;
        if split {
            let q = score_quality(prob_sum / n_sites as f64);
            rows.push(SampleSummaryRow {
                chr: cur_chr.to_string(),
                start,
                end,
                cn: cur_cn,
                quality: q,
                n_sites,
                n_hets,
            });
            cur_chr = p.chr.as_str();
            cur_cn = p.cn;
            start = p.pos;
            end = p.pos;
            n_sites = 1;
            n_hets = usize::from(p.is_het);
            prob_sum = p.probs[p.cn as usize];
            continue;
        }
        end = p.pos;
        n_sites += 1;
        n_hets += usize::from(p.is_het);
        prob_sum += p.probs[p.cn as usize];
    }

    let q = score_quality(prob_sum / n_sites as f64);
    rows.push(SampleSummaryRow {
        chr: cur_chr.to_string(),
        start,
        end,
        cn: cur_cn,
        quality: q,
        n_sites,
        n_hets,
    });

    rows
}

fn score_quality(avg_max_prob: f64) -> f64 {
    let p = (1.0 - avg_max_prob).clamp(1e-6, 1.0);
    -10.0 * p.log10()
}

fn region_summary_for(
    points: &[SamplePoint],
    chr: &str,
    start: u32,
    end: u32,
) -> Option<SampleSummaryRow> {
    let mut selected = Vec::new();
    for p in points {
        if p.chr == chr && p.pos >= start && p.pos <= end {
            selected.push(p);
        }
    }
    if selected.is_empty() {
        return None;
    }
    let mut cn_counts = [0usize; N_STATES];
    let mut n_hets = 0usize;
    let mut prob_sum = 0.0f64;
    for p in &selected {
        cn_counts[(p.cn as usize).min(N_STATES - 1)] += 1;
        prob_sum += p.probs[(p.cn as usize).min(N_STATES - 1)];
        if p.is_het {
            n_hets += 1;
        }
    }
    let mut cn = 0usize;
    for s in 1..N_STATES {
        if cn_counts[s] > cn_counts[cn] {
            cn = s;
        }
    }
    Some(SampleSummaryRow {
        chr: chr.to_string(),
        start,
        end,
        cn: cn as u8,
        quality: score_quality(prob_sum / selected.len() as f64),
        n_sites: selected.len(),
        n_hets,
    })
}

fn run_viterbi_per_chrom(points: &mut [SamplePoint], t_same: f64, baf_weight: f64, lrr_weight: f64, baf_dev: f64, lrr_dev: f64) {
    if points.is_empty() { return; }
    let mut i = 0usize;
    while i < points.len() {
        let chr = points[i].chr.clone();
        let mut j = i + 1;
        while j < points.len() && points[j].chr == chr { j += 1; }
        viterbi_chromosome(&mut points[i..j], t_same, baf_weight, lrr_weight, baf_dev, lrr_dev);
        i = j;
    }
}

#[cfg(test)]
#[path = "../../tests/unit/cnv.rs"]
mod tests;
