use crate::cli::args::GtcheckArgs;
use crate::vcf::UnifiedVcfReader;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

pub fn cmd_gtcheck(args: GtcheckArgs) -> Result<()> {
    let input = args.input.as_ref().ok_or_else(|| anyhow::anyhow!("missing input VCF/BCF"))?;
    let (samples_a, gts_a) = load_gts(input, args.homs_only)?;
    let out_path = args.output.clone();
    let mut out: Box<dyn Write> = match &out_path {
        Some(p) => Box::new(BufWriter::with_capacity(64 * 1024, File::create(p)?)),
        None => Box::new(BufWriter::with_capacity(64 * 1024, std::io::stdout())),
    };

    writeln!(out, "# This file was produced by kira_bt gtcheck (compatible with bcftools gtcheck output)")?;
    writeln!(out, "# DC, discordance: number of discordant sites, count, error rate, sample-a, sample-b")?;
    writeln!(out, "# DC\t[2]Query Sample\t[3]Genotyped Sample\t[4]Discordance\t[5]-log P(HWE)\t[6]Number of sites compared\t[7]Number of matching sites")?;

    if let Some(g) = &args.genotypes {
        let (samples_b, gts_b) = load_gts(g, args.homs_only)?;
        let query_samples = filter_samples(&samples_a, args.samples.as_deref(), args.samples_file.as_deref())?;
        for qi in &query_samples {
            let qname = &samples_a[*qi];
            let mut rows: Vec<(String, u64, f64, u64, u64)> = Vec::new();
            for (j, name) in samples_b.iter().enumerate() {
                let (ncmp, nsame) = cmp_vectors(&gts_a[*qi], &gts_b[j]);
                let mismatch = ncmp.saturating_sub(nsame);
                let score = if ncmp == 0 { 0.0 } else { mismatch as f64 / ncmp as f64 };
                rows.push((name.clone(), mismatch, score, ncmp, nsame));
            }
            rows.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));
            if args.n_matches > 0 { rows.truncate(args.n_matches); }
            for (name, mm, sc, nc, ns) in rows {
                writeln!(out, "DC\t{}\t{}\t{:.6e}\t{:.6e}\t{}\t{}", qname, name, mm, sc, nc, ns)?;
            }
        }
        return Ok(());
    }

    let pairs = pick_pairs(&samples_a, args.pairs.as_deref(), args.pairs_file.as_deref())?;
    for (a, b) in &pairs {
        let ia = samples_a.iter().position(|s| s == a).ok_or_else(|| anyhow::anyhow!("unknown sample {a}"))?;
        let ib = samples_a.iter().position(|s| s == b).ok_or_else(|| anyhow::anyhow!("unknown sample {b}"))?;
        let (ncmp, nsame) = cmp_vectors(&gts_a[ia], &gts_a[ib]);
        let mismatch = ncmp.saturating_sub(nsame);
        let score = if ncmp == 0 { 0.0 } else { mismatch as f64 / ncmp as f64 };
        writeln!(out, "DC\t{}\t{}\t{}\t{:.6e}\t{}\t{}", a, b, mismatch, score, ncmp, nsame)?;
    }
    out.flush()?;
    Ok(())
}

fn load_gts(path: &Path, homs_only: bool) -> Result<(Vec<String>, Vec<Vec<u8>>)> {
    let mut reader = UnifiedVcfReader::open(path).with_context(|| format!("open {:?}", path))?;
    let headers = reader.header()?;
    let samples = extract_samples(&headers);
    let mut gts: Vec<Vec<u8>> = vec![Vec::new(); samples.len()];

    while let Some(line) = reader.read_line()? {
        if line.is_empty() || line.as_bytes()[0] == b'#' { continue; }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 10 { continue; }
        let fmt = cols[8];
        let Some(gt_idx) = fmt.split(':').position(|k| k == "GT") else { continue; };
        for (i, raw) in cols[9..].iter().enumerate() {
            if i >= gts.len() { break; }
            let gt = raw.split(':').nth(gt_idx).unwrap_or(".");
            let cls = classify_gt(gt);
            if homs_only && cls == 1 { gts[i].push(255); continue; }
            gts[i].push(cls);
        }
    }
    Ok((samples, gts))
}

fn classify_gt(gt: &str) -> u8 {
    let alleles: Vec<&str> = gt.split(|c| c == '/' || c == '|').collect();
    if alleles.iter().any(|a| *a == "." || a.is_empty()) { return 255; }
    let first = alleles[0];
    if alleles.iter().all(|a| a == &first) {
        if first == "0" { 0 } else { 2 }
    } else { 1 }
}

fn cmp_vectors(a: &[u8], b: &[u8]) -> (u64, u64) {
    let len = a.len().min(b.len());
    let mut ncmp = 0u64; let mut nsame = 0u64;
    for i in 0..len {
        if a[i] == 255 || b[i] == 255 { continue; }
        ncmp += 1;
        if a[i] == b[i] { nsame += 1; }
    }
    (ncmp, nsame)
}

fn extract_samples(h: &[String]) -> Vec<String> {
    for line in h {
        if line.starts_with("#CHROM") {
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() > 9 { return cols[9..].iter().map(|s| s.to_string()).collect(); }
        }
    }
    Vec::new()
}

fn filter_samples(all: &[String], cli: Option<&str>, file: Option<&Path>) -> Result<Vec<usize>> {
    let mut names: Vec<String> = Vec::new();
    if let Some(s) = cli {
        names.extend(s.split(',').map(|t| t.trim().to_string()));
    }
    if let Some(p) = file {
        for line in BufReader::new(File::open(p)?).lines() {
            let l = line?; let t = l.trim();
            if !t.is_empty() && !t.starts_with('#') { names.push(t.to_string()); }
        }
    }
    if names.is_empty() { return Ok((0..all.len()).collect()); }
    let mut out = Vec::new();
    for n in &names {
        if let Some(i) = all.iter().position(|s| s == n) { out.push(i); }
    }
    Ok(out)
}

fn pick_pairs(samples: &[String], cli: Option<&str>, file: Option<&Path>) -> Result<Vec<(String, String)>> {
    let mut pairs: Vec<(String, String)> = Vec::new();
    if let Some(s) = cli {
        for tok in s.split(',') {
            let (a, b) = tok.split_once(':').ok_or_else(|| anyhow::anyhow!("--pairs: expected A:B"))?;
            pairs.push((a.trim().to_string(), b.trim().to_string()));
        }
    }
    if let Some(p) = file {
        for line in BufReader::new(File::open(p)?).lines() {
            let l = line?; let t = l.trim();
            if t.is_empty() || t.starts_with('#') { continue; }
            let parts: Vec<&str> = t.split(|c: char| c == '\t' || c == ',').collect();
            if parts.len() >= 2 { pairs.push((parts[0].trim().to_string(), parts[1].trim().to_string())); }
        }
    }
    if pairs.is_empty() {
        for i in 0..samples.len() {
            for j in i+1..samples.len() { pairs.push((samples[i].clone(), samples[j].clone())); }
        }
    }
    let _ = HashMap::<String, ()>::new();
    let _ = PathBuf::new();
    Ok(pairs)
}
