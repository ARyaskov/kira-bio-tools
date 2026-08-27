use crate::cli::args::GtcheckArgs;
use crate::vcf::UnifiedVcfReader;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

pub fn cmd_gtcheck(args: GtcheckArgs) -> Result<()> {
    let input = args.input.as_ref().ok_or_else(|| anyhow::anyhow!("missing input VCF/BCF"))?;
    let (samples_a, gts_a, af_a) = load_gts(input, args.homs_only)?;
    let out_path = args.output.clone();
    let mut out: Box<dyn Write> = match &out_path {
        Some(p) => Box::new(BufWriter::with_capacity(64 * 1024, File::create(p)?)),
        None => Box::new(BufWriter::with_capacity(64 * 1024, std::io::stdout())),
    };

    writeln!(out, "# This file was produced by kira_bt gtcheck (compatible with bcftools gtcheck output)")?;
    writeln!(out, "# DCv2\t[2]Query Sample\t[3]Genotyped Sample\t[4]Discordance\t[5]-log P(HWE)\t[6]Number of sites compared\t[7]Number of matching sites")?;
    writeln!(out, "# CN, contamination estimate (excess heterozygosity at the best match's homozygous sites)")?;
    writeln!(out, "# CN\t[2]Query Sample\t[3]Best Match\t[4]Contamination\t[5]Number of homozygous sites")?;

    if let Some(g) = &args.genotypes {
        let (samples_b, gts_b, af_b) = load_gts(g, args.homs_only)?;
        let query_samples = filter_samples(&samples_a, args.samples.as_deref(), args.samples_file.as_deref())?;
        for qi in &query_samples {
            let qname = &samples_a[*qi];
            // (genotyped name, j, discordance, -logP(HWE), ncmp, nsame)
            let mut rows: Vec<(String, usize, u64, f64, u64, u64)> = Vec::new();
            for (j, name) in samples_b.iter().enumerate() {
                let (ncmp, nsame, hwe) = cmp_full(&gts_a[*qi], &gts_b[j], &af_b);
                let mismatch = ncmp.saturating_sub(nsame);
                rows.push((name.clone(), j, mismatch, hwe, ncmp, nsame));
            }
            rows.sort_by(|a, b| a.2.cmp(&b.2));
            // Contamination is estimated against the best (lowest-discordance) match.
            if let Some((best_name, best_j, _, _, _, _)) = rows.first() {
                let (contam, n_hom) = estimate_contamination(&gts_a[*qi], &gts_b[*best_j]);
                writeln!(out, "# CN\t{}\t{}\t{:.6}\t{}", qname, best_name, contam, n_hom)?;
            }
            if args.n_matches > 0 { rows.truncate(args.n_matches); }
            for (name, _, mm, hwe, nc, ns) in rows {
                writeln!(out, "DCv2\t{}\t{}\t{}\t{}\t{}\t{}", qname, name, sci6(mm as f64), sci6(hwe), nc, ns)?;
            }
        }
        out.flush()?;
        return Ok(());
    }

    let pairs = pick_pairs(&samples_a, args.pairs.as_deref(), args.pairs_file.as_deref())?;
    for (a, b) in &pairs {
        let ia = samples_a.iter().position(|s| s == a).ok_or_else(|| anyhow::anyhow!("unknown sample {a}"))?;
        let ib = samples_a.iter().position(|s| s == b).ok_or_else(|| anyhow::anyhow!("unknown sample {b}"))?;
        let (ncmp, nsame, hwe) = cmp_full(&gts_a[ia], &gts_a[ib], &af_a);
        let mismatch = ncmp.saturating_sub(nsame);
        writeln!(out, "DCv2\t{}\t{}\t{}\t{}\t{}\t{}", a, b, sci6(mismatch as f64), sci6(hwe), ncmp, nsame)?;
    }
    out.flush()?;
    Ok(())
}

/// C `printf("%.6e")`-style: signed, zero-padded two-digit exponent
/// (e.g. `0.000000e+00`, `6.931472e-01`) to match bcftools output.
fn sci6(x: f64) -> String {
    let s = format!("{:.6e}", x);
    match s.find('e') {
        Some(epos) => {
            let (mant, exp) = s.split_at(epos);
            let exp = &exp[1..];
            let (sign, digits) = match exp.strip_prefix('-') {
                Some(d) => ("-", d),
                None => ("+", exp.strip_prefix('+').unwrap_or(exp)),
            };
            format!("{}e{}{:0>2}", mant, sign, digits)
        }
        None => s,
    }
}

/// `-ln P(genotype class | Hardy-Weinberg at alt-allele frequency `af`)`.
fn hwe_neglog(class: u8, af: f64) -> f64 {
    let af = af.clamp(1e-6, 1.0 - 1e-6);
    let p = match class {
        0 => (1.0 - af) * (1.0 - af),
        1 => 2.0 * af * (1.0 - af),
        2 => af * af,
        _ => return 0.0,
    };
    -p.max(1e-30).ln()
}

/// Compare two per-site genotype-class vectors. Returns `(sites compared,
/// matching sites, summed -logP(HWE) of the query genotypes)`.
fn cmp_full(a: &[u8], b: &[u8], af: &[f64]) -> (u64, u64, f64) {
    let len = a.len().min(b.len());
    let mut ncmp = 0u64;
    let mut nsame = 0u64;
    let mut hwe = 0.0f64;
    for i in 0..len {
        if a[i] == 255 || b[i] == 255 {
            continue;
        }
        ncmp += 1;
        if a[i] == b[i] {
            nsame += 1;
        }
        hwe += hwe_neglog(a[i], af.get(i).copied().unwrap_or(0.0));
    }
    (ncmp, nsame, hwe)
}

/// Rough contamination estimate: the fraction of the best match's homozygous
/// sites at which the query is heterozygous. Cross-sample contamination injects
/// foreign alleles, turning true homozygous sites into spurious hets, so this
/// excess-heterozygosity rate is a first-order contamination proxy.
fn estimate_contamination(query: &[u8], reference: &[u8]) -> (f64, u64) {
    let len = query.len().min(reference.len());
    let mut n_hom = 0u64;
    let mut het_at_hom = 0u64;
    for i in 0..len {
        let (q, r) = (query[i], reference[i]);
        if q == 255 || r == 255 {
            continue;
        }
        if r == 0 || r == 2 {
            n_hom += 1;
            if q == 1 {
                het_at_hom += 1;
            }
        }
    }
    let contam = if n_hom == 0 {
        0.0
    } else {
        het_at_hom as f64 / n_hom as f64
    };
    (contam, n_hom)
}

fn load_gts(path: &Path, homs_only: bool) -> Result<(Vec<String>, Vec<Vec<u8>>, Vec<f64>)> {
    let mut reader = UnifiedVcfReader::open(path).with_context(|| format!("open {:?}", path))?;
    let headers = reader.header()?;
    let samples = extract_samples(&headers);
    let mut gts: Vec<Vec<u8>> = vec![Vec::new(); samples.len()];
    let mut site_af: Vec<f64> = Vec::new();

    while let Some(line) = reader.read_line()? {
        if line.is_empty() || line.as_bytes()[0] == b'#' { continue; }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 10 { continue; }
        let fmt = cols[8];
        let Some(gt_idx) = fmt.split(':').position(|k| k == "GT") else { continue; };
        // Per-site alt-allele frequency, estimated from this file's genotypes
        // (drives the HWE term). Computed from the true class before homs_only
        // masking so the frequency is not biased by dropped hets.
        let mut alt = 0u32;
        let mut tot = 0u32;
        for (i, raw) in cols[9..].iter().enumerate() {
            if i >= gts.len() { break; }
            let gt = raw.split(':').nth(gt_idx).unwrap_or(".");
            let cls = classify_gt(gt);
            if cls <= 2 {
                alt += cls as u32;
                tot += 2;
            }
            if homs_only && cls == 1 { gts[i].push(255); continue; }
            gts[i].push(cls);
        }
        site_af.push(if tot > 0 { alt as f64 / tot as f64 } else { 0.0 });
    }
    Ok((samples, gts, site_af))
}

fn classify_gt(gt: &str) -> u8 {
    let alleles: Vec<&str> = gt.split(|c| c == '/' || c == '|').collect();
    if alleles.iter().any(|a| *a == "." || a.is_empty()) { return 255; }
    let first = alleles[0];
    if alleles.iter().all(|a| a == &first) {
        if first == "0" { 0 } else { 2 }
    } else { 1 }
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
