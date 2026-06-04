use anyhow::Result;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

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
        if cfg.collect_bafs {
            // Only retain het-like BAFs (drop hom calls) for GMM analysis.
            if (0.05..=0.95).contains(&af) {
                st.bafs.push(af);
            }
        }

        match gt_class(&rec, sample_idx) {
            GtClass::HomRef => st.hom_ref += 1,
            GtClass::Het => st.het += 1,
            GtClass::HomAlt => st.hom_alt += 1,
            GtClass::Missing => st.missing += 1,
        }
    }

    println!("#CHROM\tN_SITES\tHET\tHOM_REF\tHOM_ALT\tMEAN_AF\tABN_FRAC\tHET_RATE\tSCORE\tSTATUS\tCN_PRED\tPEAKS");
    if let Some(dir) = &cfg.plot_dir {
        let _ = std::fs::create_dir_all(dir);
    }
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

        let (cn_pred, peaks) = if cfg.collect_bafs && st.bafs.len() >= 30 {
            let fit = fit_gmm_baf(&st.bafs);
            let cn = infer_cn_from_peaks(&fit.peaks);
            let peak_str = fit
                .peaks
                .iter()
                .map(|(mu, w)| format!("{:.3}:{:.2}", mu, w))
                .collect::<Vec<_>>()
                .join(",");
            if let Some(dir) = &cfg.plot_dir {
                let svg_path = dir.join(format!("baf.{}.svg", chrom));
                let _ = write_baf_svg(&svg_path, &chrom, &st.bafs, &fit);
            }
            (cn, peak_str)
        } else {
            (0u8, "-".to_string())
        };

        println!(
            "{chrom}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.3}\t{}\t{}\t{}",
            st.n_sites, st.het, st.hom_ref, st.hom_alt, mean, abn_frac, het_rate, score, status,
            cn_pred, peaks,
        );
    }

    Ok(())
}

#[derive(Clone, Debug)]
struct GmmFit {
    peaks: Vec<(f64, f64)>, // (mean, weight)
    sigma: f64,
    log_lik: f64,
}

/// Fit a 1-vs-2 component Gaussian Mixture to the BAF distribution
/// (restricted to het sites, BAF ∈ (0,1)). Returns whichever model
/// has the higher BIC score.
fn fit_gmm_baf(bafs: &[f64]) -> GmmFit {
    let one = em_gmm(bafs, 1, 50);
    let two = em_gmm(bafs, 2, 80);
    let n = bafs.len() as f64;
    let bic1 = -2.0 * one.log_lik + 2.0 * n.ln();
    let bic2 = -2.0 * two.log_lik + 5.0 * n.ln();
    if bic2 + 6.0 < bic1 { two } else { one }
}

fn em_gmm(data: &[f64], k: usize, iters: usize) -> GmmFit {
    let mut means: Vec<f64> = if k == 1 {
        vec![0.5]
    } else {
        vec![0.4, 0.6]
    };
    let mut weights: Vec<f64> = vec![1.0 / k as f64; k];
    let mut sigma = 0.05f64;
    let n = data.len();
    if n == 0 {
        return GmmFit { peaks: means.iter().zip(weights.iter()).map(|(m, w)| (*m, *w)).collect(), sigma, log_lik: 0.0 };
    }

    let mut last_ll = f64::NEG_INFINITY;
    for _ in 0..iters {
        // E-step
        let mut resp = vec![vec![0.0f64; k]; n];
        for (i, x) in data.iter().enumerate() {
            let mut row = vec![0.0; k];
            let mut z = 0.0;
            for j in 0..k {
                let p = weights[j] * gaussian_pdf(*x, means[j], sigma);
                row[j] = p;
                z += p;
            }
            if z > 0.0 {
                for j in 0..k { row[j] /= z; }
            } else {
                for j in 0..k { row[j] = 1.0 / k as f64; }
            }
            resp[i] = row;
        }
        // M-step
        let mut new_means = vec![0.0f64; k];
        let mut new_weights = vec![0.0f64; k];
        let mut new_var = 0.0f64;
        let mut total = 0.0f64;
        for j in 0..k {
            let mut nk = 0.0;
            let mut s = 0.0;
            for (i, x) in data.iter().enumerate() {
                nk += resp[i][j];
                s += resp[i][j] * x;
            }
            if nk > 0.0 {
                new_means[j] = s / nk;
            } else {
                new_means[j] = means[j];
            }
            new_weights[j] = nk / n as f64;
            for (i, x) in data.iter().enumerate() {
                let d = x - new_means[j];
                new_var += resp[i][j] * d * d;
            }
            total += nk;
        }
        let new_sigma = ((new_var / total.max(1.0)).sqrt()).max(1e-4);
        // log-likelihood
        let mut ll = 0.0;
        for x in data {
            let mut p = 0.0;
            for j in 0..k {
                p += new_weights[j] * gaussian_pdf(*x, new_means[j], new_sigma);
            }
            ll += p.max(1e-300).ln();
        }
        means = new_means;
        weights = new_weights;
        sigma = new_sigma;
        if (ll - last_ll).abs() < 1e-5 { break; }
        last_ll = ll;
    }

    let mut peaks: Vec<(f64, f64)> = means.iter().zip(weights.iter()).map(|(m, w)| (*m, *w)).collect();
    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    GmmFit { peaks, sigma, log_lik: last_ll }
}

fn gaussian_pdf(x: f64, mu: f64, sigma: f64) -> f64 {
    let z = (x - mu) / sigma.max(1e-9);
    let inv = 1.0 / (sigma.max(1e-9) * (2.0 * std::f64::consts::PI).sqrt());
    inv * (-0.5 * z * z).exp()
}

/// Infer copy-number from fitted BAF peaks.
/// 1 peak near 0.5 → CN2; two peaks at ~0.33/0.67 → CN3; ~0.25/0.75 → CN4.
fn infer_cn_from_peaks(peaks: &[(f64, f64)]) -> u8 {
    match peaks.len() {
        0 => 0,
        1 => {
            let p = peaks[0].0;
            if (p - 0.5).abs() < 0.05 { 2 } else { 0 }
        }
        _ => {
            let (lo, hi) = (peaks[0].0, peaks[peaks.len() - 1].0);
            let split = (hi - lo).abs();
            if split < 0.1 { 2 }
            else if (lo - 1.0 / 3.0).abs() < 0.06 && (hi - 2.0 / 3.0).abs() < 0.06 { 3 }
            else if (lo - 0.25).abs() < 0.06 && (hi - 0.75).abs() < 0.06 { 4 }
            else if split > 0.55 { 4 }
            else { 3 }
        }
    }
}

/// Minimal SVG histogram + fitted Gaussian curves for the BAF distribution.
fn write_baf_svg(path: &std::path::Path, chrom: &str, bafs: &[f64], fit: &GmmFit) -> std::io::Result<()> {
    let w = 640usize;
    let h = 360usize;
    let pad = 30usize;
    let nbins = 50usize;
    let mut counts = vec![0u32; nbins];
    for x in bafs {
        let i = ((x * nbins as f64).floor() as usize).min(nbins - 1);
        counts[i] += 1;
    }
    let max_c = counts.iter().copied().max().unwrap_or(1) as f64;

    let mut svg = String::new();
    svg.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}\" height=\"{h}\" font-family=\"monospace\" font-size=\"11\">"
    ));
    svg.push_str(&format!(
        "<text x=\"{}\" y=\"15\">polysomy BAF — {} (n={})</text>",
        pad, chrom, bafs.len()
    ));
    let plot_w = w - 2 * pad;
    let plot_h = h - 2 * pad;
    let bin_w = plot_w as f64 / nbins as f64;
    for (i, c) in counts.iter().enumerate() {
        let bh = (*c as f64 / max_c) * plot_h as f64;
        let x = pad as f64 + i as f64 * bin_w;
        let y = (h - pad) as f64 - bh;
        svg.push_str(&format!(
            "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"#4a8\"/>",
            x, y, bin_w - 0.5, bh
        ));
    }
    // overlay fitted curves
    let n = bafs.len() as f64;
    let total_area = max_c * (bin_w / plot_w as f64) * n;
    let _ = total_area;
    let mut last_x = pad as f64;
    let mut last_y = (h - pad) as f64;
    for s in 0..=200 {
        let x_data = s as f64 / 200.0;
        let mut p = 0.0;
        for (mu, wgt) in &fit.peaks {
            p += wgt * gaussian_pdf(x_data, *mu, fit.sigma);
        }
        // Scale the pdf to histogram counts: pdf * bin_width * n_total → expected count
        let expected = p * (1.0 / nbins as f64) * n;
        let bh = (expected / max_c) * plot_h as f64;
        let xpos = pad as f64 + x_data * plot_w as f64;
        let ypos = (h - pad) as f64 - bh.max(0.0).min(plot_h as f64);
        if s > 0 {
            svg.push_str(&format!(
                "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#c33\" stroke-width=\"1.5\"/>",
                last_x, last_y, xpos, ypos
            ));
        }
        last_x = xpos;
        last_y = ypos;
    }
    // x-axis ticks
    for t in 0..=5 {
        let x = pad as f64 + (t as f64 / 5.0) * plot_w as f64;
        let lab = format!("{:.1}", t as f64 / 5.0);
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{}\">{}</text>",
            x - 8.0, h - 8, lab
        ));
    }
    // peaks annotation
    for (mu, wgt) in &fit.peaks {
        let xpos = pad as f64 + mu * plot_w as f64;
        svg.push_str(&format!(
            "<line x1=\"{:.1}\" y1=\"{}\" x2=\"{:.1}\" y2=\"{}\" stroke=\"#33c\" stroke-dasharray=\"3,2\"/>",
            xpos, pad, xpos, h - pad
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{}\" fill=\"#33c\">μ={:.2} w={:.2}</text>",
            xpos + 2.0, pad + 12, mu, wgt
        ));
    }
    svg.push_str("</svg>");
    let mut f = File::create(path)?;
    f.write_all(svg.as_bytes())?;
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
    bafs: Vec<f64>,
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
    collect_bafs: bool,
    plot_dir: Option<PathBuf>,
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
        collect_bafs: true,
        plot_dir: None,
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
            "-o" | "--output-dir" | "--plot-dir" => {
                i += 1;
                if let Some(v) = args.get(i) {
                    cfg.plot_dir = Some(PathBuf::from(v));
                }
            }
            "--no-gmm" => {
                cfg.collect_bafs = false;
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

#[cfg(test)]
#[path = "../../../tests/unit/cli_commands_polysomy.rs"]
mod tests;

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
