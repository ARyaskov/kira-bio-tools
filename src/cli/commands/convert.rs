use anyhow::{Result, anyhow};
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

use crate::VcfReader;
use crate::cli::args::ConvertArgs;

pub fn cmd_convert(args: ConvertArgs) -> Result<()> {
    let mut combined: Vec<String> = Vec::new();
    if let Some(p) = &args.input { combined.push(p.to_string_lossy().into_owned()); }
    if let Some(p) = &args.output { combined.push("-o".into()); combined.push(p.to_string_lossy().into_owned()); }
    if let Some(p) = &args.output_type { combined.push("-O".into()); combined.push(p.clone()); }
    if args.gvcf2vcf { combined.push("--gvcf2vcf".into()); }
    if let Some(p) = &args.tsv2vcf { combined.push("--tsv2vcf".into()); combined.push(p.to_string_lossy().into_owned()); }
    if let Some(s) = &args.gensample { combined.push("-G".into()); combined.push(s.clone()); }
    if let Some(s) = &args.gen2vcf { combined.push("-g".into()); combined.push(s.clone()); }
    if let Some(s) = &args.haplegendsample { combined.push("--haplegendsample".into()); combined.push(s.clone()); }
    if let Some(s) = &args.hapsample { combined.push("-H".into()); combined.push(s.clone()); }
    if args.haploid { combined.push("--haploid".into()); }
    if let Some(s) = &args.columns { combined.push("-c".into()); combined.push(s.clone()); }
    if let Some(p) = &args.fasta_ref { combined.push("-f".into()); combined.push(p.to_string_lossy().into_owned()); }
    if args.chrom { combined.push("--chrom".into()); }
    if let Some(s) = &args.samples { combined.push("-s".into()); combined.push(s.clone()); }
    if let Some(p) = &args.samples_file { combined.push("-S".into()); combined.push(p.to_string_lossy().into_owned()); }
    if let Some(s) = &args.regions { combined.push("-r".into()); combined.push(s.clone()); }
    if let Some(p) = &args.regions_file { combined.push("-R".into()); combined.push(p.to_string_lossy().into_owned()); }
    if let Some(s) = &args.include { combined.push("-i".into()); combined.push(s.clone()); }
    if let Some(s) = &args.exclude { combined.push("-e".into()); combined.push(s.clone()); }
    if let Some(s) = &args.tag { combined.push("--tag".into()); combined.push(s.clone()); }
    combined.extend(args.passthrough.iter().cloned());
    let cfg = parse_args(&combined)?;
    match cfg.mode {
        ConvertMode::Vcf2Gen | ConvertMode::Vcf2Hap => emit_table_from_vcf(&cfg),
        ConvertMode::Gen2Vcf => emit_vcf_from_gen(&cfg),
        ConvertMode::HapSample2Vcf => emit_vcf_from_hap(&cfg),
        ConvertMode::HapLegend2Vcf => emit_vcf_from_haplegend(&cfg),
        ConvertMode::Tsv2Vcf => emit_vcf_from_tsv(&cfg),
        ConvertMode::Gvcf2Vcf => emit_vcf_from_gvcf(&cfg),
        ConvertMode::Generic => emit_vcf_passthrough(&cfg),
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum ConvertMode {
    Vcf2Gen,
    Vcf2Hap,
    Gen2Vcf,
    HapSample2Vcf,
    HapLegend2Vcf,
    Tsv2Vcf,
    Gvcf2Vcf,
    Generic,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum HapFlavor {
    Hls,
    Hs,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum VcfTextOutput {
    None,
    Gen,
    SampleFromGen,
    Hap,
    Legend,
    SampleFromHap,
}

struct ConvertCfg {
    mode: ConvertMode,
    hap_flavor: HapFlavor,
    inputs: Vec<PathBuf>,
    gen_pair: Option<(PathBuf, PathBuf)>,
    hap_pair: Option<(PathBuf, PathBuf)>,
    hap_triple: Option<(PathBuf, PathBuf, PathBuf)>,
    tsv: Option<PathBuf>,
    col_spec: Option<String>,
    sample_name: Option<String>,
    filter_expr: Option<String>,
    vcf_ids: bool,
    three_n6: bool,
    tag: String,
    gensample_spec: Option<String>,
    hapsample_spec: Option<String>,
    fasta_ref: Option<PathBuf>,
}

fn parse_args(args: &[String]) -> Result<ConvertCfg> {
    let mut cfg = ConvertCfg {
        mode: ConvertMode::Generic,
        hap_flavor: HapFlavor::Hls,
        inputs: Vec::new(),
        gen_pair: None,
        hap_pair: None,
        hap_triple: None,
        tsv: None,
        col_spec: None,
        sample_name: None,
        filter_expr: None,
        vcf_ids: false,
        three_n6: false,
        tag: "GT".to_string(),
        gensample_spec: None,
        hapsample_spec: None,
        fasta_ref: None,
    };

    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "-g" => {
                cfg.mode = ConvertMode::Vcf2Gen;
                i += 1;
                if let Some(v) = args.get(i) {
                    cfg.gensample_spec = Some(v.clone());
                }
            }
            "-h" => {
                cfg.mode = ConvertMode::Vcf2Hap;
                cfg.hap_flavor = HapFlavor::Hls;
                i += 1;
                if let Some(v) = args.get(i) {
                    cfg.hapsample_spec = Some(v.clone());
                }
            }
            "--hapsample" => {
                cfg.mode = ConvertMode::Vcf2Hap;
                cfg.hap_flavor = HapFlavor::Hs;
                i += 1;
                if let Some(v) = args.get(i) {
                    cfg.hapsample_spec = Some(v.clone());
                }
            }
            "-G" => {
                cfg.mode = ConvertMode::Gen2Vcf;
                i += 1;
                if let Some(v) = args.get(i) {
                    let (a, b) = parse_pair(v)?;
                    cfg.gen_pair = Some((a, b));
                }
            }
            "--hapsample2vcf" => {
                cfg.mode = ConvertMode::HapSample2Vcf;
                i += 1;
                if let Some(v) = args.get(i) {
                    let (a, b) = parse_pair(v)?;
                    cfg.hap_pair = Some((a, b));
                }
            }
            "--haplegendsample2vcf" => {
                cfg.mode = ConvertMode::HapLegend2Vcf;
                i += 1;
                if let Some(v) = args.get(i) {
                    let parts: Vec<&str> = v.split(',').map(str::trim).collect();
                    if parts.len() >= 3 {
                        cfg.hap_triple = Some((PathBuf::from(parts[0]), PathBuf::from(parts[1]), PathBuf::from(parts[2])));
                    }
                }
            }
            "--haplegendsample" => {
                cfg.mode = ConvertMode::Vcf2Hap;
                cfg.hap_flavor = HapFlavor::Hls;
                i += 1;
                if let Some(v) = args.get(i) {
                    cfg.hapsample_spec = Some(v.clone());
                }
            }
            "--tsv2vcf" => {
                cfg.mode = ConvertMode::Tsv2Vcf;
                i += 1;
                if let Some(v) = args.get(i) {
                    cfg.tsv = Some(PathBuf::from(v));
                }
            }
            "--gvcf2vcf" => cfg.mode = ConvertMode::Gvcf2Vcf,
            "--3N6" => cfg.three_n6 = true,
            "--vcf-ids" => cfg.vcf_ids = true,
            "--tag" => {
                i += 1;
                if let Some(v) = args.get(i) {
                    cfg.tag = v.to_ascii_uppercase();
                }
            }
            "-c" => {
                i += 1;
                if let Some(v) = args.get(i) {
                    cfg.col_spec = Some(v.clone());
                }
            }
            "-s" => {
                i += 1;
                if let Some(v) = args.get(i) {
                    cfg.sample_name = Some(v.clone());
                }
            }
            "-i" => {
                i += 1;
                if let Some(v) = args.get(i) {
                    cfg.filter_expr = Some(v.clone());
                }
            }
            "-f" => {
                i += 1;
                if let Some(v) = args.get(i) {
                    cfg.fasta_ref = Some(PathBuf::from(v));
                }
            }
            "--no-version" => {}
            a => {
                if !a.starts_with('-') || a == "-" {
                    cfg.inputs.push(PathBuf::from(a));
                }
            }
        }
        i += 1;
    }

    Ok(cfg)
}

fn parse_pair(v: &str) -> Result<(PathBuf, PathBuf)> {
    let Some((a, b)) = v.split_once(',') else {
        return Err(anyhow!("expected file pair 'a,b', got {v}"));
    };
    Ok((PathBuf::from(a), PathBuf::from(b)))
}

fn emit_table_from_vcf(cfg: &ConvertCfg) -> Result<()> {
    let Some(path) = find_first_vcf_like(&cfg.inputs) else {
        return Ok(());
    };
    let mut r = VcfReader::open(path)?;
    let headers = r.header()?;
    let samples = extract_samples(&headers);
    let out_mode = select_text_output(cfg);

    match out_mode {
        VcfTextOutput::SampleFromGen => {
            println!("ID_1 ID_2 missing");
            println!("0 0 0");
            for s in &samples {
                println!("{s} {s} 0");
            }
            return Ok(());
        }
        VcfTextOutput::SampleFromHap => {
            if cfg.hap_flavor == HapFlavor::Hs {
                println!("ID_1 ID_2 missing");
                println!("0 0 0");
                for s in &samples {
                    println!("{s} {s} 0");
                }
            } else {
                println!("sample population group sex");
                for s in &samples {
                    println!("{s} {s} {s} 2");
                }
            }
            return Ok(());
        }
        _ => {}
    }

    if out_mode == VcfTextOutput::Legend {
        println!("id position a0 a1");
    }

    while let Some(rec) = r.next_record()? {
        if rec.alt.contains(',') {
            continue;
        }
        let id1 = format!("{}:{}_{}_{}", rec.chrom, rec.pos, rec.ref_allele, rec.alt);
        let id2 = if cfg.vcf_ids && rec.id != "." {
            rec.id.clone()
        } else {
            id1.clone()
        };

        match out_mode {
            VcfTextOutput::Gen => {
                let mut out = Vec::<String>::new();
                if cfg.three_n6 {
                    out.push(rec.chrom.clone());
                }
                out.push(id1.clone());
                out.push(id2);
                out.push(rec.pos.to_string());
                out.push(rec.ref_allele.clone());
                out.push(rec.alt.clone());
                for s in &rec.samples {
                    let (a, b, c) = sample_to_gen_probs(s, rec.format.as_deref(), &cfg.tag);
                    out.push(a);
                    out.push(b);
                    out.push(c);
                }
                println!("{}", out.join(" "));
            }
            VcfTextOutput::Hap => {
                let mut out = Vec::<String>::new();
                if cfg.hap_flavor == HapFlavor::Hs {
                    if cfg.vcf_ids {
                        out.push(id1.clone());
                        out.push(if rec.id != "." {
                            rec.id.clone()
                        } else {
                            id1.clone()
                        });
                    } else {
                        out.push(rec.chrom.clone());
                        out.push(id1.clone());
                    }
                    out.push(rec.pos.to_string());
                    out.push(rec.ref_allele.clone());
                    out.push(rec.alt.clone());
                }
                let gt_idx = gt_index(rec.format.as_deref());
                for s in &rec.samples {
                    let gt = gt_from_sample(s, gt_idx);
                    let (h1, h2) = gt_to_haps(gt);
                    out.push(h1);
                    out.push(h2);
                }
                println!("{}", out.join(" "));
            }
            VcfTextOutput::Legend => {
                if !cfg.vcf_ids || rec.id == "." {
                    println!("{id1} {} {} {}", rec.pos, rec.ref_allele, rec.alt);
                } else {
                    println!("{} {} {} {}", rec.id, rec.pos, rec.ref_allele, rec.alt);
                }
            }
            VcfTextOutput::None => {}
            VcfTextOutput::SampleFromGen | VcfTextOutput::SampleFromHap => {}
        }
    }
    Ok(())
}

fn select_text_output(cfg: &ConvertCfg) -> VcfTextOutput {
    match cfg.mode {
        ConvertMode::Vcf2Gen => {
            let spec = cfg.gensample_spec.as_deref().unwrap_or("-,.");
            let cols = spec.split(',').map(str::trim).collect::<Vec<_>>();
            if cols.get(1).copied() == Some("-") {
                VcfTextOutput::SampleFromGen
            } else if cols.first().copied() == Some("-") {
                VcfTextOutput::Gen
            } else {
                VcfTextOutput::None
            }
        }
        ConvertMode::Vcf2Hap => {
            let spec = cfg.hapsample_spec.as_deref().unwrap_or("-,.,.");
            let cols = spec.split(',').map(str::trim).collect::<Vec<_>>();
            if cfg.hap_flavor == HapFlavor::Hs {
                if cols.get(1).copied() == Some("-") {
                    VcfTextOutput::SampleFromHap
                } else if cols.first().copied() == Some("-") {
                    VcfTextOutput::Hap
                } else {
                    VcfTextOutput::None
                }
            } else if cols.get(2).copied() == Some("-") {
                VcfTextOutput::SampleFromHap
            } else if cols.get(1).copied() == Some("-") {
                VcfTextOutput::Legend
            } else if cols.first().copied() == Some("-") {
                VcfTextOutput::Hap
            } else {
                VcfTextOutput::None
            }
        }
        _ => VcfTextOutput::None,
    }
}

fn extract_samples(headers: &[String]) -> Vec<String> {
    let Some(line) = headers.iter().find(|h| h.starts_with("#CHROM\t")) else {
        return Vec::new();
    };
    let cols = line.split('\t').collect::<Vec<_>>();
    if cols.len() <= 9 {
        return Vec::new();
    }
    cols[9..].iter().map(|s| (*s).to_string()).collect()
}

fn gt_index(fmt: Option<&str>) -> Option<usize> {
    fmt.and_then(|f| f.split(':').position(|k| k == "GT"))
}

fn tag_index(fmt: Option<&str>, tag: &str) -> Option<usize> {
    fmt.and_then(|f| f.split(':').position(|k| k.eq_ignore_ascii_case(tag)))
}

fn gt_from_sample(sample: &str, gt_idx: Option<usize>) -> &str {
    let Some(idx) = gt_idx else {
        return "./.";
    };
    sample.split(':').nth(idx).unwrap_or("./.")
}

fn sample_to_gen_probs(sample: &str, fmt: Option<&str>, tag: &str) -> (String, String, String) {
    let gt = gt_from_sample(sample, gt_index(fmt));
    let haploid = gt.split(['/', '|']).count() == 1;
    if matches!(tag, "PL" | "GL" | "GP")
        && let Some(idx) = tag_index(fmt, tag)
        && let Some(raw) = sample.split(':').nth(idx)
    {
        let vals = raw
            .split(',')
            .map(|v| v.parse::<f64>().ok())
            .collect::<Vec<_>>();
        if (vals.len() == 2 || vals.len() == 3) && vals.iter().all(Option::is_some) {
            let mut p = if vals.len() == 2 {
                if haploid {
                    vec![vals[0].unwrap_or(0.0), vals[1].unwrap_or(0.0)]
                } else {
                    vec![vals[0].unwrap_or(0.0), 0.0, vals[1].unwrap_or(0.0)]
                }
            } else {
                vec![
                    vals[0].unwrap_or(0.0),
                    vals[1].unwrap_or(0.0),
                    vals[2].unwrap_or(0.0),
                ]
            };
            if tag == "PL" {
                for v in &mut p {
                    *v = 10f64.powf(-(*v) / 10.0);
                }
            } else if tag == "GL" {
                for v in &mut p {
                    *v = 10f64.powf(*v);
                }
            }
            let sum = p.iter().sum::<f64>();
            if sum > 0.0 {
                if vals.len() == 2 && haploid {
                    return (
                        format!("{:.6}", p[0] / sum),
                        "0".to_string(),
                        format!("{:.6}", p[1] / sum),
                    );
                }
                return (
                    format!("{:.6}", p[0] / sum),
                    format!("{:.6}", p[1] / sum),
                    format!("{:.6}", p[2] / sum),
                );
            }
        }
    }
    gt_to_gen_probs(gt)
}

fn gt_to_gen_probs(gt: &str) -> (String, String, String) {
    let parts = gt.split(['/', '|']).collect::<Vec<_>>();
    if parts.is_empty() {
        return ("0.33".to_string(), "0.33".to_string(), "0.33".to_string());
    }
    if parts.len() == 1 {
        return match parts[0] {
            "." => ("0.5".to_string(), "0.0".to_string(), "0.5".to_string()),
            "0" => ("1".to_string(), "0".to_string(), "0".to_string()),
            "1" => ("0".to_string(), "0".to_string(), "1".to_string()),
            _ => ("0".to_string(), "0".to_string(), "0".to_string()),
        };
    }
    let (a, b) = (parts[0], parts[1]);
    if a == "." || b == "." {
        return ("0.33".to_string(), "0.33".to_string(), "0.33".to_string());
    }
    match (a, b) {
        ("0", "0") => ("1".to_string(), "0".to_string(), "0".to_string()),
        ("1", "1") => ("0".to_string(), "0".to_string(), "1".to_string()),
        ("0", "1") | ("1", "0") => ("0".to_string(), "1".to_string(), "0".to_string()),
        _ => ("0".to_string(), "0".to_string(), "0".to_string()),
    }
}

fn gt_to_haps(gt: &str) -> (String, String) {
    let phased = gt.contains('|');
    let parts = gt.split(['/', '|']).collect::<Vec<_>>();
    if parts.is_empty() {
        return ("?".to_string(), "-".to_string());
    }
    if parts.len() == 1 {
        let a = normalize_hap_allele(parts[0]).to_string();
        return (if a == "." { "?".to_string() } else { a }, "-".to_string());
    }
    if parts[0] == "." || parts[1] == "." {
        let a = if parts[0] == "." {
            "?"
        } else {
            normalize_hap_allele(parts[0])
        };
        let b = if parts[1] == "." {
            "?"
        } else {
            normalize_hap_allele(parts[1])
        };
        return (a.to_string(), b.to_string());
    }
    let mut a = normalize_hap_allele(parts[0]).to_string();
    let mut b = normalize_hap_allele(parts[1]).to_string();
    if !phased {
        if a == "0" || a == "1" {
            a.push('*');
        }
        if b == "0" || b == "1" {
            b.push('*');
        }
    }
    (a, b)
}

fn normalize_hap_allele(v: &str) -> &'static str {
    match v {
        "0" | "0*" => "0",
        "1" | "1*" => "1",
        "?" => ".",
        "." => ".",
        "-" => "-",
        _ => ".",
    }
}

fn emit_vcf_from_gen(cfg: &ConvertCfg) -> Result<()> {
    let (gen_path, sample_path) = cfg
        .gen_pair
        .clone()
        .ok_or_else(|| anyhow!("-G requires gen,sample"))?;
    let samples = read_sample_names(&sample_path)?;

    println!("##fileformat=VCFv4.2");
    println!("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    println!("##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Estimated Genotype Probability\">");
    println!(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
        samples.join("\t")
    );

    let text = fs::read_to_string(gen_path)?;
    for line in text.lines() {
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        let cols = t.split_whitespace().collect::<Vec<_>>();
        if cols.len() < 6 {
            continue;
        }
        let id1 = cols[0];
        let id2 = cols[1];
        let chrom = id1.split(':').next().unwrap_or("1");
        let pos = cols[2];
        let r = cols[3];
        let a = cols[4];
        let id_out = if cfg.vcf_ids { id2 } else { "." };

        let probs = &cols[5..];
        let mut out_samples = Vec::<String>::new();
        let mut i = 0usize;
        while i + 2 < probs.len() {
            let p0s = normalize_num_token(probs[i]);
            let p1s = normalize_num_token(probs[i + 1]);
            let p2s = normalize_num_token(probs[i + 2]);
            let p0 = p0s.parse::<f64>().unwrap_or(0.0);
            let p1 = p1s.parse::<f64>().unwrap_or(0.0);
            let p2 = p2s.parse::<f64>().unwrap_or(0.0);
            let gt = if p1 > p0 && p1 >= p2 {
                "0/1"
            } else if p2 > p0 && p2 > p1 {
                "1/1"
            } else {
                "0/0"
            };
            out_samples.push(format!("{gt}:{p0s},{p1s},{p2s}"));
            i += 3;
        }
        while out_samples.len() < samples.len() {
            out_samples.push("./.:0.33,0.33,0.33".to_string());
        }

        println!(
            "{chrom}\t{pos}\t{id_out}\t{r}\t{a}\t.\t.\t.\tGT:GP\t{}",
            out_samples.join("\t")
        );
    }
    Ok(())
}

fn emit_vcf_from_hap(cfg: &ConvertCfg) -> Result<()> {
    let (hap_path, sample_path) = cfg
        .hap_pair
        .clone()
        .ok_or_else(|| anyhow!("--hapsample2vcf requires hap,sample"))?;
    let samples = read_sample_names(&sample_path)?;

    println!("##fileformat=VCFv4.2");
    println!("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    println!(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
        samples.join("\t")
    );

    let text = fs::read_to_string(hap_path)?;
    for line in text.lines() {
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        let cols = t.split_whitespace().collect::<Vec<_>>();
        if cols.len() < 6 {
            continue;
        }

        let id1 = if cols[0].contains(':') && cols[0].contains('_') {
            cols[0]
        } else {
            cols[1]
        };
        let chrom = id1.split(':').next().unwrap_or(cols[0]);
        let pos = cols[2];
        let r = cols[3];
        let a = cols[4];
        let id = if cfg.vcf_ids && cols[1] != id1 {
            cols[1]
        } else {
            "."
        };

        let mut info = ".".to_string();
        if a.starts_with('<')
            && let Some(end) = parse_end_from_id(id1)
        {
            info = format!("END={end}");
        }

        let haps = &cols[5..];
        let mut out_samples = Vec::<String>::new();
        let mut i = 0usize;
        while i + 1 < haps.len() {
            let h1 = normalize_hap_allele(haps[i]);
            let h2 = normalize_hap_allele(haps[i + 1]);
            let unphased = haps[i].contains('*') || haps[i + 1].contains('*');
            let gt = if h2 == "-" {
                h1.to_string()
            } else if h1 == "." && h2 == "." {
                ".|.".to_string()
            } else if unphased {
                format!("{h1}/{h2}")
            } else {
                format!("{h1}|{h2}")
            };
            out_samples.push(gt);
            i += 2;
        }
        while out_samples.len() < samples.len() {
            out_samples.push(".|.".to_string());
        }

        println!(
            "{chrom}\t{pos}\t{id}\t{r}\t{a}\t.\t.\t{info}\tGT\t{}",
            out_samples.join("\t")
        );
    }
    Ok(())
}

fn emit_vcf_from_haplegend(cfg: &ConvertCfg) -> Result<()> {
    let (hap_path, legend_path, sample_path) = cfg
        .hap_triple
        .clone()
        .ok_or_else(|| anyhow!("--haplegendsample2vcf requires hap,legend,sample"))?;
    let samples = read_sample_names(&sample_path)?;

    println!("##fileformat=VCFv4.2");
    println!("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    println!(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
        samples.join("\t")
    );

    let legend_text = fs::read_to_string(&legend_path)?;
    let hap_text = fs::read_to_string(&hap_path)?;
    let chrom_default = legend_path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("1")
        .to_string();

    let mut leg_iter = legend_text.lines();
    let header_line = leg_iter.next().unwrap_or("");
    let header_cols: Vec<&str> = header_line.split_whitespace().collect();
    let pos_idx = header_cols.iter().position(|c| c.eq_ignore_ascii_case("position")).unwrap_or(1);
    let a0_idx = header_cols.iter().position(|c| c.eq_ignore_ascii_case("a0") || c.eq_ignore_ascii_case("allele0")).unwrap_or(2);
    let a1_idx = header_cols.iter().position(|c| c.eq_ignore_ascii_case("a1") || c.eq_ignore_ascii_case("allele1")).unwrap_or(3);

    for (leg_line, hap_line) in leg_iter.zip(hap_text.lines()) {
        let lt = leg_line.trim();
        if lt.is_empty() { continue; }
        let lc: Vec<&str> = lt.split_whitespace().collect();
        if lc.len() < 4 { continue; }
        let id = lc.first().copied().unwrap_or(".");
        let pos = lc.get(pos_idx).copied().unwrap_or("0");
        let r = lc.get(a0_idx).copied().unwrap_or(".");
        let a = lc.get(a1_idx).copied().unwrap_or(".");
        let chrom = if id.contains(':') {
            id.split(':').next().unwrap_or(chrom_default.as_str())
        } else { chrom_default.as_str() };

        let haps: Vec<&str> = hap_line.split_whitespace().collect();
        let mut gts: Vec<String> = Vec::with_capacity(samples.len());
        let mut i = 0usize;
        while i + 1 < haps.len() {
            let h1 = haps[i];
            let h2 = haps[i + 1];
            let gt = if (h1 == "0" || h1 == "1") && (h2 == "0" || h2 == "1") {
                format!("{}|{}", h1, h2)
            } else if h1 == "?" || h2 == "?" {
                ".|.".to_string()
            } else {
                format!("{}|{}", h1, h2)
            };
            gts.push(gt);
            i += 2;
        }
        while gts.len() < samples.len() { gts.push(".|.".to_string()); }
        println!("{chrom}\t{pos}\t{id}\t{r}\t{a}\t.\t.\t.\tGT\t{}", gts.join("\t"));
    }
    Ok(())
}

fn parse_end_from_id(id1: &str) -> Option<String> {
    let (_, right) = id1.split_once(':')?;
    let parts = right.split('_').collect::<Vec<_>>();
    if parts.len() >= 4 {
        return Some(parts[3].to_string());
    }
    None
}

fn normalize_num_token(v: &str) -> String {
    if let Some(x) = v.strip_suffix(".0") {
        return x.to_string();
    }
    v.to_string()
}

fn emit_vcf_from_tsv(cfg: &ConvertCfg) -> Result<()> {
    let path = cfg
        .tsv
        .clone()
        .ok_or_else(|| anyhow!("--tsv2vcf requires input file"))?;
    let spec = cfg
        .col_spec
        .clone()
        .unwrap_or_else(|| "ID,CHROM,POS,REF,ALT".to_string());
    let map = spec
        .split(',')
        .map(|x| x.trim().to_string())
        .collect::<Vec<_>>();

    let (contigs, fasta) = load_fasta(cfg.fasta_ref.as_ref())?;

    println!("##fileformat=VCFv4.2");
    println!("##FILTER=<ID=PASS,Description=\"All filters passed\">");
    for (ctg, len) in &contigs {
        println!("##contig=<ID={ctg},length={len}>");
    }
    println!("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    if cfg.sample_name.is_some() {
        println!(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
            cfg.sample_name.clone().unwrap_or_default()
        );
    } else {
        println!("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    }

    let text = fs::read_to_string(path)?;
    for line in text.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let cols = t.split_whitespace().collect::<Vec<_>>();
        if cols.is_empty() {
            continue;
        }

        let mut id = ".".to_string();
        let mut chrom = "1".to_string();
        let mut pos = "1".to_string();
        let mut ref_col: Option<String> = None;
        let mut alt_col: Option<String> = None;
        let mut aa_col: Option<String> = None;

        for (i, key) in map.iter().enumerate() {
            let Some(val) = cols.get(i).copied() else {
                continue;
            };
            match key.as_str() {
                "ID" => id = val.to_string(),
                "CHROM" => chrom = val.to_string(),
                "POS" => pos = val.to_string(),
                "REF" => ref_col = Some(val.to_string()),
                "ALT" => alt_col = Some(val.to_string()),
                "AA" => aa_col = Some(val.replace(' ', "")),
                "-" => {}
                _ => {}
            }
        }

        let pos_u32 = pos.parse::<u32>().ok().unwrap_or(1);
        let ref_base_from_fa = fasta
            .get(&chrom)
            .and_then(|s| s.chars().nth((pos_u32.saturating_sub(1)) as usize))
            .map(|c| c.to_ascii_uppercase().to_string());
        let ref_base = ref_base_from_fa
            .clone()
            .or(ref_col.clone())
            .unwrap_or_else(|| "N".to_string());

        if cfg.sample_name.is_some() {
            let aa = aa_col
                .clone()
                .or_else(|| {
                    if let (Some(r), Some(a)) = (ref_col.clone(), alt_col.clone()) {
                        Some(format!("{r}{a}"))
                    } else {
                        None
                    }
                })
                .unwrap_or_else(|| "--".to_string());
            let (alt, gt) = aa_to_alt_gt(&aa, &ref_base);
            println!("{chrom}\t{pos}\t{id}\t{ref_base}\t{alt}\t.\t.\t.\tGT\t{gt}");
        } else {
            let ref_no_sample = ref_col
                .clone()
                .or(ref_base_from_fa.clone())
                .unwrap_or_else(|| ref_base.clone());
            if let Some(aa) = aa_col.clone()
                && aa.contains('-')
            {
                continue;
            }
            let mut alt = alt_col.clone().unwrap_or_else(|| ".".to_string());
            if let Some(aa) = aa_col {
                let a = aa.chars().collect::<Vec<_>>();
                if a.len() >= 2 {
                    let r = a[0].to_ascii_uppercase().to_string();
                    let x = a[1].to_ascii_uppercase().to_string();
                    let ref_final = if let Some(rc) = ref_col.clone() {
                        rc
                    } else if let Some(rf) = ref_base_from_fa.clone() {
                        rf
                    } else {
                        r
                    };
                    alt = if x == ref_final { ".".to_string() } else { x };
                    println!("{chrom}\t{pos}\t{id}\t{ref_final}\t{alt}\t.\t.\t.");
                    continue;
                }
            }
            if alt == ref_no_sample {
                alt = ".".to_string();
            }
            println!("{chrom}\t{pos}\t{id}\t{ref_no_sample}\t{alt}\t.\t.\t.");
        }
    }
    Ok(())
}

fn aa_to_alt_gt(aa: &str, ref_base: &str) -> (String, String) {
    let a = aa
        .chars()
        .map(|c| c.to_ascii_uppercase())
        .filter(|c| !c.is_whitespace())
        .collect::<Vec<_>>();
    if a.is_empty() || a.iter().all(|c| *c == '-' || *c == '.') {
        return (".".to_string(), "./.".to_string());
    }
    if a.len() == 1 {
        let x = a[0].to_string();
        if x == ref_base {
            return (".".to_string(), "0".to_string());
        }
        return (x, "1".to_string());
    }
    let a1 = if a[0] == '-' || a[0] == '.' {
        None
    } else {
        Some(a[0].to_string())
    };
    let a2 = if a[1] == '-' || a[1] == '.' {
        None
    } else {
        Some(a[1].to_string())
    };
    if a1.is_none() || a2.is_none() {
        return (".".to_string(), "./.".to_string());
    }
    let a1 = a1.unwrap_or_default();
    let a2 = a2.unwrap_or_default();
    let mut alts = Vec::<String>::new();
    if a1 != ref_base {
        alts.push(a1.clone());
    }
    if a2 != ref_base && !alts.contains(&a2) {
        alts.push(a2.clone());
    }
    if alts.is_empty() {
        return (".".to_string(), "0/0".to_string());
    }
    let i1 = if a1 == ref_base {
        0
    } else {
        alts.iter().position(|x| x == &a1).unwrap_or(0) + 1
    };
    let i2 = if a2 == ref_base {
        0
    } else {
        alts.iter().position(|x| x == &a2).unwrap_or(0) + 1
    };
    (alts.join(","), format!("{i1}/{i2}"))
}

fn load_fasta(path: Option<&PathBuf>) -> Result<(Vec<(String, usize)>, HashMap<String, String>)> {
    let mut order = Vec::<(String, usize)>::new();
    let mut out = HashMap::<String, String>::new();
    let Some(p) = path else {
        return Ok((order, out));
    };
    let text = fs::read_to_string(p)?;
    let mut name = String::new();
    let mut seq = String::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if !name.is_empty() {
                order.push((name.clone(), seq.len()));
                out.insert(name.clone(), seq.clone());
                seq.clear();
            }
            name = rest.split_whitespace().next().unwrap_or(rest).to_string();
            continue;
        }
        seq.push_str(line.trim());
    }
    if !name.is_empty() {
        order.push((name.clone(), seq.len()));
        out.insert(name, seq);
    }
    Ok((order, out))
}

fn emit_vcf_from_gvcf(cfg: &ConvertCfg) -> Result<()> {
    let Some(path) = find_first_vcf_like(&cfg.inputs) else {
        return Err(anyhow!("--gvcf2vcf requires VCF input"));
    };
    let mut r = VcfReader::open(path)?;
    let headers = r.header()?;
    let (_ctg_order, fasta) = load_fasta(cfg.fasta_ref.as_ref())?;

    for h in &headers {
        if h.starts_with("##bcftools") {
            continue;
        }
        println!("{h}");
    }

    let filter_target = parse_filter_target(cfg.filter_expr.as_deref());

    while let Some(rec) = r.next_record()? {
        if let Some(target) = filter_target.as_deref()
            && rec.filter != target
        {
            if rec.alt != "." {
                continue;
            }
        }
        if rec.alt != "." {
            print_record(&rec);
            continue;
        }

        let end = info_int(&rec.info, "END")
            .and_then(|v| u32::try_from(v).ok())
            .unwrap_or(rec.pos);
        if end <= rec.pos || rec.filter != "PASS" {
            print_record(&rec);
            continue;
        }
        let info_no_end = strip_info_end(&rec.info);
        for pos in rec.pos..=end {
            let base = fasta
                .get(&rec.chrom)
                .and_then(|s| s.chars().nth((pos.saturating_sub(1)) as usize))
                .map(|c| c.to_ascii_uppercase().to_string())
                .unwrap_or_else(|| rec.ref_allele.clone());
            print!(
                "{}\t{}\t{}\t{}\t.\t{}\t{}\t{}",
                rec.chrom, pos, rec.id, base, rec.qual, rec.filter, info_no_end
            );
            if let Some(fmt) = &rec.format {
                print!("\t{fmt}");
                for s in &rec.samples {
                    print!("\t{s}");
                }
            }
            println!();
        }
    }
    Ok(())
}

fn parse_filter_target(expr: Option<&str>) -> Option<String> {
    let e = expr?.trim();
    if let Some(v) = e.strip_prefix("FILTER=") {
        return Some(v.trim().trim_matches('"').trim_matches('\'').to_string());
    }
    None
}

fn info_int(info: &str, key: &str) -> Option<i32> {
    for kv in info.split(';') {
        if let Some((k, v)) = kv.split_once('=')
            && k == key
        {
            return v.parse::<i32>().ok();
        }
    }
    None
}

fn strip_info_end(info: &str) -> String {
    let out = info
        .split(';')
        .filter(|kv| !kv.starts_with("END=") && !kv.is_empty() && *kv != ".")
        .collect::<Vec<_>>()
        .join(";");
    if out.is_empty() { ".".to_string() } else { out }
}

fn emit_vcf_passthrough(cfg: &ConvertCfg) -> Result<()> {
    if let Some(path) = find_first_vcf_like(&cfg.inputs) {
        let mut r = VcfReader::open(path)?;
        let headers = r.header()?;
        for h in &headers {
            println!("{h}");
        }
        while let Some(rec) = r.next_record()? {
            print_record(&rec);
        }
    }
    Ok(())
}

fn print_record(rec: &crate::vcf::structs::VcfRecord) {
    print!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        rec.chrom, rec.pos, rec.id, rec.ref_allele, rec.alt, rec.qual, rec.filter, rec.info
    );
    if let Some(fmt) = &rec.format {
        print!("\t{fmt}");
        for s in &rec.samples {
            print!("\t{s}");
        }
    }
    println!();
}

fn read_sample_names(path: &PathBuf) -> Result<Vec<String>> {
    let text = fs::read_to_string(path)?;
    let mut names = Vec::<String>::new();
    for line in text.lines() {
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        let cols = t.split_whitespace().collect::<Vec<_>>();
        if cols.len() >= 3 && cols[2] == "missing" {
            continue;
        }
        if cols.len() >= 3 && cols[0] == "0" && cols[1] == "0" {
            continue;
        }
        if cols.len() >= 2 {
            names.push(cols[1].to_string());
        }
    }
    Ok(names)
}

fn find_first_vcf_like(inputs: &[PathBuf]) -> Option<&PathBuf> {
    inputs.iter().find(|p| {
        let s = p.to_string_lossy();
        s.ends_with(".vcf")
            || s.ends_with(".vcf.gz")
            || s.ends_with(".gvcf")
            || s.ends_with(".gvcf.gz")
            || s.ends_with(".bcf")
            || s.ends_with(".bcf.gz")
    })
}
