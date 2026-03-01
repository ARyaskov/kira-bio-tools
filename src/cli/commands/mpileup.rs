use anyhow::{Result, anyhow};
use std::collections::{BTreeSet, HashSet};
use std::env;
use std::fs;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::process::Command;

use crate::cli::args::MpileupArgs;

pub fn cmd_mpileup(args: MpileupArgs) -> Result<()> {
    let out_path = args
        .output
        .clone()
        .unwrap_or_else(|| PathBuf::from("out.mpileup.vcf"));
    let cfg = parse_args(&args.bcftools_args)?;
    let sample_names = derive_input_sample_names(&args.inputs)?;
    let selected = select_samples_for_output(
        &sample_names,
        cfg.sample_arg.as_deref(),
        cfg.sample_file.as_ref(),
    )?;

    let (prepared_inputs, temp_paths) = prepare_inputs(&args.inputs)?;
    let pileup_lines = run_samtools_mpileup(&cfg, &prepared_inputs)?;

    let out = File::create(&out_path)?;
    let mut w = BufWriter::new(out);

    write_headers(&mut w, &cfg)?;

    let mut contigs = BTreeSet::<String>::new();
    for line in &pileup_lines {
        if let Some(chrom) = line.split('\t').next() {
            contigs.insert(chrom.to_string());
        }
    }
    for chrom in contigs {
        writeln!(w, "##contig=<ID={chrom}>")?;
    }

    if selected.is_empty() {
        writeln!(w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;
    } else {
        let names = selected
            .iter()
            .filter_map(|&i| sample_names.get(i).cloned())
            .collect::<Vec<_>>()
            .join("\t");
        writeln!(
            w,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{names}"
        )?;
    }

    for line in pileup_lines {
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 4 {
            continue;
        }

        let chrom = fields[0];
        let pos = fields[1];
        let ref_base = normalize_base_char(fields[2].chars().next().unwrap_or('N'));
        let mut per_sample = Vec::<SampleCounts>::new();

        let mut i = 3usize;
        while i + 2 < fields.len() {
            let bases = fields[i + 1];
            per_sample.push(parse_sample_bases(ref_base, bases));
            i += 3;
        }

        if per_sample.len() < args.inputs.len() {
            per_sample.resize(args.inputs.len(), SampleCounts::default());
        }

        let mut agg = [0u32; 5];
        for c in &per_sample {
            for j in 0..5 {
                agg[j] += c.total[j];
            }
        }
        let total_dp: u32 = agg.iter().sum();
        let alt = pick_alt_base(ref_base, &agg);
        let ref_idx = base_idx(ref_base);
        let alt_idx = if alt == '.' { ref_idx } else { base_idx(alt) };

        let mut ref_fwd = 0u32;
        let mut ref_rev = 0u32;
        let mut alt_fwd = 0u32;
        let mut alt_rev = 0u32;
        for &si in &selected {
            let c = per_sample.get(si).cloned().unwrap_or_default();
            ref_fwd += c.fwd[ref_idx];
            ref_rev += c.rev[ref_idx];
            alt_fwd += c.fwd[alt_idx];
            alt_rev += c.rev[alt_idx];
        }
        let site_ref = ref_fwd + ref_rev;
        let site_dv = alt_fwd + alt_rev;
        let site_scr = if total_dp > 0 {
            site_dv as f64 / total_dp as f64 * 100.0
        } else {
            0.0
        };

        let info = build_info(
            &cfg.anno, total_dp, site_ref, site_dv, ref_fwd, ref_rev, alt_fwd, alt_rev, site_scr,
        );
        write!(w, "{chrom}\t{pos}\t.\t{ref_base}\t{alt}\t.\tPASS\t{info}")?;

        if !selected.is_empty() {
            let fmt_keys = build_format_keys(&cfg.anno);
            write!(w, "\t{}", fmt_keys.join(":"))?;

            for &si in &selected {
                let c = per_sample.get(si).cloned().unwrap_or_default();
                let sdp: u32 = c.total.iter().sum();
                let ad_ref = c.total[ref_idx];
                let ad_alt = if alt == '.' { 0 } else { c.total[alt_idx] };
                let adf_ref = c.fwd[ref_idx];
                let adf_alt = if alt == '.' { 0 } else { c.fwd[alt_idx] };
                let adr_ref = c.rev[ref_idx];
                let adr_alt = if alt == '.' { 0 } else { c.rev[alt_idx] };
                let (gt, gq, pl) = infer_gt_gq_pl(ad_ref, ad_alt, alt == '.');

                let mut vals = vec![gt.to_string()];
                if cfg.anno.fmt_dp {
                    vals.push(sdp.to_string());
                }
                if cfg.anno.fmt_ad {
                    vals.push(format!("{ad_ref},{ad_alt}"));
                }
                if cfg.anno.fmt_adf {
                    vals.push(format!("{adf_ref},{adf_alt}"));
                }
                if cfg.anno.fmt_adr {
                    vals.push(format!("{adr_ref},{adr_alt}"));
                }
                if cfg.anno.fmt_dv {
                    vals.push(ad_alt.to_string());
                }
                if cfg.anno.fmt_dpr {
                    vals.push(format!("{ad_ref},{ad_alt}"));
                }
                if cfg.anno.fmt_scr {
                    let scr = if sdp > 0 {
                        ad_alt as f64 / sdp as f64 * 100.0
                    } else {
                        0.0
                    };
                    vals.push(format!("{scr:.3}"));
                }
                if cfg.anno.fmt_gq {
                    vals.push(gq.to_string());
                }
                if cfg.anno.fmt_pl {
                    vals.push(format!("{},{},{}", pl[0], pl[1], pl[2]));
                }
                write!(w, "\t{}", vals.join(":"))?;
            }
        }

        writeln!(w)?;
    }

    w.flush()?;
    cleanup_temp_paths(&temp_paths);
    Ok(())
}

fn write_headers<W: Write>(w: &mut W, cfg: &MpileupCfg) -> Result<()> {
    writeln!(w, "##fileformat=VCFv4.2")?;
    writeln!(w, "##source=kira-bt mpileup(native)")?;

    writeln!(
        w,
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=DV,Number=1,Type=Integer,Description=\"Total non-reference depth\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=DPR,Number=2,Type=Integer,Description=\"Ref/Alt depth\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Ref/Alt fwd/rev counts\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths on forward strand\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths on reverse strand\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=SP,Number=1,Type=Integer,Description=\"Strand bias metric\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=SCR,Number=1,Type=Float,Description=\"Soft clipping ratio proxy\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=NMBZ,Number=1,Type=Float,Description=\"Mismatch bias Z-score proxy\">"
    )?;

    writeln!(
        w,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Sample depth\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths on forward strand\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths on reverse strand\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Sample non-reference depth\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=DPR,Number=2,Type=Integer,Description=\"Sample ref/alt depth\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=SCR,Number=1,Type=Float,Description=\"Soft clipping ratio proxy\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">"
    )?;

    let _ = cfg;
    Ok(())
}

fn build_info(
    anno: &AnnotationCfg,
    total_dp: u32,
    site_ref: u32,
    site_dv: u32,
    ref_fwd: u32,
    ref_rev: u32,
    alt_fwd: u32,
    alt_rev: u32,
    site_scr: f64,
) -> String {
    let mut parts = vec![format!("DP={total_dp}")];
    if anno.info_dv {
        parts.push(format!("DV={site_dv}"));
    }
    if anno.info_dpr {
        parts.push(format!("DPR={site_ref},{site_dv}"));
    }
    if anno.info_dp4 {
        parts.push(format!("DP4={ref_fwd},{ref_rev},{alt_fwd},{alt_rev}"));
    }
    if anno.info_ad {
        parts.push(format!("AD={site_ref},{site_dv}"));
    }
    if anno.info_adf {
        parts.push(format!("ADF={ref_fwd},{alt_fwd}"));
    }
    if anno.info_adr {
        parts.push(format!("ADR={ref_rev},{alt_rev}"));
    }
    if anno.info_sp {
        parts.push("SP=0".to_string());
    }
    if anno.info_scr {
        parts.push(format!("SCR={site_scr:.3}"));
    }
    if anno.info_nmbz {
        parts.push("NMBZ=0.000".to_string());
    }
    parts.join(";")
}

fn build_format_keys(anno: &AnnotationCfg) -> Vec<String> {
    let mut keys = vec!["GT".to_string()];
    if anno.fmt_dp {
        keys.push("DP".to_string());
    }
    if anno.fmt_ad {
        keys.push("AD".to_string());
    }
    if anno.fmt_adf {
        keys.push("ADF".to_string());
    }
    if anno.fmt_adr {
        keys.push("ADR".to_string());
    }
    if anno.fmt_dv {
        keys.push("DV".to_string());
    }
    if anno.fmt_dpr {
        keys.push("DPR".to_string());
    }
    if anno.fmt_scr {
        keys.push("SCR".to_string());
    }
    if anno.fmt_gq {
        keys.push("GQ".to_string());
    }
    if anno.fmt_pl {
        keys.push("PL".to_string());
    }
    keys
}

struct MpileupCfg {
    ref_path: PathBuf,
    regions: Vec<String>,
    sample_arg: Option<String>,
    sample_file: Option<PathBuf>,
    anno: AnnotationCfg,
}

#[derive(Clone)]
struct AnnotationCfg {
    fmt_dp: bool,
    fmt_ad: bool,
    fmt_adf: bool,
    fmt_adr: bool,
    fmt_dv: bool,
    fmt_dpr: bool,
    fmt_scr: bool,
    fmt_gq: bool,
    fmt_pl: bool,
    info_dv: bool,
    info_dpr: bool,
    info_dp4: bool,
    info_ad: bool,
    info_adf: bool,
    info_adr: bool,
    info_sp: bool,
    info_scr: bool,
    info_nmbz: bool,
}

fn parse_args(args: &[String]) -> Result<MpileupCfg> {
    let mut ref_path = None::<PathBuf>;
    let mut regions = Vec::<String>::new();
    let mut sample_arg = None::<String>;
    let mut sample_file = None::<PathBuf>;
    let mut anno_tokens = Vec::<String>::new();

    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "-f" => {
                i += 1;
                if let Some(v) = args.get(i) {
                    ref_path = Some(PathBuf::from(v));
                }
            }
            "-r" | "-t" => {
                i += 1;
                if let Some(v) = args.get(i) {
                    regions.push(v.clone());
                }
            }
            "-s" => {
                i += 1;
                if let Some(v) = args.get(i) {
                    sample_arg = Some(v.clone());
                }
            }
            "-S" | "-G" => {
                i += 1;
                if let Some(v) = args.get(i) {
                    sample_file = Some(PathBuf::from(v));
                }
            }
            "-a" => {
                i += 1;
                if let Some(v) = args.get(i) {
                    anno_tokens.extend(v.split(',').map(|x| x.trim().to_string()));
                }
            }
            _ => {}
        }
        i += 1;
    }

    let ref_path = ref_path.ok_or_else(|| anyhow!("missing -f <ref.fa>"))?;
    Ok(MpileupCfg {
        ref_path,
        regions,
        sample_arg,
        sample_file,
        anno: parse_annotations(&anno_tokens),
    })
}

fn parse_annotations(tokens: &[String]) -> AnnotationCfg {
    let mut include = HashSet::<String>::new();
    let mut exclude = HashSet::<String>::new();
    for t in tokens {
        let x = t.trim();
        if x.is_empty() {
            continue;
        }
        if let Some(rest) = x.strip_prefix('-') {
            if !rest.is_empty() {
                exclude.insert(rest.to_ascii_uppercase());
            }
        } else {
            include.insert(x.to_ascii_uppercase());
        }
    }

    let req = |k: &str| include.contains(k);
    let block = |k: &str| exclude.contains(k);

    AnnotationCfg {
        fmt_dp: !block("DP"),
        fmt_ad: !block("AD"),
        fmt_adf: req("ADF") && !block("ADF"),
        fmt_adr: req("ADR") && !block("ADR"),
        fmt_dv: req("DV") && !block("DV"),
        fmt_dpr: req("DPR") && !block("DPR"),
        fmt_scr: req("FMT/SCR") && !block("FMT/SCR"),
        fmt_gq: req("GQ") && !block("GQ"),
        fmt_pl: req("PL") && !block("PL"),
        info_dv: req("DV") && !block("DV"),
        info_dpr: req("INFO/DPR") && !block("INFO/DPR"),
        info_dp4: req("DP4") && !block("DP4"),
        info_ad: req("INFO/AD") && !block("INFO/AD"),
        info_adf: req("INFO/ADF") && !block("INFO/ADF"),
        info_adr: req("INFO/ADR") && !block("INFO/ADR"),
        info_sp: req("SP") && !block("SP"),
        info_scr: req("INFO/SCR") && !block("INFO/SCR"),
        info_nmbz: req("INFO/NMBZ") && !block("INFO/NMBZ"),
    }
}

fn run_samtools_mpileup(cfg: &MpileupCfg, inputs: &[PathBuf]) -> Result<Vec<String>> {
    let mut cmd = Command::new("samtools");
    cmd.arg("mpileup");
    cmd.arg("-f").arg(&cfg.ref_path);

    for r in &cfg.regions {
        cmd.arg("-r").arg(r);
    }
    for input in inputs {
        cmd.arg(input);
    }

    let out = cmd.output()?;
    if !out.status.success() {
        return Err(anyhow!(
            "samtools mpileup failed: {}",
            String::from_utf8_lossy(&out.stderr)
        ));
    }

    let text = String::from_utf8(out.stdout)?;
    Ok(text.lines().map(|s| s.to_string()).collect())
}

fn prepare_inputs(inputs: &[PathBuf]) -> Result<(Vec<PathBuf>, Vec<PathBuf>)> {
    let mut prepared = Vec::with_capacity(inputs.len());
    let mut temp_paths = Vec::new();
    let pid = std::process::id();

    for (i, input) in inputs.iter().enumerate() {
        let is_sam = input
            .extension()
            .and_then(|x| x.to_str())
            .map(|x| x.eq_ignore_ascii_case("sam"))
            .unwrap_or(false);
        if !is_sam {
            prepared.push(input.clone());
            continue;
        }

        let mut tmp = env::temp_dir();
        tmp.push(format!("kira-mpileup-{pid}-{i}.bam"));

        let mut view = Command::new("samtools");
        view.arg("view").arg("-b").arg("-o").arg(&tmp).arg(input);
        let view_out = view.output()?;
        if !view_out.status.success() {
            return Err(anyhow!(
                "samtools view failed for {}: {}",
                input.display(),
                String::from_utf8_lossy(&view_out.stderr)
            ));
        }

        let mut index = Command::new("samtools");
        index.arg("index").arg(&tmp);
        let idx_out = index.output()?;
        if !idx_out.status.success() {
            return Err(anyhow!(
                "samtools index failed for {}: {}",
                tmp.display(),
                String::from_utf8_lossy(&idx_out.stderr)
            ));
        }

        temp_paths.push(tmp.clone());
        temp_paths.push(PathBuf::from(format!("{}.bai", tmp.display())));
        prepared.push(tmp);
    }

    Ok((prepared, temp_paths))
}

fn cleanup_temp_paths(paths: &[PathBuf]) {
    for p in paths {
        let _ = fs::remove_file(p);
    }
}

fn derive_input_sample_names(inputs: &[PathBuf]) -> Result<Vec<String>> {
    let mut out = Vec::new();
    for (i, p) in inputs.iter().enumerate() {
        let mut cmd = Command::new("samtools");
        cmd.arg("view").arg("-H").arg(p);
        let hdr = cmd.output()?;
        let mut sm = None::<String>;
        if hdr.status.success() {
            let txt = String::from_utf8_lossy(&hdr.stdout);
            for line in txt.lines() {
                if !line.starts_with("@RG\t") {
                    continue;
                }
                for f in line.split('\t') {
                    if let Some(v) = f.strip_prefix("SM:") {
                        sm = Some(v.to_string());
                        break;
                    }
                }
                if sm.is_some() {
                    break;
                }
            }
        }
        out.push(sm.unwrap_or_else(|| {
            format!(
                "{}_{}",
                p.file_stem().and_then(|x| x.to_str()).unwrap_or("sample"),
                i + 1
            )
        }));
    }
    Ok(out)
}

fn select_samples_for_output(
    names: &[String],
    sample_arg: Option<&str>,
    sample_file: Option<&PathBuf>,
) -> Result<Vec<usize>> {
    let mut selected = (0..names.len()).collect::<Vec<_>>();

    if let Some(s) = sample_arg {
        let invert = s.starts_with('^');
        let set = s
            .trim_start_matches('^')
            .split(',')
            .map(str::trim)
            .filter(|x| !x.is_empty())
            .map(|x| x.to_string())
            .collect::<BTreeSet<_>>();
        if !set.is_empty() {
            if invert {
                selected.retain(|i| !set.contains(&names[*i]));
            } else {
                selected.retain(|i| set.contains(&names[*i]));
            }
        }
    }

    if let Some(path) = sample_file
        && path.exists()
    {
        let txt = fs::read_to_string(path)?;
        let set = txt
            .lines()
            .map(str::trim)
            .filter(|x| !x.is_empty() && !x.starts_with('#'))
            .map(|x| x.split_whitespace().next().unwrap_or("").to_string())
            .filter(|x| !x.is_empty())
            .collect::<BTreeSet<_>>();
        if !set.is_empty() {
            selected.retain(|i| set.contains(&names[*i]));
        }
    }

    Ok(selected)
}

fn normalize_base_char(c: char) -> char {
    match c.to_ascii_uppercase() {
        'A' | 'C' | 'G' | 'T' => c.to_ascii_uppercase(),
        _ => 'N',
    }
}

#[derive(Clone, Default)]
struct SampleCounts {
    total: [u32; 5],
    fwd: [u32; 5],
    rev: [u32; 5],
}

fn parse_sample_bases(ref_base: char, bases: &str) -> SampleCounts {
    let mut c = SampleCounts::default();
    let ref_idx = base_idx(ref_base);
    let b = bases.as_bytes();
    let mut i = 0usize;

    while i < b.len() {
        match b[i] as char {
            '^' => i = (i + 2).min(b.len()),
            '$' => i += 1,
            '.' => {
                c.total[ref_idx] += 1;
                c.fwd[ref_idx] += 1;
                i += 1;
            }
            ',' => {
                c.total[ref_idx] += 1;
                c.rev[ref_idx] += 1;
                i += 1;
            }
            'A' | 'a' => {
                add_base(&mut c, 0, b[i] as char);
                i += 1;
            }
            'C' | 'c' => {
                add_base(&mut c, 1, b[i] as char);
                i += 1;
            }
            'G' | 'g' => {
                add_base(&mut c, 2, b[i] as char);
                i += 1;
            }
            'T' | 't' => {
                add_base(&mut c, 3, b[i] as char);
                i += 1;
            }
            'N' | 'n' => {
                add_base(&mut c, 4, b[i] as char);
                i += 1;
            }
            '+' | '-' => {
                i += 1;
                let mut n = 0usize;
                while i < b.len() && (b[i] as char).is_ascii_digit() {
                    n = n * 10 + (b[i] - b'0') as usize;
                    i += 1;
                }
                i = (i + n).min(b.len());
            }
            '*' | '#' | '<' | '>' => i += 1,
            _ => i += 1,
        }
    }

    c
}

fn add_base(c: &mut SampleCounts, idx: usize, ch: char) {
    c.total[idx] += 1;
    if ch.is_ascii_uppercase() {
        c.fwd[idx] += 1;
    } else {
        c.rev[idx] += 1;
    }
}

fn base_idx(b: char) -> usize {
    match b.to_ascii_uppercase() {
        'A' => 0,
        'C' => 1,
        'G' => 2,
        'T' => 3,
        _ => 4,
    }
}

fn idx_base(i: usize) -> char {
    match i {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        _ => 'N',
    }
}

fn pick_alt_base(ref_base: char, agg: &[u32; 5]) -> char {
    let r = base_idx(ref_base);
    let mut best_i = r;
    let mut best = 0u32;
    for (i, v) in agg.iter().enumerate().take(4) {
        if i == r {
            continue;
        }
        if *v > best {
            best = *v;
            best_i = i;
        }
    }
    if best == 0 { '.' } else { idx_base(best_i) }
}

fn infer_gt_gq_pl(ad_ref: u32, ad_alt: u32, no_alt: bool) -> (&'static str, u32, [u32; 3]) {
    if no_alt {
        return ("0/0", 99, [0, 255, 255]);
    }
    let sum = ad_ref + ad_alt;
    if sum == 0 {
        return ("./.", 0, [0, 0, 0]);
    }
    let pl = compute_pl_triplet(ad_ref, ad_alt);
    let (best_i, best_v) = min_pl_idx(&pl);
    let second_v = pl
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != best_i)
        .map(|(_, v)| *v)
        .min()
        .unwrap_or(best_v);
    let gq = second_v.saturating_sub(best_v).min(99);
    let gt = match best_i {
        0 => "0/0",
        1 => "0/1",
        2 => "1/1",
        _ => "./.",
    };
    (gt, gq, pl)
}

fn compute_pl_triplet(ad_ref: u32, ad_alt: u32) -> [u32; 3] {
    let eps = 0.01f64;
    let p00 = log10_like(ad_ref, ad_alt, 1.0 - eps, eps);
    let p01 = log10_like(ad_ref, ad_alt, 0.5, 0.5);
    let p11 = log10_like(ad_ref, ad_alt, eps, 1.0 - eps);
    let max_like = p00.max(p01).max(p11);
    let to_pl = |x: f64| -> u32 { ((-10.0 * (x - max_like)).round() as i64).max(0) as u32 };
    [to_pl(p00), to_pl(p01), to_pl(p11)]
}

fn log10_like(ad_ref: u32, ad_alt: u32, p_ref: f64, p_alt: f64) -> f64 {
    (ad_ref as f64) * p_ref.log10() + (ad_alt as f64) * p_alt.log10()
}

fn min_pl_idx(pl: &[u32; 3]) -> (usize, u32) {
    let mut i = 0usize;
    let mut v = pl[0];
    for (k, x) in pl.iter().enumerate().skip(1) {
        if *x < v {
            i = k;
            v = *x;
        }
    }
    (i, v)
}
