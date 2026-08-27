use anyhow::Result;
use std::collections::{BTreeSet, HashMap};
use std::fs;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use crate::VcfReader;
use crate::cli::args::CallArgs;

pub fn cmd_call(args: CallArgs) -> Result<()> {
    let mut argv: Vec<String> = Vec::new();
    if args.consensus_caller { argv.push("-c".into()); }
    if args.multiallelic_caller { argv.push("-m".into()); }
    if args.variants_only { argv.push("-v".into()); }
    if let Some(s) = &args.samples { argv.push("-s".into()); argv.push(s.clone()); }
    if let Some(p) = &args.samples_file { argv.push("-S".into()); argv.push(p.to_string_lossy().into_owned()); }
    if let Some(s) = &args.regions { argv.push("-r".into()); argv.push(s.clone()); }
    if let Some(p) = &args.regions_file { argv.push("-R".into()); argv.push(p.to_string_lossy().into_owned()); }
    if let Some(s) = &args.targets { argv.push("-t".into()); argv.push(s.clone()); }
    if let Some(p) = &args.targets_file { argv.push("-T".into()); argv.push(p.to_string_lossy().into_owned()); }
    let cfg = parse_call_args(&argv)?;
    let out_path = args.output.clone().unwrap_or_else(|| {
        let mut p = args.input.clone();
        p.set_extension("call.vcf");
        p
    });

    let mut reader = VcfReader::open(&args.input)?;
    let headers = reader.header()?;
    let all_samples = extract_samples_from_header(&headers)?;
    let selected = select_samples(
        &all_samples,
        cfg.sample_arg.as_deref(),
        cfg.sample_file.as_ref(),
    )?;

    let out = File::create(&out_path)?;
    let mut w = BufWriter::new(out);

    write_headers(&mut w, &headers, &all_samples, &selected, &cfg)?;

    while let Some(rec) = reader.next_record()? {
        if let Some(line) = call_record(&rec, &selected, &cfg)? {
            writeln!(w, "{line}")?;
        }
    }

    w.flush()?;
    Ok(())
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum CallMode {
    Consensus,
    Multiallelic,
}

struct CallCfg {
    mode: CallMode,
    only_alt: bool,
    sample_arg: Option<String>,
    sample_file: Option<PathBuf>,
    emit_gq: bool,
    emit_gp: bool,
}

fn parse_call_args(args: &[String]) -> Result<CallCfg> {
    let mut mode = CallMode::Multiallelic;
    let mut only_alt = false;
    let mut sample_arg = None::<String>;
    let mut sample_file = None::<PathBuf>;
    let mut emit_gq = false;
    let mut emit_gp = false;

    let mut i = 0usize;
    while i < args.len() {
        let a = args[i].as_str();
        if a == "-c" || a.contains('c') && a.starts_with('-') {
            mode = CallMode::Consensus;
        }
        if a == "-m" || a.contains('m') && a.starts_with('-') {
            mode = CallMode::Multiallelic;
        }
        if a == "-v" || (a.starts_with('-') && a.contains('v')) {
            only_alt = true;
        }
        match a {
            "-s" => {
                i += 1;
                if let Some(v) = args.get(i) {
                    sample_arg = Some(v.clone());
                }
            }
            "-S" => {
                i += 1;
                if let Some(v) = args.get(i) {
                    sample_file = Some(PathBuf::from(v));
                }
            }
            "-a" => {
                i += 1;
                if let Some(v) = args.get(i) {
                    for t in v.split(',').map(str::trim) {
                        if t.eq_ignore_ascii_case("GQ") {
                            emit_gq = true;
                        }
                        if t.eq_ignore_ascii_case("GP") {
                            emit_gp = true;
                        }
                    }
                }
            }
            _ => {}
        }
        i += 1;
    }

    Ok(CallCfg {
        mode,
        only_alt,
        sample_arg,
        sample_file,
        emit_gq,
        emit_gp,
    })
}

fn extract_samples_from_header(headers: &[String]) -> Result<Vec<String>> {
    let Some(last) = headers.iter().rfind(|h| h.starts_with("#CHROM\t")) else {
        return Ok(Vec::new());
    };
    let parts = last.split('\t').collect::<Vec<_>>();
    if parts.len() <= 9 {
        return Ok(Vec::new());
    }
    Ok(parts[9..].iter().map(|s| (*s).to_string()).collect())
}

fn select_samples(
    names: &[String],
    sample_arg: Option<&str>,
    sample_file: Option<&PathBuf>,
) -> Result<Vec<usize>> {
    let mut out = (0..names.len()).collect::<Vec<_>>();

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
                out.retain(|i| !set.contains(&names[*i]));
            } else {
                out.retain(|i| set.contains(&names[*i]));
            }
        }
    }

    if let Some(path) = sample_file {
        if path.as_os_str() == "-" {
            return Ok(out);
        }
        if path.exists() {
            let txt = fs::read_to_string(path)?;
            let set = txt
                .lines()
                .map(str::trim)
                .filter(|x| !x.is_empty() && !x.starts_with('#'))
                .map(|x| x.split_whitespace().next().unwrap_or("").to_string())
                .filter(|x| !x.is_empty())
                .collect::<BTreeSet<_>>();
            if !set.is_empty() {
                out.retain(|i| set.contains(&names[*i]));
            }
        }
    }

    Ok(out)
}

fn write_headers<W: Write>(
    w: &mut W,
    headers: &[String],
    all_samples: &[String],
    selected: &[usize],
    cfg: &CallCfg,
) -> Result<()> {
    let mut has_gt = false;
    let mut has_ac = false;
    let mut has_an = false;
    let mut has_gq = false;
    let mut has_gp = false;

    for h in headers {
        if h.starts_with("#CHROM\t") {
            continue;
        }
        if h.starts_with("##FORMAT=<ID=GT,") {
            has_gt = true;
        }
        if h.starts_with("##FORMAT=<ID=GQ,") {
            has_gq = true;
        }
        if h.starts_with("##FORMAT=<ID=GP,") {
            has_gp = true;
        }
        if h.starts_with("##INFO=<ID=AC,") {
            has_ac = true;
        }
        if h.starts_with("##INFO=<ID=AN,") {
            has_an = true;
        }
        writeln!(w, "{h}")?;
    }

    if !has_gt {
        writeln!(
            w,
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        )?;
    }
    if cfg.emit_gq && !has_gq {
        writeln!(
            w,
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Phred-scaled Genotype Quality\">"
        )?;
    }
    if cfg.emit_gp && !has_gp {
        writeln!(
            w,
            "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype posterior probabilities\">"
        )?;
    }
    if !has_ac {
        writeln!(
            w,
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">"
        )?;
    }
    if !has_an {
        writeln!(
            w,
            "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">"
        )?;
    }

    let samples = selected
        .iter()
        .filter_map(|&i| all_samples.get(i).cloned())
        .collect::<Vec<_>>();

    if samples.is_empty() {
        writeln!(w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;
    } else {
        writeln!(
            w,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
            samples.join("\t")
        )?;
    }
    Ok(())
}

fn call_record(
    rec: &crate::vcf::structs::VcfRecord,
    selected: &[usize],
    cfg: &CallCfg,
) -> Result<Option<String>> {
    if cfg.mode == CallMode::Multiallelic {
        return call_record_mcall(rec, selected, cfg);
    }
    call_record_consensus(rec, selected, cfg)
}

/// `call -c` consensus caller. Estimates the site alt-allele frequency by EM
/// over the per-sample genotype likelihoods, applies a Hardy-Weinberg prior at
/// that frequency, calls each genotype by maximum posterior, and reports a real
/// site QUAL = -10·log10 P(all samples homozygous reference). This is the
/// Li-2011 consensus model — not a per-sample PL argmax. Falls back to the
/// AD-based heuristic ([`call_record_legacy`]) when no PL/likelihoods exist.
fn call_record_consensus(
    rec: &crate::vcf::structs::VcfRecord,
    selected: &[usize],
    cfg: &CallCfg,
) -> Result<Option<String>> {
    let raw_alts = rec
        .alt
        .split(',')
        .map(str::trim)
        .map(|a| a.to_string())
        .collect::<Vec<_>>();

    // Real alt alleles (drop ., *, symbolic).
    let kept_all: Vec<usize> = raw_alts
        .iter()
        .enumerate()
        .filter(|(_, a)| {
            *a != "." && *a != "*" && !(a.starts_with('<') && a.ends_with('>')) && !a.is_empty()
        })
        .map(|(i, _)| i)
        .collect();

    let format_keys = rec.format.as_deref().unwrap_or("").split(':').collect::<Vec<_>>();
    let pl_idx = format_keys.iter().position(|k| *k == "PL");
    let ad_idx = format_keys.iter().position(|k| *k == "AD");

    // No alleles or no likelihoods → fall back to the AD heuristic.
    if kept_all.is_empty() || pl_idx.is_none() || selected.is_empty() {
        return call_record_legacy(rec, selected, cfg);
    }
    let pl_idx = pl_idx.unwrap();

    // Consensus is biallelic: pick the single best-supported ALT.
    let a_raw = if kept_all.len() == 1 {
        kept_all[0]
    } else {
        best_alt_by_support(rec, selected, &kept_all, ad_idx)
    };
    // Genotype-index allele number for the kept ALT (0 = ref).
    let alt_n = a_raw + 1;
    let i_00 = tri_pl_index(0, 0);
    let i_0a = tri_pl_index(0, alt_n);
    let i_aa = tri_pl_index(alt_n, alt_n);

    // Per-sample biallelic likelihoods L(00), L(0A), L(AA); None = missing.
    let mut gls: Vec<Option<[f64; 3]>> = Vec::with_capacity(selected.len());
    for &si in selected {
        let sval = rec.samples.get(si).map(|s| s.as_str()).unwrap_or(".");
        let parts: Vec<&str> = sval.split(':').collect();
        let pl = parts.get(pl_idx).copied().and_then(parse_pl_vec);
        match pl {
            Some(v) if v.len() > i_aa => {
                gls.push(Some([
                    pl_to_lik(v[i_00]),
                    pl_to_lik(v[i_0a]),
                    pl_to_lik(v[i_aa]),
                ]));
            }
            _ => gls.push(None),
        }
    }

    let f = em_alt_freq(&gls);
    let prior = [
        (1.0 - f) * (1.0 - f),
        2.0 * f * (1.0 - f),
        f * f,
    ];

    let mut ac = 0usize;
    let mut an = 0usize;
    let mut sample_out: Vec<String> = Vec::with_capacity(selected.len());
    // P(all-ref) accumulated in log10 across samples.
    let mut log10_p_allref = 0.0f64;

    for g in &gls {
        let Some(l) = g else {
            let mut fields = vec!["./.".to_string()];
            if cfg.emit_gq {
                fields.push(".".to_string());
            }
            if cfg.emit_gp {
                fields.push(".,.,.".to_string());
            }
            sample_out.push(fields.join(":"));
            continue;
        };
        let mut post = [l[0] * prior[0], l[1] * prior[1], l[2] * prior[2]];
        let z: f64 = post.iter().sum();
        if z > 0.0 {
            for p in post.iter_mut() {
                *p /= z;
            }
        } else {
            post = [1.0, 0.0, 0.0];
        }
        log10_p_allref += post[0].max(1e-30).log10();

        // argmax posterior; GQ from the gap to the runner-up (phred).
        let mut order = [0usize, 1, 2];
        order.sort_by(|&x, &y| post[y].partial_cmp(&post[x]).unwrap());
        let best = order[0];
        let second = order[1];
        let gq = phred_gap(post[best], post[second]);

        let (a0, a1) = match best {
            0 => (0, 0),
            1 => (0, 1),
            _ => (1, 1),
        };
        an += 2;
        ac += (a0 == 1) as usize + (a1 == 1) as usize;

        let mut fields = vec![format!("{}/{}", a0, a1)];
        if cfg.emit_gq {
            fields.push(gq.to_string());
        }
        if cfg.emit_gp {
            fields.push(format!("{:.6},{:.6},{:.6}", post[0], post[1], post[2]));
        }
        sample_out.push(fields.join(":"));
    }

    if cfg.only_alt && ac == 0 {
        return Ok(None);
    }

    let qual = (-10.0 * log10_p_allref).max(0.0);
    let alt_field = raw_alts.get(a_raw).cloned().unwrap_or_else(|| ".".into());

    let mut info_map = parse_info_map(&rec.info);
    info_map.insert("AN".to_string(), an.to_string());
    if ac == 0 {
        info_map.remove("AC");
    } else {
        info_map.insert("AC".to_string(), ac.to_string());
    }
    let info = render_info_map(&info_map);

    let mut fmt_keys = vec!["GT".to_string()];
    if cfg.emit_gq {
        fmt_keys.push("GQ".to_string());
    }
    if cfg.emit_gp {
        fmt_keys.push("GP".to_string());
    }

    Ok(Some(format!(
        "{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}",
        rec.chrom,
        rec.pos,
        rec.id,
        rec.ref_allele,
        alt_field,
        qual,
        rec.filter,
        info,
        fmt_keys.join(":"),
        sample_out.join("\t")
    )))
}

/// PL (phred) → relative likelihood. PLs are already normalised so the best
/// genotype is 0, hence this yields values in (0, 1].
#[inline]
fn pl_to_lik(pl: u32) -> f64 {
    10f64.powf(-(pl as f64) / 10.0)
}

/// EM estimate of the alt-allele frequency from per-sample biallelic
/// likelihoods. Missing samples are skipped.
fn em_alt_freq(gls: &[Option<[f64; 3]>]) -> f64 {
    let mut f = 0.5f64;
    for _ in 0..20 {
        let mut alt = 0.0f64;
        let mut n = 0.0f64;
        let prior = [(1.0 - f) * (1.0 - f), 2.0 * f * (1.0 - f), f * f];
        for l in gls.iter().flatten() {
            let mut post = [l[0] * prior[0], l[1] * prior[1], l[2] * prior[2]];
            let z: f64 = post.iter().sum();
            if z <= 0.0 {
                continue;
            }
            for p in post.iter_mut() {
                *p /= z;
            }
            alt += post[1] + 2.0 * post[2];
            n += 2.0;
        }
        if n > 0.0 {
            f = (alt / n).clamp(1e-6, 1.0 - 1e-6);
        }
    }
    f
}

/// GQ = phred gap between the best and second-best genotype posteriors, capped 99.
fn phred_gap(best: f64, second: f64) -> u32 {
    if second <= 0.0 {
        return 99;
    }
    let pl_best = -10.0 * best.max(1e-30).log10();
    let pl_second = -10.0 * second.max(1e-30).log10();
    (pl_second - pl_best).round().clamp(0.0, 99.0) as u32
}

/// Pick the ALT (raw index) with the most read support across selected samples,
/// using AD when present, else the first kept ALT.
fn best_alt_by_support(
    rec: &crate::vcf::structs::VcfRecord,
    selected: &[usize],
    kept: &[usize],
    ad_idx: Option<usize>,
) -> usize {
    let Some(ad_idx) = ad_idx else {
        return kept[0];
    };
    let mut best = kept[0];
    let mut best_sup = 0u32;
    for &k in kept {
        let mut sup = 0u32;
        for &si in selected {
            let sval = rec.samples.get(si).map(|s| s.as_str()).unwrap_or(".");
            if let Some(adf) = sval.split(':').nth(ad_idx) {
                let ad = parse_u32_list(adf);
                sup += ad.get(k + 1).copied().unwrap_or(0);
            }
        }
        if sup > best_sup {
            best_sup = sup;
            best = k;
        }
    }
    best
}

fn call_record_mcall(
    rec: &crate::vcf::structs::VcfRecord,
    selected: &[usize],
    cfg: &CallCfg,
) -> Result<Option<String>> {
    use crate::call::{Caller, CallerOpts};
    use crate::call::mcall::{CallSite, CallResult};

    let raw_alts: Vec<String> = rec.alt.split(',').map(str::trim).filter(|a| !a.is_empty() && *a != "." && !(a.starts_with('<') && a.ends_with('>'))).map(|a| a.to_string()).collect();
    let n_als = 1 + raw_alts.len();
    let n_gt = n_als * (n_als + 1) / 2;

    let format_keys: Vec<&str> = rec.format.as_deref().unwrap_or("").split(':').collect();
    let pl_idx = format_keys.iter().position(|k| *k == "PL");
    let n_smpl = selected.len();
    if n_smpl == 0 || pl_idx.is_none() {
        return call_record_legacy(rec, selected, cfg);
    }
    let pl_idx = pl_idx.unwrap();

    let mut pls: Vec<i32> = vec![0; n_smpl * n_gt];
    for (out_i, &si) in selected.iter().enumerate() {
        let sval = rec.samples.get(si).map(|s| s.as_str()).unwrap_or(".");
        let parts: Vec<&str> = sval.split(':').collect();
        let pl_str = parts.get(pl_idx).copied().unwrap_or(".");
        let row = &mut pls[out_i * n_gt..(out_i + 1) * n_gt];
        if pl_str == "." {
            for v in row.iter_mut() { *v = i32::MIN; }
            continue;
        }
        let parsed: Vec<i32> = pl_str.split(',').map(|s| s.parse::<i32>().unwrap_or(i32::MIN)).collect();
        for (i, v) in parsed.iter().enumerate().take(n_gt) { row[i] = *v; }
        for i in parsed.len()..n_gt { row[i] = i32::MIN + 1; }
    }

    let opts = CallerOpts {
        theta: 1.1e-3,
        keep_alts: cfg.mode == CallMode::Consensus,
        variants_only: cfg.only_alt,
        min_ac: 0,
        ploidy: 2,
        ..CallerOpts::default()
    };
    let caller = Caller::new(opts, n_smpl);
    let is_indel = rec.ref_allele.len() > 1 || raw_alts.iter().any(|a| a.len() != rec.ref_allele.len());
    let mut site = CallSite { n_samples: n_smpl, n_alleles: n_als, pls, is_indel, depths: None };

    let result = caller.call_site(&mut site);
    match result {
        CallResult::Skip => Ok(None),
        CallResult::Called { alleles_kept, qual, gts, gqs, pls: _, ac, an } => {
            let out_alts: Vec<String> = alleles_kept.iter().skip(1)
                .filter_map(|&i| raw_alts.get((i - 1) as usize).cloned()).collect();
            let alt_field = if out_alts.is_empty() { ".".to_string() } else { out_alts.join(",") };

            let mut info_map = parse_info_map(&rec.info);
            info_map.insert("AN".to_string(), an.to_string());
            if out_alts.is_empty() { info_map.remove("AC"); }
            else { info_map.insert("AC".to_string(), ac.iter().map(u32::to_string).collect::<Vec<_>>().join(",")); }
            let info = render_info_map(&info_map);

            let mut sample_out: Vec<String> = Vec::with_capacity(n_smpl);
            for (i, &(a, b)) in gts.iter().enumerate() {
                let pa = alleles_kept.iter().position(|x| *x == a).map(|p| p.to_string()).unwrap_or_else(|| ".".into());
                let pb = alleles_kept.iter().position(|x| *x == b).map(|p| p.to_string()).unwrap_or_else(|| ".".into());
                let gt = format!("{}/{}", pa, pb);
                let mut fields = vec![gt];
                if cfg.emit_gq { fields.push(gqs[i].to_string()); }
                sample_out.push(fields.join(":"));
            }

            let qual_s = format!("{:.2}", qual);
            let mut fmt_keys = vec!["GT".to_string()];
            if cfg.emit_gq { fmt_keys.push("GQ".to_string()); }

            Ok(Some(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                rec.chrom, rec.pos, rec.id, rec.ref_allele, alt_field, qual_s, rec.filter, info,
                fmt_keys.join(":"), sample_out.join("\t"))))
        }
    }
}

fn call_record_legacy(
    rec: &crate::vcf::structs::VcfRecord,
    selected: &[usize],
    cfg: &CallCfg,
) -> Result<Option<String>> {
    let raw_alts = rec
        .alt
        .split(',')
        .map(str::trim)
        .filter(|a| !a.is_empty())
        .map(|a| a.to_string())
        .collect::<Vec<_>>();

    let mut kept = Vec::<usize>::new();
    for (i, a) in raw_alts.iter().enumerate() {
        if a == "." || a == "*" {
            continue;
        }
        if a.starts_with('<') && a.ends_with('>') {
            continue;
        }
        kept.push(i);
    }
    if cfg.mode == CallMode::Consensus && kept.len() > 1 {
        kept.truncate(1);
    }

    let out_alts = kept
        .iter()
        .filter_map(|&i| raw_alts.get(i).cloned())
        .collect::<Vec<_>>();
    let alt_field = if out_alts.is_empty() {
        ".".to_string()
    } else {
        out_alts.join(",")
    };

    let format_keys = rec
        .format
        .as_deref()
        .unwrap_or("")
        .split(':')
        .map(|s| s.to_string())
        .collect::<Vec<_>>();
    let key_pos = format_keys
        .iter()
        .enumerate()
        .map(|(i, k)| (k.as_str(), i))
        .collect::<HashMap<_, _>>();

    let pl_idx = key_pos.get("PL").copied();
    let ad_idx = key_pos.get("AD").copied();
    let dp_idx = key_pos.get("DP").copied();
    let passthrough_keys = format_keys
        .iter()
        .filter(|k| k.as_str() != "GT")
        .cloned()
        .collect::<Vec<_>>();

    let mut ac = vec![0usize; out_alts.len()];
    let mut an = 0usize;
    let mut sample_out = Vec::<String>::new();
    let mut gq_sum = 0f64;
    let mut gq_n = 0usize;

    for &si in selected {
        let sval = rec.samples.get(si).map(|s| s.as_str()).unwrap_or(".");
        let parts = sval.split(':').collect::<Vec<_>>();
        let pl = pl_idx
            .and_then(|i| parts.get(i).copied())
            .and_then(parse_pl_vec);
        let ad = ad_idx
            .and_then(|i| parts.get(i).copied())
            .map(parse_u32_list)
            .unwrap_or_default();
        let dp = dp_idx
            .and_then(|i| parts.get(i).copied())
            .and_then(|x| x.parse::<u32>().ok())
            .unwrap_or_else(|| ad.iter().sum());

        let call = infer_genotype(
            pl.as_deref(),
            &ad,
            dp,
            raw_alts.len(),
            &kept,
            out_alts.is_empty(),
        );

        for a in [call.allele0, call.allele1] {
            if let Some(x) = a {
                an += 1;
                if x > 0 {
                    let idx = (x - 1) as usize;
                    if let Some(slot) = ac.get_mut(idx) {
                        *slot += 1;
                    }
                }
            }
        }

        if let Some(v) = call.gq {
            gq_sum += v as f64;
            gq_n += 1;
        }

        let mut sample_fmt = vec![call.gt.clone()];
        for k in &passthrough_keys {
            let v = key_pos
                .get(k.as_str())
                .and_then(|idx| parts.get(*idx).copied())
                .unwrap_or(".");
            sample_fmt.push(v.to_string());
        }
        if cfg.emit_gq && !passthrough_keys.iter().any(|k| k == "GQ") {
            sample_fmt.push(
                call.gq
                    .map(|x| x.to_string())
                    .unwrap_or_else(|| ".".to_string()),
            );
        }
        if cfg.emit_gp && !passthrough_keys.iter().any(|k| k == "GP") {
            sample_fmt.push(call.gp.clone().unwrap_or_else(|| ".,.,.".to_string()));
        }
        sample_out.push(sample_fmt.join(":"));
    }

    if cfg.only_alt && ac.iter().all(|&x| x == 0) {
        return Ok(None);
    }

    if out_alts.is_empty() && cfg.only_alt {
        return Ok(None);
    }

    let mut info_map = parse_info_map(&rec.info);
    info_map.insert("AN".to_string(), an.to_string());
    if out_alts.is_empty() {
        info_map.remove("AC");
    } else {
        info_map.insert(
            "AC".to_string(),
            ac.iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(","),
        );
    }

    let info = render_info_map(&info_map);
    let qual = if gq_n > 0 {
        format!("{:.3}", gq_sum / gq_n as f64)
    } else {
        rec.qual.clone()
    };

    if sample_out.is_empty() {
        return Ok(Some(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            rec.chrom, rec.pos, rec.id, rec.ref_allele, alt_field, qual, rec.filter, info
        )));
    }

    let mut fmt_keys = vec!["GT".to_string()];
    fmt_keys.extend(passthrough_keys.clone());
    if cfg.emit_gq && !fmt_keys.iter().any(|k| k == "GQ") {
        fmt_keys.push("GQ".to_string());
    }
    if cfg.emit_gp && !fmt_keys.iter().any(|k| k == "GP") {
        fmt_keys.push("GP".to_string());
    }

    Ok(Some(format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        rec.chrom,
        rec.pos,
        rec.id,
        rec.ref_allele,
        alt_field,
        qual,
        rec.filter,
        info,
        fmt_keys.join(":"),
        sample_out.join("\t")
    )))
}

#[derive(Clone)]
struct GenotypeCall {
    gt: String,
    allele0: Option<u8>,
    allele1: Option<u8>,
    gq: Option<u32>,
    gp: Option<String>,
}

fn infer_genotype(
    pl: Option<&[u32]>,
    ad: &[u32],
    dp: u32,
    _raw_alt_count: usize,
    kept_alt_idx: &[usize],
    no_alt: bool,
) -> GenotypeCall {
    let out_alt_count = kept_alt_idx.len();
    if dp == 0 {
        return GenotypeCall {
            gt: "./.".to_string(),
            allele0: None,
            allele1: None,
            gq: None,
            gp: None,
        };
    }

    if no_alt {
        return GenotypeCall {
            gt: "0/0".to_string(),
            allele0: Some(0),
            allele1: Some(0),
            gq: Some(99),
            gp: Some("1".to_string()),
        };
    }

    if let Some(pls) = pl {
        let mut raw_by_out = Vec::<usize>::with_capacity(out_alt_count + 1);
        raw_by_out.push(0);
        for &k in kept_alt_idx {
            raw_by_out.push(k + 1);
        }

        let pairs = genotype_pairs(out_alt_count + 1);
        let mut best_pair = (0usize, 0usize);
        let mut best_idx = 0usize;
        let mut best_pl = u32::MAX;
        let mut second = u32::MAX;
        let mut probs = Vec::<f64>::with_capacity(pairs.len());

        for (i, (a_out, b_out)) in pairs.iter().copied().enumerate() {
            let a_raw = raw_by_out[a_out];
            let b_raw = raw_by_out[b_out];
            let idx = tri_pl_index(a_raw, b_raw);
            if idx >= pls.len() {
                probs.push(0.0);
                continue;
            }
            let v = pls[idx];
            let p = 10f64.powf(-(v as f64) / 10.0);
            probs.push(p);
            if v < best_pl {
                second = best_pl;
                best_pl = v;
                best_pair = (a_out, b_out);
                best_idx = i;
            } else if v < second {
                second = v;
            }
        }

        if best_pl == u32::MAX {
            return fallback_ad_call(ad, kept_alt_idx);
        }

        let sum_p: f64 = probs.iter().sum();
        let gp = if sum_p > 0.0 {
            probs
                .iter()
                .map(|p| format!("{:.6}", p / sum_p))
                .collect::<Vec<_>>()
                .join(",")
        } else {
            hard_gp_for_best(out_alt_count + 1, best_idx)
        };

        let gq = second.saturating_sub(best_pl).min(99);
        let gt = format!("{}/{}", best_pair.0, best_pair.1);

        return GenotypeCall {
            gt,
            allele0: Some(best_pair.0 as u8),
            allele1: Some(best_pair.1 as u8),
            gq: Some(gq),
            gp: Some(gp),
        };
    }

    fallback_ad_call(ad, kept_alt_idx)
}

fn fallback_ad_call(ad: &[u32], kept_alt_idx: &[usize]) -> GenotypeCall {
    let ref_dp = ad.first().copied().unwrap_or(0);
    let mut best_alt = 0u32;
    let mut best_out_idx = 1usize;
    for (out_i, &k) in kept_alt_idx.iter().enumerate() {
        if let Some(v) = ad.get(k + 1)
            && *v > best_alt
        {
            best_alt = *v;
            best_out_idx = out_i + 1;
        }
    }
    let sum = ref_dp + best_alt;
    if sum == 0 {
        return GenotypeCall {
            gt: "./.".to_string(),
            allele0: None,
            allele1: None,
            gq: None,
            gp: None,
        };
    }

    let frac = best_alt as f64 / sum as f64;
    let (a0, a1) = if frac > 0.8 {
        (best_out_idx, best_out_idx)
    } else if frac < 0.2 {
        (0usize, 0usize)
    } else {
        (0usize, best_out_idx)
    };
    let gt = format!("{a0}/{a1}");
    let gp = hard_gp_for_best(
        kept_alt_idx.len() + 1,
        genotype_index(kept_alt_idx.len() + 1, a0, a1).unwrap_or(0),
    );

    GenotypeCall {
        gt,
        allele0: Some(a0 as u8),
        allele1: Some(a1 as u8),
        gq: Some(60),
        gp: Some(gp),
    }
}

fn tri_pl_index(a: usize, b: usize) -> usize {
    b * (b + 1) / 2 + a
}

fn genotype_pairs(n_alleles: usize) -> Vec<(usize, usize)> {
    let mut out = Vec::<(usize, usize)>::new();
    for a in 0..n_alleles {
        for b in a..n_alleles {
            out.push((a, b));
        }
    }
    out
}

fn genotype_index(n_alleles: usize, a: usize, b: usize) -> Option<usize> {
    let aa = a.min(b);
    let bb = a.max(b);
    let pairs = genotype_pairs(n_alleles);
    pairs.iter().position(|(x, y)| *x == aa && *y == bb)
}

fn hard_gp_for_best(n_alleles: usize, best_idx: usize) -> String {
    let n = n_alleles * (n_alleles + 1) / 2;
    let mut out = vec!["0".to_string(); n];
    if let Some(v) = out.get_mut(best_idx) {
        *v = "1".to_string();
    }
    out.join(",")
}

fn parse_pl_vec(s: &str) -> Option<Vec<u32>> {
    let vals = s
        .split(',')
        .map(|x| x.trim().parse::<u32>().ok())
        .collect::<Option<Vec<_>>>()?;
    Some(vals)
}

fn parse_u32_list(s: &str) -> Vec<u32> {
    s.split(',')
        .map(|x| x.trim().parse::<u32>().unwrap_or(0))
        .collect()
}

fn parse_info_map(info: &str) -> HashMap<String, String> {
    let mut out = HashMap::<String, String>::new();
    if info == "." || info.is_empty() {
        return out;
    }
    for kv in info.split(';') {
        if kv.is_empty() {
            continue;
        }
        if let Some((k, v)) = kv.split_once('=') {
            out.insert(k.to_string(), v.to_string());
        } else {
            out.insert(kv.to_string(), String::new());
        }
    }
    out
}

fn render_info_map(info: &HashMap<String, String>) -> String {
    if info.is_empty() {
        return ".".to_string();
    }
    let mut keys = info.keys().cloned().collect::<Vec<_>>();
    keys.sort();
    let mut parts = Vec::<String>::with_capacity(keys.len());
    for k in keys {
        let v = info.get(&k).map(|x| x.as_str()).unwrap_or("");
        if v.is_empty() {
            parts.push(k);
        } else {
            parts.push(format!("{k}={v}"));
        }
    }
    parts.join(";")
}
