use anyhow::Result;
use std::collections::HashSet;

use crate::VcfReader;
use crate::cli::args::StatsArgs;

pub fn cmd_stats(args: StatsArgs) -> Result<()> {
    let cfg = parse_stats_args(&args.bcftools_args)?;

    let mut out = String::new();
    out.push_str("# This file was produced by kira-bt stats (native)\n");
    out.push_str("# ID\t[2]id\t[3]file\n");

    for (idx, input) in args.inputs.iter().enumerate() {
        let mut reader = VcfReader::open(input)?;
        let headers = reader.header()?;
        let sample_names = header_samples(&headers);
        let selected = select_samples(&sample_names, cfg.samples.as_deref());
        let mut sample_stats = vec![SampleStats::default(); selected.len()];
        let selected_name_refs: Vec<&str> = selected.iter().map(|(_, n)| n.as_str()).collect();

        let mut s = FileStats::default();
        while let Some(rec) = reader.next_record()? {
            if !matches_region(&cfg.regions, &rec.chrom) {
                continue;
            }
            if !matches_filter_list(&cfg.filters, &rec.filter) {
                continue;
            }
            let vtype = record_type(&rec.ref_allele, &rec.alt, cfg.first_alt_only);
            if !matches_expr(
                &cfg.include_expr,
                &cfg.exclude_expr,
                rec.qual.as_str(),
                vtype,
            ) {
                continue;
            }

            s.records += 1;
            if rec.alt == "." || rec.alt == rec.ref_allele {
                s.no_alt += 1;
            }
            if vtype.has_snp {
                s.snps += 1;
            }
            if vtype.has_mnp {
                s.mnps += 1;
            }
            if vtype.has_indel {
                s.indels += 1;
            }
            if vtype.has_other {
                s.others += 1;
            }
            if vtype.multiallelic {
                s.multiallelic += 1;
            }
            if vtype.multiallelic && vtype.all_snp {
                s.multiallelic_snp += 1;
            }
            s.ts += vtype.ts;
            s.tv += vtype.tv;
            if rec.id != "." && !rec.id.is_empty() {
                s.known_id += 1;
            } else {
                s.novel_id += 1;
            }
            update_sample_stats(&rec.format, &rec.samples, &selected, &mut sample_stats);
        }

        out.push_str(&format!("ID\t{idx}\t{}\n", input.display()));
        out.push_str(&format!(
            "SN\t{idx}\tnumber_of_samples\t{}\n",
            selected.len()
        ));
        out.push_str(&format!("SN\t{idx}\tnumber_of_records\t{}\n", s.records));
        out.push_str(&format!("SN\t{idx}\tnumber_of_no-ALTs\t{}\n", s.no_alt));
        out.push_str(&format!("SN\t{idx}\tnumber_of_SNPs\t{}\n", s.snps));
        out.push_str(&format!("SN\t{idx}\tnumber_of_MNPs\t{}\n", s.mnps));
        out.push_str(&format!("SN\t{idx}\tnumber_of_indels\t{}\n", s.indels));
        out.push_str(&format!("SN\t{idx}\tnumber_of_others\t{}\n", s.others));
        out.push_str(&format!(
            "SN\t{idx}\tnumber_of_multiallelic_sites\t{}\n",
            s.multiallelic
        ));
        out.push_str(&format!(
            "SN\t{idx}\tnumber_of_multiallelic_snp_sites\t{}\n",
            s.multiallelic_snp
        ));

        let ratio = if s.tv == 0 {
            "inf".to_string()
        } else {
            format!("{:.4}", s.ts as f64 / s.tv as f64)
        };
        out.push_str(&format!("TSTV\t{idx}\t{}\t{}\t{}\n", s.ts, s.tv, ratio));

        if cfg.verbose {
            out.push_str(&format!("IDS\t{idx}\tknown\t{}\n", s.known_id));
            out.push_str(&format!("IDS\t{idx}\tnovel\t{}\n", s.novel_id));
        }

        for (i, (_, name)) in selected.iter().enumerate() {
            let st = &sample_stats[i];
            let avg_dp = if st.dp_n == 0 {
                0.0
            } else {
                st.dp_sum as f64 / st.dp_n as f64
            };
            out.push_str(&format!(
                "PSC\t{idx}\t{name}\t{}\t{}\t{}\t{}\t{avg_dp:.2}\n",
                st.hom_ref, st.hom_alt, st.het, st.missing
            ));
        }

        if idx + 1 != args.inputs.len() {
            out.push('\n');
        }

        let _ = selected_name_refs;
    }

    print!("{out}");
    Ok(())
}

#[derive(Default)]
struct StatsCfg {
    samples: Option<String>,
    include_expr: Option<Expr>,
    exclude_expr: Option<Expr>,
    regions: Option<HashSet<String>>,
    filters: Option<HashSet<String>>,
    first_alt_only: bool,
    verbose: bool,
}

#[derive(Clone, Copy)]
enum Expr {
    QualGt(f64),
    TypeSnp,
}

fn parse_stats_args(args: &[String]) -> Result<StatsCfg> {
    let mut cfg = StatsCfg::default();
    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "-s" => {
                i += 1;
                cfg.samples = Some(
                    args.get(i)
                        .ok_or_else(|| anyhow::anyhow!("missing value for -s"))?
                        .to_string(),
                );
            }
            "-i" => {
                i += 1;
                cfg.include_expr = Some(parse_expr(
                    args.get(i)
                        .ok_or_else(|| anyhow::anyhow!("missing value for -i"))?,
                )?);
            }
            "-e" => {
                i += 1;
                cfg.exclude_expr = Some(parse_expr(
                    args.get(i)
                        .ok_or_else(|| anyhow::anyhow!("missing value for -e"))?,
                )?);
            }
            "-r" => {
                i += 1;
                let raw = args
                    .get(i)
                    .ok_or_else(|| anyhow::anyhow!("missing value for -r"))?;
                let mut set = HashSet::new();
                for item in raw.split(',') {
                    let chrom = item.split(':').next().unwrap_or(item).to_string();
                    if !chrom.is_empty() {
                        set.insert(chrom);
                    }
                }
                cfg.regions = Some(set);
            }
            "-f" => {
                i += 1;
                let raw = args
                    .get(i)
                    .ok_or_else(|| anyhow::anyhow!("missing value for -f"))?;
                let set = raw
                    .split(',')
                    .filter(|s| !s.is_empty())
                    .map(|s| s.to_string())
                    .collect::<HashSet<_>>();
                cfg.filters = Some(set);
            }
            "-1" => cfg.first_alt_only = true,
            "-v" => cfg.verbose = true,
            _ => {}
        }
        i += 1;
    }
    Ok(cfg)
}

fn parse_expr(s: &str) -> Result<Expr> {
    let t = s.trim();
    if let Some(v) = t.strip_prefix("QUAL>") {
        let q = v.trim().parse::<f64>()?;
        return Ok(Expr::QualGt(q));
    }
    if t == "type=\"snp\"" {
        return Ok(Expr::TypeSnp);
    }
    anyhow::bail!("unsupported expression in native stats: {t}")
}

fn header_samples(headers: &[String]) -> Vec<String> {
    for h in headers {
        if h.starts_with("#CHROM\t") {
            let cols: Vec<&str> = h.split('\t').collect();
            if cols.len() > 9 {
                return cols[9..].iter().map(|s| s.to_string()).collect();
            }
        }
    }
    Vec::new()
}

fn select_samples(all: &[String], raw: Option<&str>) -> Vec<(usize, String)> {
    match raw {
        None | Some("-") => all
            .iter()
            .enumerate()
            .map(|(i, s)| (i, s.clone()))
            .collect(),
        Some(v) => {
            let wanted = v
                .split(',')
                .filter(|s| !s.is_empty())
                .collect::<HashSet<_>>();
            all.iter()
                .enumerate()
                .filter(|(_, s)| wanted.contains(s.as_str()))
                .map(|(i, s)| (i, s.clone()))
                .collect()
        }
    }
}

#[derive(Default, Clone)]
struct FileStats {
    records: u64,
    no_alt: u64,
    snps: u64,
    mnps: u64,
    indels: u64,
    others: u64,
    multiallelic: u64,
    multiallelic_snp: u64,
    ts: u64,
    tv: u64,
    known_id: u64,
    novel_id: u64,
}

#[derive(Default, Clone)]
struct SampleStats {
    hom_ref: u64,
    hom_alt: u64,
    het: u64,
    missing: u64,
    dp_sum: u64,
    dp_n: u64,
}

#[derive(Clone, Copy)]
struct RecordType {
    has_snp: bool,
    has_mnp: bool,
    has_indel: bool,
    has_other: bool,
    multiallelic: bool,
    all_snp: bool,
    ts: u64,
    tv: u64,
}

fn record_type(r: &str, alt: &str, first_alt_only: bool) -> RecordType {
    let mut alts: Vec<&str> = alt.split(',').collect();
    if first_alt_only && !alts.is_empty() {
        alts.truncate(1);
    }
    let multiallelic = alt.contains(',');
    let mut has_snp = false;
    let mut has_mnp = false;
    let mut has_indel = false;
    let mut has_other = false;
    let mut all_snp = !alts.is_empty();
    let mut ts = 0u64;
    let mut tv = 0u64;

    for a in &alts {
        if *a == "." || *a == r {
            all_snp = false;
            continue;
        }
        if is_snp(r, a) {
            has_snp = true;
            if is_transition(r.as_bytes()[0], a.as_bytes()[0]) {
                ts += 1;
            } else {
                tv += 1;
            }
        } else {
            all_snp = false;
            if r.len() == a.len() && r.len() > 1 {
                has_mnp = true;
            } else if r.len() != a.len() {
                has_indel = true;
            } else {
                has_other = true;
            }
        }
    }

    RecordType {
        has_snp,
        has_mnp,
        has_indel,
        has_other,
        multiallelic,
        all_snp,
        ts,
        tv,
    }
}

fn is_snp(r: &str, a: &str) -> bool {
    r.len() == 1
        && a.len() == 1
        && r.as_bytes()[0].is_ascii_alphabetic()
        && a.as_bytes()[0].is_ascii_alphabetic()
}

fn is_transition(r: u8, a: u8) -> bool {
    matches!(
        (r.to_ascii_uppercase(), a.to_ascii_uppercase()),
        (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C')
    )
}

fn matches_region(regions: &Option<HashSet<String>>, chrom: &str) -> bool {
    match regions {
        None => true,
        Some(set) => set.contains(chrom),
    }
}

fn matches_filter_list(filters: &Option<HashSet<String>>, filter: &str) -> bool {
    match filters {
        None => true,
        Some(set) => set.contains(filter),
    }
}

fn matches_expr(include: &Option<Expr>, exclude: &Option<Expr>, qual: &str, t: RecordType) -> bool {
    if let Some(expr) = include {
        if !eval_expr(*expr, qual, t) {
            return false;
        }
    }
    if let Some(expr) = exclude {
        if eval_expr(*expr, qual, t) {
            return false;
        }
    }
    true
}

fn eval_expr(expr: Expr, qual: &str, t: RecordType) -> bool {
    match expr {
        Expr::QualGt(v) => qual.parse::<f64>().ok().map(|q| q > v).unwrap_or(false),
        Expr::TypeSnp => t.has_snp,
    }
}

fn update_sample_stats(
    format: &Option<String>,
    samples: &[String],
    selected: &[(usize, String)],
    sample_stats: &mut [SampleStats],
) {
    let Some(fmt) = format else {
        return;
    };
    let keys: Vec<&str> = fmt.split(':').collect();
    let gt_idx = keys.iter().position(|k| *k == "GT");
    let dp_idx = keys.iter().position(|k| *k == "DP");

    for (out_idx, (src_idx, _name)) in selected.iter().enumerate() {
        let Some(sample) = samples.get(*src_idx) else {
            continue;
        };
        let vals: Vec<&str> = sample.split(':').collect();
        if let Some(i) = dp_idx {
            if let Some(v) = vals.get(i).and_then(|x| x.parse::<u64>().ok()) {
                sample_stats[out_idx].dp_sum += v;
                sample_stats[out_idx].dp_n += 1;
            }
        }
        if let Some(i) = gt_idx {
            if let Some(gt) = vals.get(i) {
                classify_gt(gt, &mut sample_stats[out_idx]);
            } else {
                sample_stats[out_idx].missing += 1;
            }
        } else {
            sample_stats[out_idx].missing += 1;
        }
    }
}

fn classify_gt(gt: &str, s: &mut SampleStats) {
    if gt == "." || gt == "./." || gt == ".|." || gt.contains('.') {
        s.missing += 1;
        return;
    }
    let sep = if gt.contains('|') { '|' } else { '/' };
    let alleles: Vec<&str> = gt.split(sep).collect();
    if alleles.is_empty() {
        s.missing += 1;
        return;
    }
    let mut parsed = Vec::with_capacity(alleles.len());
    for a in &alleles {
        if let Ok(v) = a.parse::<u32>() {
            parsed.push(v);
        } else {
            s.missing += 1;
            return;
        }
    }
    if parsed.iter().all(|v| *v == 0) {
        s.hom_ref += 1;
    } else if parsed.iter().all(|v| *v == parsed[0]) {
        s.hom_alt += 1;
    } else {
        s.het += 1;
    }
}
