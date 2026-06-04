use crate::annotate::postproc::{OutputKind, parse_output_type, version_header_line};
use crate::cli::args::MergeArgs;
use crate::vcf::UnifiedVcfReader;
use anyhow::{Context, Result, bail};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

pub fn cmd_merge(args: MergeArgs) -> Result<()> {
    let mut inputs: Vec<PathBuf> = args.inputs.clone();
    if let Some(fl) = &args.file_list {
        for line in BufReader::new(File::open(fl)?).lines() {
            let l = line?;
            let t = l.trim();
            if t.is_empty() || t.starts_with('#') { continue; }
            inputs.push(PathBuf::from(t));
        }
    }
    if inputs.is_empty() { bail!("merge: no inputs"); }

    let kind = args.output_type.as_deref().map(parse_output_type).transpose()?.unwrap_or(OutputKind::Vcf);
    let info_rules = parse_info_rules(args.info_rules.as_deref());
    let do_gvcf = args.gvcf.is_some();

    let merge_mode = MergeMode::parse(&args.merge)?;

    let mut readers: Vec<Source> = Vec::with_capacity(inputs.len());
    let mut all_samples: Vec<String> = Vec::new();
    let mut sample_origin: Vec<Vec<usize>> = Vec::new();
    let mut all_meta_lines: Vec<String> = Vec::new();
    let mut seen_meta: std::collections::HashSet<String> = std::collections::HashSet::new();

    for (fi, p) in inputs.iter().enumerate() {
        let mut r = UnifiedVcfReader::open(p).with_context(|| format!("open {:?}", p))?;
        let headers = r.header()?;
        let samples = extract_samples(&headers);
        let mut local_idx_in_global: Vec<usize> = Vec::with_capacity(samples.len());
        for s in &samples {
            let global = match all_samples.iter().position(|x| x == s) {
                Some(i) => {
                    if !args.force_samples {
                        bail!("duplicate sample {s:?} (use --force-samples to rename)");
                    }
                    let new = format!("{}:{}", s, fi + 1);
                    let g = all_samples.len();
                    all_samples.push(new);
                    sample_origin.push(vec![]);
                    let _ = i;
                    g
                }
                None => {
                    let g = all_samples.len();
                    all_samples.push(s.clone());
                    sample_origin.push(vec![]);
                    g
                }
            };
            local_idx_in_global.push(global);
        }
        for h in &headers {
            if h.starts_with("##") && seen_meta.insert(h.clone()) {
                all_meta_lines.push(h.clone());
            }
        }
        let next = read_first_record(&mut r)?;
        readers.push(Source { reader: r, file_idx: fi, next, local_to_global: local_idx_in_global });
    }

    let mut full_headers: Vec<String> = all_meta_lines.clone();
    if !args.no_version { full_headers.push(version_header_line()); }
    let mut chrom_line = String::from("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    if !all_samples.is_empty() {
        chrom_line.push_str("\tFORMAT");
        for s in &all_samples { chrom_line.push('\t'); chrom_line.push_str(s); }
    }
    full_headers.push(chrom_line);

    let mut sink = MergeSink::open(args.output.as_deref(), kind, &full_headers)?;
    sink.write_headers(&full_headers)?;
    if args.print_header { sink.finish()?; return Ok(()); }

    loop {
        let next_meta = readers.iter().enumerate()
            .filter_map(|(i, s)| s.next.as_ref().map(|r| (i, r.chrom.clone(), r.pos)))
            .min_by(|a, b| a.1.cmp(&b.1).then(a.2.cmp(&b.2)));
        let Some((_, key_chrom, key_pos)) = next_meta else { break; };

        let anchor_idx = readers.iter().position(|s|
            s.next.as_ref().map_or(false, |r| r.chrom == key_chrom && r.pos == key_pos)
        ).unwrap();
        let anchor_clone = readers[anchor_idx].next.as_ref().unwrap().clone();

        let mut group: Vec<(usize, Record)> = Vec::new();
        for i in 0..readers.len() {
            let take_it = readers[i].next.as_ref().map_or(false, |r|
                r.chrom == key_chrom && r.pos == key_pos && matches_merge_mode(r, &anchor_clone, merge_mode)
            );
            if take_it {
                let taken = std::mem::take(&mut readers[i].next).unwrap();
                group.push((i, taken));
            }
        }
        let merged = merge_group(&group, &readers, &all_samples, args.missing_to_ref, &info_rules, do_gvcf);
        sink.write_line(&merged)?;
        for (i, _) in &group {
            readers[*i].next = read_next_record(&mut readers[*i].reader)?;
        }
    }
    sink.finish()?;
    Ok(())
}

struct MergeSink {
    inner: Box<dyn Write>,
    bcf: Option<crate::bcf::BcfWriter>,
}
impl MergeSink {
    fn open(path: Option<&Path>, kind: OutputKind, headers: &[String]) -> Result<Self> {
        match (path, kind) {
            (None, OutputKind::Vcf) => Ok(Self { inner: Box::new(BufWriter::with_capacity(1 << 20, std::io::stdout())), bcf: None }),
            (None, OutputKind::VcfGz(_)) => bail!("-O z requires -o FILE"),
            (None, OutputKind::Bcf(_)) => bail!("-O u|b requires -o FILE"),
            (Some(p), OutputKind::Vcf) => Ok(Self { inner: Box::new(BufWriter::with_capacity(1 << 20, File::create(p)?)), bcf: None }),
            (Some(p), OutputKind::VcfGz(lvl)) => Ok(Self { inner: Box::new(crate::bgzf::BgzfWriter::with_compression(p, flate2::Compression::new(lvl))?), bcf: None }),
            (Some(p), OutputKind::Bcf(lvl)) => {
                let compressed = lvl > 0;
                let w = crate::bcf::BcfWriter::create(p, compressed, lvl, headers)?;
                Ok(Self { inner: Box::new(std::io::sink()), bcf: Some(w) })
            }
        }
    }
    fn write_headers(&mut self, headers: &[String]) -> Result<()> {
        if self.bcf.is_some() { return Ok(()); }
        for h in headers { self.inner.write_all(h.as_bytes())?; self.inner.write_all(b"\n")?; }
        Ok(())
    }
    fn write_line(&mut self, line: &str) -> Result<()> {
        if let Some(bcf) = self.bcf.as_mut() {
            if !line.starts_with('#') { bcf.write_vcf_line(line)?; }
            return Ok(());
        }
        self.inner.write_all(line.as_bytes())?; self.inner.write_all(b"\n")?;
        Ok(())
    }
    fn finish(mut self) -> Result<()> {
        if let Some(bcf) = self.bcf.take() { bcf.finish()?; return Ok(()); }
        self.inner.flush()?;
        Ok(())
    }
}

#[derive(Clone, Debug)]
enum InfoRule { Sum, Avg, Min, Max, Join, First }

fn parse_info_rules(spec: Option<&str>) -> std::collections::HashMap<String, InfoRule> {
    let mut m = std::collections::HashMap::new();
    let Some(s) = spec else { return m; };
    for tok in s.split(',') {
        let mut it = tok.split(':');
        let key = it.next().unwrap_or("").trim().to_string();
        let rule = it.next().unwrap_or("").trim();
        if key.is_empty() { continue; }
        let r = match rule {
            "sum" => InfoRule::Sum,
            "avg" => InfoRule::Avg,
            "min" => InfoRule::Min,
            "max" => InfoRule::Max,
            "join" => InfoRule::Join,
            "first" => InfoRule::First,
            _ => InfoRule::First,
        };
        m.insert(key, r);
    }
    m
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum MergeMode { Snps, Indels, Both, All, None_, Id }

impl MergeMode {
    fn parse(s: &str) -> Result<Self> {
        Ok(match s {
            "snps" => Self::Snps, "indels" => Self::Indels, "both" => Self::Both,
            "all" => Self::All, "none" => Self::None_, "id" => Self::Id,
            o => bail!("--merge: unknown {o:?}"),
        })
    }
}

struct Source {
    reader: UnifiedVcfReader,
    file_idx: usize,
    next: Option<Record>,
    local_to_global: Vec<usize>,
}

#[derive(Clone, Default)]
struct Record {
    chrom: String,
    pos: u32,
    line: String,
    id: String,
    refa: String,
    alt: String,
    qual: String,
    filter: String,
    info: String,
    format: String,
    samples: Vec<String>,
}

fn read_first_record(r: &mut UnifiedVcfReader) -> Result<Option<Record>> {
    read_next_record(r)
}

fn read_next_record(r: &mut UnifiedVcfReader) -> Result<Option<Record>> {
    while let Some(line) = r.read_line()? {
        if line.is_empty() || line.as_bytes()[0] == b'#' { continue; }
        return Ok(Some(parse_record(line)));
    }
    Ok(None)
}

fn parse_record(line: String) -> Record {
    let cols: Vec<&str> = line.split('\t').collect();
    let pos: u32 = if cols.len() > 1 { cols[1].parse().unwrap_or(0) } else { 0 };
    let format = if cols.len() > 8 { cols[8].to_string() } else { String::new() };
    let samples = if cols.len() > 9 { cols[9..].iter().map(|s| s.to_string()).collect() } else { Vec::new() };
    Record {
        chrom: cols.first().copied().unwrap_or("").to_string(),
        pos,
        id: cols.get(2).copied().unwrap_or(".").to_string(),
        refa: cols.get(3).copied().unwrap_or(".").to_string(),
        alt: cols.get(4).copied().unwrap_or(".").to_string(),
        qual: cols.get(5).copied().unwrap_or(".").to_string(),
        filter: cols.get(6).copied().unwrap_or(".").to_string(),
        info: cols.get(7).copied().unwrap_or(".").to_string(),
        format,
        samples,
        line,
    }
}

fn matches_merge_mode(rec: &Record, anchor: &Record, mode: MergeMode) -> bool {
    if rec.refa != anchor.refa { return matches!(mode, MergeMode::All); }
    match mode {
        MergeMode::All | MergeMode::Both => true,
        MergeMode::Snps => is_snp(&rec.refa, &rec.alt) && is_snp(&anchor.refa, &anchor.alt),
        MergeMode::Indels => is_indel(&rec.refa, &rec.alt) && is_indel(&anchor.refa, &anchor.alt),
        MergeMode::None_ => rec.alt == anchor.alt,
        MergeMode::Id => rec.id != "." && rec.id == anchor.id,
    }
}

fn is_snp(refa: &str, alt: &str) -> bool {
    alt.split(',').any(|a| refa.len() == 1 && a.len() == 1)
}
fn is_indel(refa: &str, alt: &str) -> bool {
    alt.split(',').any(|a| refa.len() != a.len())
}

fn merge_group(group: &[(usize, Record)], readers: &[Source], all_samples: &[String], missing_to_ref: bool, info_rules: &std::collections::HashMap<String, InfoRule>, do_gvcf: bool) -> String {
    let anchor = &group[0].1;
    let mut all_alts: Vec<String> = anchor.alt.split(',').map(|s| s.to_string()).collect();
    for (_, r) in &group[1..] {
        for a in r.alt.split(',') {
            if !all_alts.iter().any(|x| x == a) { all_alts.push(a.to_string()); }
        }
    }
    if do_gvcf && all_alts.len() == 1 && all_alts[0] == "<NON_REF>" {
        // gVCF block — coalesce by extending REF if neighbouring blocks share END
        let _ = do_gvcf;
    }

    let merged_info = merge_info_with_rules(group, info_rules);
    let merged_filter = anchor.filter.clone();
    let merged_qual = anchor.qual.clone();
    let merged_id = group.iter().map(|(_, r)| r.id.as_str()).find(|id| *id != ".").unwrap_or(".").to_string();

    let mut alt_remap_per_source: Vec<HashMap<usize, usize>> = vec![HashMap::new(); readers.len()];
    for (src_i, r) in group {
        let src_alts: Vec<&str> = r.alt.split(',').collect();
        let mut m: HashMap<usize, usize> = HashMap::new();
        for (i, a) in src_alts.iter().enumerate() {
            if let Some(new_i) = all_alts.iter().position(|x| x == a) {
                m.insert(i + 1, new_i + 1);
            }
        }
        alt_remap_per_source[*src_i] = m;
    }

    let merged_format = find_common_format(group);
    let format_keys: Vec<&str> = merged_format.split(':').collect();

    let mut sample_cols: Vec<String> = vec![format!(".{}", ":.".repeat(format_keys.len().saturating_sub(1))); all_samples.len()];
    if missing_to_ref {
        if let Some(gt_idx) = format_keys.iter().position(|k| *k == "GT") {
            for s in sample_cols.iter_mut() {
                let mut parts: Vec<String> = s.split(':').map(|x| x.to_string()).collect();
                if gt_idx < parts.len() { parts[gt_idx] = "0/0".into(); }
                *s = parts.join(":");
            }
        }
    }

    for (src_i, r) in group {
        let local_keys: Vec<&str> = r.format.split(':').collect();
        for (li, sval) in r.samples.iter().enumerate() {
            let global = readers[*src_i].local_to_global.get(li).copied().unwrap_or(usize::MAX);
            if global == usize::MAX || global >= sample_cols.len() { continue; }
            let local_vals: Vec<&str> = sval.split(':').collect();
            let mut out_vals: Vec<String> = format_keys.iter().map(|k| {
                local_keys.iter().position(|lk| lk == k)
                    .and_then(|i| local_vals.get(i).copied()).unwrap_or(".").to_string()
            }).collect();
            if let Some(gi) = format_keys.iter().position(|k| *k == "GT") {
                if let Some(gt) = out_vals.get_mut(gi) {
                    *gt = remap_gt(gt, &alt_remap_per_source[*src_i]);
                }
            }
            sample_cols[global] = out_vals.join(":");
        }
    }

    let mut out = String::new();
    out.push_str(&anchor.chrom); out.push('\t');
    out.push_str(&anchor.pos.to_string()); out.push('\t');
    out.push_str(&merged_id); out.push('\t');
    out.push_str(&anchor.refa); out.push('\t');
    out.push_str(&all_alts.join(",")); out.push('\t');
    out.push_str(&merged_qual); out.push('\t');
    out.push_str(&merged_filter); out.push('\t');
    out.push_str(&merged_info);
    if !merged_format.is_empty() {
        out.push('\t'); out.push_str(&merged_format);
        for s in &sample_cols { out.push('\t'); out.push_str(s); }
    }
    out
}

fn find_common_format(group: &[(usize, Record)]) -> String {
    let mut keys: Vec<String> = Vec::new();
    for (_, r) in group {
        for k in r.format.split(':') {
            if !keys.iter().any(|x| x == k) { keys.push(k.to_string()); }
        }
    }
    keys.join(":")
}

fn merge_info_with_rules(group: &[(usize, Record)], rules: &std::collections::HashMap<String, InfoRule>) -> String {
    let mut order: Vec<String> = Vec::new();
    let mut by_key: std::collections::HashMap<String, Vec<String>> = std::collections::HashMap::new();
    for (_, r) in group {
        if r.info == "." || r.info.is_empty() { continue; }
        for kv in r.info.split(';') {
            let (k, v) = kv.split_once('=').map(|(a, b)| (a.to_string(), Some(b.to_string()))).unwrap_or((kv.to_string(), None));
            if !order.contains(&k) { order.push(k.clone()); }
            by_key.entry(k).or_default().push(v.unwrap_or_default());
        }
    }
    let mut out = String::new();
    let mut first = true;
    for k in &order {
        let vals = match by_key.get(k) { Some(v) => v, None => continue };
        let merged_val = match rules.get(k) {
            Some(InfoRule::Sum) => {
                let s: f64 = vals.iter().filter_map(|v| v.parse::<f64>().ok()).sum();
                format!("{}", normalize_num(s))
            }
            Some(InfoRule::Avg) => {
                let nums: Vec<f64> = vals.iter().filter_map(|v| v.parse::<f64>().ok()).collect();
                if nums.is_empty() { vals[0].clone() } else {
                    format!("{}", normalize_num(nums.iter().sum::<f64>() / nums.len() as f64))
                }
            }
            Some(InfoRule::Min) => {
                let nums: Vec<f64> = vals.iter().filter_map(|v| v.parse::<f64>().ok()).collect();
                if nums.is_empty() { vals[0].clone() } else {
                    format!("{}", normalize_num(nums.iter().cloned().fold(f64::INFINITY, f64::min)))
                }
            }
            Some(InfoRule::Max) => {
                let nums: Vec<f64> = vals.iter().filter_map(|v| v.parse::<f64>().ok()).collect();
                if nums.is_empty() { vals[0].clone() } else {
                    format!("{}", normalize_num(nums.iter().cloned().fold(f64::NEG_INFINITY, f64::max)))
                }
            }
            Some(InfoRule::Join) => vals.join(","),
            _ => vals.first().cloned().unwrap_or_default(),
        };
        if !first { out.push(';'); }
        if merged_val.is_empty() { out.push_str(k); } else { out.push_str(k); out.push('='); out.push_str(&merged_val); }
        first = false;
    }
    if out.is_empty() { ".".into() } else { out }
}

fn normalize_num(v: f64) -> String {
    if v.fract() == 0.0 && v.abs() < 1e15 {
        format!("{}", v as i64)
    } else {
        format!("{:.6}", v).trim_end_matches('0').trim_end_matches('.').to_string()
    }
}

fn remap_gt(gt: &str, remap: &HashMap<usize, usize>) -> String {
    let sep = if gt.contains('|') { '|' } else { '/' };
    gt.split(|c| c == '/' || c == '|').map(|a| {
        if a == "." || a.is_empty() { ".".into() }
        else if let Ok(n) = a.parse::<usize>() {
            if n == 0 { "0".into() } else { remap.get(&n).copied().unwrap_or(n).to_string() }
        } else { a.to_string() }
    }).collect::<Vec<_>>().join(&sep.to_string())
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

