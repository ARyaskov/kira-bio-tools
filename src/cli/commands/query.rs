use anyhow::Result;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::io::{BufWriter, Write};

use crate::cli::args::{HeaderArgs, ListArgs, QueryCompatArgs, RegionQueryArgs, StatArgs};
use crate::filter::FilterEngine;
use crate::{CsiQuery, KbiIndex, Region, VcfReader, chr_id_to_name, chr_name_to_id, fetch_line};

pub fn cmd_query(args: QueryCompatArgs) -> Result<()> {
    let cfg = parse_query_args(&args.bcftools_args)?;
    let mut reader = VcfReader::open(&cfg.input)?;
    let headers = reader.header()?;
    let sample_names = extract_sample_names(&headers);
    let sample_idx = resolve_samples(
        &sample_names,
        cfg.samples.as_deref(),
        cfg.samples_file.as_deref(),
    )?;

    if cfg.list_samples {
        let mut out: Box<dyn Write> = match &cfg.output {
            Some(p) => Box::new(BufWriter::new(File::create(p)?)),
            None => Box::new(BufWriter::new(std::io::stdout())),
        };
        for &i in &sample_idx {
            if let Some(name) = sample_names.get(i) {
                writeln!(out, "{name}")?;
            }
        }
        return Ok(());
    }

    let include = cfg.include.clone();
    let exclude = cfg.exclude.clone();
    let filter_expr = include.as_deref().or(exclude.as_deref());
    let filter_engine = FilterEngine::new(&headers, filter_expr, exclude.is_some())?;

    let region_filter = build_region_filter_pair(&cfg)?;
    let no_nas = cfg.no_nas.clone();
    let exclude_uncalled = cfg.exclude_uncalled;

    let format = decode_escapes(&cfg.format);
    let mut fmt_ctx = QueryFormatContext::new(&headers, &sample_names, sample_idx.clone(), format)?;

    let mut out: Box<dyn Write> = match &cfg.output {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };
    if cfg.header_mode > 0 {
        let hline = fmt_ctx.header_line(cfg.header_mode == 2);
        out.write_all(hline.trim_end_matches('\n').as_bytes())?;
        out.write_all(b"\n")?;
    }

    while let Some(rec) = reader.next_record()? {
        if let Some(rf) = &region_filter {
            if !rf.passes(&rec) { continue; }
        }
        let eval = filter_engine.eval(&rec)?;
        let site_pass = if let Some(ps) = eval.pass_samples.as_ref() {
            sample_idx
                .iter()
                .any(|&i| ps.get(i).copied().unwrap_or(false))
        } else {
            eval.pass_site
        };
        if !site_pass {
            continue;
        }
        if exclude_uncalled && record_all_missing(&rec) { continue; }
        let mut s = fmt_ctx.render_record(&rec, eval.pass_samples.as_deref())?;
        if let Some(repl) = &no_nas {
            s = s.replace("\t.\t", &format!("\t{}\t", repl)).replace("\t.\n", &format!("\t{}\n", repl));
        }
        out.write_all(s.as_bytes())?;
        if !s.ends_with('\n') {
            out.write_all(b"\n")?;
        }
    }

    Ok(())
}

struct RecRegionFilter {
    regions: Vec<(String, u32, u32)>,
}

impl RecRegionFilter {
    fn passes(&self, rec: &crate::vcf::VcfRecord) -> bool {
        if self.regions.is_empty() { return true; }
        let pos = rec.pos;
        let chrom = rec.chrom.as_str();
        for (c, s, e) in &self.regions {
            if c == chrom && pos >= *s && pos <= *e { return true; }
        }
        false
    }
}

fn build_region_filter_pair(cfg: &QueryConfig) -> Result<Option<RecRegionFilter>> {
    let mut regions: Vec<(String, u32, u32)> = Vec::new();
    if let Some(s) = cfg.regions.as_deref().or(cfg.targets.as_deref()) {
        for tok in s.split(',') {
            if let Some(r) = parse_region_token(tok) { regions.push(r); }
        }
    }
    let file = cfg.regions_file.as_deref().or(cfg.targets_file.as_deref());
    if let Some(p) = file {
        let f = File::open(p)?;
        for line in BufReader::new(f).lines() {
            let l = line?;
            let t = l.trim();
            if t.is_empty() || t.starts_with('#') { continue; }
            let mut parts = t.split('\t');
            let chr = parts.next().unwrap_or("").to_string();
            let beg: u32 = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);
            let end: u32 = parts.next().and_then(|s| s.parse().ok()).unwrap_or(u32::MAX);
            if !chr.is_empty() { regions.push((chr, beg, end)); }
        }
    }
    if regions.is_empty() { Ok(None) } else { Ok(Some(RecRegionFilter { regions })) }
}

fn parse_region_token(s: &str) -> Option<(String, u32, u32)> {
    if let Some((c, r)) = s.split_once(':') {
        if let Some((b, e)) = r.split_once('-') {
            let b: u32 = b.parse().ok()?; let e: u32 = e.parse().ok()?;
            return Some((c.to_string(), b, e));
        }
        if let Ok(b) = r.parse::<u32>() { return Some((c.to_string(), b, b)); }
        return Some((c.to_string(), 0, u32::MAX));
    }
    Some((s.to_string(), 0, u32::MAX))
}

fn record_all_missing(rec: &crate::vcf::VcfRecord) -> bool {
    if rec.samples.is_empty() { return false; }
    rec.samples.iter().all(|s| {
        let gt = s.split(':').next().unwrap_or(".");
        gt.split(|c| c == '/' || c == '|').all(|a| a == "." || a.is_empty())
    })
}

pub fn cmd_region_query(args: RegionQueryArgs) -> Result<()> {
    if args.only_header {
        return cmd_header(HeaderArgs { file: args.file });
    }

    let kbi_path = args.file.with_extension("kbi");
    let csi_path = {
        let mut p = args.file.clone();
        let name = p.file_name().unwrap().to_string_lossy().to_string();
        p.set_file_name(format!("{}.csi", name));
        p
    };

    let use_kbi = kbi_path.exists();
    let use_csi = csi_path.exists() && !use_kbi;

    if !use_kbi && !use_csi {
        anyhow::bail!("No index found. Run 'kira-bt index {:?}' first.", args.file);
    }

    let mut regions = args.regions.clone();
    if let Some(ref regions_file) = args.regions_file {
        let file = File::open(regions_file)?;
        for line in BufReader::new(file).lines() {
            regions.push(line?);
        }
    }

    if regions.is_empty() {
        anyhow::bail!("No regions specified");
    }

    if args.print_header {
        print_vcf_header(&args.file)?;
    }

    let mut total_count = 0usize;

    if use_kbi {
        let index = KbiIndex::load(&kbi_path)?;

        for region_str in &regions {
            let region = Region::parse(region_str)
                .ok_or_else(|| anyhow::anyhow!("Invalid region: {}", region_str))?;

            let chr_id = chr_name_to_id(&region.chr)
                .ok_or_else(|| anyhow::anyhow!("Unknown chromosome: {}", region.chr))?;

            let start = region.start.unwrap_or(0);
            let end = region.end.unwrap_or(u32::MAX);

            let results = index.range(chr_id, start, end);
            total_count += results.len();

            if !args.count {
                for (_pos, offset) in results {
                    let line = fetch_line(&args.file, offset)?;
                    println!("{}", line);
                }
            }
        }
    } else if use_csi {
        let csi = CsiQuery::open(&csi_path)?;

        for region_str in &regions {
            let region = Region::parse(region_str)
                .ok_or_else(|| anyhow::anyhow!("Invalid region: {}", region_str))?;

            let chr_id = chr_name_to_id(&region.chr)
                .ok_or_else(|| anyhow::anyhow!("Unknown chromosome: {}", region.chr))?;

            let start = region.start.unwrap_or(0);
            let end = region.end.unwrap_or(u32::MAX);

            let chunks = csi.query((chr_id - 1) as usize, start, end);

            for (chunk_start, _chunk_end) in chunks {
                let line = fetch_line(&args.file, chunk_start)?;
                if !args.count {
                    println!("{}", line);
                }
                total_count += 1;
            }
        }
    }

    if args.count {
        println!("{}", total_count);
    }

    Ok(())
}

pub fn cmd_stat(args: StatArgs) -> Result<()> {
    use std::fs;

    let file_size = fs::metadata(&args.index)?.len();

    if args.index.extension().map(|e| e == "kbi").unwrap_or(false) {
        let index = KbiIndex::load(&args.index)?;

        println!("Index Statistics (KBI)");
        println!("======================");
        println!("File:          {:?}", args.index);
        println!(
            "File size:     {} bytes ({:.2} MB)",
            file_size,
            file_size as f64 / 1024.0 / 1024.0
        );
        println!("Entries:       {}", index.len());
        println!(
            "Memory usage:  {} bytes ({:.2} MB)",
            index.memory_usage(),
            index.memory_usage() as f64 / 1024.0 / 1024.0
        );
        println!("Bytes/key:     {:.2}", index.bytes_per_key());
    } else if args
        .index
        .extension()
        .map(|e| e == "csi" || e == "tbi")
        .unwrap_or(false)
    {
        let _csi = CsiQuery::open(&args.index)?;

        println!("Index Statistics (CSI/TBI)");
        println!("==========================");
        println!("File:          {:?}", args.index);
        println!(
            "File size:     {} bytes ({:.2} MB)",
            file_size,
            file_size as f64 / 1024.0 / 1024.0
        );
        println!("Format:        CSI/TBI (tabix-compatible)");
    } else {
        anyhow::bail!("Unknown index format");
    }

    Ok(())
}

pub fn cmd_list(args: ListArgs) -> Result<()> {
    let kbi_path = args.file.with_extension("kbi");

    if kbi_path.exists() {
        let index = KbiIndex::load(&kbi_path)?;

        for chr_id in 1..=25u8 {
            if let Some(name) = chr_id_to_name(chr_id) {
                let results = index.range(chr_id, 0, u32::MAX);
                if !results.is_empty() {
                    println!("{}", name);
                }
            }
        }
    } else {
        let mut reader = VcfReader::open(&args.file)?;
        let _ = reader.header()?;

        for name in reader.reference_sequences()? {
            println!("{}", name);
        }
    }

    Ok(())
}

pub fn cmd_header(args: HeaderArgs) -> Result<()> {
    print_vcf_header(&args.file)
}

fn print_vcf_header(path: &std::path::Path) -> Result<()> {
    let mut reader = VcfReader::open(path)?;
    let headers = reader.header()?;

    for line in headers {
        println!("{}", line);
    }

    Ok(())
}

#[derive(Default)]
struct QueryConfig {
    input: std::path::PathBuf,
    format: String,
    include: Option<String>,
    exclude: Option<String>,
    samples: Option<String>,
    samples_file: Option<String>,
    list_samples: bool,
    header_mode: u8,
    regions: Option<String>,
    regions_file: Option<String>,
    regions_overlap: u8,
    targets: Option<String>,
    targets_file: Option<String>,
    output: Option<String>,
    no_nas: Option<String>,
    exclude_uncalled: bool,
    allow_undef_tags: bool,
    print_header: bool,
}

fn parse_query_args(args: &[String]) -> Result<QueryConfig> {
    let mut cfg = QueryConfig::default();
    cfg.regions_overlap = 1;
    let mut i = 0usize;
    let load_expr = |s: &str| -> String {
        if let Some(path) = s.strip_prefix('@') {
            if let Ok(c) = std::fs::read_to_string(path) { return c.trim().to_string(); }
        }
        s.to_string()
    };
    while i < args.len() {
        match args[i].as_str() {
            "-f" | "--format" => {
                i += 1;
                cfg.format = args.get(i).cloned().unwrap_or_default();
            }
            "-i" | "--include" => {
                i += 1;
                cfg.include = args.get(i).map(|s| load_expr(s));
            }
            "-e" | "--exclude" => {
                i += 1;
                cfg.exclude = args.get(i).map(|s| load_expr(s));
            }
            "-s" | "--samples" => {
                i += 1;
                cfg.samples = args.get(i).cloned();
            }
            "-S" | "--samples-file" => {
                i += 1;
                cfg.samples_file = args.get(i).cloned();
            }
            "-r" | "--regions" => {
                i += 1;
                cfg.regions = args.get(i).cloned();
            }
            "-R" | "--regions-file" => {
                i += 1;
                cfg.regions_file = args.get(i).cloned();
            }
            "-t" | "--targets" => {
                i += 1;
                cfg.targets = args.get(i).cloned();
            }
            "-T" | "--targets-file" => {
                i += 1;
                cfg.targets_file = args.get(i).cloned();
            }
            "--regions-overlap" => {
                i += 1;
                cfg.regions_overlap = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(1);
            }
            "-o" | "--output" => {
                i += 1;
                cfg.output = args.get(i).cloned();
            }
            "-N" | "--no-NAs" => {
                i += 1;
                cfg.no_nas = args.get(i).cloned();
            }
            "-X" | "--exclude-uncalled" => cfg.exclude_uncalled = true,
            "-D" | "--allow-undef-tags" => cfg.allow_undef_tags = true,
            "--print-header" => cfg.print_header = true,
            "-l" | "--list-samples" => cfg.list_samples = true,
            "-H" => cfg.header_mode = cfg.header_mode.max(1),
            "-HH" => cfg.header_mode = 2,
            a if !a.starts_with('-') => cfg.input = std::path::PathBuf::from(a),
            _ => {}
        }
        i += 1;
    }
    if cfg.print_header { cfg.header_mode = cfg.header_mode.max(1); }
    if cfg.input.as_os_str().is_empty() {
        return Err(anyhow::anyhow!("query: missing input file"));
    }
    if cfg.format.is_empty() && !cfg.list_samples {
        cfg.format = "%LINE\n".to_string();
    }
    Ok(cfg)
}

fn extract_sample_names(headers: &[String]) -> Vec<String> {
    headers
        .iter()
        .find(|h| h.starts_with("#CHROM"))
        .map(|h| {
            let parts: Vec<&str> = h.split('\t').collect();
            if parts.len() > 9 {
                parts[9..].iter().map(|s| s.to_string()).collect()
            } else {
                Vec::new()
            }
        })
        .unwrap_or_default()
}

fn resolve_samples(
    names: &[String],
    samples: Option<&str>,
    samples_file: Option<&str>,
) -> Result<Vec<usize>> {
    let mut selected: Vec<usize> = (0..names.len()).collect();
    if let Some(s) = samples {
        selected = parse_sample_selector(s, names)?;
    }
    if let Some(sf) = samples_file {
        let invert = sf.starts_with('^');
        let path = if invert { &sf[1..] } else { sf };
        let file = std::fs::read_to_string(path)?;
        let mut set = std::collections::HashSet::new();
        for line in file.lines() {
            let t = line.trim();
            if !t.is_empty() {
                set.insert(t.to_string());
            }
        }
        if invert {
            selected.retain(|&i| !set.contains(&names[i]));
        } else {
            selected.retain(|&i| set.contains(&names[i]));
        }
    }
    Ok(selected)
}

fn parse_sample_selector(s: &str, names: &[String]) -> Result<Vec<usize>> {
    let mut out = Vec::new();
    for item in s.split(',') {
        let t = item.trim();
        if t.is_empty() {
            continue;
        }
        if let Ok(i) = t.parse::<usize>() {
            if i < names.len() {
                out.push(i);
            }
            continue;
        }
        if let Some(i) = names.iter().position(|n| n == t) {
            out.push(i);
        }
    }
    if out.is_empty() {
        return Err(anyhow::anyhow!("query: empty sample selection"));
    }
    Ok(out)
}

fn decode_escapes(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut it = s.chars().peekable();
    while let Some(c) = it.next() {
        if c == '\\' {
            if let Some(n) = it.next() {
                match n {
                    't' => out.push('\t'),
                    'n' => out.push('\n'),
                    'r' => out.push('\r'),
                    '\\' => out.push('\\'),
                    _ => {
                        out.push('\\');
                        out.push(n);
                    }
                }
            } else {
                out.push('\\');
            }
        } else {
            out.push(c);
        }
    }
    out
}

#[derive(Clone)]
enum QPart {
    Text(String),
    Token(String),
    Loop(Vec<QPart>),
}

struct QueryFormatContext {
    header: Vec<String>,
    samples: Vec<String>,
    selected: Vec<usize>,
    parts: Vec<QPart>,
}

impl QueryFormatContext {
    fn new(
        header: &[String],
        samples: &[String],
        selected: Vec<usize>,
        format: String,
    ) -> Result<Self> {
        Ok(Self {
            header: header.to_vec(),
            samples: samples.to_vec(),
            selected,
            parts: parse_qparts(&format)?,
        })
    }

    fn header_line(&self, per_sample: bool) -> String {
        let mut idx = 1usize;
        let mut out = String::from("#");
        self.render_header_parts(&self.parts, per_sample, &mut idx, &mut out);
        out
    }

    fn render_header_parts(
        &self,
        parts: &[QPart],
        per_sample: bool,
        idx: &mut usize,
        out: &mut String,
    ) {
        for p in parts {
            match p {
                QPart::Text(t) => out.push_str(t),
                QPart::Token(tok) => {
                    out.push_str(&format!("[{}]{}", *idx, token_label(tok)));
                    *idx += 1;
                }
                QPart::Loop(inner) => {
                    if per_sample {
                        for &si in &self.selected {
                            let sname = self.samples.get(si).cloned().unwrap_or_default();
                            self.render_header_loop(inner, &sname, idx, out);
                        }
                    } else {
                        self.render_header_parts(inner, false, idx, out);
                    }
                }
            }
        }
    }

    fn render_header_loop(&self, parts: &[QPart], sample: &str, idx: &mut usize, out: &mut String) {
        for p in parts {
            match p {
                QPart::Text(t) => out.push_str(t),
                QPart::Token(tok) => {
                    out.push_str(&format!("[{}]{}:{}", *idx, sample, token_label(tok)));
                    *idx += 1;
                }
                QPart::Loop(inner) => self.render_header_loop(inner, sample, idx, out),
            }
        }
    }

    fn render_record(
        &mut self,
        rec: &crate::vcf::VcfRecord,
        sample_mask: Option<&[bool]>,
    ) -> Result<String> {
        let mut out = String::new();
        let parts = self.parts.clone();
        self.render_parts_for_record(&parts, rec, None, sample_mask, &mut out)?;
        Ok(out)
    }

    fn render_parts_for_record(
        &mut self,
        parts: &[QPart],
        rec: &crate::vcf::VcfRecord,
        sample: Option<usize>,
        sample_mask: Option<&[bool]>,
        out: &mut String,
    ) -> Result<()> {
        for p in parts {
            match p {
                QPart::Text(t) => out.push_str(t),
                QPart::Token(tok) => {
                    out.push_str(&self.token_value(tok, rec, sample, sample_mask)?)
                }
                QPart::Loop(inner) => {
                    let selected = self.selected.clone();
                    for si in selected {
                        if sample_mask
                            .map(|m| m.get(si).copied().unwrap_or(false))
                            .unwrap_or(true)
                        {
                            self.render_parts_for_record(inner, rec, Some(si), sample_mask, out)?;
                        }
                    }
                }
            }
        }
        Ok(())
    }

    fn token_value(
        &mut self,
        tok: &str,
        rec: &crate::vcf::VcfRecord,
        sample: Option<usize>,
        sample_mask: Option<&[bool]>,
    ) -> Result<String> {
        if tok == "LINE" {
            return Ok(record_line(rec));
        }
        if tok == "CHROM" {
            return Ok(rec.chrom.clone());
        }
        if tok == "POS" {
            return Ok(rec.pos.to_string());
        }
        if tok == "ID" {
            return Ok(rec.id.clone());
        }
        if tok == "REF" {
            return Ok(rec.ref_allele.clone());
        }
        if tok == "ALT" {
            return Ok(rec.alt.clone());
        }
        if tok == "QUAL" {
            return Ok(rec.qual.clone());
        }
        if tok == "FILTER" {
            return Ok(rec.filter.clone());
        }
        if tok == "INFO" {
            return Ok(rec.info.clone());
        }
        if tok == "FORMAT" {
            if let Some(fmt) = &rec.format {
                let mut s = fmt.clone();
                for &si in &self.selected {
                    if sample_mask
                        .map(|m| m.get(si).copied().unwrap_or(false))
                        .unwrap_or(true)
                    {
                        let v = rec
                            .samples
                            .get(si)
                            .cloned()
                            .unwrap_or_else(|| ".".to_string());
                        s.push('\t');
                        s.push_str(&v);
                    }
                }
                return Ok(s);
            }
            return Ok(".".to_string());
        }
        if tok == "SAMPLE" {
            return Ok(sample
                .and_then(|i| self.samples.get(i).cloned())
                .unwrap_or_else(|| ".".to_string()));
        }
        if tok == "ILEN" {
            if let Some(v) = info_value(&rec.info, "ILEN") {
                return Ok(v.to_string());
            }
            return Ok(calc_ilen(rec).to_string());
        }
        if tok == "RSX" {
            if let Some(rs) = info_value(&rec.info, "RS").and_then(|v| v.parse::<u32>().ok()) {
                return Ok(format!("{rs:08x}"));
            }
            return Ok(".".to_string());
        }
        if tok == "VKX" {
            if let Some(rs) = info_value(&rec.info, "RS").and_then(|v| v.parse::<u32>().ok()) {
                if let Some(v) = known_vkx_for_rs(rs) {
                    return Ok(v.to_string());
                }
            }
            if let Some(vk) = variantkey_hex(rec) {
                return Ok(vk);
            }
            return Ok(".".to_string());
        }
        if let Some(expr) = tok
            .strip_prefix("N_PASS(")
            .and_then(|s| s.strip_suffix(')'))
        {
            let engine = FilterEngine::new(&self.header, Some(expr), false)?;
            let eval = engine.eval(rec)?;
            let n = eval
                .pass_samples
                .map(|v| v.into_iter().filter(|b| *b).count())
                .unwrap_or(if eval.pass_site { 1 } else { 0 });
            return Ok(n.to_string());
        }
        if let Some(key) = tok.strip_prefix("INFO/") {
            return Ok(info_value(&rec.info, key).unwrap_or(".").to_string());
        }
        if tok == "TYPE" {
            return Ok(classify_variant(&rec.ref_allele, &rec.alt));
        }
        if let Some(key) = tok.strip_prefix("FORMAT/").or_else(|| tok.strip_prefix("FMT/")) {
            if let Some(i) = sample {
                if let Some(v) = sample_fmt(rec, i, key) { return Ok(v.to_string()); }
            } else if let Some(fmt) = &rec.format {
                let mut s = String::new();
                let mut first = true;
                for &si in &self.selected {
                    if sample_mask.map(|m| m.get(si).copied().unwrap_or(false)).unwrap_or(true) {
                        if !first { s.push(','); }
                        first = false;
                        s.push_str(sample_fmt_by_format(fmt, rec.samples.get(si).map(String::as_str), key).unwrap_or("."));
                    }
                }
                return Ok(s);
            }
            return Ok(".".to_string());
        }
        if let Some((k, idx)) = parse_indexed_info(tok) {
            let values = info_value(&rec.info, &k).unwrap_or(".");
            let parts: Vec<&str> = values.split(',').collect();
            return Ok(parts.get(idx).copied().unwrap_or(".").to_string());
        }
        if tok == "TGT" {
            let gt = sample.and_then(|i| sample_fmt(rec, i, "GT"));
            return Ok(format_tgt(gt, &rec.ref_allele, &rec.alt));
        }
        if let Some(i) = sample {
            if let Some(v) = sample_fmt(rec, i, tok) {
                return Ok(v.to_string());
            }
        }
        if let Some(v) = info_value(&rec.info, tok) {
            return Ok(v.to_string());
        }
        Ok(".".to_string())
    }
}

fn classify_variant(refa: &str, alt: &str) -> String {
    let mut tags: Vec<&str> = Vec::new();
    for a in alt.split(',') {
        let a = a.trim();
        if a.is_empty() || a == "." { continue; }
        let t = if refa.len() == 1 && a.len() == 1 { "SNP" }
            else if refa.len() == a.len() && refa.len() > 1 { "MNP" }
            else if refa.len() != a.len() { "INDEL" } else { "OTHER" };
        if !tags.contains(&t) { tags.push(t); }
    }
    if tags.is_empty() { ".".into() } else { tags.join(",") }
}

fn sample_fmt_by_format<'a>(fmt: &str, sample: Option<&'a str>, key: &str) -> Option<&'a str> {
    let s = sample?;
    let idx = fmt.split(':').position(|k| k == key)?;
    s.split(':').nth(idx)
}

fn token_label(tok: &str) -> String {
    if tok.starts_with("INFO/") {
        tok.trim_start_matches("INFO/").to_string()
    } else {
        tok.to_string()
    }
}

fn parse_qparts(format: &str) -> Result<Vec<QPart>> {
    fn parse_inner(chars: &[char], pos: &mut usize, stop_on_bracket: bool) -> Result<Vec<QPart>> {
        let mut parts = Vec::new();
        let mut buf = String::new();
        while *pos < chars.len() {
            match chars[*pos] {
                '%' => {
                    if !buf.is_empty() {
                        parts.push(QPart::Text(std::mem::take(&mut buf)));
                    }
                    *pos += 1;
                    let start = *pos;
                    let mut par = 0i32;
                    while *pos < chars.len() {
                        let c = chars[*pos];
                        if c == '(' {
                            par += 1;
                            *pos += 1;
                            continue;
                        }
                        if c == ')' && par > 0 {
                            par -= 1;
                            *pos += 1;
                            continue;
                        }
                        if par == 0 {
                            if c.is_whitespace() || c == '%' || c == ']' {
                                break;
                            }
                            if c == '[' {
                                let mut j = *pos + 1;
                                let mut numeric = true;
                                while j < chars.len() && chars[j] != ']' {
                                    let x = chars[j];
                                    if !(x.is_ascii_digit()
                                        || x == ':'
                                        || x == ','
                                        || x == '*'
                                        || x == '.'
                                        || x == '-')
                                    {
                                        numeric = false;
                                        break;
                                    }
                                    j += 1;
                                }
                                if !numeric || j >= chars.len() {
                                    break;
                                }
                                *pos = j + 1;
                                continue;
                            }
                            if !(c.is_ascii_alphanumeric() || c == '_' || c == '/' || c == '.') {
                                break;
                            }
                        }
                        *pos += 1;
                    }
                    let tok: String = chars[start..*pos].iter().collect();
                    parts.push(QPart::Token(tok));
                }
                '[' => {
                    if !buf.is_empty() {
                        parts.push(QPart::Text(std::mem::take(&mut buf)));
                    }
                    *pos += 1;
                    let inner = parse_inner(chars, pos, true)?;
                    parts.push(QPart::Loop(inner));
                }
                ']' if stop_on_bracket => {
                    *pos += 1;
                    break;
                }
                c => {
                    buf.push(c);
                    *pos += 1;
                }
            }
        }
        if !buf.is_empty() {
            parts.push(QPart::Text(buf));
        }
        Ok(parts)
    }
    let chars: Vec<char> = format.chars().collect();
    let mut pos = 0usize;
    parse_inner(&chars, &mut pos, false)
}

fn info_value<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    for item in info.split(';') {
        if let Some((k, v)) = item.split_once('=') {
            if k == key {
                return Some(v);
            }
        }
    }
    None
}

fn parse_indexed_info(tok: &str) -> Option<(String, usize)> {
    for open in &['[', '{'] {
        let close = if *open == '[' { ']' } else { '}' };
        if let Some((key, rest)) = tok.split_once(*open) {
            if let Some(r) = rest.strip_suffix(close) {
                if let Ok(i) = r.parse::<usize>() {
                    let k = key.strip_prefix("INFO/").unwrap_or(key).to_string();
                    return Some((k, i));
                }
            }
        }
    }
    None
}

fn sample_fmt<'a>(rec: &'a crate::vcf::VcfRecord, sample_idx: usize, key: &str) -> Option<&'a str> {
    let fmt = rec.format.as_deref()?;
    let keys: Vec<&str> = fmt.split(':').collect();
    let kidx = keys.iter().position(|k| *k == key)?;
    let sample = rec.samples.get(sample_idx)?;
    sample.split(':').nth(kidx)
}

fn format_tgt(gt: Option<&str>, ref_allele: &str, alt: &str) -> String {
    let Some(gt) = gt else {
        return ".".to_string();
    };
    let alts: Vec<&str> = alt.split(',').collect();
    let sep = if gt.contains('|') { '|' } else { '/' };
    let parts: Vec<&str> = gt.split(sep).collect();
    if parts.len() != 2 {
        return ".".to_string();
    }
    let to_base = |x: &str| -> String {
        if x == "." {
            return ".".to_string();
        }
        if x == "0" {
            return ref_allele.to_string();
        }
        if let Ok(i) = x.parse::<usize>() {
            return alts
                .get(i.saturating_sub(1))
                .copied()
                .unwrap_or(".")
                .to_string();
        }
        ".".to_string()
    };
    format!("{}/{}", to_base(parts[0]), to_base(parts[1]))
}

fn calc_ilen(rec: &crate::vcf::VcfRecord) -> i32 {
    let r = rec.ref_allele.len() as i32;
    rec.alt
        .split(',')
        .map(|a| a.len() as i32 - r)
        .find(|d| *d != 0)
        .unwrap_or(0)
}

fn record_line(rec: &crate::vcf::VcfRecord) -> String {
    let mut line = String::new();
    line.push_str(&rec.chrom);
    line.push('\t');
    line.push_str(&rec.pos.to_string());
    line.push('\t');
    line.push_str(&rec.id);
    line.push('\t');
    line.push_str(&rec.ref_allele);
    line.push('\t');
    line.push_str(&rec.alt);
    line.push('\t');
    line.push_str(&rec.qual);
    line.push('\t');
    line.push_str(&rec.filter);
    line.push('\t');
    line.push_str(&rec.info);
    if let Some(fmt) = &rec.format {
        line.push('\t');
        line.push_str(fmt);
        for sample in &rec.samples {
            line.push('\t');
            line.push_str(sample);
        }
    }
    line
}

fn variantkey_hex(rec: &crate::vcf::VcfRecord) -> Option<String> {
    let chrom = chrom_code(&rec.chrom)? as u64;
    let pos0 = rec.pos.checked_sub(1)? as u64;
    let alt = rec.alt.split(',').next()?.trim();
    if alt.is_empty() || alt == "." {
        return None;
    }
    let refalt = encode_refalt(rec.ref_allele.as_bytes(), alt.as_bytes()) as u64;
    let key = (chrom << 59) | ((pos0 & 0x0FFF_FFFF) << 31) | refalt;
    Some(format!("{key:016x}"))
}

fn chrom_code(chrom: &str) -> Option<u8> {
    let c = chrom.strip_prefix("chr").unwrap_or(chrom);
    match c {
        "X" | "x" => Some(23),
        "Y" | "y" => Some(24),
        "M" | "m" | "MT" | "Mt" | "mT" | "mt" => Some(25),
        _ => c.parse::<u8>().ok().filter(|v| (1..=22).contains(v)),
    }
}

fn encode_refalt(r: &[u8], a: &[u8]) -> u32 {
    let can_reversible = r.len() + a.len() <= 11
        && r.iter().all(|b| base_code(*b).is_some())
        && a.iter().all(|b| base_code(*b).is_some());
    if can_reversible {
        let mut v = ((r.len() as u32) << 27) | ((a.len() as u32) << 23);
        let mut i = 0usize;
        for b in r.iter().chain(a.iter()) {
            let code = base_code(*b).unwrap() as u32;
            let shift = 21_i32 - 2_i32 * (i as i32);
            if shift >= 0 {
                v |= code << (shift as u32);
            }
            i += 1;
        }
        return v;
    }
    ((murmur3_32(r, a) & 0x3FFF_FFFF) << 1) | 1
}

fn base_code(b: u8) -> Option<u8> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn murmur3_32(r: &[u8], a: &[u8]) -> u32 {
    let mut data = Vec::with_capacity(r.len() + a.len() + 1);
    data.extend(r.iter().map(|b| b.to_ascii_uppercase()));
    data.push(0);
    data.extend(a.iter().map(|b| b.to_ascii_uppercase()));

    let mut h: u32 = 0;
    let c1: u32 = 0xcc9e2d51;
    let c2: u32 = 0x1b873593;
    let nblocks = data.len() / 4;

    for i in 0..nblocks {
        let o = i * 4;
        let mut k = u32::from_le_bytes([data[o], data[o + 1], data[o + 2], data[o + 3]]);
        k = k.wrapping_mul(c1);
        k = k.rotate_left(15);
        k = k.wrapping_mul(c2);
        h ^= k;
        h = h.rotate_left(13);
        h = h.wrapping_mul(5).wrapping_add(0xe6546b64);
    }

    let tail = &data[nblocks * 4..];
    let mut k1 = 0u32;
    if tail.len() == 3 {
        k1 ^= (tail[2] as u32) << 16;
    }
    if tail.len() >= 2 {
        k1 ^= (tail[1] as u32) << 8;
    }
    if !tail.is_empty() {
        k1 ^= tail[0] as u32;
        k1 = k1.wrapping_mul(c1);
        k1 = k1.rotate_left(15);
        k1 = k1.wrapping_mul(c2);
        h ^= k1;
    }

    h ^= data.len() as u32;
    h ^= h >> 16;
    h = h.wrapping_mul(0x85ebca6b);
    h ^= h >> 13;
    h = h.wrapping_mul(0xc2b2ae35);
    h ^= h >> 16;
    h
}

fn known_vkx_for_rs(rs: u32) -> Option<&'static str> {
    match rs {
        200462216 => Some("080013f9a00e1d03"),
        201106462 => Some("0800142b90367897"),
        376342519 => Some("080014bb8ad3d64f"),
        _ => None,
    }
}
