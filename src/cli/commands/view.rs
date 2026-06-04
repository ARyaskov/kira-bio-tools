use crate::annotate::postproc::{
    HeaderOptions, LineAction, OutputKind, PostProcessor, Predicate, RegionFilter,
    apply_to_header, parse_output_type, parse_samples_cli, process_record_line,
    read_samples_file, resolve_samples_keep, version_header_line,
};
use crate::cli::args::ViewArgs;
use crate::filter::FilterEngine;
use crate::vcf::UnifiedVcfReader;
use anyhow::{Context, Result, bail};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn cmd_view(args: ViewArgs) -> Result<()> {
    let region_filter = if let Some(s) = &args.regions {
        Some(RegionFilter::from_cli(s)?)
    } else if let Some(p) = &args.regions_file {
        Some(RegionFilter::from_file(p)?)
    } else { None };

    let target_filter = if let Some(s) = &args.targets {
        Some(RegionFilter::from_cli(s.trim_start_matches('^'))?)
    } else if let Some(p) = &args.targets_file {
        Some(RegionFilter::from_file(p)?)
    } else { None };
    let targets_inverse = args.targets.as_deref().map_or(false, |s| s.starts_with('^'));

    let mut reader = UnifiedVcfReader::open(&args.input).context("open input")?;
    let headers = reader.header()?;
    let input_samples = extract_samples(&headers);

    let samples_keep = resolve_samples(&args, &input_samples)?;
    let mut pp = PostProcessor::default();
    pp.no_version = args.no_version;
    pp.samples_keep = samples_keep.clone();

    if let Some(expr) = &args.include {
        let engine = FilterEngine::new(&headers, Some(expr.as_str()), false)
            .context("-i/--include expression")?;
        pp.include = Some(Predicate { engine });
    }
    if let Some(expr) = &args.exclude {
        let engine = FilterEngine::new(&headers, Some(expr.as_str()), false)
            .context("-e/--exclude expression")?;
        pp.exclude = Some(Predicate { engine });
    }

    let type_filter = TypeFilter::parse(args.types.as_deref(), args.exclude_types.as_deref())?;
    let ac_filter = AcAfFilter::build(&args)?;
    let gt_filter = GenotypeFilter::parse(args.genotype.as_deref())?;
    let apply_filters = parse_apply_filters(args.apply_filters.as_deref());

    let kind = args.output_type.as_deref().map(parse_output_type).transpose()?
        .unwrap_or(OutputKind::Vcf);
    let mut sink = OutputSink::open(args.output.as_deref(), kind, args.compression_level, &headers)?;

    let version = version_header_line();
    let opts = HeaderOptions {
        no_version: args.no_version,
        extra_header_lines: &[],
        remove: None,
        rename_chrs: None,
        rename_annots: None,
        mark_sites: None,
        set_id: false,
        samples_keep: samples_keep.as_deref(),
        version_line: Some(&version),
    };
    let mut out_headers = apply_to_header(headers, &opts);
    if args.drop_genotypes {
        out_headers = drop_genotypes_in_header(out_headers);
    }
    if !args.no_header {
        for h in &out_headers {
            sink.write_line(h)?;
        }
    }
    if args.header_only { return sink.finish(); }

    while let Some(line) = reader.read_line()? {
        if line.is_empty() || line.as_bytes()[0] == b'#' { continue; }

        if let Some(rf) = &region_filter { if !rf.line_passes_mode(&line, args.regions_overlap) { continue; } }
        if let Some(tf) = &target_filter {
            let pass = tf.line_passes(&line);
            if pass == targets_inverse { continue; }
        }

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 8 { continue; }

        if let Some(af) = &apply_filters {
            if !filter_passes_apply(cols[6], af) { continue; }
        }

        if args.known && cols[2] == "." { continue; }
        if args.novel && cols[2] != "." { continue; }

        if let Some(t) = &type_filter {
            if !t.passes(cols[3], cols[4]) { continue; }
        }

        let n_alts = cols[4].split(',').filter(|s| !s.is_empty() && *s != "<*>" && *s != "<NON_REF>").count() as u32;
        if let Some(min) = args.min_alleles { if 1 + n_alts < min { continue; } }
        if let Some(max) = args.max_alleles { if 1 + n_alts > max { continue; } }

        let samples_slice: &[&str] = if cols.len() > 9 { &cols[9..] } else { &[] };
        let format_str: &str = if cols.len() > 8 { cols[8] } else { "" };

        if let Some(f) = &gt_filter { if !f.passes(format_str, samples_slice, samples_keep.as_deref()) { continue; } }
        if args.phased && !any_phased(format_str, samples_slice) { continue; }
        if args.exclude_phased && any_phased(format_str, samples_slice) { continue; }
        if args.uncalled && !any_uncalled(format_str, samples_slice) { continue; }
        if args.exclude_uncalled && any_uncalled(format_str, samples_slice) { continue; }

        if let Some(f) = &ac_filter {
            if !f.passes(format_str, samples_slice, cols[7], n_alts as usize) { continue; }
        }

        if args.private && !is_private(format_str, samples_slice) { continue; }
        if args.exclude_private && is_private(format_str, samples_slice) { continue; }

        let mut out_line = line.clone();
        if args.drop_genotypes && cols.len() > 8 {
            out_line = cols[..8].join("\t");
        } else {
            match process_record_line(&line, &pp, false) {
                LineAction::Replace(s) => out_line = s,
                LineAction::Drop => continue,
                LineAction::Keep => {}
            }
        }

        if args.trim_alt_alleles {
            out_line = trim_unused_alts(&out_line);
        }

        sink.write_line(&out_line)?;
    }
    sink.finish()?;
    if let (Some(kind), Some(out)) = (args.write_index.as_deref(), args.output.as_deref()) {
        let _ = match crate::csi::build_csi_index(out, &std::path::PathBuf::from(format!("{}.{}", out.display(), if kind == "tbi" { "tbi" } else { "csi" }))) {
            Ok(_) => (), Err(e) => eprintln!("[view] -W: {}", e),
        };
    }
    Ok(())
}

struct OutputSink {
    inner: Box<dyn Write>,
    bgzf: bool,
    bcf: Option<crate::bcf::BcfWriter>,
}

impl OutputSink {
    fn open(path: Option<&Path>, kind: OutputKind, lvl: i32, headers: &[String]) -> Result<Self> {
        let level = match kind { OutputKind::VcfGz(l) | OutputKind::Bcf(l) => l, _ => if lvl >= 0 { lvl as u32 } else { 6 } };
        match (path, kind) {
            (None, OutputKind::Vcf) => Ok(Self { inner: Box::new(BufWriter::with_capacity(1 << 20, std::io::stdout())), bgzf: false, bcf: None }),
            (None, OutputKind::VcfGz(_)) => bail!("-O z without -o FILE not supported (BGZF needs seekable file)"),
            (None, OutputKind::Bcf(_)) => bail!("-O u|b without -o FILE not supported"),
            (Some(p), OutputKind::Vcf) => Ok(Self { inner: Box::new(BufWriter::with_capacity(1 << 20, File::create(p)?)), bgzf: false, bcf: None }),
            (Some(p), OutputKind::VcfGz(_)) => {
                let w = crate::bgzf::BgzfWriter::with_compression(p, flate2::Compression::new(level))?;
                Ok(Self { inner: Box::new(w), bgzf: true, bcf: None })
            }
            (Some(p), OutputKind::Bcf(_)) => {
                let compressed = matches!(kind, OutputKind::Bcf(l) if l > 0);
                let w = crate::bcf::BcfWriter::create(p, compressed, level, headers)?;
                Ok(Self { inner: Box::new(std::io::sink()), bgzf: false, bcf: Some(w) })
            }
        }
    }
    fn write_line(&mut self, line: &str) -> Result<()> {
        if let Some(bcf) = self.bcf.as_mut() {
            if !line.starts_with('#') { bcf.write_vcf_line(line)?; }
            return Ok(());
        }
        self.inner.write_all(line.as_bytes())?;
        self.inner.write_all(b"\n")?;
        Ok(())
    }
    fn finish(mut self) -> Result<()> {
        if let Some(bcf) = self.bcf.take() { bcf.finish()?; return Ok(()); }
        self.inner.flush()?;
        if self.bgzf {}
        Ok(())
    }
}

fn extract_samples(headers: &[String]) -> Vec<String> {
    for h in headers {
        if h.starts_with("#CHROM") {
            let cols: Vec<&str> = h.split('\t').collect();
            if cols.len() > 9 { return cols[9..].iter().map(|s| s.to_string()).collect(); }
        }
    }
    Vec::new()
}

fn resolve_samples(args: &ViewArgs, input: &[String]) -> Result<Option<Vec<usize>>> {
    if let Some(s) = &args.samples {
        let (names, inv) = parse_samples_cli(s);
        if !args.force_samples {
            for n in &names {
                if !input.iter().any(|s| s == n) { bail!("sample {n:?} not found (use --force-samples to ignore)"); }
            }
        }
        return Ok(Some(resolve_samples_keep(input, &names, inv)));
    }
    if let Some(p) = &args.samples_file {
        let (names, inv) = read_samples_file(p)?;
        if !args.force_samples {
            for n in &names {
                if !input.iter().any(|s| s == n) { bail!("sample {n:?} not found (use --force-samples to ignore)"); }
            }
        }
        return Ok(Some(resolve_samples_keep(input, &names, inv)));
    }
    Ok(None)
}

fn drop_genotypes_in_header(headers: Vec<String>) -> Vec<String> {
    headers.into_iter().filter_map(|h| {
        if h.starts_with("##FORMAT=") { return None; }
        if h.starts_with("#CHROM") {
            let cols: Vec<&str> = h.split('\t').collect();
            if cols.len() > 8 { return Some(cols[..8].join("\t")); }
        }
        Some(h)
    }).collect()
}

struct TypeFilter { incl: Option<u8>, excl: Option<u8> }
const T_SNP: u8 = 1; const T_INDEL: u8 = 2; const T_MNP: u8 = 4; const T_OTHER: u8 = 8;

impl TypeFilter {
    fn parse(inc: Option<&str>, exc: Option<&str>) -> Result<Option<Self>> {
        if inc.is_none() && exc.is_none() { return Ok(None); }
        let p = |s: &str| -> Result<u8> {
            let mut m = 0u8;
            for t in s.split(',') {
                m |= match t.trim() {
                    "snps" | "snp" => T_SNP,
                    "indels" | "indel" => T_INDEL,
                    "mnps" | "mnp" => T_MNP,
                    "other" => T_OTHER,
                    "ref" => 0,
                    o => bail!("unknown variant type {o:?}"),
                };
            }
            Ok(m)
        };
        Ok(Some(Self { incl: inc.map(p).transpose()?, excl: exc.map(p).transpose()? }))
    }
    fn classify(refa: &str, alt: &str) -> u8 {
        let mut m = 0u8;
        for a in alt.split(',') {
            let a = a.trim();
            if a.is_empty() || a == "." { continue; }
            if refa.len() == 1 && a.len() == 1 { m |= T_SNP; }
            else if refa.len() == a.len() && refa.len() > 1 { m |= T_MNP; }
            else if refa.len() != a.len() { m |= T_INDEL; }
            else { m |= T_OTHER; }
        }
        m
    }
    fn passes(&self, refa: &str, alt: &str) -> bool {
        let m = Self::classify(refa, alt);
        if let Some(i) = self.incl { if m & i == 0 { return false; } }
        if let Some(e) = self.excl { if m & e != 0 { return false; } }
        true
    }
}

struct GenotypeFilter { mode: GtMode }
enum GtMode { Hom, Het, Miss, NoMiss, NeedRef }

impl GenotypeFilter {
    fn parse(s: Option<&str>) -> Result<Option<Self>> {
        let Some(s) = s else { return Ok(None); };
        let mode = match s {
            "hom" => GtMode::Hom, "het" => GtMode::Het,
            "miss" => GtMode::Miss, "^miss" => GtMode::NoMiss,
            "^ref" => GtMode::NeedRef,
            o => bail!("--genotype: unknown mode {o:?}"),
        };
        Ok(Some(Self { mode }))
    }
    fn passes(&self, format: &str, samples: &[&str], keep: Option<&[usize]>) -> bool {
        let Some(gt_idx) = format.split(':').position(|k| k == "GT") else { return matches!(self.mode, GtMode::Miss); };
        let iter: Box<dyn Iterator<Item = &&str>> = match keep {
            Some(k) => Box::new(k.iter().filter_map(|&i| samples.get(i))),
            None => Box::new(samples.iter()),
        };
        for s in iter {
            let gt = s.split(':').nth(gt_idx).unwrap_or(".");
            let cls = gt_class(gt);
            let pass = matches!((&self.mode, cls),
                (GtMode::Hom, GtClass::Hom) |
                (GtMode::Het, GtClass::Het) |
                (GtMode::Miss, GtClass::Miss) |
                (GtMode::NoMiss, GtClass::Hom | GtClass::Het) |
                (GtMode::NeedRef, GtClass::Het | GtClass::Hom));
            if pass { return true; }
        }
        false
    }
}

enum GtClass { Hom, Het, Miss }
fn gt_class(gt: &str) -> GtClass {
    let alleles: Vec<&str> = gt.split(|c| c == '/' || c == '|').collect();
    if alleles.is_empty() || alleles.iter().any(|a| *a == "." || a.is_empty()) { return GtClass::Miss; }
    let first = alleles[0];
    if alleles.iter().all(|a| a == &first) { GtClass::Hom } else { GtClass::Het }
}

fn any_phased(format: &str, samples: &[&str]) -> bool {
    let Some(gt_idx) = format.split(':').position(|k| k == "GT") else { return false; };
    samples.iter().any(|s| s.split(':').nth(gt_idx).unwrap_or("").contains('|'))
}

fn any_uncalled(format: &str, samples: &[&str]) -> bool {
    let Some(gt_idx) = format.split(':').position(|k| k == "GT") else { return false; };
    samples.iter().any(|s| matches!(gt_class(s.split(':').nth(gt_idx).unwrap_or(".")), GtClass::Miss))
}

fn is_private(format: &str, samples: &[&str]) -> bool {
    let Some(gt_idx) = format.split(':').position(|k| k == "GT") else { return false; };
    let mut non_ref = 0usize;
    for s in samples {
        let gt = s.split(':').nth(gt_idx).unwrap_or(".");
        for a in gt.split(|c| c == '/' || c == '|') {
            if let Ok(n) = a.parse::<u32>() { if n > 0 { non_ref += 1; break; } }
        }
    }
    non_ref == 1
}

struct AcAfFilter { min_ac: Option<u32>, max_ac: Option<u32>, min_af: Option<f64>, max_af: Option<f64> }

impl AcAfFilter {
    fn build(a: &ViewArgs) -> Result<Option<Self>> {
        let p_u32 = |s: &str| -> Result<u32> { Ok(s.split(':').next().unwrap_or(s).parse()?) };
        let p_f64 = |s: &str| -> Result<f64> { Ok(s.split(':').next().unwrap_or(s).parse()?) };
        let min_ac = a.min_ac.as_deref().map(p_u32).transpose()?;
        let max_ac = a.max_ac.as_deref().map(p_u32).transpose()?;
        let min_af = a.min_af.as_deref().map(p_f64).transpose()?;
        let max_af = a.max_af.as_deref().map(p_f64).transpose()?;
        if min_ac.is_none() && max_ac.is_none() && min_af.is_none() && max_af.is_none() { return Ok(None); }
        Ok(Some(Self { min_ac, max_ac, min_af, max_af }))
    }
    fn passes(&self, format: &str, samples: &[&str], info: &str, n_alts: usize) -> bool {
        let (ac, an) = if !samples.is_empty() && !format.is_empty() {
            compute_ac_an_from_gts(format, samples, n_alts)
        } else {
            info_ac_an(info)
        };
        let af: f64 = if an > 0 { ac as f64 / an as f64 } else { 0.0 };
        if let Some(v) = self.min_ac { if ac < v { return false; } }
        if let Some(v) = self.max_ac { if ac > v { return false; } }
        if let Some(v) = self.min_af { if af < v { return false; } }
        if let Some(v) = self.max_af { if af > v { return false; } }
        true
    }
}

fn compute_ac_an_from_gts(format: &str, samples: &[&str], _n_alts: usize) -> (u32, u32) {
    let Some(gt_idx) = format.split(':').position(|k| k == "GT") else { return (0, 0); };
    let (mut ac, mut an) = (0u32, 0u32);
    for s in samples {
        let gt = s.split(':').nth(gt_idx).unwrap_or(".");
        for a in gt.split(|c| c == '/' || c == '|') {
            if a == "." || a.is_empty() { continue; }
            an += 1;
            if let Ok(n) = a.parse::<u32>() { if n > 0 { ac += 1; } }
        }
    }
    (ac, an)
}

fn info_ac_an(info: &str) -> (u32, u32) {
    let mut ac = 0u32; let mut an = 0u32;
    for kv in info.split(';') {
        if let Some(v) = kv.strip_prefix("AC=") {
            ac = v.split(',').filter_map(|s| s.parse::<u32>().ok()).sum();
        } else if let Some(v) = kv.strip_prefix("AN=") {
            an = v.parse().unwrap_or(0);
        }
    }
    (ac, an)
}

fn parse_apply_filters(s: Option<&str>) -> Option<Vec<String>> {
    s.map(|x| x.split(',').map(|t| t.trim().to_string()).collect())
}

fn filter_passes_apply(filter_col: &str, allowed: &[String]) -> bool {
    if filter_col == "." || filter_col.is_empty() { return false; }
    filter_col.split(';').any(|t| allowed.iter().any(|a| a == t || (a == "PASS" && t == "PASS")))
}

fn trim_unused_alts(line: &str) -> String {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 8 { return line.to_string(); }
    let alts: Vec<&str> = cols[4].split(',').collect();
    if alts.len() <= 1 { return line.to_string(); }
    let format_str = if cols.len() > 8 { cols[8] } else { "" };
    let samples = if cols.len() > 9 { &cols[9..] } else { &[][..] };
    let Some(gt_idx) = format_str.split(':').position(|k| k == "GT") else { return line.to_string(); };
    let mut used = vec![false; alts.len() + 1];
    for s in samples {
        let gt = s.split(':').nth(gt_idx).unwrap_or(".");
        for a in gt.split(|c| c == '/' || c == '|') {
            if let Ok(n) = a.parse::<usize>() { if n < used.len() { used[n] = true; } }
        }
    }
    let new_alts: Vec<&str> = alts.iter().enumerate().filter(|(i, _)| used[i + 1]).map(|(_, a)| *a).collect();
    if new_alts.len() == alts.len() { return line.to_string(); }
    let mut cols = cols;
    let joined = new_alts.join(",");
    cols[4] = &joined;
    cols.join("\t")
}
