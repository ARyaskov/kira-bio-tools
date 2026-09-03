use crate::annotate::postproc::{
    HeaderOptions, LineAction, OutputKind, PostProcessor, Predicate, RegionFilter,
    apply_to_header, parse_output_type, parse_samples_cli, read_samples_file,
    resolve_samples_keep, version_header_line,
};
use crate::cli::args::ViewArgs;
use crate::csi::{IndexKind, build_index, find_index_for};
use crate::filter::FilterEngine;
use crate::vcf::alleles::{gt_alleles, remap_info, remap_samples};
use crate::vcf::header::{HeaderInfo, extract_samples};
use crate::vcf::variant_type::{VT_REF, allele_types, parse_type_mask, record_type};
use crate::vcf::{UnifiedVcfReader, VcfSink};
use anyhow::{Context, Result, bail};
use std::path::Path;

pub fn cmd_view(args: ViewArgs) -> Result<()> {
    let region_filter = if let Some(s) = &args.regions {
        Some(RegionFilter::from_cli(s)?)
    } else if let Some(p) = &args.regions_file {
        Some(RegionFilter::from_file(p)?)
    } else {
        None
    };

    let target_filter = if let Some(s) = &args.targets {
        Some(RegionFilter::from_cli(s.trim_start_matches('^'))?)
    } else if let Some(p) = &args.targets_file {
        Some(RegionFilter::from_file(p)?)
    } else {
        None
    };
    let targets_inverse = args.targets.as_deref().is_some_and(|s| s.starts_with('^'));

    let mut reader = UnifiedVcfReader::open(&args.input).context("open input")?;
    let headers = reader.header()?;
    let hdr_info = HeaderInfo::parse(&headers);
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

    let kind = args
        .output_type
        .as_deref()
        .map(parse_output_type)
        .transpose()?
        .unwrap_or(OutputKind::Vcf);
    let level = if args.compression_level >= 0 { Some(args.compression_level as u32) } else { None };

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
    let mut sink = VcfSink::open_with_level(args.output.as_deref(), kind, level, &out_headers)?;
    if !args.no_header {
        sink.write_header(&out_headers)?;
    }
    if args.header_only {
        return sink.finish();
    }

    let update_ac_an = samples_keep.is_some() && !args.no_update && !args.drop_genotypes;
    let ctx = RecordCtx {
        args: &args,
        pp: &pp,
        hdr: &hdr_info,
        target_filter: target_filter.as_ref(),
        targets_inverse,
        type_filter: type_filter.as_ref(),
        ac_filter: ac_filter.as_ref(),
        gt_filter: gt_filter.as_ref(),
        apply_filters: apply_filters.as_deref(),
        samples_keep: samples_keep.as_deref(),
        update_ac_an,
    };

    // `-r` uses the index when there is one; otherwise every record is scanned.
    let indexed = region_filter
        .as_ref()
        .filter(|_| args.input != Path::new("-") && find_index_for(&args.input).is_some());
    if let Some(rf) = indexed {
        rf.stream_with_index(&args.input, args.regions_overlap, |line| {
            if let Some(out_line) = process_view_record(line, &ctx) {
                sink.write_line(&out_line)?;
            }
            Ok(())
        })?;
    } else {
        while let Some(line) = reader.read_line()? {
            if line.is_empty() || line.as_bytes()[0] == b'#' {
                continue;
            }
            if let Some(rf) = &region_filter {
                if !rf.line_passes_mode(&line, args.regions_overlap) {
                    continue;
                }
            }
            if let Some(out_line) = process_view_record(&line, &ctx) {
                sink.write_line(&out_line)?;
            }
        }
    }
    sink.finish()?;

    if let (Some(kind_s), Some(out)) = (args.write_index.as_deref(), args.output.as_deref()) {
        if matches!(kind, OutputKind::VcfGz(_) | OutputKind::Bcf(_)) && out != Path::new("-") {
            let (ikind, ext) = if kind_s == "tbi" { (IndexKind::Tbi, "tbi") } else { (IndexKind::Csi, "csi") };
            let idx = std::path::PathBuf::from(format!("{}.{}", out.display(), ext));
            build_index(out, &idx, ikind, None).with_context(|| format!("-W: write index {}", idx.display()))?;
        } else {
            eprintln!("[view] -W: an index needs BGZF or BCF output to a file; skipped");
        }
    }
    Ok(())
}

struct RecordCtx<'a> {
    args: &'a ViewArgs,
    pp: &'a PostProcessor,
    hdr: &'a HeaderInfo,
    target_filter: Option<&'a RegionFilter>,
    targets_inverse: bool,
    type_filter: Option<&'a TypeFilter>,
    ac_filter: Option<&'a AcAfFilter>,
    gt_filter: Option<&'a GenotypeFilter>,
    apply_filters: Option<&'a [String]>,
    samples_keep: Option<&'a [usize]>,
    update_ac_an: bool,
}

/// All per-record filtering and rewriting; `None` drops the record.
fn process_view_record(line: &str, ctx: &RecordCtx<'_>) -> Option<String> {
    let args = ctx.args;
    if let Some(tf) = ctx.target_filter {
        let pass = tf.line_passes_mode(line, args.targets_overlap);
        if pass == ctx.targets_inverse {
            return None;
        }
    }

    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 8 {
        return None;
    }

    if let Some(af) = ctx.apply_filters {
        if !filter_passes_apply(cols[6], af) {
            return None;
        }
    }
    if args.known && cols[2] == "." {
        return None;
    }
    if args.novel && cols[2] != "." {
        return None;
    }
    if let Some(t) = ctx.type_filter {
        if !t.passes(cols[3], cols[4]) {
            return None;
        }
    }

    let n_alleles = 1 + if cols[4] == "." || cols[4].is_empty() { 0 } else { cols[4].split(',').count() as u32 };
    if let Some(min) = args.min_alleles {
        if n_alleles < min {
            return None;
        }
    }
    if let Some(max) = args.max_alleles {
        if n_alleles > max {
            return None;
        }
    }

    let samples_slice: &[&str] = if cols.len() > 9 { &cols[9..] } else { &[] };
    let format_str: &str = if cols.len() > 8 { cols[8] } else { "" };

    if let Some(f) = ctx.gt_filter {
        if !f.passes(format_str, samples_slice, ctx.samples_keep) {
            return None;
        }
    }
    if args.phased && !any_phased(format_str, samples_slice) {
        return None;
    }
    if args.exclude_phased && any_phased(format_str, samples_slice) {
        return None;
    }
    if args.uncalled && !any_uncalled(format_str, samples_slice) {
        return None;
    }
    if args.exclude_uncalled && any_uncalled(format_str, samples_slice) {
        return None;
    }
    if let Some(f) = ctx.ac_filter {
        if !f.passes(format_str, samples_slice, ctx.samples_keep, cols[7], n_alleles as usize) {
            return None;
        }
    }
    if args.private && !is_private(format_str, samples_slice) {
        return None;
    }
    if args.exclude_private && is_private(format_str, samples_slice) {
        return None;
    }

    let mut out_line = line.to_string();
    if args.drop_genotypes && cols.len() > 8 {
        out_line = cols[..8].join("\t");
    } else {
        match crate::annotate::postproc::process_record_cols(line, &cols, ctx.pp, false) {
            LineAction::Replace(s) => out_line = s,
            LineAction::Drop => return None,
            LineAction::Keep => {}
        }
        if ctx.update_ac_an {
            out_line = update_ac_an(&out_line, ctx.hdr);
        }
    }

    if args.trim_alt_alleles || args.trim_unseen_allele {
        out_line = trim_alt_alleles(&out_line, ctx.hdr, args.trim_alt_alleles, args.trim_unseen_allele);
    }
    Some(out_line)
}

fn resolve_samples(args: &ViewArgs, input: &[String]) -> Result<Option<Vec<usize>>> {
    if let Some(s) = &args.samples {
        let (names, inv) = parse_samples_cli(s);
        if !args.force_samples {
            for n in &names {
                if !input.iter().any(|s| s == n) {
                    bail!("sample {n:?} not found (use --force-samples to ignore)");
                }
            }
        }
        return Ok(Some(resolve_samples_keep(input, &names, inv)));
    }
    if let Some(p) = &args.samples_file {
        let (names, inv) = read_samples_file(p)?;
        if !args.force_samples {
            for n in &names {
                if !input.iter().any(|s| s == n) {
                    bail!("sample {n:?} not found (use --force-samples to ignore)");
                }
            }
        }
        return Ok(Some(resolve_samples_keep(input, &names, inv)));
    }
    Ok(None)
}

fn drop_genotypes_in_header(headers: Vec<String>) -> Vec<String> {
    headers
        .into_iter()
        .filter_map(|h| {
            if h.starts_with("##FORMAT=") {
                return None;
            }
            if h.starts_with("#CHROM") {
                let cols: Vec<&str> = h.split('\t').collect();
                if cols.len() > 8 {
                    return Some(cols[..8].join("\t"));
                }
            }
            Some(h)
        })
        .collect()
}

struct TypeFilter {
    incl: Option<u32>,
    excl: Option<u32>,
}

impl TypeFilter {
    fn parse(inc: Option<&str>, exc: Option<&str>) -> Result<Option<Self>> {
        if inc.is_none() && exc.is_none() {
            return Ok(None);
        }
        let p = |s: &str| -> Result<u32> {
            parse_type_mask(s).ok_or_else(|| anyhow::anyhow!("unknown variant type in {s:?}"))
        };
        Ok(Some(Self { incl: inc.map(p).transpose()?, excl: exc.map(p).transpose()? }))
    }

    fn matches(mask: u32, ty: u32) -> bool {
        if mask == VT_REF {
            ty == VT_REF
        } else {
            ty & mask != 0
        }
    }

    fn passes(&self, refa: &str, alt: &str) -> bool {
        let ty = record_type(refa, alt);
        if let Some(i) = self.incl {
            if !Self::matches(i, ty) {
                return false;
            }
        }
        if let Some(e) = self.excl {
            if Self::matches(e, ty) {
                return false;
            }
        }
        true
    }
}

struct GenotypeFilter {
    mode: GtMode,
    negate: bool,
}

#[derive(Clone, Copy)]
enum GtMode {
    Hom,
    Het,
    Miss,
}

impl GenotypeFilter {
    fn parse(s: Option<&str>) -> Result<Option<Self>> {
        let Some(s) = s else { return Ok(None) };
        let (negate, body) = match s.strip_prefix('^') {
            Some(b) => (true, b),
            None => (false, s),
        };
        let mode = match body {
            "hom" => GtMode::Hom,
            "het" => GtMode::Het,
            "miss" => GtMode::Miss,
            o => bail!("--genotype: unknown mode {o:?} (expected hom|het|miss with optional ^)"),
        };
        Ok(Some(Self { mode, negate }))
    }

    /// bcftools: `-g hom` keeps sites with at least one hom genotype; `^hom`
    /// drops sites with any hom genotype (same for het and miss).
    fn passes(&self, format: &str, samples: &[&str], keep: Option<&[usize]>) -> bool {
        let gt_idx = format.split(':').position(|k| k == "GT");
        let mut any = false;
        let iter: Box<dyn Iterator<Item = &&str>> = match keep {
            Some(k) => Box::new(k.iter().filter_map(|&i| samples.get(i))),
            None => Box::new(samples.iter()),
        };
        for s in iter {
            let gt = gt_idx.and_then(|i| s.split(':').nth(i)).unwrap_or(".");
            let cls = gt_class(gt);
            let hit = matches!(
                (self.mode, cls),
                (GtMode::Hom, GtClass::Hom) | (GtMode::Het, GtClass::Het) | (GtMode::Miss, GtClass::Miss)
            );
            if hit {
                any = true;
                break;
            }
        }
        if self.negate { !any } else { any }
    }
}

enum GtClass {
    Hom,
    Het,
    Miss,
}

fn gt_class(gt: &str) -> GtClass {
    let alleles = gt_alleles(gt);
    if alleles.is_empty() || alleles.iter().any(|a| a.is_none()) {
        return GtClass::Miss;
    }
    let first = alleles[0];
    if alleles.iter().all(|a| *a == first) { GtClass::Hom } else { GtClass::Het }
}

fn any_phased(format: &str, samples: &[&str]) -> bool {
    let Some(gt_idx) = format.split(':').position(|k| k == "GT") else { return false };
    samples.iter().any(|s| s.split(':').nth(gt_idx).unwrap_or("").contains('|'))
}

fn any_uncalled(format: &str, samples: &[&str]) -> bool {
    let Some(gt_idx) = format.split(':').position(|k| k == "GT") else { return false };
    samples.iter().any(|s| matches!(gt_class(s.split(':').nth(gt_idx).unwrap_or(".")), GtClass::Miss))
}

fn is_private(format: &str, samples: &[&str]) -> bool {
    let Some(gt_idx) = format.split(':').position(|k| k == "GT") else { return false };
    let mut non_ref = 0usize;
    for s in samples {
        let gt = s.split(':').nth(gt_idx).unwrap_or(".");
        if gt_alleles(gt).iter().any(|a| a.is_some_and(|n| n > 0)) {
            non_ref += 1;
        }
    }
    non_ref == 1
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum AcMode {
    Nref,
    Alt1,
    Minor,
    Major,
    Nonmajor,
}

impl AcMode {
    fn parse(s: &str) -> Result<Self> {
        Ok(match s {
            "" | "nref" => Self::Nref,
            "alt1" => Self::Alt1,
            "minor" => Self::Minor,
            "major" => Self::Major,
            "nonmajor" => Self::Nonmajor,
            o => bail!("unknown allele-count mode {o:?} (expected nref|alt1|minor|major|nonmajor)"),
        })
    }
}

struct AcAfFilter {
    min_ac: Option<(u32, AcMode)>,
    max_ac: Option<(u32, AcMode)>,
    min_af: Option<(f64, AcMode)>,
    max_af: Option<(f64, AcMode)>,
}

impl AcAfFilter {
    fn build(a: &ViewArgs) -> Result<Option<Self>> {
        fn split(s: &str) -> (&str, &str) {
            match s.split_once(':') {
                Some((v, m)) => (v, m),
                None => (s, ""),
            }
        }
        let p_u32 = |s: &str| -> Result<(u32, AcMode)> {
            let (v, m) = split(s);
            Ok((v.parse().with_context(|| format!("bad allele count {v:?}"))?, AcMode::parse(m)?))
        };
        let p_f64 = |s: &str| -> Result<(f64, AcMode)> {
            let (v, m) = split(s);
            Ok((v.parse().with_context(|| format!("bad allele frequency {v:?}"))?, AcMode::parse(m)?))
        };
        let min_ac = a.min_ac.as_deref().map(p_u32).transpose()?;
        let max_ac = a.max_ac.as_deref().map(p_u32).transpose()?;
        let min_af = a.min_af.as_deref().map(p_f64).transpose()?;
        let max_af = a.max_af.as_deref().map(p_f64).transpose()?;
        if min_ac.is_none() && max_ac.is_none() && min_af.is_none() && max_af.is_none() {
            return Ok(None);
        }
        Ok(Some(Self { min_ac, max_ac, min_af, max_af }))
    }

    /// Allele counts per allele (index 0 = REF) and AN, from genotypes when
    /// present, else from INFO/AC and INFO/AN.
    fn counts(format: &str, samples: &[&str], keep: Option<&[usize]>, info: &str, n_alleles: usize) -> (Vec<u32>, u32) {
        let mut ac = vec![0u32; n_alleles.max(1)];
        let mut an = 0u32;
        let gt_idx = format.split(':').position(|k| k == "GT");
        if let (Some(gi), false) = (gt_idx, samples.is_empty()) {
            let iter: Box<dyn Iterator<Item = &&str>> = match keep {
                Some(k) => Box::new(k.iter().filter_map(|&i| samples.get(i))),
                None => Box::new(samples.iter()),
            };
            for s in iter {
                let gt = s.split(':').nth(gi).unwrap_or(".");
                for a in gt_alleles(gt).into_iter().flatten() {
                    an += 1;
                    if a < ac.len() {
                        ac[a] += 1;
                    }
                }
            }
            return (ac, an);
        }
        let mut alt_counts: Vec<u32> = Vec::new();
        for kv in info.split(';') {
            if let Some(v) = kv.strip_prefix("AC=") {
                alt_counts = v.split(',').map(|s| s.parse::<u32>().unwrap_or(0)).collect();
            } else if let Some(v) = kv.strip_prefix("AN=") {
                an = v.parse().unwrap_or(0);
            }
        }
        for (i, c) in alt_counts.iter().enumerate() {
            if i + 1 < ac.len() {
                ac[i + 1] = *c;
            }
        }
        let alt_sum: u32 = alt_counts.iter().sum();
        ac[0] = an.saturating_sub(alt_sum);
        (ac, an)
    }

    fn select(ac: &[u32], mode: AcMode) -> u32 {
        let nref: u32 = ac.iter().skip(1).sum();
        match mode {
            AcMode::Nref => nref,
            AcMode::Alt1 => ac.get(1).copied().unwrap_or(0),
            AcMode::Minor => {
                // bcftools: min over alleles of the biallelic view REF vs non-REF.
                ac[0].min(nref)
            }
            AcMode::Major => ac.iter().copied().max().unwrap_or(0),
            AcMode::Nonmajor => {
                let total: u32 = ac.iter().sum();
                total - ac.iter().copied().max().unwrap_or(0)
            }
        }
    }

    fn passes(&self, format: &str, samples: &[&str], keep: Option<&[usize]>, info: &str, n_alleles: usize) -> bool {
        let (ac, an) = Self::counts(format, samples, keep, info, n_alleles);
        let af = |mode: AcMode| -> f64 {
            if an > 0 { Self::select(&ac, mode) as f64 / an as f64 } else { 0.0 }
        };
        if let Some((v, m)) = self.min_ac {
            if Self::select(&ac, m) < v {
                return false;
            }
        }
        if let Some((v, m)) = self.max_ac {
            if Self::select(&ac, m) > v {
                return false;
            }
        }
        if let Some((v, m)) = self.min_af {
            if af(m) < v {
                return false;
            }
        }
        if let Some((v, m)) = self.max_af {
            if af(m) > v {
                return false;
            }
        }
        true
    }
}

fn parse_apply_filters(s: Option<&str>) -> Option<Vec<String>> {
    s.map(|x| x.split(',').map(|t| t.trim().to_string()).collect())
}

fn filter_passes_apply(filter_col: &str, allowed: &[String]) -> bool {
    if filter_col == "." || filter_col.is_empty() {
        return allowed.iter().any(|a| a == ".");
    }
    filter_col.split(';').any(|t| allowed.iter().any(|a| a == t))
}

/// Recompute INFO/AC and INFO/AN from the genotypes present on the line
/// (after sample subsetting), when those tags are declared in the header.
fn update_ac_an(line: &str, hdr: &HeaderInfo) -> String {
    if !hdr.info.contains_key("AC") && !hdr.info.contains_key("AN") {
        return line.to_string();
    }
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 10 {
        return line.to_string();
    }
    let Some(gi) = cols[8].split(':').position(|k| k == "GT") else { return line.to_string() };
    let n_alt = if cols[4] == "." { 0 } else { cols[4].split(',').count() };
    let mut ac = vec![0u32; n_alt + 1];
    let mut an = 0u32;
    for s in &cols[9..] {
        let gt = s.split(':').nth(gi).unwrap_or(".");
        for a in gt_alleles(gt).into_iter().flatten() {
            an += 1;
            if a < ac.len() {
                ac[a] += 1;
            }
        }
    }
    let mut items: Vec<(String, Option<String>)> = crate::vcf::alleles::split_info(cols[7])
        .into_iter()
        .map(|(k, v)| (k.to_string(), v.map(|s| s.to_string())))
        .collect();
    let has_ac = items.iter().any(|(k, _)| k == "AC");
    let has_an = items.iter().any(|(k, _)| k == "AN");
    let ac_val = ac[1..].iter().map(u32::to_string).collect::<Vec<_>>().join(",");
    for (k, v) in items.iter_mut() {
        if k == "AC" && n_alt > 0 {
            *v = Some(ac_val.clone());
        } else if k == "AN" {
            *v = Some(an.to_string());
        }
    }
    if !has_ac && hdr.info.contains_key("AC") && n_alt > 0 {
        items.push(("AC".into(), Some(ac_val)));
    }
    if !has_an && hdr.info.contains_key("AN") {
        items.push(("AN".into(), Some(an.to_string())));
    }
    let mut out = cols[..7].join("\t");
    out.push('\t');
    out.push_str(&crate::vcf::alleles::join_info(&items));
    for c in &cols[8..] {
        out.push('\t');
        out.push_str(c);
    }
    out
}

/// `-a`: drop ALT alleles not observed in any genotype; `--trim-unseen-allele`:
/// drop `<*>`/`<NON_REF>` when unobserved. INFO/FORMAT Number=A/R/G fields and
/// GT indices are rewritten for the new allele set.
fn trim_alt_alleles(line: &str, hdr: &HeaderInfo, trim_all: bool, trim_unseen: bool) -> String {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 8 || cols[4] == "." {
        return line.to_string();
    }
    let alts: Vec<&str> = cols[4].split(',').collect();
    let format_str = if cols.len() > 8 { cols[8] } else { "" };
    let samples: Vec<&str> = if cols.len() > 9 { cols[9..].to_vec() } else { Vec::new() };
    let gt_idx = format_str.split(':').position(|k| k == "GT");
    let mut used = vec![false; alts.len() + 1];
    used[0] = true;
    if let Some(gi) = gt_idx {
        for s in &samples {
            let gt = s.split(':').nth(gi).unwrap_or(".");
            for a in gt_alleles(gt).into_iter().flatten() {
                if a < used.len() {
                    used[a] = true;
                }
            }
        }
    }
    let types = allele_types(cols[3], cols[4]);
    let mut keep_alt = vec![true; alts.len()];
    for (i, a) in alts.iter().enumerate() {
        let unseen_symbolic = *a == "<*>" || *a == "<NON_REF>";
        let observed = used[i + 1];
        if trim_all && !observed && gt_idx.is_some() {
            keep_alt[i] = false;
        }
        if trim_unseen && unseen_symbolic && !observed {
            keep_alt[i] = false;
        }
        let _ = types.get(i);
    }
    if keep_alt.iter().all(|k| *k) {
        return line.to_string();
    }
    let n_old = alts.len() + 1;
    let mut map: Vec<Option<usize>> = vec![Some(0)];
    let mut new_alts: Vec<&str> = Vec::new();
    for (i, a) in alts.iter().enumerate() {
        if keep_alt[i] {
            new_alts.push(a);
            map.push(Some(new_alts.len()));
        } else {
            map.push(None);
        }
    }
    let n_new = new_alts.len() + 1;
    let info = remap_info(cols[7], hdr, n_old, n_new, &map);
    let alt_field = if new_alts.is_empty() { ".".to_string() } else { new_alts.join(",") };
    let mut out = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", cols[0], cols[1], cols[2], cols[3], alt_field, cols[5], cols[6], info);
    if !format_str.is_empty() {
        out.push('\t');
        out.push_str(format_str);
        for s in remap_samples(format_str, &samples, hdr, n_old, n_new, &map) {
            out.push('\t');
            out.push_str(&s);
        }
    }
    out
}

#[cfg(test)]
#[path = "../../../tests/unit/cli_commands_view.rs"]
mod tests;
