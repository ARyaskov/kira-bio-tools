use crate::annotate::postproc::{RegionFilter, version_header_line};
use crate::cli::args::NormArgs;
use crate::vcf::alleles::{gt_alleles, remap_info, remap_samples, remap_value};
use crate::vcf::header::{FieldNumber, HeaderInfo};
use crate::vcf::sink::{OutputKind, parse_output_type};
use crate::vcf::variant_type::{VT_INDEL, VT_SNP, record_type};
use crate::vcf::{UnifiedVcfReader, VcfSink};
use anyhow::{Context, Result, bail};
use std::path::Path;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum MultiMode { None, SplitAll, SplitSnps, SplitIndels, JoinAll, JoinSnps, JoinIndels, JoinBoth }

impl MultiMode {
    fn parse(s: Option<&str>) -> Result<Self> {
        let Some(s) = s else { return Ok(Self::None); };
        Ok(match s {
            "-" | "-any" | "-both" => Self::SplitAll,
            "-snps" => Self::SplitSnps,
            "-indels" => Self::SplitIndels,
            "+" | "+any" => Self::JoinAll,
            "+both" => Self::JoinBoth,
            "+snps" => Self::JoinSnps,
            "+indels" => Self::JoinIndels,
            other => bail!("--multiallelics: unknown {other:?}"),
        })
    }

    fn is_join(self) -> bool {
        matches!(self, Self::JoinAll | Self::JoinSnps | Self::JoinIndels | Self::JoinBoth)
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum CheckRef { Exit, Warn, Exclude, Set }

impl CheckRef {
    fn parse(s: &str) -> Result<Self> {
        Ok(match s {
            "e" => Self::Exit, "w" => Self::Warn, "x" => Self::Exclude, "s" => Self::Set,
            o => bail!("-c, --check-ref: unknown {o:?} (expected e|w|x|s)"),
        })
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum RmDup { None, All, Exact, Snps, Indels, Both, AnyId }

impl RmDup {
    fn parse(s: Option<&str>) -> Result<Self> {
        let Some(s) = s else { return Ok(Self::None); };
        Ok(match s {
            "none" => Self::None,
            "exact" => Self::Exact,
            "snps" => Self::Snps,
            "indels" => Self::Indels,
            "both" => Self::Both,
            "any" | "all" => Self::All,
            "id" => Self::AnyId,
            o => bail!("--rm-dup: unknown {o:?} (none|exact|snps|indels|both|all|id)"),
        })
    }
}

#[derive(Clone, Debug)]
struct SplitOpts {
    missing_for_overlap: bool,
    keep_sum_keys: Vec<String>,
    old_rec_tag: Option<String>,
}

pub fn cmd_norm(args: NormArgs) -> Result<()> {
    let multi = MultiMode::parse(args.multiallelics.as_deref())?;
    let check_ref = CheckRef::parse(&args.check_ref)?;
    let rm_dup = RmDup::parse(args.rm_dup.as_deref())?;
    let split_opts = SplitOpts {
        missing_for_overlap: args.multi_overlaps == ".",
        keep_sum_keys: args
            .keep_sum
            .as_deref()
            .map(|s| s.split(',').map(|t| t.trim().to_string()).filter(|t| !t.is_empty()).collect())
            .unwrap_or_default(),
        old_rec_tag: args.old_rec_tag.clone(),
    };

    let mut fasta = args.fasta_ref.as_ref().map(|p| load_fasta(p)).transpose()?;

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

    let out_path = args.output.clone();
    let kind = args.output_type.as_deref().map(parse_output_type).transpose()?.unwrap_or(OutputKind::Vcf);

    let mut reader = UnifiedVcfReader::open(&args.input).context("open input")?;
    let headers = reader.header()?;
    let hdr = HeaderInfo::parse(&headers);

    let mut out_headers: Vec<String> = Vec::with_capacity(headers.len() + 2);
    let version = version_header_line();
    for h in &headers {
        if h.starts_with("#CHROM") {
            if let Some(tag) = &split_opts.old_rec_tag {
                if !hdr.info.contains_key(tag) {
                    out_headers.push(format!(
                        "##INFO=<ID={tag},Number=1,Type=String,Description=\"Original record before normalization (CHROM|POS|REF|ALT|USED_ALT_IDX)\">"
                    ));
                }
            }
            if !args.no_version {
                out_headers.push(version.clone());
            }
        }
        out_headers.push(h.clone());
    }
    let mut sink = VcfSink::open(out_path.as_deref(), kind, &out_headers)?;
    sink.write_header(&out_headers)?;

    let mut dedup = DedupWindow::default();
    let mut joiner = Joiner::new(multi);
    let mut outbuf = OutBuf::new(args.site_win);

    while let Some(line) = reader.read_line()? {
        if line.is_empty() || line.as_bytes()[0] == b'#' { continue; }
        if let Some(rf) = &region_filter {
            if !rf.line_passes_mode(&line, args.regions_overlap) { continue; }
        }
        if let Some(tf) = &target_filter {
            if tf.line_passes_mode(&line, 0) == targets_inverse { continue; }
        }
        let mut fields = line.splitn(3, '\t');
        let in_chrom = fields.next().unwrap_or("");
        let in_pos: u32 = fields
            .next()
            .and_then(|p| p.parse().ok())
            .ok_or_else(|| anyhow::anyhow!("norm: malformed POS in record {line:?}"))?;
        if let Some(fa) = fasta.as_mut() {
            // Chromosome-sorted input keeps one contig in memory at a time.
            if outbuf.chrom != in_chrom {
                fa.evict_except(in_chrom);
            }
        }
        let records = expand_record(&line, multi, args.atomize, &split_opts, &hdr)?;
        for rec in records {
            let mut cols: Vec<String> = rec.split('\t').map(|s| s.to_string()).collect();
            if cols.len() < 8 { continue; }

            if let Some(fa) = &fasta {
                match verify_ref(&cols, fa, check_ref) {
                    RefAction::Keep => {}
                    RefAction::Skip => continue,
                    RefAction::Fix(new_ref) => fix_ref(&mut cols, &new_ref, &hdr),
                    RefAction::Fail(msg) => bail!("-c e: {msg}"),
                    RefAction::Warn(msg) => eprintln!("[norm] warn: {msg}"),
                }
                if !args.do_not_normalize && fa.has(&cols[0]) {
                    if let Some(aligned) = left_align(&cols, fa) {
                        cols = aligned;
                    }
                    // A reference-only record keeps just its first base (bcftools right-trims it).
                    if (cols[4] == "." || cols[4].is_empty()) && cols[3].len() > 1 {
                        cols[3].truncate(1);
                    }
                }
            }

            for out_rec in joiner.push(cols, &hdr) {
                outbuf.push(out_rec, in_pos, &mut sink, &mut dedup, rm_dup)?;
            }
        }
    }
    for out_rec in joiner.finish(&hdr) {
        outbuf.push(out_rec, u32::MAX, &mut sink, &mut dedup, rm_dup)?;
    }
    outbuf.flush_all(&mut sink, &mut dedup, rm_dup)?;
    sink.finish()?;

    if let (Some(kind_s), Some(out)) = (args.write_index.as_deref(), out_path.as_deref()) {
        if matches!(kind, OutputKind::VcfGz(_) | OutputKind::Bcf(_)) && out != Path::new("-") {
            let (ik, ext) = if kind_s == "tbi" { (crate::csi::IndexKind::Tbi, "tbi") } else { (crate::csi::IndexKind::Csi, "csi") };
            let idx = std::path::PathBuf::from(format!("{}.{}", out.display(), ext));
            crate::csi::build_index(out, &idx, ik, None).with_context(|| format!("-W: write {}", idx.display()))?;
        }
    }
    Ok(())
}

/// Split a multiallelic record (`-m -`) and/or atomize MNPs.
fn expand_record(line: &str, multi: MultiMode, atomize: bool, opts: &SplitOpts, hdr: &HeaderInfo) -> Result<Vec<String>> {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 8 { return Ok(vec![line.to_string()]); }
    let refa = cols[3];
    let alt = cols[4];
    let alts: Vec<&str> = alt.split(',').collect();
    let vt = record_type(refa, alt);

    let split_now = alts.len() > 1
        && match multi {
            MultiMode::SplitAll => true,
            MultiMode::SplitSnps => vt & VT_SNP != 0 && vt & VT_INDEL == 0,
            MultiMode::SplitIndels => vt & VT_INDEL != 0,
            _ => false,
        };

    let mut result: Vec<String> = if split_now {
        split_multiallelic(&cols, &alts, opts, hdr)
    } else {
        vec![line.to_string()]
    };

    if atomize {
        let mut atomic = Vec::new();
        for rec in result {
            atomic.extend(atomize_record(&rec, opts));
        }
        result = atomic;
    }

    Ok(result)
}

/// `-m -`: one biallelic record per ALT. Number=A/R/G fields are subset to
/// `[REF, ALT_i]`; other ALTs in GT become REF (or missing for `*` with
/// `--multi-overlaps .`); `--keep-sum` folds the other alleles' counts into REF.
fn split_multiallelic(cols: &[&str], alts: &[&str], opts: &SplitOpts, hdr: &HeaderInfo) -> Vec<String> {
    let n_old = alts.len() + 1;
    let info = cols[7];
    let format = if cols.len() > 8 { cols[8] } else { "" };
    let samples: Vec<&str> = if cols.len() > 9 { cols[9..].to_vec() } else { Vec::new() };
    let mut out = Vec::with_capacity(alts.len());
    for (i, alt_i) in alts.iter().enumerate() {
        let mut map: Vec<Option<usize>> = vec![None; n_old];
        map[0] = Some(0);
        map[i + 1] = Some(1);

        let mut new_info = remap_info_keep_sum(info, hdr, n_old, i + 1, &map, &opts.keep_sum_keys);
        if let Some(tag) = &opts.old_rec_tag {
            let extra = format!("{}={}|{}|{}|{}|{}", tag, cols[0], cols[1], cols[3], cols[4], i + 1);
            if new_info == "." || new_info.is_empty() { new_info = extra; } else { new_info.push(';'); new_info.push_str(&extra); }
        }

        let mut s = String::with_capacity(cols.iter().map(|c| c.len() + 1).sum());
        for (idx, col) in cols.iter().take(8).enumerate() {
            if idx > 0 { s.push('\t'); }
            match idx {
                4 => s.push_str(alt_i),
                7 => s.push_str(&new_info),
                _ => s.push_str(col),
            }
        }
        if !format.is_empty() {
            s.push('\t');
            s.push_str(format);
            // GT: other ALTs collapse to REF (or '.' for '*'), everything else is subset.
            let mut gt_map: Vec<Option<usize>> = vec![Some(0); n_old];
            gt_map[i + 1] = Some(1);
            for (j, a) in alts.iter().enumerate() {
                if j != i && *a == "*" && opts.missing_for_overlap {
                    gt_map[j + 1] = None;
                }
            }
            let keys: Vec<&str> = format.split(':').collect();
            for smp in &samples {
                let parts: Vec<&str> = smp.split(':').collect();
                let mut new_parts: Vec<String> = Vec::with_capacity(parts.len());
                for (k, v) in parts.iter().enumerate() {
                    let key = keys.get(k).copied().unwrap_or("");
                    if key == "GT" {
                        new_parts.push(crate::vcf::alleles::remap_gt(v, &gt_map));
                    } else {
                        let num = hdr.format_number(key);
                        if num.is_per_allele() {
                            let nv = if opts.keep_sum_keys.iter().any(|t| t == key) && num == FieldNumber::R {
                                keep_sum_r(v, i + 1)
                            } else {
                                remap_value(v, num, n_old, 2, &map)
                            };
                            new_parts.push(nv.unwrap_or_else(|| v.to_string()));
                        } else {
                            new_parts.push(v.to_string());
                        }
                    }
                }
                s.push('\t');
                s.push_str(&new_parts.join(":"));
            }
        }
        out.push(s);
    }
    out
}

/// Number=R value for allele `alt_idx` with the other alleles' counts folded into REF.
fn keep_sum_r(value: &str, alt_idx: usize) -> Option<String> {
    if value == "." { return Some(value.to_string()); }
    let vals: Vec<&str> = value.split(',').collect();
    if alt_idx >= vals.len() { return None; }
    let mut rest = 0f64;
    let mut all_int = true;
    for (j, v) in vals.iter().enumerate() {
        if j == alt_idx { continue; }
        if let Ok(n) = v.parse::<f64>() {
            rest += n;
            if v.contains('.') { all_int = false; }
        }
    }
    let rest_s = if all_int { format!("{}", rest as i64) } else { format!("{}", rest) };
    Some(format!("{},{}", rest_s, vals[alt_idx]))
}

fn remap_info_keep_sum(info: &str, hdr: &HeaderInfo, n_old: usize, alt_idx: usize, map: &[Option<usize>], keep_sum: &[String]) -> String {
    if keep_sum.is_empty() {
        return remap_info(info, hdr, n_old, 2, map);
    }
    if info == "." || info.is_empty() { return info.to_string(); }
    let mut out = String::with_capacity(info.len());
    let mut first = true;
    for kv in info.split(';') {
        if !first { out.push(';'); }
        first = false;
        match kv.split_once('=') {
            Some((k, v)) => {
                let num = hdr.info_number(k);
                let nv = if keep_sum.iter().any(|t| t == k) && num == FieldNumber::R {
                    keep_sum_r(v, alt_idx)
                } else if num.is_per_allele() {
                    remap_value(v, num, n_old, 2, map)
                } else {
                    None
                };
                out.push_str(k);
                out.push('=');
                out.push_str(nv.as_deref().unwrap_or(v));
            }
            None => out.push_str(kv),
        }
    }
    out
}

/// `--atomize`: MNPs become one SNP per differing base (complex records are
/// left as they are).
fn atomize_record(line: &str, opts: &SplitOpts) -> Vec<String> {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 8 { return vec![line.to_string()]; }
    let refa = cols[3].as_bytes();
    let alt = cols[4];
    if alt.contains(',') { return vec![line.to_string()]; }
    let alta = alt.as_bytes();
    let Ok(pos) = cols[1].parse::<u32>() else { return vec![line.to_string()] };
    if refa.len() < 2 || alta.len() < 2 || refa.len() != alta.len() { return vec![line.to_string()]; }

    let mut out = Vec::new();
    for (k, (&r, &a)) in refa.iter().zip(alta.iter()).enumerate() {
        if r.eq_ignore_ascii_case(&a) { continue; }
        let mut c: Vec<String> = cols.iter().map(|s| s.to_string()).collect();
        c[1] = (pos + k as u32).to_string();
        c[3] = (r as char).to_string();
        c[4] = (a as char).to_string();
        if let Some(tag) = &opts.old_rec_tag {
            let extra = format!("{}={}|{}|{}|{}|1", tag, cols[0], cols[1], cols[3], cols[4]);
            if c[7] == "." || c[7].is_empty() { c[7] = extra; } else { c[7].push(';'); c[7].push_str(&extra); }
        }
        out.push(c.join("\t"));
    }
    if out.is_empty() { out.push(line.to_string()); }
    out
}

/// `-m +`: join records at the same position into one multiallelic record.
struct Joiner {
    mode: MultiMode,
    group: Vec<Vec<String>>,
}

impl Joiner {
    fn new(mode: MultiMode) -> Self {
        Self { mode, group: Vec::new() }
    }

    fn same_site(&self, cols: &[String]) -> bool {
        match self.group.first() {
            Some(first) => first[0] == cols[0] && first[1] == cols[1],
            None => false,
        }
    }

    fn eligible(&self, cols: &[String]) -> bool {
        if cols[4] == "." { return false; }
        let vt = record_type(&cols[3], &cols[4]);
        match self.mode {
            MultiMode::JoinAll => true,
            MultiMode::JoinSnps => vt & VT_SNP != 0 && vt & VT_INDEL == 0,
            MultiMode::JoinIndels => vt & VT_INDEL != 0,
            MultiMode::JoinBoth => vt & VT_SNP != 0 || vt & VT_INDEL != 0,
            _ => false,
        }
    }

    /// Feed a record; returns records that are ready to be written.
    fn push(&mut self, cols: Vec<String>, hdr: &HeaderInfo) -> Vec<Vec<String>> {
        if !self.mode.is_join() {
            return vec![cols];
        }
        let mut out = Vec::new();
        if !self.same_site(&cols) {
            out.extend(self.flush(hdr));
        }
        if self.eligible(&cols) {
            self.group.push(cols);
        } else {
            out.extend(self.flush(hdr));
            out.push(cols);
        }
        out
    }

    fn finish(&mut self, hdr: &HeaderInfo) -> Vec<Vec<String>> {
        self.flush(hdr)
    }

    fn flush(&mut self, hdr: &HeaderInfo) -> Vec<Vec<String>> {
        let group = std::mem::take(&mut self.group);
        if group.is_empty() {
            return Vec::new();
        }
        if self.mode == MultiMode::JoinBoth {
            // SNPs and indels are joined separately.
            let (snps, rest): (Vec<_>, Vec<_>) = group.into_iter().partition(|c| {
                let vt = record_type(&c[3], &c[4]);
                vt & VT_SNP != 0 && vt & VT_INDEL == 0
            });
            let mut out = Vec::new();
            for g in [snps, rest] {
                if !g.is_empty() {
                    out.push(join_group(g, hdr));
                }
            }
            return out;
        }
        vec![join_group(group, hdr)]
    }
}

/// Join biallelic records at one site (bcftools `merge_alleles`): the longest
/// REF wins, shorter alleles are padded with its tail, A/R/G fields are
/// expanded to the union, GT alleles are combined per sample.
fn join_group(group: Vec<Vec<String>>, hdr: &HeaderInfo) -> Vec<String> {
    if group.len() == 1 {
        return group.into_iter().next().unwrap();
    }
    let longest = group.iter().map(|c| c[3].len()).max().unwrap_or(0);
    let long_ref: String = group.iter().find(|c| c[3].len() == longest).map(|c| c[3].clone()).unwrap_or_default();
    // REFs must be prefixes of the longest; otherwise leave the records alone.
    let compatible = group.iter().all(|c| long_ref.as_bytes()[..c[3].len()].eq_ignore_ascii_case(c[3].as_bytes()));
    if !compatible {
        eprintln!("[norm] warn: incompatible REF alleles at {}:{}, not joined", group[0][0], group[0][1]);
        // Emit them individually by returning the first and pushing back the rest is not possible here;
        // concatenate lines with a newline separator so the caller writes them verbatim.
        let mut merged = group[0].clone();
        for g in &group[1..] {
            merged.last_mut().unwrap().push('\n');
            merged.last_mut().unwrap().push_str(&g.join("\t"));
        }
        return merged;
    }

    let mut new_alts: Vec<String> = Vec::new();
    let mut maps: Vec<Vec<Option<usize>>> = Vec::with_capacity(group.len());
    for c in &group {
        let tail = &long_ref[c[3].len()..];
        let mut map = vec![Some(0)];
        for a in c[4].split(',') {
            let padded = if a.starts_with('<') || a == "*" || a.starts_with('.') { a.to_string() } else { format!("{a}{tail}") };
            let idx = match new_alts.iter().position(|x| x.eq_ignore_ascii_case(&padded)) {
                Some(i) => i,
                None => {
                    new_alts.push(padded);
                    new_alts.len() - 1
                }
            };
            map.push(Some(idx + 1));
        }
        maps.push(map);
    }
    let n_new = new_alts.len() + 1;

    // ID: unique ids joined; QUAL: max; FILTER: union (PASS only if all PASS).
    let mut ids: Vec<&str> = Vec::new();
    for c in &group {
        for id in c[2].split(';') {
            if id != "." && !ids.contains(&id) { ids.push(id); }
        }
    }
    let id = if ids.is_empty() { ".".to_string() } else { ids.join(";") };
    let qual = group
        .iter()
        .filter_map(|c| c[5].parse::<f64>().ok())
        .fold(None, |m: Option<f64>, q| Some(m.map_or(q, |v| v.max(q))))
        .map(|q| if q.fract() == 0.0 { format!("{}", q as i64) } else { q.to_string() })
        .unwrap_or_else(|| ".".to_string());
    let filter = merge_filters(group.iter().map(|c| c[6].as_str()));

    // INFO: A/R/G expanded and combined; other tags from the first record that has them.
    let mut info_keys: Vec<String> = Vec::new();
    let mut info_vals: Vec<Option<String>> = Vec::new();
    for (gi, c) in group.iter().enumerate() {
        let n_old = maps[gi].len();
        for (k, v) in crate::vcf::alleles::split_info(&c[7]) {
            let num = hdr.info_number(k);
            let expanded = match v {
                Some(v) if num.is_per_allele() => remap_value(v, num, n_old, n_new, &maps[gi]),
                Some(v) => Some(v.to_string()),
                None => None,
            };
            match info_keys.iter().position(|x| x == k) {
                Some(p) => {
                    if num.is_per_allele() {
                        if let (Some(cur), Some(new)) = (info_vals[p].clone(), expanded) {
                            info_vals[p] = Some(merge_missing(&cur, &new));
                        }
                    }
                }
                None => {
                    info_keys.push(k.to_string());
                    info_vals.push(expanded);
                }
            }
        }
    }
    let info_items: Vec<(String, Option<String>)> = info_keys.into_iter().zip(info_vals).collect();
    let info = crate::vcf::alleles::join_info(&info_items);

    let mut out = vec![
        group[0][0].clone(),
        group[0][1].clone(),
        id,
        long_ref.clone(),
        new_alts.join(","),
        qual,
        filter,
        info,
    ];

    if group[0].len() > 8 {
        // FORMAT union, then per-sample combination.
        let mut keys: Vec<String> = Vec::new();
        for c in &group {
            for k in c[8].split(':') {
                if !keys.iter().any(|x| x == k) { keys.push(k.to_string()); }
            }
        }
        let n_samples = group[0].len().saturating_sub(9);
        let remapped: Vec<Vec<String>> = group
            .iter()
            .enumerate()
            .map(|(gi, c)| {
                let samples: Vec<&str> = c[9..].iter().map(String::as_str).collect();
                remap_samples(&c[8], &samples, hdr, maps[gi].len(), n_new, &maps[gi])
            })
            .collect();
        out.push(keys.join(":"));
        for si in 0..n_samples {
            let mut vals: Vec<Option<String>> = vec![None; keys.len()];
            for (gi, c) in group.iter().enumerate() {
                let local_keys: Vec<&str> = c[8].split(':').collect();
                let Some(smp) = remapped[gi].get(si) else { continue };
                let parts: Vec<&str> = smp.split(':').collect();
                for (li, lk) in local_keys.iter().enumerate() {
                    let Some(v) = parts.get(li) else { continue };
                    let ki = keys.iter().position(|k| k == lk).unwrap();
                    if *lk == "GT" {
                        vals[ki] = Some(match &vals[ki] {
                            Some(cur) => combine_gt(cur, v),
                            None => v.to_string(),
                        });
                    } else if hdr.format_number(lk).is_per_allele() {
                        vals[ki] = Some(match &vals[ki] {
                            Some(cur) => merge_missing(cur, v),
                            None => v.to_string(),
                        });
                    } else if vals[ki].is_none() {
                        vals[ki] = Some(v.to_string());
                    }
                }
            }
            out.push(vals.into_iter().map(|v| v.unwrap_or_else(|| ".".into())).collect::<Vec<_>>().join(":"));
        }
    }
    out
}

/// Fill missing entries of `cur` from `new` (comma vectors of equal length).
fn merge_missing(cur: &str, new: &str) -> String {
    let a: Vec<&str> = cur.split(',').collect();
    let b: Vec<&str> = new.split(',').collect();
    if a.len() != b.len() { return cur.to_string(); }
    a.iter().zip(b.iter()).map(|(x, y)| if *x == "." { *y } else { *x }).collect::<Vec<_>>().join(",")
}

/// Combine two remapped GTs of one sample. Non-reference alleles of the later
/// record take over reference (or missing) slots: phased genotypes keep their
/// haplotype, unphased ones use any free slot and are reported sorted. When no
/// slot is free the earlier record wins.
fn combine_gt(cur: &str, new: &str) -> String {
    let a = gt_alleles(cur);
    let b = gt_alleles(new);
    if a.len() != b.len() { return cur.to_string(); }
    let phased = cur.contains('|') && new.contains('|');
    let mut out: Vec<Option<usize>> = a.clone();
    for i in 0..b.len() {
        let Some(nb) = b[i] else { continue };
        if nb == 0 {
            if out[i].is_none() { out[i] = Some(0); }
            continue;
        }
        let free = |x: &Option<usize>| x.is_none_or(|v| v == 0);
        if free(&out[i]) {
            out[i] = Some(nb);
        } else if !phased {
            if let Some(slot) = out.iter().position(free) {
                out[slot] = Some(nb);
            }
        }
    }
    if !phased {
        out.sort_by_key(|x| x.map(|v| v as i64).unwrap_or(i64::MAX));
    }
    let sep = if phased { "|" } else { "/" };
    out.iter().map(|x| x.map(|v| v.to_string()).unwrap_or_else(|| ".".into())).collect::<Vec<_>>().join(sep)
}

fn merge_filters<'a, I: Iterator<Item = &'a str>>(filters: I) -> String {
    let mut all_pass = true;
    let mut set: Vec<&str> = Vec::new();
    let mut any = false;
    for f in filters {
        any = true;
        if f == "PASS" { continue; }
        all_pass = false;
        if f == "." || f.is_empty() { continue; }
        for t in f.split(';') {
            if !set.contains(&t) { set.push(t); }
        }
    }
    if !any || (all_pass && set.is_empty()) {
        return if any { "PASS".into() } else { ".".into() };
    }
    if set.is_empty() { ".".into() } else { set.join(";") }
}

/// `-d`: duplicates within the current position, classified with the shared
/// variant typing.
#[derive(Default)]
struct DedupWindow {
    pos_key: (String, String),
    recs: Vec<Vec<String>>,
}

impl DedupWindow {
    fn push(&mut self, rec: Vec<String>) {
        let key = (rec[0].clone(), rec[1].clone());
        if key != self.pos_key {
            self.pos_key = key;
            self.recs.clear();
        }
        self.recs.push(rec);
    }

    fn is_dup(&self, rec: &[String], mode: RmDup) -> bool {
        if mode == RmDup::None { return false; }
        if (rec[0].as_str(), rec[1].as_str()) != (self.pos_key.0.as_str(), self.pos_key.1.as_str()) {
            return false;
        }
        let vt = record_type(&rec[3], &rec[4]);
        for prev in &self.recs {
            let pvt = record_type(&prev[3], &prev[4]);
            let hit = match mode {
                RmDup::Exact => rec[3].eq_ignore_ascii_case(&prev[3]) && rec[4].eq_ignore_ascii_case(&prev[4]),
                RmDup::All => true,
                RmDup::Snps => vt & VT_SNP != 0 && pvt & VT_SNP != 0,
                RmDup::Indels => vt & VT_INDEL != 0 && pvt & VT_INDEL != 0,
                RmDup::Both => (vt & VT_SNP != 0 && pvt & VT_SNP != 0) || (vt & VT_INDEL != 0 && pvt & VT_INDEL != 0),
                RmDup::AnyId => rec[2] != "." && rec[2] == prev[2],
                RmDup::None => false,
            };
            if hit { return true; }
        }
        false
    }
}

/// Output records ordered by position: normalization moves records, so a
/// record is held back until the input is `win` bases past it (bcftools
/// `--site-win`). Ties keep input order; `-d` is applied on the sorted stream.
struct OutBuf {
    win: u32,
    chrom: String,
    pending: std::collections::VecDeque<(u32, Vec<String>)>,
}

impl OutBuf {
    fn new(win: u32) -> Self {
        Self { win, chrom: String::new(), pending: std::collections::VecDeque::new() }
    }

    fn push(&mut self, rec: Vec<String>, in_pos: u32, sink: &mut VcfSink, dedup: &mut DedupWindow, rm_dup: RmDup) -> Result<()> {
        if rec[0] != self.chrom {
            self.flush_all(sink, dedup, rm_dup)?;
            self.chrom = rec[0].clone();
        }
        let pos: u32 = rec[1].parse().map_err(|_| anyhow::anyhow!("norm: invalid POS {:?}", rec[1]))?;
        let at = self.pending.partition_point(|(p, _)| *p <= pos);
        self.pending.insert(at, (pos, rec));
        let limit = in_pos.saturating_sub(self.win);
        while self.pending.front().is_some_and(|(p, _)| *p < limit) {
            let (_, r) = self.pending.pop_front().unwrap();
            Self::emit(r, sink, dedup, rm_dup)?;
        }
        Ok(())
    }

    fn flush_all(&mut self, sink: &mut VcfSink, dedup: &mut DedupWindow, rm_dup: RmDup) -> Result<()> {
        while let Some((_, r)) = self.pending.pop_front() {
            Self::emit(r, sink, dedup, rm_dup)?;
        }
        Ok(())
    }

    fn emit(r: Vec<String>, sink: &mut VcfSink, dedup: &mut DedupWindow, rm_dup: RmDup) -> Result<()> {
        if dedup.is_dup(&r, rm_dup) {
            return Ok(());
        }
        sink.write_line(&r.join("\t"))?;
        dedup.push(r);
        Ok(())
    }
}

enum RefAction { Keep, Skip, Fix(String), Fail(String), Warn(String) }

/// Reference contigs are loaded on demand through the `.fai` index.
pub(crate) type Fasta = crate::fasta::IndexedFasta;

pub(crate) fn load_fasta(p: &Path) -> Result<Fasta> {
    Fasta::open(p)
}

fn verify_ref(cols: &[String], fa: &Fasta, mode: CheckRef) -> RefAction {
    let pos: u32 = match cols[1].parse() { Ok(v) => v, Err(_) => return RefAction::Keep };
    let refa = cols[3].as_bytes();
    if cols[3].starts_with('<') || cols[3] == "N" && refa.len() == 1 && cols[4].starts_with('<') {
        return RefAction::Keep;
    }
    let Some(seq) = fa.slice(&cols[0], pos, refa.len()) else {
        let msg = format!("REF mismatch at {}:{} (no fasta sequence)", cols[0], pos);
        return match mode { CheckRef::Exit => RefAction::Fail(msg), CheckRef::Warn => RefAction::Warn(msg), CheckRef::Exclude => RefAction::Skip, CheckRef::Set => RefAction::Keep };
    };
    if seq.eq_ignore_ascii_case(refa) { return RefAction::Keep; }
    let new_ref = std::str::from_utf8(seq).unwrap_or("N").to_string();
    let msg = format!("REF mismatch at {}:{} (vcf={}, fasta={})", cols[0], pos, cols[3], new_ref);
    match mode {
        CheckRef::Exit => RefAction::Fail(msg),
        CheckRef::Warn => RefAction::Warn(msg),
        CheckRef::Exclude => RefAction::Skip,
        CheckRef::Set => RefAction::Fix(new_ref),
    }
}

/// `-c s`: set REF from the reference. When the reference base is one of the
/// ALTs, REF and that ALT are swapped and genotypes and A/R/G fields follow.
fn fix_ref(cols: &mut Vec<String>, new_ref: &str, hdr: &HeaderInfo) {
    let alts: Vec<String> = cols[4].split(',').map(|s| s.to_string()).collect();
    let swap_idx = alts.iter().position(|a| a.eq_ignore_ascii_case(new_ref));
    let Some(k) = swap_idx else {
        cols[3] = new_ref.to_string();
        return;
    };
    // Permutation: old allele k+1 -> 0, old 0 -> k+1.
    let n = alts.len() + 1;
    let mut map: Vec<Option<usize>> = (0..n).map(Some).collect();
    map[0] = Some(k + 1);
    map[k + 1] = Some(0);
    let mut new_alts = alts.clone();
    new_alts[k] = cols[3].clone();
    cols[3] = new_ref.to_string();
    cols[4] = new_alts.join(",");
    cols[7] = remap_info(&cols[7], hdr, n, n, &map);
    if cols.len() > 9 {
        let format = cols[8].clone();
        let samples: Vec<&str> = cols[9..].iter().map(String::as_str).collect();
        let new_samples = remap_samples(&format, &samples, hdr, n, n, &map);
        for (i, s) in new_samples.into_iter().enumerate() {
            cols[9 + i] = s;
        }
    }
}

fn is_symbolic_alt(a: &str) -> bool {
    a.starts_with('<') || a == "*" || a == "."
}

/// Largest symbolic SV span left-alignable without materialising the whole
/// reference allele. Beyond this the record passes through unshifted.
const MAX_SV_SPAN: u32 = 1_000_000;

fn left_align(cols: &[String], fa: &Fasta) -> Option<Vec<String>> {
    let chr = cols[0].as_str();
    let mut pos: u32 = cols[1].parse().ok()?;
    let alt_raw: Vec<&str> = cols[4].split(',').collect();
    if alt_raw.is_empty() {
        return None;
    }

    if alt_raw.iter().any(|a| is_symbolic_alt(a)) {
        return left_align_symbolic(cols, chr, pos, &alt_raw, fa);
    }

    let mut r: Vec<u8> = cols[3].as_bytes().to_ascii_uppercase();
    let mut alts: Vec<Vec<u8>> = alt_raw.iter().map(|a| a.as_bytes().to_ascii_uppercase()).collect();
    if r.is_empty() || alts.iter().any(|a| a.is_empty()) {
        return None;
    }

    left_shift(&mut r, &mut alts, &mut pos, chr, fa)?;
    // Left-trim a shared leading base while every allele keeps length >= 2.
    while r.len() >= 2 && alts.iter().all(|a| a.len() >= 2 && a[0] == r[0]) {
        r.remove(0);
        for a in alts.iter_mut() { a.remove(0); }
        pos += 1;
    }

    let new_ref = String::from_utf8(r).ok()?;
    let new_alt = alts
        .into_iter()
        .map(|a| String::from_utf8(a).unwrap_or_default())
        .collect::<Vec<_>>()
        .join(",");
    let mut out: Vec<String> = cols.to_vec();
    out[1] = pos.to_string();
    out[3] = new_ref;
    out[4] = new_alt;
    Some(out)
}

/// Canonical left-alignment loop (Tan/Abecasis/Kang 2015; bcftools/vt `norm`).
fn left_shift(r: &mut Vec<u8>, alts: &mut [Vec<u8>], pos: &mut u32, chr: &str, fa: &Fasta) -> Option<()> {
    loop {
        let last = *r.last().unwrap();
        if !alts.iter().all(|a| *a.last().unwrap() == last) {
            break;
        }
        if r.len() == 1 || alts.iter().any(|a| a.len() == 1) {
            if *pos <= 1 { break; }
            *pos -= 1;
            let prev = fa.base(chr, *pos)?.to_ascii_uppercase();
            r.insert(0, prev);
            for a in alts.iter_mut() { a.insert(0, prev); }
        } else {
            r.pop();
            for a in alts.iter_mut() { a.pop(); }
        }
    }
    Some(())
}

fn left_align_symbolic(
    cols: &[String],
    chr: &str,
    orig_pos: u32,
    alt_raw: &[&str],
    fa: &Fasta,
) -> Option<Vec<String>> {
    // Only symbolic deletions are realigned (bcftools leaves <DUP> etc. alone).
    if !alt_raw.iter().any(|a| a.eq_ignore_ascii_case("<DEL>") || a.starts_with("<DEL:")) {
        return None;
    }
    // Without INFO/END the span is the anchor base alone (rlen 1), as in htslib.
    let has_end = info_get_u32(&cols[7], "END").is_some();
    let end = info_get_u32(&cols[7], "END").unwrap_or(orig_pos);
    if end < orig_pos || end - orig_pos > MAX_SV_SPAN {
        return None;
    }
    let seq_idx: Vec<usize> = alt_raw
        .iter()
        .enumerate()
        .filter(|(_, a)| !is_symbolic_alt(a))
        .map(|(i, _)| i)
        .collect();

    let mut pos = orig_pos;
    let symbolic_only = seq_idx.is_empty();
    let (mut r, mut seq_alts): (Vec<u8>, Vec<Vec<u8>>) = if symbolic_only {
        let mut refseq = Vec::with_capacity((end - orig_pos + 1) as usize);
        for p in orig_pos..=end {
            refseq.push(fa.base(chr, p)?.to_ascii_uppercase());
        }
        let anchor = vec![refseq[0]];
        (refseq, vec![anchor])
    } else {
        let r = cols[3].as_bytes().to_ascii_uppercase();
        let alts = seq_idx
            .iter()
            .map(|&i| alt_raw[i].as_bytes().to_ascii_uppercase())
            .collect();
        (r, alts)
    };
    if r.is_empty() || seq_alts.iter().any(|a| a.is_empty()) {
        return None;
    }

    left_shift(&mut r, &mut seq_alts, &mut pos, chr, fa)?;

    let delta = orig_pos - pos;
    if delta == 0 {
        return None;
    }
    let new_end = end - delta;

    let new_ref = if symbolic_only {
        (r[0] as char).to_string()
    } else {
        String::from_utf8(r).ok()?
    };
    let mut seq_iter = seq_alts.into_iter();
    let new_alt = alt_raw
        .iter()
        .map(|a| {
            if is_symbolic_alt(a) {
                a.to_string()
            } else {
                seq_iter
                    .next()
                    .map(|v| String::from_utf8(v).unwrap_or_default())
                    .unwrap_or_default()
            }
        })
        .collect::<Vec<_>>()
        .join(",");

    let mut out: Vec<String> = cols.to_vec();
    out[1] = pos.to_string();
    out[3] = new_ref;
    out[4] = new_alt;
    if has_end {
        out[7] = info_set_u32(&cols[7], "END", new_end);
    }
    Some(out)
}

fn info_get_u32(info: &str, key: &str) -> Option<u32> {
    for f in info.split(';') {
        if let Some((k, v)) = f.split_once('=') {
            if k == key {
                return v.parse().ok();
            }
        }
    }
    None
}

fn info_set_u32(info: &str, key: &str, val: u32) -> String {
    info.split(';')
        .map(|f| match f.split_once('=') {
            Some((k, _)) if k == key => format!("{key}={val}"),
            _ => f.to_string(),
        })
        .collect::<Vec<_>>()
        .join(";")
}

#[cfg(test)]
#[path = "../../../tests/unit/cli_commands_norm.rs"]
mod tests;
