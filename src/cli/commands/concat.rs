use crate::annotate::postproc::{RegionFilter, version_header_line};
use crate::bgzf::{BGZF_EOF, is_bgzf};
use crate::cli::args::ConcatArgs;
use crate::vcf::alleles::gt_alleles;
use crate::vcf::header::{ContigDict, extract_samples};
use crate::vcf::sink::{OutputKind, parse_output_type};
use crate::vcf::variant_type::{VT_INDEL, VT_SNP, record_type};
use crate::vcf::{UnifiedVcfReader, VcfSink};
use anyhow::{Context, Result, bail};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RmDups { None, Exact, Snps, Indels, Both, All }

impl RmDups {
    fn parse(s: Option<&str>, flag_d: bool) -> Result<Self> {
        if let Some(s) = s {
            return Ok(match s {
                "none" => Self::None,
                "exact" => Self::Exact,
                "snps" => Self::Snps,
                "indels" => Self::Indels,
                "both" => Self::Both,
                "all" | "any" => Self::All,
                o => bail!("--rm-dups: unknown {o:?} (none|snps|indels|both|all|exact)"),
            });
        }
        Ok(if flag_d { Self::Exact } else { Self::None })
    }
}

pub fn cmd_concat(args: ConcatArgs) -> Result<()> {
    let mut inputs: Vec<PathBuf> = args.inputs.clone();
    if let Some(fl) = &args.file_list {
        for line in BufReader::new(File::open(fl)?).lines() {
            let l = line?;
            let t = l.trim();
            if t.is_empty() || t.starts_with('#') { continue; }
            inputs.push(PathBuf::from(t));
        }
    }
    if inputs.is_empty() { bail!("concat: no inputs"); }

    let kind = args.output_type.as_deref().map(parse_output_type).transpose()?.unwrap_or(OutputKind::Vcf);
    let rm_dups = RmDups::parse(args.rm_dups.as_deref(), args.remove_duplicates)?;
    let region = if let Some(s) = &args.regions {
        Some(RegionFilter::from_cli(s)?)
    } else if let Some(p) = &args.regions_file {
        Some(RegionFilter::from_file(p)?)
    } else {
        None
    };

    if args.naive || args.naive_force {
        return concat_naive(&inputs, &args, kind);
    }

    // Header: the first file's lines plus new meta lines from the others
    // (bcftools merges the headers); sample sets must match.
    let (headers, samples) = {
        let r = UnifiedVcfReader::open(&inputs[0]).context("open first input")?;
        let h = r.header()?;
        let s = extract_samples(&h);
        (h, s)
    };
    let mut meta: Vec<String> = headers.iter().filter(|h| h.starts_with("##")).cloned().collect();
    let chrom_line = headers.iter().find(|h| h.starts_with("#CHROM")).cloned().unwrap_or_default();
    for p in &inputs[1..] {
        let r = UnifiedVcfReader::open(p).with_context(|| format!("open {:?}", p))?;
        let h = r.header()?;
        if extract_samples(&h) != samples {
            bail!("concat: sample names differ between {} and {}", inputs[0].display(), p.display());
        }
        for l in h.iter().filter(|l| l.starts_with("##") && !l.starts_with("##fileformat")) {
            if !meta.contains(l) { meta.push(l.clone()); }
        }
    }
    let ligate = args.ligate || args.ligate_force || args.ligate_warn;
    if ligate {
        if !meta.iter().any(|l| l.starts_with("##FORMAT=<ID=PQ,")) {
            meta.push("##FORMAT=<ID=PQ,Number=1,Type=Integer,Description=\"Phasing Quality (bigger is better)\">".into());
        }
        if !meta.iter().any(|l| l.starts_with("##FORMAT=<ID=PS,")) {
            meta.push("##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Set\">".into());
        }
    }
    let mut out_headers: Vec<String> = Vec::with_capacity(meta.len() + 2);
    for h in &meta {
        if args.drop_genotypes && h.starts_with("##FORMAT=") { continue; }
        out_headers.push(h.clone());
    }
    if !args.no_version { out_headers.push(version_header_line()); }
    out_headers.push(if args.drop_genotypes { trim_chrom_samples(&chrom_line) } else { chrom_line.clone() });
    let mut sink = VcfSink::open(args.output.as_deref(), kind, &out_headers)?;
    sink.write_header(&out_headers)?;

    let ctx = Ctx { drop_genotypes: args.drop_genotypes, region: region.as_ref(), regions_overlap: args.regions_overlap, rm_dups };

    if ligate {
        concat_ligate(&inputs, &args, &samples, &out_headers, &ctx, &mut sink)?;
    } else if args.allow_overlaps {
        concat_sorted(&inputs, &ctx, &mut sink)?;
    } else {
        concat_streaming(&inputs, &out_headers, &ctx, &mut sink)?;
    }
    sink.finish()?;

    if let (Some(kind_s), Some(out)) = (args.write_index.as_deref(), args.output.as_deref()) {
        if matches!(kind, OutputKind::VcfGz(_) | OutputKind::Bcf(_)) && out != Path::new("-") {
            let (ik, ext) = if kind_s == "tbi" { (crate::csi::IndexKind::Tbi, "tbi") } else { (crate::csi::IndexKind::Csi, "csi") };
            let idx = PathBuf::from(format!("{}.{}", out.display(), ext));
            crate::csi::build_index(out, &idx, ik, None).with_context(|| format!("-W: write {}", idx.display()))?;
        }
    }
    Ok(())
}

struct Ctx<'a> {
    drop_genotypes: bool,
    region: Option<&'a RegionFilter>,
    regions_overlap: u8,
    rm_dups: RmDups,
}

/// Duplicate detection across consecutive records at one position.
#[derive(Default)]
struct DupWindow {
    key: (String, String),
    recs: Vec<(String, String, u32)>,
}

impl DupWindow {
    fn check(&mut self, line: &str, mode: RmDups) -> bool {
        if mode == RmDups::None { return false; }
        let mut it = line.splitn(6, '\t');
        let chrom = it.next().unwrap_or("").to_string();
        let pos = it.next().unwrap_or("").to_string();
        let _id = it.next();
        let refa = it.next().unwrap_or("").to_string();
        let alt = it.next().unwrap_or("").to_string();
        let vt = record_type(&refa, &alt);
        if (chrom.as_str(), pos.as_str()) != (self.key.0.as_str(), self.key.1.as_str()) {
            self.key = (chrom, pos);
            self.recs.clear();
        }
        let dup = self.recs.iter().any(|(r, a, pvt)| same_site(mode, r, a, *pvt, &refa, &alt, vt));
        if !dup {
            self.recs.push((refa, alt, vt));
        }
        dup
    }
}

/// Synced-reader pairing rule for a `-d` mode (`Exact`/`None` = identical alleles).
fn same_site(mode: RmDups, r1: &str, a1: &str, t1: u32, r2: &str, a2: &str, t2: u32) -> bool {
    match mode {
        RmDups::All => true,
        RmDups::Snps => (t1 & VT_SNP != 0 && t2 & VT_SNP != 0) || (r1.eq_ignore_ascii_case(r2) && a1.eq_ignore_ascii_case(a2)),
        RmDups::Indels => (t1 & VT_INDEL != 0 && t2 & VT_INDEL != 0) || (r1.eq_ignore_ascii_case(r2) && a1.eq_ignore_ascii_case(a2)),
        RmDups::Both => {
            (t1 & VT_SNP != 0 && t2 & VT_SNP != 0)
                || (t1 & VT_INDEL != 0 && t2 & VT_INDEL != 0)
                || (r1.eq_ignore_ascii_case(r2) && a1.eq_ignore_ascii_case(a2))
        }
        RmDups::Exact | RmDups::None => r1.eq_ignore_ascii_case(r2) && a1.eq_ignore_ascii_case(a2),
    }
}

fn record_key(line: &str) -> Option<(String, u32)> {
    let mut it = line.splitn(3, '\t');
    let chrom = it.next()?.to_string();
    let pos: u32 = it.next()?.parse().ok()?;
    Some((chrom, pos))
}

/// CHROM, POS, REF, ALT of a data line.
fn line_site(line: &str) -> Option<(String, u32, String, String)> {
    let mut it = line.splitn(6, '\t');
    let chrom = it.next()?.to_string();
    let pos: u32 = it.next()?.parse().ok()?;
    it.next()?;
    let refa = it.next()?.to_string();
    let alt = it.next()?.to_string();
    Some((chrom, pos, refa, alt))
}

fn concat_streaming(inputs: &[PathBuf], headers: &[String], ctx: &Ctx<'_>, sink: &mut VcfSink) -> Result<()> {
    let mut contigs = ContigDict::from_header_lines(headers.iter().map(String::as_str));
    let mut dups = DupWindow::default();
    let mut last: Option<(usize, u32)> = None;
    for p in inputs {
        let mut r = UnifiedVcfReader::open(p).with_context(|| format!("open {:?}", p))?;
        let _ = r.header()?;
        let mut first_in_file = true;
        while let Some(line) = r.read_line()? {
            if line.is_empty() || line.as_bytes()[0] == b'#' { continue; }
            if let Some(rf) = ctx.region {
                if !rf.line_passes_mode(&line, ctx.regions_overlap) { continue; }
            }
            if first_in_file {
                first_in_file = false;
                if let (Some((c, pos)), Some((lr, lp))) = (record_key(&line), last) {
                    let rank = contigs.insert(&c) as usize;
                    if (rank, pos) < (lr, lp) {
                        bail!(
                            "concat: {} starts before the end of the previous file ({}:{} after contig #{lr}:{lp}); use -a to sort overlapping inputs",
                            p.display(), c, pos
                        );
                    }
                }
            }
            if let Some((c, pos)) = record_key(&line) {
                last = Some((contigs.insert(&c) as usize, pos));
            }
            if dups.check(&line, ctx.rm_dups) { continue; }
            let out_line = if ctx.drop_genotypes { drop_genotypes_line(&line) } else { line };
            sink.write_line(&out_line)?;
        }
    }
    Ok(())
}

/// Chromosomes of a file in order of first appearance: from the index names
/// when the index carries them, otherwise from a scan of the records.
fn chrom_order(p: &Path) -> Result<Vec<String>> {
    let is_bcf = is_bcf_file(p).unwrap_or(false);
    let index = crate::csi::find_index_for(p).and_then(|ip| crate::csi::BinIndex::load(&ip).ok());
    if !is_bcf {
        if let Some(idx) = &index {
            let names: Vec<String> = idx.names().to_vec();
            if !names.is_empty() { return Ok(names); }
        }
    }
    let mut r = UnifiedVcfReader::open(p).with_context(|| format!("open {:?}", p))?;
    let headers = r.header()?;
    if is_bcf {
        // BCF sequences are numbered by the header, so bcftools walks them in header order.
        let contigs = ContigDict::from_header_lines(headers.iter().map(String::as_str));
        let names: Vec<String> = contigs.names().to_vec();
        if let Some(idx) = &index {
            return Ok(names.into_iter().enumerate().filter(|(rid, _)| idx.n_records(*rid).unwrap_or(0) > 0).map(|(_, n)| n).collect());
        }
        let mut present: std::collections::HashSet<String> = std::collections::HashSet::new();
        while let Some(line) = r.read_line()? {
            if line.is_empty() || line.as_bytes()[0] == b'#' { continue; }
            present.insert(line.split('\t').next().unwrap_or("").to_string());
        }
        return Ok(names.into_iter().filter(|n| present.contains(n)).collect());
    }
    let mut order: Vec<String> = Vec::new();
    while let Some(line) = r.read_line()? {
        if line.is_empty() || line.as_bytes()[0] == b'#' { continue; }
        let chrom = line.split('\t').next().unwrap_or("");
        if order.last().map(String::as_str) != Some(chrom) && !order.iter().any(|c| c == chrom) {
            order.push(chrom.to_string());
        }
    }
    Ok(order)
}

/// `-a`: synced walk of overlapping inputs. Chromosomes follow their first
/// appearance across the inputs; at one position records with identical
/// alleles (or the `-d` class) are grouped, each group in input order.
fn concat_sorted(inputs: &[PathBuf], ctx: &Ctx<'_>, sink: &mut VcfSink) -> Result<()> {
    let mut order: Vec<String> = Vec::new();
    for p in inputs {
        for c in chrom_order(p)? {
            if !order.contains(&c) { order.push(c); }
        }
    }
    // Inputs are not sorted by the same chromosome order, so each chromosome
    // is fetched from every input (index query, or a filtered scan).
    let fetch = |p: &Path, chrom: &str| -> Result<Vec<(u32, String, String, u32, String)>> {
        let mut out = Vec::new();
        let mut push = |line: &str| {
            if let Some(rf) = ctx.region {
                if !rf.line_passes_mode(line, ctx.regions_overlap) { return; }
            }
            if let Some((c, pos, refa, alt)) = line_site(line) {
                if c == chrom {
                    let vt = record_type(&refa, &alt);
                    out.push((pos, refa, alt, vt, line.to_string()));
                }
            }
        };
        if crate::csi::find_index_for(p).is_some() {
            if let Ok(mut ir) = crate::csi::IndexedVcfReader::open(p) {
                ir.query(chrom, 1, u32::MAX, |line| { push(line); Ok(true) })?;
                return Ok(out);
            }
        }
        let mut r = UnifiedVcfReader::open(p).with_context(|| format!("open {:?}", p))?;
        let _ = r.header()?;
        while let Some(line) = r.read_line()? {
            if line.is_empty() || line.as_bytes()[0] == b'#' { continue; }
            push(&line);
        }
        Ok(out)
    };

    for chrom in &order {
        let mut lists: Vec<std::collections::VecDeque<(u32, String, String, u32, String)>> = Vec::with_capacity(inputs.len());
        for p in inputs {
            lists.push(fetch(p, chrom)?.into());
        }
        loop {
            let Some(pos) = lists.iter().filter_map(|l| l.front().map(|r| r.0)).min() else { break };
            // All records at this position, per input.
            let mut recs: Vec<Vec<(String, String, u32, String)>> = vec![Vec::new(); lists.len()];
            for (i, l) in lists.iter_mut().enumerate() {
                while l.front().is_some_and(|r| r.0 == pos) {
                    let (_, refa, alt, vt, line) = l.pop_front().unwrap();
                    recs[i].push((refa, alt, vt, line));
                }
            }
            let mut used: Vec<Vec<bool>> = recs.iter().map(|v| vec![false; v.len()]).collect();
            for i in 0..recs.len() {
                for j in 0..recs[i].len() {
                    if used[i][j] { continue; }
                    used[i][j] = true;
                    let (r1, a1, t1, _) = &recs[i][j];
                    let mut set: Vec<(usize, usize)> = vec![(i, j)];
                    for k in i + 1..recs.len() {
                        if let Some(m) = (0..recs[k].len()).find(|&m| !used[k][m] && same_site(ctx.rm_dups, r1, a1, *t1, &recs[k][m].0, &recs[k][m].1, recs[k][m].2)) {
                            used[k][m] = true;
                            set.push((k, m));
                        }
                    }
                    for (idx, (a, b)) in set.iter().enumerate() {
                        if ctx.rm_dups != RmDups::None && idx > 0 { break; }
                        let line = &recs[*a][*b].3;
                        let out_line = if ctx.drop_genotypes { drop_genotypes_line(line) } else { line.clone() };
                        sink.write_line(&out_line)?;
                    }
                }
            }
        }
    }
    Ok(())
}

/// One buffered data line of a ligation input.
#[derive(Clone)]
struct LLine {
    chrom: String,
    pos: u32,
    refa: String,
    alt: String,
    text: String,
}

struct LReader {
    reader: UnifiedVcfReader,
    peek: Option<LLine>,
    done: bool,
}

impl LReader {
    fn open(p: &Path, ctx: &Ctx<'_>) -> Result<Self> {
        let reader = UnifiedVcfReader::open(p).with_context(|| format!("open {:?}", p))?;
        let _ = reader.header()?;
        let mut r = Self { reader, peek: None, done: false };
        r.fill(ctx)?;
        Ok(r)
    }

    fn fill(&mut self, ctx: &Ctx<'_>) -> Result<()> {
        if self.peek.is_some() || self.done { return Ok(()); }
        while let Some(line) = self.reader.read_line()? {
            if line.is_empty() || line.as_bytes()[0] == b'#' { continue; }
            if let Some(rf) = ctx.region {
                if !rf.line_passes_mode(&line, ctx.regions_overlap) { continue; }
            }
            let Some((chrom, pos, refa, alt)) = line_site(&line) else { continue };
            self.peek = Some(LLine { chrom, pos, refa, alt, text: line });
            return Ok(());
        }
        self.done = true;
        Ok(())
    }

    /// Skip records before `pos` on `chrom` (an index seek in bcftools).
    fn seek(&mut self, chrom: &str, pos: u32, ctx: &Ctx<'_>) -> Result<()> {
        loop {
            self.fill(ctx)?;
            match &self.peek {
                Some(l) if l.chrom == chrom && l.pos < pos => { self.peek = None; }
                _ => return Ok(()),
            }
        }
    }
}

/// Phase bookkeeping of `bcftools concat -l`.
struct Ligator<'a> {
    n: usize,
    swap: Vec<bool>,
    nswap: usize,
    phase_set: Vec<Option<u32>>,
    phase_set_changed: bool,
    phase_qual: Vec<i64>,
    nmatch: Vec<u32>,
    nmism: Vec<u32>,
    prev_chr: Option<String>,
    seen: Vec<String>,
    /// Overlap buffer: (record of the earlier file, record of the later file).
    buf: Vec<(Option<LLine>, Option<LLine>)>,
    min_pq: i64,
    compact_ps: bool,
    drop_genotypes: bool,
    sink: &'a mut VcfSink,
}

impl<'a> Ligator<'a> {
    fn write(&mut self, l: &LLine, with_pq: bool) -> Result<()> {
        let mut text = if self.nswap > 0 { phase_update(&l.text, &self.swap) } else { l.text.clone() };
        if with_pq {
            let vals: Vec<String> = self.phase_qual.iter().map(|q| q.to_string()).collect();
            text = add_format_field(&text, "PQ", &vals);
        }
        if !self.compact_ps || self.phase_set_changed {
            let vals: Vec<String> = self.phase_set.iter().map(|p| p.map(|v| v.to_string()).unwrap_or_else(|| ".".into())).collect();
            text = add_format_field(&text, "PS", &vals);
            self.phase_set_changed = false;
        }
        let out = if self.drop_genotypes { drop_genotypes_line(&text) } else { text };
        self.sink.write_line(&out)
    }

    fn push(&mut self, a: Option<LLine>, b: Option<LLine>, is_overlap: bool) -> Result<()> {
        let first = a.as_ref().or(b.as_ref()).expect("push needs a record");
        if self.prev_chr.as_deref() != Some(first.chrom.as_str()) {
            if self.prev_chr.is_some() { self.flush()?; }
            let ps = first.pos;
            for p in self.phase_set.iter_mut() { *p = Some(ps); }
            self.phase_set_changed = true;
            if self.seen.contains(&first.chrom) {
                bail!("The chromosome block {} is not contiguous", first.chrom);
            }
            self.seen.push(first.chrom.clone());
            self.prev_chr = Some(first.chrom.clone());
        }
        if !is_overlap {
            let l = a.or(b).expect("non-overlap record");
            return self.write(&l, false);
        }
        self.buf.push((a, b));
        Ok(())
    }

    fn flush(&mut self) -> Result<()> {
        if self.buf.is_empty() { return Ok(()); }
        let n = self.n;
        // Phase agreement over the overlap, relative to the current swap state.
        for (a, b) in &self.buf {
            let (Some(a), Some(b)) = (a, b) else { continue };
            let (Some(ga), Some(gb)) = (diploid_gts(&a.text, n), diploid_gts(&b.text, n)) else { continue };
            for j in 0..n {
                let (Some((a0, a1, pa)), Some((b0, b1, pb))) = (ga[j], gb[j]) else { continue };
                if !pa || !pb || a0 == a1 || b0 == b1 { continue; }
                if a0 == b0 && a1 == b1 {
                    if self.swap[j] { self.nmism[j] += 1 } else { self.nmatch[j] += 1 }
                }
                if a0 == b1 && a1 == b0 {
                    if self.swap[j] { self.nmatch[j] += 1 } else { self.nmism[j] += 1 }
                }
            }
        }
        let buf = std::mem::take(&mut self.buf);
        let npairs = buf.len();
        let half = npairs.div_ceil(2);
        // First half from the earlier file under the current swap state.
        for (a, b) in &buf[..half] {
            match (a, b) {
                (Some(a), _) => self.write(a, false)?,
                (None, Some(b)) => {
                    // Only the earlier file's records are re-phased here.
                    let saved = self.nswap;
                    self.nswap = 0;
                    self.write(b, false)?;
                    self.nswap = saved;
                }
                (None, None) => {}
            }
        }
        // New swap state and phasing quality per sample.
        self.nswap = 0;
        for j in 0..n {
            if self.nmatch[j] >= self.nmism[j] {
                self.swap[j] = false;
            } else {
                self.swap[j] = true;
                self.nswap += 1;
            }
            self.phase_qual[j] = if self.nmatch[j] > 0 && self.nmism[j] > 0 {
                let f = self.nmatch[j] as f64 / (self.nmatch[j] + self.nmism[j]) as f64;
                (99.0 * (0.7 + f * f.ln() + (1.0 - f) * (1.0 - f).ln()) / 0.7) as i64
            } else {
                99
            };
            self.nmatch[j] = 0;
            self.nmism[j] = 0;
        }
        // Second half from the later file; the first shared site carries PQ.
        let mut pq_printed = false;
        for (a, b) in &buf[half..] {
            let rec = match (a, b) {
                (_, Some(b)) => b,
                (Some(a), None) => a,
                (None, None) => continue,
            };
            let with_pq = !pq_printed && a.is_some() && b.is_some();
            if with_pq {
                pq_printed = true;
                for j in 0..n {
                    if self.phase_qual[j] < self.min_pq {
                        self.phase_set[j] = Some(rec.pos);
                        self.phase_set_changed = true;
                    } else if self.compact_ps {
                        self.phase_set[j] = None;
                    }
                }
            }
            self.write(rec, with_pq)?;
        }
        Ok(())
    }
}

/// Diploid genotype per sample as (allele, allele, phased); `None` for
/// missing or non-diploid entries, or when GT is absent / ploidy is not two.
fn diploid_gts(line: &str, n: usize) -> Option<Vec<Option<(usize, usize, bool)>>> {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 10 { return None; }
    let gi = cols[8].split(':').position(|k| k == "GT")?;
    let mut out = Vec::with_capacity(n);
    let mut max_ploidy = 0;
    for s in cols[9..].iter().take(n) {
        let gt = s.split(':').nth(gi).unwrap_or(".");
        let al = gt_alleles(gt);
        max_ploidy = max_ploidy.max(al.len());
        out.push(match (al.first().copied().flatten(), al.get(1).copied().flatten()) {
            (Some(a), Some(b)) if al.len() == 2 => Some((a, b, gt.contains('|'))),
            _ => None,
        });
    }
    while out.len() < n { out.push(None); }
    if max_ploidy != 2 { return None; }
    Some(out)
}

/// Swap the haplotypes of phased diploid genotypes for the flagged samples.
fn phase_update(line: &str, swap: &[bool]) -> String {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 10 { return line.to_string(); }
    let fmt: Vec<&str> = cols[8].split(':').collect();
    let Some(gi) = fmt.iter().position(|k| *k == "GT") else { return line.to_string() };
    let mut out: Vec<String> = cols[..9].iter().map(|s| s.to_string()).collect();
    for (i, samp) in cols[9..].iter().enumerate() {
        if !swap.get(i).copied().unwrap_or(false) {
            out.push(samp.to_string());
            continue;
        }
        let mut parts: Vec<String> = samp.split(':').map(|s| s.to_string()).collect();
        if let Some(gt) = parts.get(gi) {
            let al = gt_alleles(gt);
            if al.len() == 2 && gt.contains('|') && al[0].is_some() {
                let a = al[0].map(|v| v.to_string()).unwrap_or_else(|| ".".into());
                let b = al[1].map(|v| v.to_string()).unwrap_or_else(|| ".".into());
                parts[gi] = format!("{b}|{a}");
            }
        }
        out.push(parts.join(":"));
    }
    out.join("\t")
}

/// Set a FORMAT field on every sample (appending it when absent); shorter
/// sample columns are padded with `.` first.
fn add_format_field(line: &str, key: &str, vals: &[String]) -> String {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 9 { return line.to_string(); }
    let mut fmt: Vec<&str> = if cols[8].is_empty() || cols[8] == "." { Vec::new() } else { cols[8].split(':').collect() };
    let idx = match fmt.iter().position(|k| *k == key) {
        Some(i) => i,
        None => { fmt.push(key); fmt.len() - 1 }
    };
    let mut out: Vec<String> = cols[..8].iter().map(|s| s.to_string()).collect();
    out.push(fmt.join(":"));
    for (i, samp) in cols[9..].iter().enumerate() {
        let mut parts: Vec<String> = if samp.is_empty() { Vec::new() } else { samp.split(':').map(|s| s.to_string()).collect() };
        while parts.len() < fmt.len() { parts.push(".".into()); }
        parts[idx] = vals.get(i).cloned().unwrap_or_else(|| ".".into());
        out.push(parts.join(":"));
    }
    out.join("\t")
}

/// `-l`: ligate phased chunks the way bcftools does: two files are open at a
/// time, the overlap is buffered, the phase is compared per sample, the first
/// half of the overlap comes from the earlier file and the second half from
/// the later one (with PQ at the switch and PS phase sets).
fn concat_ligate(inputs_all: &[PathBuf], args: &ConcatArgs, samples: &[String], headers: &[String], ctx: &Ctx<'_>, sink: &mut VcfSink) -> Result<()> {
    let n = samples.len();
    let contigs = ContigDict::from_header_lines(headers.iter().map(String::as_str));
    // First record of every non-empty file; `None` marks a file starting a new chromosome.
    let mut inputs: Vec<PathBuf> = Vec::with_capacity(inputs_all.len());
    let mut start_pos: Vec<Option<u32>> = Vec::with_capacity(inputs_all.len());
    let mut prev_chr: Option<String> = None;
    for p in inputs_all {
        let r = LReader::open(p, ctx)?;
        if let Some(l) = &r.peek {
            start_pos.push(if prev_chr.as_deref() == Some(l.chrom.as_str()) { Some(l.pos) } else { None });
            prev_chr = Some(l.chrom.clone());
            inputs.push(p.clone());
        }
    }
    let inputs = &inputs[..];

    let mut lig = Ligator {
        n,
        swap: vec![false; n],
        nswap: 0,
        phase_set: vec![None; n],
        phase_set_changed: false,
        phase_qual: vec![99; n],
        nmatch: vec![0; n],
        nmism: vec![0; n],
        prev_chr: None,
        seen: Vec::new(),
        buf: Vec::new(),
        min_pq: args.min_pq as i64,
        compact_ps: args.compact_ps,
        drop_genotypes: ctx.drop_genotypes,
        sink,
    };
    let mut warned = false;
    let mut readers: Vec<LReader> = Vec::new();
    let mut ifname = 0usize;
    let rank = |c: &str| contigs.id(c).map(|v| v as usize).unwrap_or(usize::MAX);

    while ifname < inputs.len() {
        while readers.len() < 2 && ifname < inputs.len() {
            readers.push(LReader::open(&inputs[ifname], ctx)?);
            ifname += 1;
            if start_pos[ifname - 1].is_none() { break; }
            if ifname < inputs.len() && start_pos[ifname].is_none() { break; }
        }
        // Continue from the position the previous round stopped at.
        let mut seek: Option<(String, u32)> = None;
        if let Some(l) = readers[0].peek.clone() {
            for r in readers.iter_mut() { r.seek(&l.chrom, l.pos, ctx)?; }
            seek = Some((l.chrom, l.pos));
        }

        loop {
            for r in readers.iter_mut() { r.fill(ctx)?; }
            // Next position across the open readers.
            let Some((chrom, pos)) = readers
                .iter()
                .filter_map(|r| r.peek.as_ref().map(|l| (rank(&l.chrom), l.pos, l.chrom.clone())))
                .min()
                .map(|(_, p, c)| (c, p))
            else { break };
            if let Some((sc, sp)) = &seek {
                if *sc == chrom && *sp > pos {
                    // Records starting before the seek position are skipped (bcftools seek semantics).
                    for r in readers.iter_mut() {
                        if r.peek.as_ref().is_some_and(|l| l.chrom == chrom && l.pos == pos) { r.peek = None; }
                    }
                    continue;
                }
            }
            seek = None;
            // Open the next file once its first position is reached; it joins this position.
            while ifname < inputs.len() && start_pos[ifname].is_some_and(|s| pos >= s) {
                let mut r = LReader::open(&inputs[ifname], ctx)?;
                r.seek(&chrom, pos, ctx)?;
                readers.push(r);
                ifname += 1;
            }
            // Records of each reader at this position, paired by identical alleles.
            let mut here: Vec<Vec<LLine>> = Vec::new();
            for r in readers.iter_mut() {
                let mut v = Vec::new();
                loop {
                    r.fill(ctx)?;
                    match &r.peek {
                        Some(l) if l.chrom == chrom && l.pos == pos => v.push(r.peek.take().unwrap()),
                        _ => break,
                    }
                }
                here.push(v);
            }
            let mut sets: Vec<(Option<LLine>, Option<LLine>)> = Vec::new();
            if here.len() == 1 {
                for l in here[0].drain(..) { sets.push((Some(l), None)); }
            } else {
                let second: Vec<LLine> = std::mem::take(&mut here[1]);
                let first: Vec<LLine> = std::mem::take(&mut here[0]);
                let mut used = vec![false; second.len()];
                for a in first {
                    match (0..second.len()).find(|&m| !used[m] && second[m].refa == a.refa && second[m].alt == a.alt) {
                        Some(m) => { used[m] = true; sets.push((Some(a), Some(second[m].clone()))); }
                        None => sets.push((Some(a), None)),
                    }
                }
                for (m, b) in second.into_iter().enumerate() {
                    if !used[m] { sets.push((None, Some(b))); }
                }
            }

            let mut sets: std::collections::VecDeque<(Option<LLine>, Option<LLine>)> = sets.into();
            while let Some((a, b)) = sets.pop_front() {
                let mut is_overlap = readers.len() >= 2;
                if a.is_none() {
                    if readers[0].done && readers[0].peek.is_none() {
                        lig.flush()?;
                        readers.remove(0);
                        is_overlap = false;
                        // The remaining records of this position now belong to the first reader.
                        for s in sets.iter_mut() {
                            if s.0.is_none() { s.0 = s.1.take(); }
                        }
                    } else if args.ligate_warn {
                        if !warned {
                            let l = b.as_ref().unwrap();
                            eprintln!("Warning: Dropping the site {}:{}. The --ligate option is intended for VCFs with perfect\n         overlap, sites in overlapping regions present in one but missing in other are dropped.\n         This warning is printed only once.", l.chrom, l.pos);
                            warned = true;
                        }
                        continue;
                    } else if !args.ligate_force {
                        let l = b.as_ref().unwrap();
                        bail!("Error: The --ligate option is intended for VCFs with perfect overlap, the site {}:{} breaks the assumption", l.chrom, l.pos);
                    }
                }
                let line = a.as_ref().or(b.as_ref()).unwrap().clone();
                if readers.len() > 1 && b.is_none() && !(readers[1].done && readers[1].peek.is_none()) && !args.ligate_force {
                    if args.ligate_warn {
                        if !warned {
                            eprintln!("Warning: Dropping the site {}:{}. The --ligate option is intended for VCFs with perfect\n         overlap, sites in overlapping regions present in one but missing in other are dropped.\n         This warning is printed only once.", line.chrom, line.pos);
                            warned = true;
                        }
                    } else {
                        bail!("Error: The --ligate option is intended for VCFs with perfect overlap, the site {}:{} breaks the assumption", line.chrom, line.pos);
                    }
                    continue;
                }
                let (a, b) = if readers.len() == 1 || !is_overlap { (a.or(b), None) } else { (a, b) };
                lig.push(a, b, is_overlap)?;
            }
        }
        if !readers.is_empty() {
            lig.flush()?;
            readers.clear();
        }
    }
    Ok(())
}

/// `--naive`: concatenate without decoding. BGZF inputs are copied block by
/// block (the block holding the header/data boundary is re-compressed);
/// anything else is streamed through.
fn concat_naive(inputs: &[PathBuf], args: &ConcatArgs, kind: OutputKind) -> Result<()> {
    let all_bgzf = inputs.iter().all(|p| is_bgzf(p).unwrap_or(false));
    let all_vcfgz = all_bgzf && inputs.iter().all(|p| !is_bcf_file(p).unwrap_or(false));
    // Block copy is possible when every input is BGZF VCF and the output is too.
    let can_copy = all_vcfgz && matches!(kind, OutputKind::VcfGz(_));

    // Header compatibility: sample names must match unless --naive-force.
    let mut first_samples: Option<Vec<String>> = None;
    for p in inputs {
        let r = UnifiedVcfReader::open(p).with_context(|| format!("open {:?}", p))?;
        let s = extract_samples(&r.header()?);
        match &first_samples {
            None => first_samples = Some(s),
            Some(f) if *f != s && !args.naive_force => bail!("concat --naive: sample names differ in {} (use --naive-force)", p.display()),
            _ => {}
        }
    }

    let out_path = args.output.as_deref();
    if can_copy {
        let mut out: Box<dyn Write> = match out_path {
            Some(p) if p != Path::new("-") => Box::new(std::io::BufWriter::with_capacity(1 << 20, File::create(p)?)),
            _ => Box::new(std::io::BufWriter::with_capacity(1 << 20, std::io::stdout())),
        };
        for (i, p) in inputs.iter().enumerate() {
            copy_bgzf_data(p, i == 0, &mut out)?;
        }
        out.write_all(&BGZF_EOF)?;
        out.flush()?;
        return Ok(());
    }

    // Fallback: decode and re-encode in file order without any checks.
    let mut first = UnifiedVcfReader::open(&inputs[0])?;
    let headers = first.header()?;
    let mut sink = VcfSink::open(out_path, kind, &headers)?;
    sink.write_header(&headers)?;
    while let Some(line) = first.read_line()? {
        if !line.starts_with('#') { sink.write_line(&line)?; }
    }
    for p in &inputs[1..] {
        let mut r = UnifiedVcfReader::open(p)?;
        let _ = r.header()?;
        while let Some(line) = r.read_line()? {
            if !line.starts_with('#') { sink.write_line(&line)?; }
        }
    }
    sink.finish()
}

fn is_bcf_file(p: &Path) -> Result<bool> {
    let mut r = noodles_bgzf::io::Reader::new(File::open(p)?);
    let mut magic = [0u8; 5];
    let n = r.read(&mut magic)?;
    Ok(n == 5 && magic == crate::bcf::BCF_MAGIC)
}

/// Copy the data blocks of a BGZF VCF. For the first file the header blocks
/// are included; for later files copying starts at the first data record,
/// re-compressing the partial block it lives in. Trailing EOF blocks are dropped.
fn copy_bgzf_data(path: &Path, with_header: bool, out: &mut dyn Write) -> Result<()> {
    let mut r = noodles_bgzf::io::Reader::new(File::open(path)?);
    let mut line = String::new();
    let first_data = loop {
        let vpos = r.virtual_position();
        line.clear();
        let n = r.read_line(&mut line)?;
        if n == 0 || !line.starts_with('#') { break vpos; }
    };
    let start = if with_header { noodles_bgzf::VirtualPosition::from(0) } else { first_data };
    let block_off = start.compressed();
    let in_block = start.uncompressed() as usize;

    let mut f = File::open(path)?;
    let len = f.metadata()?.len();
    // Strip a trailing EOF marker.
    let mut end = len;
    if len >= 28 {
        f.seek(SeekFrom::Start(len - 28))?;
        let mut tail = [0u8; 28];
        f.read_exact(&mut tail)?;
        if tail == BGZF_EOF { end = len - 28; }
    }
    let mut copy_from = block_off;
    if in_block > 0 {
        // Re-compress the remainder of the boundary block.
        r.seek(start)?;
        let rest: Vec<u8> = {
            let mut buf = Vec::new();
            loop {
                let cur = r.virtual_position();
                if cur.compressed() != block_off { break; }
                let avail = r.fill_buf()?;
                if avail.is_empty() { break; }
                let n = avail.len();
                buf.extend_from_slice(avail);
                r.consume(n);
            }
            buf
        };
        let mut w = noodles_bgzf::io::Writer::new(Vec::new());
        w.write_all(&rest)?;
        let bytes = w.finish()?;
        // Drop the EOF block the noodles writer appends.
        let cut = if bytes.ends_with(&BGZF_EOF) { bytes.len() - 28 } else { bytes.len() };
        out.write_all(&bytes[..cut])?;
        copy_from = r.virtual_position().compressed();
    }
    if copy_from < end {
        f.seek(SeekFrom::Start(copy_from))?;
        let mut remaining = end - copy_from;
        let mut buf = vec![0u8; 1 << 20];
        while remaining > 0 {
            let want = buf.len().min(remaining as usize);
            let n = f.read(&mut buf[..want])?;
            if n == 0 { break; }
            out.write_all(&buf[..n])?;
            remaining -= n as u64;
        }
    }
    Ok(())
}

fn trim_chrom_samples(h: &str) -> String {
    let cols: Vec<&str> = h.split('\t').collect();
    cols.iter().take(8).copied().collect::<Vec<_>>().join("\t")
}

fn drop_genotypes_line(line: &str) -> String {
    let cols: Vec<&str> = line.split('\t').collect();
    cols.iter().take(8).copied().collect::<Vec<_>>().join("\t")
}

#[cfg(test)]
#[path = "../../../tests/unit/cli_commands_concat.rs"]
mod tests;
