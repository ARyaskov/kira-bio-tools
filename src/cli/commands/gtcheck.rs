//! `bcftools gtcheck` port: DCv2 discordance with the same scoring, HWE
//! term, site pairing and report layout (cross-check and `-g` panel modes).

use crate::annotate::postproc::RegionFilter;
use crate::cli::args::GtcheckArgs;
use crate::filter::FilterEngine;
use crate::vcf::alleles::gt_alleles;
use crate::vcf::header::extract_samples;
use crate::vcf::variant_type::{VT_REF, record_type};
use crate::vcf::{UnifiedVcfReader, parse_vcf_line};
use anyhow::{Context, Result, bail};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Dosage bitmask: bit 0 hom-ref, bit 1 het, bit 2 hom-alt; 0 = missing.
type Dsg = u8;

/// C `printf("%.6e")`-style: signed, zero-padded two-digit exponent.
pub(crate) fn sci6(x: f64) -> String {
    let s = format!("{:.6e}", x);
    match s.find('e') {
        Some(epos) => {
            let (mant, exp) = s.split_at(epos);
            let exp = &exp[1..];
            let (sign, digits) = match exp.strip_prefix('-') {
                Some(d) => ("-", d),
                None => ("+", exp.strip_prefix('+').unwrap_or(exp)),
            };
            format!("{}e{}{:0>2}", mant, sign, digits)
        }
        None => s,
    }
}

/// htslib `hts_lrand48`: the 48-bit LCG behind the distinctive-sites shuffle.
struct Lrand48(u64);

impl Lrand48 {
    fn new(seed: u32) -> Self {
        Self(((seed as u64) << 16) | 0x330E)
    }

    fn next(&mut self) -> u32 {
        self.0 = (self.0.wrapping_mul(0x5DEECE66D).wrapping_add(0xB)) & ((1u64 << 48) - 1);
        (self.0 >> 17) as u32
    }
}

/// One record with the fields the comparison needs.
struct Rec {
    chrom: String,
    pos: u32,
    refa: String,
    alt: String,
    n_allele: usize,
    is_ref_only: bool,
    format: Vec<String>,
    samples: Vec<String>,
    info: String,
    line: String,
}

impl Rec {
    fn parse(line: String) -> Option<Self> {
        let mut it = line.split('\t');
        let chrom = it.next()?.to_string();
        let pos: u32 = it.next()?.parse().ok()?;
        it.next()?;
        let refa = it.next()?.to_string();
        let alt = it.next()?.to_string();
        it.next()?;
        it.next()?;
        let info = it.next()?.to_string();
        let format: Vec<String> = it.next().map(|f| f.split(':').map(|s| s.to_string()).collect()).unwrap_or_default();
        let samples: Vec<String> = it.map(|s| s.to_string()).collect();
        let n_allele = if alt == "." || alt.is_empty() { 1 } else { alt.split(',').count() + 1 };
        let is_ref_only = record_type(&refa, &alt) == VT_REF;
        Some(Self { chrom, pos, refa, alt, n_allele, is_ref_only, format, samples, info, line })
    }

    fn fmt_idx(&self, key: &str) -> Option<usize> {
        self.format.iter().position(|k| k == key)
    }

    /// Diploid genotypes (haploid samples are padded and read as missing);
    /// `None` when GT is absent, `Some(None)` when the ploidy is not two.
    fn genotypes(&self) -> Option<Option<Vec<Option<[usize; 2]>>>> {
        let gi = self.fmt_idx("GT")?;
        let mut max_ploidy = 0usize;
        let mut out = Vec::with_capacity(self.samples.len());
        for s in &self.samples {
            let gt = s.split(':').nth(gi).unwrap_or(".");
            let al = gt_alleles(gt);
            max_ploidy = max_ploidy.max(al.len());
            out.push(match (al.first().copied().flatten(), al.get(1).copied().flatten()) {
                (Some(a), Some(b)) if al.len() == 2 => Some([a, b]),
                _ => None,
            });
        }
        if max_ploidy != 2 {
            return Some(None);
        }
        Some(Some(out))
    }

    /// Three PL values per sample; `None` when PL is absent, `Some(None)` when
    /// the vectors are not diploid.
    fn pls(&self) -> Option<Option<Vec<Option<[i32; 3]>>>> {
        let pi = self.fmt_idx("PL")?;
        let mut max_len = 0usize;
        let mut out = Vec::with_capacity(self.samples.len());
        for s in &self.samples {
            let v = s.split(':').nth(pi).unwrap_or(".");
            let vals: Vec<Option<i32>> = v.split(',').map(|x| x.parse::<i32>().ok()).collect();
            max_len = max_len.max(if v == "." { 1 } else { vals.len() });
            out.push(match (vals.first().copied().flatten(), vals.get(1).copied().flatten(), vals.get(2).copied().flatten()) {
                (Some(a), Some(b), Some(c)) => Some([a, b, c]),
                _ => None,
            });
        }
        if max_len != 3 {
            return Some(None);
        }
        Some(Some(out))
    }

    /// `bcf_calc_ac`: INFO AN/AC when both are present, else all genotypes.
    fn allele_counts(&self) -> Option<(u32, u32)> {
        let mut an: Option<u32> = None;
        let mut ac: Option<u32> = None;
        for kv in self.info.split(';') {
            if let Some(v) = kv.strip_prefix("AN=") {
                an = v.parse().ok();
            } else if let Some(v) = kv.strip_prefix("AC=") {
                ac = v.split(',').next().and_then(|x| x.parse().ok());
            }
        }
        if let (Some(an), Some(ac)) = (an, ac) {
            return Some((an.saturating_sub(ac), ac));
        }
        let gi = self.fmt_idx("GT")?;
        let (mut n0, mut n1) = (0u32, 0u32);
        for s in &self.samples {
            for a in gt_alleles(s.split(':').nth(gi).unwrap_or(".")).into_iter().flatten() {
                if a == 0 { n0 += 1 } else { n1 += 1 }
            }
        }
        Some((n0, n1))
    }
}

/// Synced-reader pairing: the same ALT alleles in any order.
fn same_alt_set(a: &str, b: &str) -> bool {
    if a == b { return true; }
    let mut x: Vec<&str> = a.split(',').collect();
    let mut y: Vec<&str> = b.split(',').collect();
    x.sort_unstable();
    y.sort_unstable();
    x == y
}

fn gt_to_dsg(g: Option<[usize; 2]>) -> Dsg {
    match g {
        Some([a, b]) => 1 << ((a > 0) as u8 + (b > 0) as u8),
        None => 0,
    }
}

fn pl_to_dsg(p: Option<[i32; 3]>) -> Dsg {
    let Some(p) = p else { return 0 };
    let min = p[0].min(p[1]).min(p[2]);
    let mut d = 0;
    if p[0] == min { d |= 1 }
    if p[1] == min { d |= 2 }
    if p[2] == min { d |= 4 }
    d
}

/// `-log P(genotype | dosage)` under the single-allele error model.
fn gt_to_prob(dsg: Dsg, neg_log_e: f64) -> [f64; 3] {
    match dsg {
        1 => [0.0, neg_log_e, 2.0 * neg_log_e],
        2 => [neg_log_e, 0.0, neg_log_e],
        4 => [2.0 * neg_log_e, neg_log_e, 0.0],
        _ => [f64::INFINITY; 3],
    }
}

fn pl_to_prob(p: [i32; 3]) -> [f64; 3] {
    let f = |x: i32| 10f64.powf(-0.1 * (x.clamp(0, 255) as f64));
    let mut v = [f(p[0]), f(p[1]), f(p[2])];
    let sum: f64 = v.iter().sum();
    for x in v.iter_mut() {
        *x = -(*x / sum).ln();
    }
    v
}

#[derive(Clone, Copy, Default)]
struct PairAcc {
    ncnt: u32,
    ndiff: u32,
    pdiff: f64,
    hwe: f64,
    nmatch: u32,
}

struct Side {
    samples: Vec<String>,
    /// Selected sample indices, sorted.
    sel: Vec<usize>,
    use_gt: bool,
    include: Option<FilterEngine>,
    exclude: Option<FilterEngine>,
}

struct Counters {
    ncmp: u32,
    no_match: u32,
    not_ba: u32,
    mono: u32,
    no_data: u32,
    dip_gt: u32,
    dip_pl: u32,
    filter: u32,
    used: [[u32; 2]; 2],
}

struct Gtcheck {
    use_pls: u32,
    neg_log_e: f64,
    hom_only: bool,
    calc_hwe: bool,
    pairs: Option<Vec<(usize, usize)>>,
    cross: bool,
    acc: Vec<PairAcc>,
    cnt: Counters,
    /// Distinctive sites: (ndiff, chrom, pos, rand, differing pair indices).
    ds: Option<Vec<(u32, String, u32, u32, Vec<usize>)>>,
    rng: Lrand48,
}

/// Split `[qry:|gt:]VALUE` specs into the query and panel parts.
/// `-g` panel streamed alongside the query. Both inputs are position-sorted
/// (contig order from the query header, then first appearance), so only the
/// records at the current site are buffered; earlier records are unmatched.
struct PanelStream {
    reader: UnifiedVcfReader,
    peek: Option<(usize, Rec)>,
    buf: Vec<(Rec, bool)>,
    key: Option<(usize, u32)>,
    unmatched: u32,
}

impl PanelStream {
    fn new(reader: UnifiedVcfReader) -> Self {
        Self { reader, peek: None, buf: Vec::new(), key: None, unmatched: 0 }
    }

    fn fill(&mut self, contigs: &mut crate::vcf::ContigDict, keep: &dyn Fn(&str) -> bool) -> Result<()> {
        while self.peek.is_none() {
            let Some(line) = self.reader.read_line()? else { return Ok(()) };
            if line.is_empty() || line.as_bytes()[0] == b'#' || !keep(&line) { continue; }
            if let Some(rec) = Rec::parse(line) {
                let rank = contigs.insert(&rec.chrom) as usize;
                self.peek = Some((rank, rec));
            }
        }
        Ok(())
    }

    /// Panel records at `key` (contig rank, position).
    fn at(&mut self, key: (usize, u32), contigs: &mut crate::vcf::ContigDict, keep: &dyn Fn(&str) -> bool) -> Result<&mut Vec<(Rec, bool)>> {
        if self.key != Some(key) {
            self.unmatched += self.buf.iter().filter(|(_, used)| !used).count() as u32;
            self.buf.clear();
            self.key = Some(key);
            loop {
                self.fill(contigs, keep)?;
                let Some((rank, rec)) = self.peek.take() else { break };
                let k = (rank, rec.pos);
                if k < key {
                    self.unmatched += 1;
                    continue;
                }
                if k == key {
                    self.buf.push((rec, false));
                    continue;
                }
                self.peek = Some((rank, rec));
                break;
            }
        }
        Ok(&mut self.buf)
    }

    /// Records never paired with a query site, including the rest of the file.
    fn finish(&mut self, contigs: &mut crate::vcf::ContigDict, keep: &dyn Fn(&str) -> bool) -> Result<u32> {
        let mut n = self.unmatched + self.buf.iter().filter(|(_, used)| !used).count() as u32;
        self.buf.clear();
        self.unmatched = 0;
        if self.peek.take().is_some() { n += 1; }
        loop {
            self.fill(contigs, keep)?;
            if self.peek.take().is_none() { break; }
            n += 1;
        }
        Ok(n)
    }
}

fn split_side(specs: &[String]) -> (Vec<String>, Vec<String>) {
    let mut q = Vec::new();
    let mut g = Vec::new();
    for s in specs {
        if let Some(v) = s.strip_prefix("gt:") {
            g.push(v.to_string());
        } else if let Some(v) = s.strip_prefix("qry:") {
            q.push(v.to_string());
        } else {
            q.push(s.clone());
        }
    }
    (q, g)
}

fn init_samples(list: Option<&str>, file: Option<&str>, all: &[String], fname: &Path) -> Result<Vec<usize>> {
    let mut names: Vec<String> = Vec::new();
    if let Some(l) = list {
        if l == "-" {
            return Ok((0..all.len()).collect());
        }
        names.extend(l.split(',').map(|s| s.trim().to_string()).filter(|s| !s.is_empty()));
    }
    if let Some(f) = file {
        for line in BufReader::new(std::fs::File::open(f).with_context(|| format!("open {f}"))?).lines() {
            let l = line?;
            let t = l.trim();
            if !t.is_empty() && !t.starts_with('#') {
                names.push(t.split_whitespace().next().unwrap_or("").to_string());
            }
        }
    }
    if list.is_none() && file.is_none() {
        return Ok((0..all.len()).collect());
    }
    let mut idx: Vec<usize> = Vec::new();
    for n in &names {
        match all.iter().position(|s| s == n) {
            Some(i) => idx.push(i),
            None => bail!("No such sample in {}: [{}]", fname.display(), n),
        }
    }
    idx.sort_unstable();
    if idx.windows(2).any(|w| w[0] == w[1]) {
        bail!("Error: a sample is listed twice in the sample list");
    }
    Ok(idx)
}

fn make_filter(headers: &[String], specs: &[String], exclude: bool) -> Result<Option<FilterEngine>> {
    match specs.first() {
        Some(e) => Ok(Some(FilterEngine::new(headers, Some(e), exclude).context("filter expression")?)),
        None => Ok(None),
    }
}

fn passes(side: &Side, rec: &Rec) -> bool {
    if side.include.is_none() && side.exclude.is_none() {
        return true;
    }
    let Some(v) = parse_vcf_line(&rec.line) else { return true };
    if let Some(e) = &side.include {
        if !e.eval(&v).map(|r| r.pass_site).unwrap_or(true) {
            return false;
        }
    }
    if let Some(e) = &side.exclude {
        if e.eval(&v).map(|r| r.pass_site).unwrap_or(false) {
            return false;
        }
    }
    true
}

/// Dosage (and -log probabilities) for every sample of a record, following
/// `set_data`: fall back to the other tag when the requested one is absent.
fn site_data(rec: &Rec, use_gt: bool, use_pls: u32, neg_log_e: f64, cnt: &mut Counters) -> Option<(bool, Vec<(Dsg, [f64; 3])>)> {
    let mut want_gt = use_gt;
    for attempt in 0..2 {
        if want_gt {
            match rec.genotypes() {
                None => {
                    if attempt == 0 { want_gt = false; continue; }
                    cnt.no_data += 1;
                    return None;
                }
                Some(None) => { cnt.dip_gt += 1; return None; }
                Some(Some(g)) => {
                    let v = g.into_iter().map(|x| {
                        let d = gt_to_dsg(x);
                        (d, if use_pls > 0 { gt_to_prob(d, neg_log_e) } else { [0.0; 3] })
                    }).collect();
                    return Some((true, v));
                }
            }
        } else {
            match rec.pls() {
                None => {
                    if attempt == 0 { want_gt = true; continue; }
                    cnt.no_data += 1;
                    return None;
                }
                Some(None) => { cnt.dip_pl += 1; return None; }
                Some(Some(p)) => {
                    let v = p.into_iter().map(|x| {
                        let d = pl_to_dsg(x);
                        (d, match (x, use_pls > 0) { (Some(pl), true) => pl_to_prob(pl), _ => [0.0; 3] })
                    }).collect();
                    return Some((false, v));
                }
            }
        }
    }
    None
}

impl Gtcheck {
    fn score(&mut self, idx: usize, q: &(Dsg, [f64; 3]), g: &(Dsg, [f64; 3]), hwe_dsg: &[f64; 8], diff: &mut Vec<usize>) {
        let a = &mut self.acc[idx];
        let m = q.0 & g.0;
        if self.use_pls > 0 {
            let mut min = q.1[0] + g.1[0];
            min = min.min(q.1[1] + g.1[1]);
            min = min.min(q.1[2] + g.1[2]);
            a.pdiff += min;
            if self.calc_hwe {
                a.hwe += hwe_dsg[m as usize];
                a.nmatch += 1;
            }
        } else if m == 0 {
            a.ndiff += 1;
            diff.push(idx);
        } else if self.calc_hwe {
            a.hwe += hwe_dsg[m as usize];
            a.nmatch += 1;
        }
        a.ncnt += 1;
    }

    fn process(&mut self, qry: &Rec, gt: Option<&Rec>, qside: &Side, gside: &Side) {
        let Some((q_is_gt, qdata)) = site_data(qry, qside.use_gt, self.use_pls, self.neg_log_e, &mut self.cnt) else { return };
        let (g_is_gt, gdata) = match gt {
            Some(g) => match site_data(g, gside.use_gt, self.use_pls, self.neg_log_e, &mut self.cnt) {
                Some(d) => d,
                None => return,
            },
            None => (q_is_gt, qdata.clone()),
        };
        self.cnt.ncmp += 1;
        self.cnt.used[q_is_gt as usize][g_is_gt as usize] += 1;

        let mut hwe_dsg = [0.0f64; 8];
        if self.calc_hwe {
            let (n0, n1) = gt.unwrap_or(qry).allele_counts().unwrap_or((0, 0));
            let af = if n0 + n1 > 0 { n1 as f64 / (n0 + n1) as f64 } else { 1e-6 };
            let hwe = [-((1.0 - af) * (1.0 - af)).ln(), -(2.0 * af * (1.0 - af)).ln(), -(af * af).ln()];
            for (i, slot) in hwe_dsg.iter_mut().enumerate().skip(1) {
                let mut v = f64::INFINITY;
                for (k, h) in hwe.iter().enumerate() {
                    if (i >> k) & 1 == 1 && v > *h { v = *h; }
                }
                *slot = v;
            }
        }

        let mut diff: Vec<usize> = Vec::new();
        if let Some(pairs) = self.pairs.clone() {
            for (i, (iq, ig)) in pairs.iter().enumerate() {
                let g = gdata[*ig];
                if g.0 == 0 { continue; }
                if self.hom_only && g.0 & 5 == 0 { continue; }
                let q = qdata[*iq];
                if q.0 == 0 { continue; }
                self.score(i, &q, &g, &hwe_dsg, &mut diff);
            }
            if let Some(ds) = self.ds.as_mut() {
                if !diff.is_empty() {
                    let r = self.rng.next();
                    ds.push((diff.len() as u32, qry.chrom.clone(), qry.pos, r, diff));
                }
            }
            return;
        }
        let qsel = &qside.sel;
        let gsel = &gside.sel;
        let mut idx = 0usize;
        for (i, &iq) in qsel.iter().enumerate() {
            let ngt = if self.cross { i } else { gsel.len() };
            let q = qdata[iq];
            if q.0 == 0 {
                idx += ngt;
                continue;
            }
            for &ig in gsel.iter().take(ngt) {
                let mut g = gdata[ig];
                if self.hom_only && g.0 & 5 == 0 { g.0 = 0; }
                if g.0 == 0 {
                    idx += 1;
                    continue;
                }
                self.score(idx, &q, &g, &hwe_dsg, &mut diff);
                idx += 1;
            }
        }
    }
}

pub fn cmd_gtcheck(args: GtcheckArgs) -> Result<()> {
    let qry_path = args.input.clone().ok_or_else(|| anyhow::anyhow!("missing query VCF/BCF"))?;

    // The old `-e INT` form still selects the error probability.
    let mut use_pls: u32 = args.error_probability.unwrap_or(40);
    let mut excludes: Vec<String> = Vec::new();
    for e in &args.exclude {
        match e.parse::<u32>() {
            Ok(v) if !e.contains(':') => {
                if args.error_probability.is_none() { use_pls = v; }
                eprintln!("[warning] auto-detected the old format --error-probability option, please switch from -e to -E.");
            }
            _ => excludes.push(e.clone()),
        }
    }
    let neg_log_e = -(10f64.powf(-0.1 * use_pls as f64)).ln();

    let (use_q, use_g) = {
        let toks: Vec<&str> = args.use_tag.split(',').map(|s| s.trim()).collect();
        let parse = |t: &str| -> Result<Option<bool>> {
            match t.to_ascii_uppercase().as_str() {
                "GT" => Ok(Some(true)),
                "PL" => Ok(Some(false)),
                _ => bail!("Failed to parse --use {}; only GT and PL are supported", args.use_tag),
            }
        };
        (parse(toks[0])?, if toks.len() > 1 { parse(toks[1])? } else { None })
    };

    let region = if let Some(s) = &args.regions {
        Some(RegionFilter::from_cli(s)?)
    } else if let Some(p) = &args.regions_file {
        Some(RegionFilter::from_file(p)?)
    } else {
        None
    };
    let target = if let Some(s) = &args.targets {
        Some(RegionFilter::from_cli(s)?)
    } else if let Some(p) = &args.targets_file {
        Some(RegionFilter::from_file(p)?)
    } else {
        None
    };
    let keep = |line: &str| -> bool {
        if let Some(r) = &region { if !r.line_passes_mode(line, args.regions_overlap) { return false; } }
        if let Some(t) = &target { if !t.line_passes_mode(line, args.targets_overlap) { return false; } }
        true
    };

    let (inc_q, inc_g) = split_side(&args.include);
    let (exc_q, exc_g) = split_side(&excludes);
    let (smp_q, smp_g) = split_side(&args.samples);
    let (sfile_q, sfile_g) = split_side(&args.samples_file);

    // Query header and sides.
    let mut qreader = UnifiedVcfReader::open(&qry_path).with_context(|| format!("open {}", qry_path.display()))?;
    let qheaders = qreader.header()?;
    let qsamples = extract_samples(&qheaders);
    if qsamples.is_empty() { bail!("No samples in {}?", qry_path.display()); }
    let has_tag = |h: &[String], tag: &str| h.iter().any(|l| l.starts_with(&format!("##FORMAT=<ID={tag},")));
    let q_use_gt = match use_q {
        Some(v) => v,
        None => {
            if has_tag(&qheaders, "PL") { false } else if has_tag(&qheaders, "GT") { true } else { bail!("Neither PL nor GT tag is present in the header of {}", qry_path.display()) }
        }
    };
    let mut qside = Side {
        sel: init_samples(smp_q.first().map(String::as_str), sfile_q.first().map(String::as_str), &qsamples, &qry_path)?,
        samples: qsamples,
        use_gt: q_use_gt,
        include: make_filter(&qheaders, &inc_q, false)?,
        exclude: make_filter(&qheaders, &exc_q, true)?,
    };

    // The -g panel is streamed in lock-step with the query: both inputs are
    // position-sorted, so only the records at the current site are held.
    let mut contigs = crate::vcf::ContigDict::from_header_lines(qheaders.iter().map(String::as_str));
    let mut panel: Option<PanelStream> = None;
    let mut gside: Option<Side> = None;
    if let Some(gp) = &args.genotypes {
        let r = UnifiedVcfReader::open(gp).with_context(|| format!("open {}", gp.display()))?;
        let gheaders = r.header()?;
        let gsamples = extract_samples(&gheaders);
        if gsamples.is_empty() { bail!("No samples in {}?", gp.display()); }
        let g_use_gt = match use_g {
            Some(v) => v,
            None => {
                if has_tag(&gheaders, "GT") { true } else if has_tag(&gheaders, "PL") { false } else { bail!("Neither PL nor GT tag is present in the header of {}", gp.display()) }
            }
        };
        gside = Some(Side {
            sel: init_samples(smp_g.first().map(String::as_str), sfile_g.first().map(String::as_str), &gsamples, gp)?,
            samples: gsamples,
            use_gt: g_use_gt,
            include: make_filter(&gheaders, &inc_g, false)?,
            exclude: make_filter(&gheaders, &exc_g, true)?,
        });
        panel = Some(PanelStream::new(r));
    } else if !smp_g.is_empty() || !sfile_g.is_empty() {
        // `gt:` lists without -g select the second side of the same file.
        let sel = init_samples(smp_g.first().map(String::as_str), sfile_g.first().map(String::as_str), &qside.samples, &qry_path)?;
        gside = Some(Side { samples: qside.samples.clone(), sel, use_gt: q_use_gt, include: None, exclude: None });
    }

    // Explicit pairs (-p/-P), sorted by sample indices like bcftools.
    let pairs: Option<Vec<(usize, usize)>> = {
        let gnames: &[String] = gside.as_ref().map(|g| g.samples.as_slice()).unwrap_or(&qside.samples);
        let find = |names: &[String], n: &str, what: &Path| -> Result<usize> {
            names.iter().position(|s| s == n).ok_or_else(|| anyhow::anyhow!("No such sample in {}: [{}]", what.display(), n))
        };
        let gpath: &Path = args.genotypes.as_deref().unwrap_or(&qry_path);
        let mut v: Vec<(usize, usize)> = Vec::new();
        if let Some(p) = &args.pairs {
            let toks: Vec<&str> = p.split(',').map(|t| t.trim()).filter(|t| !t.is_empty()).collect();
            if toks.len() % 2 != 0 { bail!("Expected even number of comma-delimited samples with -p"); }
            for c in toks.chunks(2) {
                v.push((find(&qside.samples, c[0], &qry_path)?, find(gnames, c[1], gpath)?));
            }
        }
        if let Some(f) = &args.pairs_file {
            for line in BufReader::new(std::fs::File::open(f)?).lines() {
                let l = line?;
                let t = l.trim();
                if t.is_empty() || t.starts_with('#') { continue; }
                let mut it = t.split(|c: char| c == '\t' || c == ',' || c == ' ').filter(|x| !x.is_empty());
                let (Some(a), Some(b)) = (it.next(), it.next()) else { bail!("Could not parse {}: {}", f.display(), t) };
                v.push((find(&qside.samples, a, &qry_path)?, find(gnames, b, gpath)?));
            }
        }
        if v.is_empty() { None } else { v.sort_unstable(); Some(v) }
    };

    let cross = gside.is_none();
    if cross {
        gside = Some(Side { samples: qside.samples.clone(), sel: qside.sel.clone(), use_gt: q_use_gt, include: None, exclude: None });
    }
    let gside = gside.unwrap();
    if pairs.is_some() && !args.genotypes.is_some() && (!smp_g.is_empty() || !sfile_g.is_empty()) {
        // pairs win over sample lists
    }
    let npairs = match &pairs {
        Some(p) => p.len(),
        None if cross => qside.sel.len() * (qside.sel.len() + 1) / 2,
        None => qside.sel.len() * gside.sel.len(),
    };
    let ds_spec: Option<f64> = args.distinctive_sites.as_deref().map(|s| s.split(',').next().unwrap_or("").parse::<f64>()).transpose().context("--distinctive-sites")?;
    let ds_min: Option<usize> = ds_spec.map(|d| {
        let n = if d <= 1.0 { (npairs as f64 * d) as usize } else { d as usize };
        n.min(npairs).max(1)
    });

    let mut gc = Gtcheck {
        use_pls,
        neg_log_e,
        hom_only: args.homs_only,
        calc_hwe: !args.no_hwe_prob,
        pairs: pairs.clone(),
        cross,
        acc: vec![PairAcc::default(); npairs.max(1)],
        cnt: Counters { ncmp: 0, no_match: 0, not_ba: 0, mono: 0, no_data: 0, dip_gt: 0, dip_pl: 0, filter: 0, used: [[0; 2]; 2] },
        ds: if ds_spec.is_some() && pairs.is_some() { Some(Vec::new()) } else { None },
        rng: Lrand48::new(0),
    };

    // Walk the query; -g sites are looked up by position and paired by alleles.
    while let Some(line) = qreader.read_line()? {
        if line.is_empty() || line.as_bytes()[0] == b'#' || !keep(&line) { continue; }
        let Some(q) = Rec::parse(line) else { continue };
        if let Some(ps) = panel.as_mut() {
            let rank = contigs.insert(&q.chrom) as usize;
            let cands = ps.at((rank, q.pos), &mut contigs, &keep)?;
            let Some(k) = cands.iter().position(|(g, used)| !used && g.refa == q.refa && same_alt_set(&g.alt, &q.alt)) else { gc.cnt.no_match += 1; continue };
            cands[k].1 = true;
            let g = &cands[k].0;
            if q.n_allele > 2 || g.n_allele > 2 { gc.cnt.not_ba += 1; continue; }
            if q.is_ref_only || g.is_ref_only { gc.cnt.mono += 1; continue; }
            if !passes(&qside, &q) || !passes(&gside, g) { gc.cnt.filter += 1; continue; }
            if args.dry_run { break; }
            gc.process(&q, Some(g), &qside, &gside);
        } else {
            if q.n_allele > 2 { gc.cnt.not_ba += 1; continue; }
            if q.is_ref_only { gc.cnt.mono += 1; continue; }
            if !passes(&qside, &q) { gc.cnt.filter += 1; continue; }
            if args.dry_run { break; }
            gc.process(&q, None, &qside, &gside);
        }
    }
    if let Some(ps) = panel.as_mut() {
        gc.cnt.no_match += ps.finish(&mut contigs, &keep)?;
    }
    qside.sel.shrink_to_fit();

    // Report.
    let mut out: Box<dyn Write> = match &args.output {
        Some(p) if p.as_os_str() != "-" => {
            if args.output_type.as_deref().is_some_and(|t| t.starts_with('z')) {
                Box::new(crate::bgzf::BgzfWriter::create(p)?)
            } else {
                Box::new(BufWriter::with_capacity(64 * 1024, std::fs::File::create(p)?))
            }
        }
        _ => Box::new(BufWriter::with_capacity(64 * 1024, std::io::stdout())),
    };
    let argv: Vec<String> = std::env::args().skip(2).filter(|a| a != "--").collect();
    writeln!(out, "# This file was produced by kira-bt gtcheck ({}), the command line was:", env!("CARGO_PKG_VERSION"))?;
    writeln!(out, "#\tkira-bt gtcheck {}", argv.join(" "))?;
    writeln!(out, "#")?;
    let c = &gc.cnt;
    writeln!(out, "INFO\tsites-compared\t{}", c.ncmp)?;
    writeln!(out, "INFO\tsites-skipped-no-match\t{}", c.no_match)?;
    writeln!(out, "INFO\tsites-skipped-multiallelic\t{}", c.not_ba)?;
    writeln!(out, "INFO\tsites-skipped-monoallelic\t{}", c.mono)?;
    writeln!(out, "INFO\tsites-skipped-no-data\t{}", c.no_data)?;
    writeln!(out, "INFO\tsites-skipped-GT-not-diploid\t{}", c.dip_gt)?;
    writeln!(out, "INFO\tsites-skipped-PL-not-diploid\t{}", c.dip_pl)?;
    writeln!(out, "INFO\tsites-skipped-filtering-expression\t{}", c.filter)?;
    writeln!(out, "INFO\tsites-used-PL-vs-PL\t{}", c.used[0][0])?;
    writeln!(out, "INFO\tsites-used-PL-vs-GT\t{}", c.used[0][1])?;
    writeln!(out, "INFO\tsites-used-GT-vs-PL\t{}", c.used[1][0])?;
    writeln!(out, "INFO\tsites-used-GT-vs-GT\t{}", c.used[1][1])?;
    writeln!(out, "# DCv2, discordance version 2:")?;
    writeln!(out, "#     - Query sample")?;
    writeln!(out, "#     - Genotyped sample")?;
    writeln!(out, "#     - Discordance, given either as an abstract score or number of mismatches, see the options -E/-u")?;
    writeln!(out, "#       in man page for details. Note that samples with high missingness have fewer sites compared,")?;
    writeln!(out, "#       which results in lower overall discordance. Therefore it is advisable to use the average score")?;
    writeln!(out, "#       per site rather than the absolute value, i.e. divide the value by the number of sites compared")?;
    writeln!(out, "#       (smaller value = better match)")?;
    writeln!(out, "#     - Average negative log of HWE probability at matching sites, attempts to quantify the following")?;
    writeln!(out, "#       intuition: rare genotype matches are more informative than common genotype matches, hence two")?;
    writeln!(out, "#       samples with similar discordance can be further stratified by the HWE score (bigger value = better")?;
    writeln!(out, "#       match, the observed concordance was less likely to occur by chance)")?;
    writeln!(out, "#     - Number of sites compared for this pair of samples (bigger = more informative)")?;
    writeln!(out, "#     - Number of matching genotypes")?;
    writeln!(out, "#DCv2\t[2]Query Sample\t[3]Genotyped Sample\t[4]Discordance\t[5]Average -log P(HWE)\t[6]Number of sites compared\t[6]Number of matching genotypes")?;

    let row = |out: &mut dyn Write, qn: &str, gn: &str, a: &PairAcc| -> Result<()> {
        let hwe = if gc.calc_hwe && a.nmatch > 0 { a.hwe / a.nmatch as f64 } else { 0.0 };
        if gc.use_pls == 0 {
            writeln!(out, "DCv2\t{}\t{}\t{}\t{}\t{}\t{}", qn, gn, a.ndiff, sci6(hwe), a.ncnt, a.nmatch)?;
        } else {
            writeln!(out, "DCv2\t{}\t{}\t{}\t{}\t{}\t{}", qn, gn, sci6(a.pdiff), sci6(hwe), a.ncnt, a.nmatch)?;
        }
        Ok(())
    };
    let sort_by_hwe = args.n_matches < 0;
    let ntop = args.n_matches.unsigned_abs() as usize;
    let val = |a: &PairAcc| -> f64 {
        if sort_by_hwe {
            if a.nmatch > 0 { -a.hwe / a.nmatch as f64 } else { 0.0 }
        } else if gc.use_pls == 0 {
            if a.ncnt > 0 { a.ndiff as f64 / a.ncnt as f64 } else { 0.0 }
        } else if a.ncnt > 0 {
            a.pdiff / a.ncnt as f64
        } else {
            0.0
        }
    };
    let mut trim = ntop;
    if pairs.is_none() && gside.sel.len() <= ntop { trim = 0; }

    if let Some(p) = &pairs {
        for (i, (iq, ig)) in p.iter().enumerate() {
            row(&mut out, &qside.samples[*iq], &gside.samples[*ig], &gc.acc[i])?;
        }
    } else if trim == 0 {
        let mut idx = 0usize;
        for (i, &iq) in qside.sel.iter().enumerate() {
            let ngt = if cross { i } else { gside.sel.len() };
            for &ig in gside.sel.iter().take(ngt) {
                row(&mut out, &qside.samples[iq], &gside.samples[ig], &gc.acc[idx])?;
                idx += 1;
            }
        }
    } else if !cross {
        let ngt = gside.sel.len();
        for (i, &iq) in qside.sel.iter().enumerate() {
            let mut arr: Vec<(f64, usize, usize)> = (0..ngt).map(|j| (val(&gc.acc[i * ngt + j]), j, i * ngt + j)).collect();
            arr.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            for (_, j, idx) in arr.iter().take(ntop) {
                row(&mut out, &qside.samples[iq], &gside.samples[gside.sel[*j]], &gc.acc[*idx])?;
            }
        }
    } else {
        let n = qside.sel.len();
        for i in 0..n {
            let mut arr: Vec<(f64, usize, usize)> = Vec::with_capacity(n.saturating_sub(1));
            for j in 0..i {
                let idx = i * (i.saturating_sub(1)) / 2 + j;
                arr.push((val(&gc.acc[idx]), j, idx));
            }
            for j in i..n.saturating_sub(1) {
                let idx = j * (j + 1) / 2 + i;
                arr.push((val(&gc.acc[idx]), j + 1, idx));
            }
            arr.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            for (_, ism, idx) in arr.iter().take(ntop) {
                if i <= *ism { continue; }
                row(&mut out, &qside.samples[qside.sel[i]], &qside.samples[qside.sel[*ism]], &gc.acc[*idx])?;
            }
        }
    }

    if let (Some(mut sites), Some(ds_min)) = (gc.ds.take(), ds_min) {
        sites.sort_by(|a, b| b.0.cmp(&a.0).then(a.3.cmp(&b.3)));
        writeln!(out, "# DS, distinctive sites:")?;
        writeln!(out, "#     - chromosome")?;
        writeln!(out, "#     - position")?;
        writeln!(out, "#     - cumulative number of pairs distinguished by this block")?;
        writeln!(out, "#     - block id")?;
        writeln!(out, "#DS\t[2]Chromosome\t[3]Position\t[4]Cumulative number of distinct pairs\t[5]Block id")?;
        let mut blk = vec![false; npairs];
        let (mut tot, mut iblock) = (0usize, 0usize);
        for (_, chrom, pos, _, diff) in &sites {
            let mut new = 0usize;
            for &p in diff {
                if !blk[p] { blk[p] = true; new += 1; }
            }
            if new == 0 { continue; }
            tot += new;
            writeln!(out, "DS\t{}\t{}\t{}\t{}", chrom, pos, tot, iblock)?;
            if tot < ds_min { continue; }
            iblock += 1;
            tot = 0;
            blk.iter_mut().for_each(|b| *b = false);
        }
    }
    out.flush()?;
    Ok(())
}

#[cfg(test)]
#[path = "../../../tests/unit/cli_commands_gtcheck.rs"]
mod tests;
