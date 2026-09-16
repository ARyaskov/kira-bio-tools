//! Runs of autozygosity: a port of bcftools `vcfroh.c` on the shared
//! [`crate::hmm`] kernel. Two states, HW (Hardy-Weinberg) and AZ
//! (autozygous); emissions come from genotype likelihoods (or genotypes with
//! `-G`) and the allele frequency; transitions can follow a genetic map or a
//! constant recombination rate.

use std::io::{BufRead, Write};
use std::path::Path;

use anyhow::{Result, bail};

use crate::af_file::AfFile;
use crate::hmm::{Hmm, Snapshot, TransitionFn};
use crate::util::{fmt_e, phred_score};
use crate::vcf::VcfRecord;

pub const STATE_HW: usize = 0;
pub const STATE_AZ: usize = 1;

/// IMPUTE2 genetic map of one chromosome: 0-based position and cumulative
/// distance scaled from cM to a probability.
pub struct GenMap {
    map: Vec<(i64, f64)>,
    igenmap: usize,
}

impl GenMap {
    /// `load_genmap`: `{CHROM}` in the file name is replaced by `chr`;
    /// `Ok(None)` when the file does not exist (the chromosome is skipped).
    pub fn load(fname: &str, chr: &str) -> Result<Option<Self>> {
        let fname = fname.replace("{CHROM}", chr);
        let path = Path::new(&fname);
        if !path.exists() {
            return Ok(None);
        }
        let mut rd = crate::util::open_text_maybe_gz(path)?;
        let mut line = String::new();
        rd.read_line(&mut line)?;
        let head = line.trim_end_matches(['\n', '\r']);
        if head != "position COMBINED_rate(cM/Mb) Genetic_Map(cM)" {
            bail!(
                "Unexpected header in {fname}, found:\n\t[{head}], but expected:\n\t[position COMBINED_rate(cM/Mb) Genetic_Map(cM)]"
            );
        }
        let mut map = Vec::new();
        loop {
            line.clear();
            if rd.read_line(&mut line)? == 0 {
                break;
            }
            let l = line.trim_end_matches(['\n', '\r']);
            if l.is_empty() {
                continue;
            }
            let mut cols = l.split_whitespace();
            let pos = cols.next().and_then(|s| s.parse::<i64>().ok());
            let _skip = cols.next();
            let rate = cols.next().and_then(|s| s.parse::<f64>().ok());
            let (Some(pos), Some(rate)) = (pos, rate) else {
                bail!("Could not parse {fname}: {l}")
            };
            map.push((pos - 1, rate * 0.01));
        }
        if map.is_empty() {
            bail!("Genetic map empty?");
        }
        Ok(Some(Self { map, igenmap: 0 }))
    }

    /// `get_genmap_rate`: recombination probability between two positions.
    pub fn rate(&mut self, start: i64, end: i64) -> f64 {
        let n = self.map.len();
        let mut i = self.igenmap;
        if self.map[i].0 > start {
            while i > 0 && self.map[i].0 > start {
                i -= 1;
            }
        } else {
            while i + 1 < n && self.map[i + 1].0 < start {
                i += 1;
            }
        }
        let mut j = i;
        while j + 1 < n && self.map[j].0 < end {
            j += 1;
        }
        if i == j {
            self.igenmap = i;
            return 0.0;
        }
        let start = start.max(self.map[i].0);
        let end = end.min(self.map[j].0);
        let rate = (self.map[j].1 - self.map[i].1) / (self.map[j].0 - self.map[i].0) as f64
            * (end - start) as f64;
        self.igenmap = j;
        rate
    }
}

/// Transition adjustment: none, constant rate (`-M`) or genetic map (`-m`).
pub enum RohTransition {
    None,
    RecRate(f64),
    GenMap { map: GenMap, rec_rate: f64 },
}

fn scale_tprob(tprob: &mut [f64], mut ci: f64) {
    if ci > 1.0 {
        ci = 1.0;
    }
    tprob[2 * STATE_HW + STATE_AZ] *= ci;
    tprob[2 * STATE_AZ + STATE_HW] *= ci;
    tprob[2 * STATE_AZ + STATE_AZ] = 1.0 - tprob[2 * STATE_HW + STATE_AZ];
    tprob[2 * STATE_HW + STATE_HW] = 1.0 - tprob[2 * STATE_AZ + STATE_HW];
}

impl TransitionFn for RohTransition {
    fn apply(&mut self, prev_pos: u32, pos: u32, tprob: &mut [f64]) {
        match self {
            Self::None => {}
            Self::RecRate(r) => scale_tprob(tprob, (pos.wrapping_sub(prev_pos)) as f64 * *r),
            Self::GenMap { map, rec_rate } => {
                let mut ci = map.rate(prev_pos as i64, pos as i64);
                if *rec_rate != 0.0 {
                    ci *= *rec_rate;
                }
                scale_tprob(tprob, ci);
            }
        }
    }
}

#[derive(Clone, Debug)]
pub struct RohOpts {
    /// P(AZ|HW) and P(HW|AZ).
    pub t2az: f64,
    pub t2hw: f64,
    pub rec_rate: f64,
    pub genmap_fname: Option<String>,
    pub af_tag: Option<String>,
    pub af_file: Option<String>,
    pub dflt_af: f64,
    /// `[GT,|PL,]<file|->`.
    pub estimate_af: Option<String>,
    /// `-G`: phred of the two unlikely genotypes.
    pub gts_only: Option<f64>,
    pub snps_only: bool,
    pub skip_homref: bool,
    pub include_noalt: bool,
    /// `-V` convergence threshold.
    pub vi_training: Option<f64>,
    /// `-b int[,int]`.
    pub buffer_size: Option<String>,
    pub output_st: bool,
    pub output_rg: bool,
    pub cmdline: String,
}

impl Default for RohOpts {
    fn default() -> Self {
        Self {
            t2az: 6.7e-8,
            t2hw: 5e-9,
            rec_rate: 0.0,
            genmap_fname: None,
            af_tag: None,
            af_file: None,
            dflt_af: 0.0,
            estimate_af: None,
            gts_only: None,
            snps_only: false,
            skip_homref: false,
            include_noalt: false,
            vi_training: None,
            buffer_size: None,
            output_st: true,
            output_rg: true,
            cmdline: String::new(),
        }
    }
}

#[derive(Default)]
struct Rg {
    beg: u32,
    end: u32,
    nqual: u32,
    qual: f64,
    rid: String,
    state: u8,
}

#[derive(Default)]
struct Smpl {
    /// `[HW, AZ]` per site.
    eprob: Vec<f64>,
    sites: Vec<u32>,
    igenmap: usize,
    nused: usize,
    /// Chromosome blocks kept for Viterbi training.
    rid: Vec<String>,
    rid_off: Vec<usize>,
    snapshot: Option<Snapshot>,
    rg: Rg,
}

/// Counters printed at the end (`Number of lines ...`).
#[derive(Default, Debug, Clone)]
pub struct RohStats {
    pub ntot: usize,
    pub nfiltered: usize,
    pub nno_af: usize,
    pub nno_alt: usize,
    pub nmultiallelic: usize,
    pub ndup: usize,
}

pub struct RohCaller {
    opts: RohOpts,
    samples: Vec<String>,
    /// Header indexes of the analysed samples.
    roh_smpl: Vec<usize>,
    /// Samples for `--estimate-AF`; `None` = all.
    af_smpl: Option<Vec<usize>>,
    af_from_pl: bool,
    estimate_af: bool,
    unseen_pl: f64,
    fake_pls: bool,
    af_file: Option<AfFile>,
    pl2p: [f64; 256],
    nbuf_max: usize,
    nbuf_olap: usize,
    smpl: Vec<Smpl>,
    hmm: Hmm<RohTransition>,
    prev_rid: Option<String>,
    prev_pos: u32,
    skip_rid: Option<String>,
    pub stats: RohStats,
    out: Box<dyn Write>,
}

/// One genotype allele as htslib sees it: `None` missing, `Some(-65)` the
/// vector-end padding of a haploid call in a diploid column.
fn parse_gt_pairs(rec: &VcfRecord) -> Option<Vec<[Option<i32>; 2]>> {
    let fmt = rec.format.as_deref()?;
    let igt = fmt.split(':').position(|k| k == "GT")?;
    let mut out = Vec::with_capacity(rec.samples.len());
    let mut maxn = 0;
    for s in &rec.samples {
        let f = s.split(':').nth(igt).unwrap_or(".");
        let mut v: Vec<Option<i32>> = Vec::with_capacity(2);
        for tok in f.split(['/', '|']) {
            v.push(if tok == "." || tok.is_empty() {
                None
            } else {
                tok.parse::<i32>().ok()
            });
        }
        maxn = maxn.max(v.len());
        out.push(v);
    }
    if maxn != 2 {
        return None;
    }
    Some(
        out.into_iter()
            .map(|v| {
                [
                    v.first().copied().unwrap_or(Some(-65)),
                    v.get(1).copied().unwrap_or(Some(-65)),
                ]
            })
            .collect(),
    )
}

/// PL values per sample; `-1` for `.`, vector-end (`-2`) for short vectors.
/// Returns `(values, n)` with `n` the maximum count.
fn parse_pls(rec: &VcfRecord) -> Option<(Vec<Vec<i64>>, usize)> {
    let fmt = rec.format.as_deref()?;
    let ipl = fmt.split(':').position(|k| k == "PL")?;
    let mut out = Vec::with_capacity(rec.samples.len());
    let mut n = 0;
    for s in &rec.samples {
        let f = s.split(':').nth(ipl).unwrap_or(".");
        let v: Vec<i64> = if f == "." || f.is_empty() {
            vec![-1]
        } else {
            f.split(',')
                .map(|t| {
                    if t == "." {
                        -1
                    } else {
                        t.parse::<i64>().unwrap_or(-1)
                    }
                })
                .collect()
        };
        n = n.max(v.len());
        out.push(v);
    }
    for v in &mut out {
        v.resize(n, -2);
    }
    Some((out, n))
}

fn info_value<'a>(info: &'a str, tag: &str) -> Option<&'a str> {
    info.split(';').find_map(|kv| {
        let (k, v) = kv.split_once('=')?;
        (k == tag).then_some(v)
    })
}

fn bcf_is_snp(alleles: &[String]) -> bool {
    alleles.iter().all(|a| {
        let b = a.as_bytes();
        (b.len() == 1 && b[0] != b'*') || a == "<X>" || a == "<*>"
    })
}

fn alleles2gt(a: usize, b: usize) -> usize {
    let (a, b) = if a > b { (b, a) } else { (a, b) };
    b * (b + 1) / 2 + a
}

/// Sample list of `-s/-S`/`--estimate-AF`: names or a file, `^` negates;
/// unknown names are reported and skipped (`SMPL_VERBOSE`).
pub fn sample_list(samples: &[String], spec: &str, is_file: bool) -> Result<Vec<usize>> {
    let (negate, body) = match spec.strip_prefix('^') {
        Some(r) => (true, r),
        None => (false, spec),
    };
    let names: Vec<String> = if is_file {
        std::fs::read_to_string(body)
            .map_err(|e| anyhow::anyhow!("Could not read {body}: {e}"))?
            .lines()
            .map(|l| l.split_whitespace().next().unwrap_or("").to_string())
            .filter(|s| !s.is_empty())
            .collect()
    } else {
        body.split(',').map(str::to_string).collect()
    };
    let mut keep = vec![false; samples.len()];
    for n in &names {
        match samples.iter().position(|s| s == n) {
            Some(i) => keep[i] = true,
            None => eprintln!("No such sample: \"{n}\""),
        }
    }
    Ok((0..samples.len()).filter(|&i| keep[i] != negate).collect())
}

impl RohCaller {
    /// `init_data`: `headers` are the VCF header lines, `samples` the header
    /// sample names, `roh_smpl` the analysed samples.
    pub fn new(
        opts: RohOpts,
        headers: &[String],
        samples: Vec<String>,
        roh_smpl: Vec<usize>,
        out: Box<dyn Write>,
    ) -> Result<Self> {
        if samples.is_empty() {
            bail!("No samples in the VCF?");
        }
        let has_fmt = |tag: &str| {
            headers
                .iter()
                .any(|h| h.starts_with(&format!("##FORMAT=<ID={tag},")))
        };
        let fake_pls = opts.gts_only.is_some();
        if !fake_pls {
            let Some(pl) = headers.iter().find(|h| h.starts_with("##FORMAT=<ID=PL,")) else {
                bail!("Error: The FORMAT/PL tag not found in the header, consider running with -G");
            };
            if !pl.contains("Type=Integer") {
                bail!("Error: The FORMAT/PL tag not defined as Integer in the header");
            }
        }
        let mut af_from_pl = false;
        let mut af_smpl = None;
        let mut estimate_af = false;
        if let Some(spec) = &opts.estimate_af {
            estimate_af = true;
            let mut s = spec.as_str();
            if let Some(r) = s.strip_prefix("GT,") {
                s = r;
            } else if let Some(r) = s.strip_prefix("PL,") {
                s = r;
                af_from_pl = true;
            }
            if s != "-" {
                let list = sample_list(&samples, s, true)?;
                if list.len() != samples.len() {
                    af_smpl = Some(list);
                }
            }
        }
        if estimate_af || fake_pls {
            if af_from_pl {
                if !has_fmt("PL") {
                    bail!("Error: The FORMAT/PL tag not found in the header");
                }
            } else if !has_fmt("GT") {
                bail!("Error: The FORMAT/GT tag not found in the header");
            }
        }
        let af_file = match &opts.af_file {
            Some(p) => Some(AfFile::load(Path::new(p))?),
            None => None,
        };
        let (mut nbuf_max, mut nbuf_olap) = (0usize, 0usize);
        if let Some(b) = &opts.buffer_size {
            let (first, olap) = match b.split_once(',') {
                Some((a, o)) => (a, Some(o)),
                None => (b.as_str(), None),
            };
            let Ok(tmp) = first.parse::<f64>() else {
                bail!("Could not parse: --buffer-size {b}")
            };
            let mut olap_v: i64 = -1;
            if let Some(o) = olap {
                let Ok(v) = o.parse::<i64>() else {
                    bail!("Could not parse: --bufer-size {b}")
                };
                if v < 0 {
                    bail!("Could not parse: --bufer-size {b}");
                }
                olap_v = v;
            }
            nbuf_max = if tmp < 0.0 {
                (tmp.abs() * 1e6 / (4.0 + 8.0 * 2.0) / roh_smpl.len() as f64) as usize
            } else {
                tmp as usize
            };
            nbuf_olap = if olap_v < 0 {
                (nbuf_max as f64 * 0.01) as usize
            } else {
                olap_v as usize
            };
        }
        eprintln!("Number of target samples: {}", roh_smpl.len());
        eprintln!(
            "Number of --estimate-AF samples: {}",
            match &af_smpl {
                Some(l) => l.len(),
                None if estimate_af => samples.len(),
                None => 0,
            }
        );
        if nbuf_max > 0 {
            eprintln!("Number of sites in the buffer/overlap: {nbuf_max}/{nbuf_olap}");
        } else {
            eprintln!("Number of sites in the buffer/overlap: unlimited");
        }
        let mut pl2p = [0.0f64; 256];
        for (i, v) in pl2p.iter_mut().enumerate() {
            *v = 10f64.powf(-(i as f64) / 10.0);
        }
        let tprob = [1.0 - opts.t2az, opts.t2hw, opts.t2az, 1.0 - opts.t2hw];
        let tfn = if opts.genmap_fname.is_some() {
            // Replaced per chromosome by `load_genmap`.
            RohTransition::None
        } else if opts.rec_rate > 0.0 {
            RohTransition::RecRate(opts.rec_rate)
        } else {
            RohTransition::None
        };
        let hmm = Hmm::new(2, &tprob, 10000, tfn);
        let unseen_pl = opts.gts_only.map(|g| 10f64.powf(-g / 10.0)).unwrap_or(0.0);
        let mut c = Self {
            smpl: (0..roh_smpl.len()).map(|_| Smpl::default()).collect(),
            opts,
            samples,
            roh_smpl,
            af_smpl,
            af_from_pl,
            estimate_af,
            unseen_pl,
            fake_pls,
            af_file,
            pl2p,
            nbuf_max,
            nbuf_olap,
            hmm,
            prev_rid: None,
            prev_pos: 0,
            skip_rid: None,
            stats: RohStats::default(),
            out,
        };
        c.write_header()?;
        Ok(c)
    }

    fn write_header(&mut self) -> Result<()> {
        let mut s = String::new();
        s.push_str(&format!(
            "# This file was produced by: kira-bt roh({})\n",
            env!("CARGO_PKG_VERSION")
        ));
        s.push_str(&format!(
            "# The command line was:\tkira-bt {}\n#\n",
            self.opts.cmdline
        ));
        if self.opts.output_rg {
            s.push_str("# RG\t[2]Sample\t[3]Chromosome\t[4]Start\t[5]End\t[6]Length (bp)\t[7]Number of markers\t[8]Quality (average fwd-bwd phred score)\n");
        }
        if self.opts.output_st {
            s.push_str("# ST\t[2]Sample\t[3]Chromosome\t[4]Position\t[5]State (0:HW, 1:AZ)\t[6]Quality (fwd-bwd phred score)\n");
        }
        if self.opts.vi_training.is_some() {
            s.push_str("# VT, Viterbi Training\t[2]Sample\t[3]Iteration\t[4]dAZ\t[5]dHW\t[6]1 - P(HW|HW)\t[7]P(AZ|HW)\t[8]1 - P(AZ|AZ)\t[9]P(HW|AZ)\n");
        }
        self.out.write_all(s.as_bytes())?;
        Ok(())
    }

    /// Targets semantics of `--AF-file`: the record must be listed.
    pub fn af_file_has(&self, chrom: &str, pos1: u32) -> bool {
        self.af_file
            .as_ref()
            .is_none_or(|f| f.has_site(chrom, pos1))
    }

    fn load_genmap(&mut self, chr: &str) -> Result<bool> {
        let Some(fname) = self.opts.genmap_fname.clone() else {
            return Ok(false);
        };
        match GenMap::load(&fname, chr)? {
            Some(map) => {
                self.hmm.tfn = RohTransition::GenMap {
                    map,
                    rec_rate: self.opts.rec_rate,
                };
                Ok(false)
            }
            None => Ok(true),
        }
    }

    fn genmap_index(&self) -> usize {
        match &self.hmm.tfn {
            RohTransition::GenMap { map, .. } => map.igenmap,
            _ => 0,
        }
    }

    fn set_genmap_index(&mut self, i: usize) {
        if let RohTransition::GenMap { map, .. } = &mut self.hmm.tfn {
            map.igenmap = i;
        }
    }

    /// `vcfroh`: one record; `pass` is the `-i/-e` result.
    pub fn process(&mut self, rec: &VcfRecord, pass: bool) -> Result<()> {
        self.stats.ntot += 1;
        if !pass {
            self.stats.nfiltered += 1;
            return Ok(());
        }
        if self.skip_rid.as_deref() == Some(rec.chrom.as_str()) {
            return Ok(());
        }
        let mut alleles = vec![rec.ref_allele.clone()];
        if rec.alt != "." && !rec.alt.is_empty() {
            alleles.extend(rec.alt.split(',').map(str::to_string));
        }
        let mut ial = 0usize;
        let mut nalt = alleles.len() as i64 - 1;
        for (i, a) in alleles.iter().enumerate().skip(1) {
            if a == "<*>" || a == "<NON_REF>" {
                nalt -= 1;
            } else if ial == 0 {
                ial = i;
            }
        }
        if nalt == 0 {
            self.stats.nno_alt += 1;
            if !self.opts.include_noalt {
                return Ok(());
            }
        } else if nalt > 1 {
            self.stats.nmultiallelic += 1;
            return Ok(());
        }
        if self.opts.snps_only && !bcf_is_snp(&alleles) {
            return Ok(());
        }
        let pos = rec.pos.saturating_sub(1);
        let mut skip = false;
        if self.prev_rid.is_none() {
            skip = self.load_genmap(&rec.chrom)?;
        }
        if self.prev_rid.as_deref() != Some(rec.chrom.as_str()) {
            if self.opts.vi_training.is_none() {
                for i in 0..self.smpl.len() {
                    self.flush_viterbi(i)?;
                    let snap = self.smpl[i].snapshot.as_mut();
                    self.hmm.reset(snap);
                }
            }
            self.prev_rid = Some(rec.chrom.clone());
            self.prev_pos = pos;
            skip = self.load_genmap(&rec.chrom)?;
        } else if self.prev_pos == pos {
            self.stats.ndup += 1;
            return Ok(());
        }
        if skip {
            eprintln!("Skipping the sequence, no genmap for {}", rec.chrom);
            self.skip_rid = Some(rec.chrom.clone());
            return Ok(());
        }
        if self.prev_pos > pos {
            bail!("The file is not sorted?!");
        }
        self.prev_rid = Some(rec.chrom.clone());
        self.prev_pos = pos;
        self.process_line(rec, &alleles, ial)?;
        Ok(())
    }

    fn estimate_af_from_gt(&self, gts: &[[Option<i32>; 2]]) -> Option<f64> {
        let (mut nalt, mut nref) = (0usize, 0usize);
        let idx: Vec<usize> = match &self.af_smpl {
            Some(l) => l.clone(),
            None => (0..self.samples.len()).collect(),
        };
        for i in idx {
            let g = gts[i];
            let (Some(a), Some(b)) = (g[0], g[1]) else {
                continue;
            };
            if a != 0 {
                nalt += 1
            } else {
                nref += 1
            }
            if b != 0 { nalt += 1 } else { nref += 1 }
        }
        if nalt == 0 && nref == 0 {
            return None;
        }
        Some(nalt as f64 / (nalt + nref) as f64)
    }

    fn estimate_af_from_pl(&self, pls: &[Vec<i64>], n: usize, ial: usize) -> Option<f64> {
        let (irr, ira, iaa) = (alleles2gt(0, 0), alleles2gt(0, ial), alleles2gt(ial, ial));
        if iaa >= n {
            return None;
        }
        let idx: Vec<usize> = match &self.af_smpl {
            Some(l) => l.clone(),
            None => (0..self.samples.len()).collect(),
        };
        let mut af = 0.0;
        let mut naf = 0usize;
        for i in idx {
            let p = &pls[i];
            if p[irr] < 0 || p[ira] < 0 || p[iaa] < 0 {
                continue;
            }
            if p[irr] == p[ira] && p[irr] == p[iaa] {
                continue;
            }
            let mut prob = [
                self.pl2p[p[irr].min(255) as usize],
                self.pl2p[p[ira].min(255) as usize],
                self.pl2p[p[iaa].min(255) as usize],
            ];
            let norm: f64 = prob.iter().sum();
            for v in &mut prob {
                *v /= norm;
            }
            af += 0.5 * prob[1] + prob[2];
            naf += 1;
        }
        if naf == 0 {
            return None;
        }
        Some(af / naf as f64)
    }

    /// `process_line`: allele frequency, per-sample P(D|G) and emissions.
    fn process_line(&mut self, rec: &VcfRecord, alleles: &[String], ial: usize) -> Result<()> {
        let mut alt_freq = 0.0f64;
        let mut ret: i32 = 0;
        let mut gts: Option<Vec<[Option<i32>; 2]>> = None;
        let mut pls: Option<(Vec<Vec<i64>>, usize)> = None;
        if let Some(tag) = &self.opts.af_tag {
            match info_value(&rec.info, tag) {
                Some(v) => {
                    let vals: Vec<&str> = v.split(',').collect();
                    let s = vals.get(ial - 1).copied().unwrap_or(".");
                    alt_freq = if s == "." {
                        f64::NAN
                    } else {
                        s.parse::<f64>().map_err(|_| {
                            anyhow::anyhow!(
                                "Type mismatch for INFO/{} tag at {}:{}",
                                tag,
                                rec.chrom,
                                rec.pos
                            )
                        })?
                    };
                }
                None => ret = -1,
            }
        } else if let Some(af) = &self.af_file {
            match af.read_af(&rec.chrom, rec.pos, alleles)? {
                Some(v) => alt_freq = v,
                None => ret = -1,
            }
        } else if self.opts.dflt_af > 0.0 {
            alt_freq = self.opts.dflt_af;
        } else if self.estimate_af {
            if self.af_from_pl {
                let Some(p) = parse_pls(rec) else {
                    return Ok(());
                };
                match self.estimate_af_from_pl(&p.0, p.1, ial) {
                    Some(v) => alt_freq = v,
                    None => ret = -1,
                }
                pls = Some(p);
            } else {
                let Some(g) = parse_gt_pairs(rec) else {
                    return Ok(());
                };
                match self.estimate_af_from_gt(&g) {
                    Some(v) => alt_freq = v,
                    None => ret = -1,
                }
                gts = Some(g);
            }
        } else {
            let an =
                info_value(&rec.info, "AN").and_then(|v| v.split(',').next()?.parse::<i64>().ok());
            let ac =
                info_value(&rec.info, "AC").and_then(|v| v.split(',').next()?.parse::<i64>().ok());
            match (an, ac) {
                (Some(an), Some(ac)) if an > 0 && ac >= 0 => alt_freq = ac as f64 / an as f64,
                _ => ret = -1,
            }
        }
        if self.opts.dflt_af > 0.0 && (ret < 0 || alt_freq == 0.0) {
            alt_freq = self.opts.dflt_af;
        } else if ret < 0 {
            self.stats.nno_af += 1;
            return Ok(());
        } else if alt_freq == 0.0 {
            self.stats.nno_af += 1;
            return Ok(());
        }
        let (irr, ira, iaa) = (alleles2gt(0, 0), alleles2gt(0, ial), alleles2gt(ial, ial));
        if self.fake_pls {
            if gts.is_none() {
                gts = parse_gt_pairs(rec);
            }
            if gts.is_none() {
                return Ok(());
            }
        } else {
            if pls.is_none() {
                pls = parse_pls(rec);
            }
            let Some((_, n)) = &pls else { return Ok(()) };
            if iaa >= *n {
                return Ok(());
            }
        }
        let u = self.unseen_pl;
        for i in 0..self.roh_smpl.len() {
            let ismpl = self.roh_smpl[i];
            let mut pdg = [0.0f64; 3];
            if self.fake_pls {
                let g = gts.as_ref().unwrap()[ismpl];
                let (Some(a), Some(b)) = (g[0], g[1]) else {
                    continue;
                };
                if a != b {
                    pdg = [u, 1.0 - 2.0 * u, u];
                } else if a == 0 {
                    pdg = [1.0 - u - u * u, u, u * u];
                } else {
                    pdg = [u * u, u, 1.0 - u - u * u];
                }
            } else {
                let p = &pls.as_ref().unwrap().0[ismpl];
                if p[irr] < 0 || p[ira] < 0 || p[iaa] < 0 {
                    continue;
                }
                if p[irr] == p[ira] && p[irr] == p[iaa] {
                    continue;
                }
                pdg[0] = self.pl2p[p[irr].min(255) as usize];
                pdg[1] = self.pl2p[p[ira].min(255) as usize];
                pdg[2] = self.pl2p[p[iaa].min(255) as usize];
            }
            let sum = pdg[0] + pdg[1] + pdg[2];
            if sum == 0.0 {
                continue;
            }
            for v in &mut pdg {
                *v /= sum;
            }
            if self.opts.skip_homref && pdg[0] > 0.99 {
                continue;
            }
            let s = &mut self.smpl[i];
            s.nused += 1;
            let e_az = pdg[0] * (1.0 - alt_freq) + pdg[2] * alt_freq;
            let e_hw = pdg[0] * (1.0 - alt_freq) * (1.0 - alt_freq)
                + 2.0 * pdg[1] * (1.0 - alt_freq) * alt_freq
                + pdg[2] * alt_freq * alt_freq;
            s.eprob.push(e_hw);
            s.eprob.push(e_az);
            s.sites.push(rec.pos.saturating_sub(1));
            if self.opts.vi_training.is_some() {
                if s.rid.last().map(String::as_str) != Some(rec.chrom.as_str()) {
                    s.rid.push(rec.chrom.clone());
                    s.rid_off.push(s.sites.len() - 1);
                }
            } else if self.nbuf_max > 0 && s.sites.len() >= self.nbuf_max {
                self.flush_viterbi(i)?;
            }
        }
        Ok(())
    }

    fn write_rg(out: &mut dyn Write, name: &str, rg: &Rg) -> Result<()> {
        writeln!(
            out,
            "RG\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}",
            name,
            rg.rid,
            rg.beg + 1,
            rg.end + 1,
            rg.end - rg.beg + 1,
            rg.nqual,
            rg.qual / rg.nqual as f64
        )?;
        Ok(())
    }

    /// `flush_viterbi`: run the HMM on the buffered sites of one sample and
    /// print ST/RG (or VT/ROH with `-V`) lines.
    fn flush_viterbi(&mut self, i: usize) -> Result<()> {
        if self.smpl[i].sites.is_empty() {
            return Ok(());
        }
        let name = self.samples[self.roh_smpl[i]].clone();
        if self.opts.vi_training.is_none() {
            let snap = self.smpl[i].snapshot.as_ref();
            self.hmm.restore(snap);
            let nsites = self.smpl[i].sites.len();
            let end = if self.nbuf_max > 0 && nsites >= self.nbuf_max && nsites > self.nbuf_olap {
                nsites - self.nbuf_olap
            } else {
                nsites
            };
            let requested = end < nsites;
            if requested {
                let pos = self.smpl[i].sites[nsites - self.nbuf_olap - 1];
                let prev = self.smpl[i].snapshot.take();
                self.hmm.request_snapshot(prev, pos);
            }
            let ig = self.smpl[i].igenmap;
            self.set_genmap_index(ig);
            let (eprob, sites) = {
                let s = &self.smpl[i];
                (s.eprob.clone(), s.sites.clone())
            };
            self.hmm.run_viterbi(&eprob, &sites);
            self.hmm.run_fwd_bwd(&eprob, &sites);
            if requested {
                self.smpl[i].snapshot = self.hmm.take_snapshot();
            }
            let chr = self.prev_rid.clone().unwrap_or_default();
            let mut lines = String::new();
            for k in 0..end {
                let state = if self.hmm.vstate(k) == STATE_AZ as u8 {
                    1usize
                } else {
                    0
                };
                // bcftools reads row k, the posterior of the previous site.
                let qual = phred_score(1.0 - self.hmm.fwd_bwd_prob()[k * 2 + state]);
                if self.opts.output_st {
                    lines.push_str(&format!(
                        "ST\t{}\t{}\t{}\t{}\t{:.1}\n",
                        name,
                        chr,
                        sites[k] + 1,
                        state,
                        qual
                    ));
                }
                if self.opts.output_rg {
                    let s = &mut self.smpl[i];
                    if state as u8 != s.rg.state {
                        if state == 0 {
                            self.out.write_all(lines.as_bytes())?;
                            lines.clear();
                            Self::write_rg(&mut *self.out, &name, &s.rg)?;
                            s.rg.state = 0;
                        } else {
                            s.rg.state = 1;
                            s.rg.beg = sites[k];
                            s.rg.end = sites[k];
                            s.rg.rid = chr.clone();
                            s.rg.qual = qual;
                            s.rg.nqual = 1;
                        }
                    } else if state == 1 {
                        s.rg.nqual += 1;
                        s.rg.qual += qual;
                        s.rg.end = sites[k];
                    }
                }
            }
            self.out.write_all(lines.as_bytes())?;
            let ig = self.genmap_index();
            let s = &mut self.smpl[i];
            if end < nsites {
                let keep = nsites - self.nbuf_olap;
                s.sites.drain(..keep);
                s.eprob.drain(..keep * 2);
                s.igenmap = ig;
            } else {
                s.sites.clear();
                s.eprob.clear();
                s.igenmap = 0;
                if s.rg.state != 0 {
                    Self::write_rg(&mut *self.out, &name, &s.rg)?;
                    s.rg.state = 0;
                }
            }
            return Ok(());
        }
        // Viterbi training over all chromosomes of the sample.
        let th = self.opts.vi_training.unwrap();
        let mut tprob = self.hmm.tprob().to_vec();
        tprob[2 * STATE_HW + STATE_HW] = 1.0 - self.opts.t2az;
        tprob[2 * STATE_HW + STATE_AZ] = self.opts.t2hw;
        tprob[2 * STATE_AZ + STATE_HW] = self.opts.t2az;
        tprob[2 * STATE_AZ + STATE_AZ] = 1.0 - self.opts.t2hw;
        self.hmm.set_tprob(&tprob, 10000);
        let nrid = self.smpl[i].rid.len();
        let mut niter = 0;
        loop {
            let cur = self.hmm.tprob().to_vec();
            let t2az_prev = cur[2 * STATE_AZ + STATE_HW];
            let t2hw_prev = cur[2 * STATE_HW + STATE_AZ];
            let mut tprob_new = [0.0f64; 4];
            for r in 0..nrid {
                let ioff = self.smpl[i].rid_off[r];
                let nsites = if r + 1 == nrid {
                    self.smpl[i].sites.len()
                } else {
                    self.smpl[i].rid_off[r + 1]
                } - ioff;
                self.set_genmap_index(0);
                let (eprob, sites) = {
                    let s = &self.smpl[i];
                    (
                        s.eprob[ioff * 2..(ioff + nsites) * 2].to_vec(),
                        s.sites[ioff..ioff + nsites].to_vec(),
                    )
                };
                let t = self.hmm.run_baum_welch(&eprob, &sites);
                for k in 0..4 {
                    tprob_new[k] += t[k];
                }
            }
            for v in &mut tprob_new {
                *v /= nrid as f64;
            }
            self.hmm.set_tprob(&tprob_new, 10000);
            let deltaz = (tprob_new[2 * STATE_AZ + STATE_HW] - t2az_prev).abs();
            let delthw = (tprob_new[2 * STATE_HW + STATE_AZ] - t2hw_prev).abs();
            niter += 1;
            writeln!(
                self.out,
                "VT\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                name,
                niter,
                fmt_e(deltaz, 6),
                fmt_e(delthw, 6),
                fmt_e(1.0 - tprob_new[2 * STATE_HW + STATE_HW], 6),
                fmt_e(tprob_new[2 * STATE_AZ + STATE_HW], 6),
                fmt_e(1.0 - tprob_new[2 * STATE_AZ + STATE_AZ], 6),
                fmt_e(tprob_new[2 * STATE_HW + STATE_AZ], 6)
            )?;
            if !(deltaz > th || delthw > th) {
                break;
            }
        }
        for r in 0..nrid {
            let ioff = self.smpl[i].rid_off[r];
            let nsites = if r + 1 == nrid {
                self.smpl[i].sites.len()
            } else {
                self.smpl[i].rid_off[r + 1]
            } - ioff;
            self.set_genmap_index(0);
            let (eprob, sites, chr) = {
                let s = &self.smpl[i];
                (
                    s.eprob[ioff * 2..(ioff + nsites) * 2].to_vec(),
                    s.sites[ioff..ioff + nsites].to_vec(),
                    s.rid[r].clone(),
                )
            };
            self.hmm.run_viterbi(&eprob, &sites);
            self.hmm.run_fwd_bwd(&eprob, &sites);
            for j in 0..nsites {
                let state = if self.hmm.vstate(j) == STATE_AZ as u8 {
                    1usize
                } else {
                    0
                };
                let pval = self.hmm.fwd_bwd_prob()[j * 2 + state];
                writeln!(
                    self.out,
                    "ROH\t{}\t{}\t{}\t{}\t{:.1}",
                    name,
                    chr,
                    sites[j] + 1,
                    state,
                    phred_score(1.0 - pval)
                )?;
            }
        }
        Ok(())
    }

    /// Flush every sample and print the counters.
    pub fn finish(&mut self) -> Result<()> {
        for i in 0..self.smpl.len() {
            self.flush_viterbi(i)?;
        }
        self.out.flush()?;
        let nmin = self.smpl.iter().map(|s| s.nused).min().unwrap_or(0);
        if self.af_file.is_some() {
            eprintln!(
                "Number of lines overlapping with --AF-file/processed: {}/{}",
                self.stats.ntot, nmin
            );
        } else {
            eprintln!(
                "Number of lines total/processed: {}/{}",
                self.stats.ntot, nmin
            );
        }
        eprintln!(
            "Number of lines filtered/no AF/no alt/multiallelic/dup: {}/{}/{}/{}/{}",
            self.stats.nfiltered,
            self.stats.nno_af,
            self.stats.nno_alt,
            self.stats.nmultiallelic,
            self.stats.ndup
        );
        if nmin == 0 {
            eprintln!("No usable sites were found.");
            if self.opts.af_tag.is_none()
                && self.opts.af_file.is_none()
                && !self.estimate_af
                && self.opts.dflt_af == 0.0
            {
                eprintln!(" Consider using one of the AF options.");
            }
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "../tests/unit/roh.rs"]
mod tests;
