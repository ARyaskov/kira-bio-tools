//! Copy-number variation from Illumina BAF/LRR intensities: a port of
//! bcftools `vcfcnv.c` on the shared [`crate::hmm`] kernel. Four states
//! (CN0..CN3) for one sample, 16 paired states with a control sample; BAF
//! emissions are mixtures of Gaussians per genotype, LRR emissions
//! Gaussians around the expected shifts. `f32` is kept where the C code
//! stores floats so printed probabilities agree.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};

use crate::af_file::AfFile;
use crate::cli::getopt::{OptSpec, opt, parse_f64, parse_i64};
use crate::hmm::Hmm;
use crate::regions::RegionSet;
use crate::util::{erfc, phred_score};
use crate::vcf::header::extract_samples;
use crate::vcf::{ContigDict, UnifiedVcfReader, VcfRecord, parse_vcf_record};

pub const N_STATES: usize = 4;
const CN0: usize = 0;
const CN1: usize = 1;
const CN2: usize = 2;
const CN3: usize = 3;

#[derive(Clone, Copy, Default, Debug)]
struct GaussParam {
    mean: f32,
    dev2: f32,
    norm: f32,
}

#[derive(Default)]
struct Sample {
    name: Option<String>,
    idx: Option<usize>,
    lrr: Vec<f32>,
    baf: Vec<f32>,
    baf_dev2: f32,
    baf_dev2_dflt: f32,
    lrr_dev2: f32,
    cell_frac: f32,
    cell_frac_dflt: f32,
    gauss: [GaussParam; 18],
    pobs: [f64; N_STATES],
    dat_fh: Option<BufWriter<File>>,
    cn_fh: Option<BufWriter<File>>,
    summary_fh: Option<BufWriter<File>>,
    dat_fname: String,
    cn_fname: String,
}

const CN1_PK_R: usize = 0;
const CN1_PK_A: usize = 1;
const CN2_PK_RR: usize = 2;
const CN2_PK_RA: usize = 3;
const CN2_PK_AA: usize = 4;
const CN3_PK_RRR: usize = 5;
const CN3_PK_RRA: usize = 6;
const CN3_PK_RAA: usize = 7;
const CN3_PK_AAA: usize = 8;

/// Options of `bcftools cnv`.
#[derive(Clone, Debug)]
pub struct CnvConfig {
    pub input: PathBuf,
    pub output_dir: PathBuf,
    pub query_sample: Option<String>,
    pub control_sample: Option<String>,
    pub af_file: Option<PathBuf>,
    pub regions: Option<String>,
    pub regions_is_file: bool,
    pub regions_overlap: u8,
    pub targets: Option<String>,
    pub targets_is_file: bool,
    pub targets_overlap: u8,
    pub plot_th: f32,
    pub lrr_smooth_win: i64,
    pub lrr_bias: f64,
    pub baf_bias: f64,
    pub err_prob: f64,
    pub ij_prob: f64,
    pub same_prob: f64,
    pub optimize_frac: f64,
    pub baum_welch_th: f64,
    /// Squared BAF deviation of query and control.
    pub baf_dev2: (f32, f32),
    pub lrr_dev2: (f32, f32),
    pub cell_frac: (f32, f32),
    pub cmdline: String,
}

impl Default for CnvConfig {
    fn default() -> Self {
        Self {
            input: PathBuf::from("-"),
            output_dir: PathBuf::new(),
            query_sample: None,
            control_sample: None,
            af_file: None,
            regions: None,
            regions_is_file: false,
            regions_overlap: 1,
            targets: None,
            targets_is_file: false,
            targets_overlap: 0,
            plot_th: 1e9,
            lrr_smooth_win: 10,
            lrr_bias: 0.2,
            baf_bias: 1.0,
            err_prob: 1e-4,
            ij_prob: 1e-9,
            same_prob: 0.5,
            optimize_frac: 0.0,
            baum_welch_th: 0.0,
            baf_dev2: (0.04 * 0.04, 0.04 * 0.04),
            lrr_dev2: (0.2 * 0.2, 0.2 * 0.2),
            cell_frac: (1.0, 1.0),
            cmdline: String::new(),
        }
    }
}

const SPECS: &[OptSpec] = &[
    opt(Some('d'), "BAF-dev", true),
    opt(Some('k'), "LRR-dev", true),
    opt(Some('L'), "LRR-smooth-win", true),
    opt(Some('f'), "AF-file", true),
    opt(Some('W'), "baum-welch", true),
    opt(Some('O'), "optimize", true),
    opt(Some('a'), "aberrant", true),
    opt(Some('e'), "err-prob", true),
    opt(Some('b'), "BAF-weight", true),
    opt(Some('l'), "LRR-weight", true),
    opt(Some('P'), "same-prob", true),
    opt(Some('x'), "xy-prob", true),
    opt(Some('s'), "query-sample", true),
    opt(Some('c'), "control-sample", true),
    opt(Some('t'), "targets", true),
    opt(Some('T'), "targets-file", true),
    opt(None, "targets-overlap", true),
    opt(Some('r'), "regions", true),
    opt(Some('R'), "regions-file", true),
    opt(None, "regions-overlap", true),
    opt(Some('p'), "plot-threshold", true),
    opt(Some('o'), "output-dir", true),
];

fn parse_pair(v: &str, what: &str) -> Result<(f64, Option<f64>)> {
    match v.split_once(',') {
        Some((a, b)) => Ok((parse_f64(a, what)?, Some(parse_f64(b, what)?))),
        None => Ok((parse_f64(v, what)?, None)),
    }
}

fn overlap_opt(v: &str, what: &str) -> Result<u8> {
    match v {
        "0" | "pos" => Ok(0),
        "1" | "record" => Ok(1),
        "2" | "variant" => Ok(2),
        _ => bail!("Could not parse: {what} {v}"),
    }
}

/// Parse a bcftools `cnv` command line (without the program name).
pub fn parse_args(args: &[String]) -> Result<CnvConfig> {
    let mut cfg = CnvConfig::default();
    let p = crate::cli::getopt::parse(args, SPECS)?;
    for (name, val) in p.opts {
        let v = val.unwrap_or_default();
        match name {
            "LRR-smooth-win" => cfg.lrr_smooth_win = parse_i64(&v, "--LRR-smooth-win")?,
            "AF-file" => cfg.af_file = Some(PathBuf::from(v)),
            "optimize" => cfg.optimize_frac = parse_f64(&v, "-O")?,
            "BAF-dev" => {
                let (a, b) = parse_pair(&v, "-d")?;
                let b = b.unwrap_or(a);
                cfg.baf_dev2 = ((a as f32) * (a as f32), (b as f32) * (b as f32));
            }
            "LRR-dev" => {
                let (a, b) = parse_pair(&v, "-k")?;
                let b = b.unwrap_or(a);
                cfg.lrr_dev2 = ((a as f32) * (a as f32), (b as f32) * (b as f32));
            }
            "aberrant" => {
                let (a, b) = parse_pair(&v, "-a")?;
                cfg.cell_frac = (a as f32, b.map(|x| x as f32).unwrap_or(cfg.cell_frac.1));
            }
            "baum-welch" => cfg.baum_welch_th = parse_f64(&v, "-W")?,
            "err-prob" => cfg.err_prob = parse_f64(&v, "-e")?,
            "BAF-weight" => cfg.baf_bias = parse_f64(&v, "-b")?,
            "xy-prob" => cfg.ij_prob = parse_f64(&v, "-x")?,
            "same-prob" => cfg.same_prob = parse_f64(&v, "-P")?,
            "LRR-weight" => cfg.lrr_bias = parse_f64(&v, "-l")?,
            "plot-threshold" => cfg.plot_th = parse_f64(&v, "-p")? as f32,
            "output-dir" => cfg.output_dir = PathBuf::from(v),
            "query-sample" => cfg.query_sample = Some(v),
            "control-sample" => cfg.control_sample = Some(v),
            "targets" => {
                cfg.targets = Some(v);
                cfg.targets_is_file = false;
            }
            "targets-file" => {
                cfg.targets = Some(v);
                cfg.targets_is_file = true;
            }
            "regions" => {
                cfg.regions = Some(v);
                cfg.regions_is_file = false;
            }
            "regions-file" => {
                cfg.regions = Some(v);
                cfg.regions_is_file = true;
            }
            "regions-overlap" => cfg.regions_overlap = overlap_opt(&v, "--regions-overlap")?,
            "targets-overlap" => cfg.targets_overlap = overlap_opt(&v, "--targets-overlap")?,
            _ => {}
        }
    }
    match p.positional.first() {
        Some(f) => cfg.input = PathBuf::from(f),
        None => bail!(
            "Usage:   kira-bt cnv [OPTIONS] FILE.vcf\n         requires -s/--query-sample (multi-sample files) and -o/--output-dir"
        ),
    }
    if cfg.output_dir.as_os_str().is_empty() {
        bail!("Expected -o option");
    }
    cfg.cmdline = args.join(" ");
    Ok(cfg)
}

pub fn run_from_args(args: &[String]) -> Result<()> {
    let cfg = parse_args(args)?;
    run(cfg)
}

/// `init_tprob_matrix`: single-sample or paired transition matrix.
fn init_tprob_matrix(ndim: usize, ij_prob: f64, same_prob: f64) -> Result<Vec<f64>> {
    let mut mat = vec![0.0f64; ndim * ndim];
    if ndim == N_STATES {
        let pii = 1.0 - ij_prob * (N_STATES as f64 - 1.0);
        if pii < ij_prob {
            bail!(
                "Error: -x set a bit too high, P(x|x) < P(x|y): {} vs {}",
                crate::util::fmt_e(pii, 6),
                crate::util::fmt_e(ij_prob, 6)
            );
        }
        for j in 0..ndim {
            for i in 0..ndim {
                mat[ndim * i + j] = if i == j { pii } else { ij_prob };
            }
        }
    } else {
        let pii = 1.0 - ij_prob * (N_STATES as f64 - 1.0);
        let ij_prob = (1.0 - pii) / (ndim as f64 - 1.0);
        for j in 0..ndim {
            let (ja, jb) = (j / N_STATES, j % N_STATES);
            let mut sum = 0.0;
            for i in 0..ndim {
                let (ia, ib) = (i / N_STATES, i % N_STATES);
                let pa = if ja == ia { pii } else { ij_prob };
                let pb = if jb == ib { pii } else { ij_prob };
                mat[ndim * i + j] = if ia == ib && ja == jb {
                    pa * pb - pa * pb * same_prob + (pa * pb).sqrt() * same_prob
                } else if ia == ib {
                    pa * pb
                } else {
                    pa * pb * (1.0 - same_prob)
                };
                sum += mat[ndim * i + j];
            }
            for i in 0..ndim {
                mat[ndim * i + j] /= sum;
            }
        }
    }
    Ok(mat)
}

/// `init_iprobs`: prior on CN2.
fn init_iprobs(ndim: usize, same_prob: f64) -> Vec<f64> {
    let mut probs = vec![0.0f64; ndim];
    if ndim == N_STATES {
        for (i, p) in probs.iter_mut().enumerate() {
            *p = if i == CN2 { 0.5 } else { 0.5 / 3.0 };
        }
    } else {
        let mut norm = 0.0;
        for (i, p) in probs.iter_mut().enumerate() {
            let (ia, ib) = (i / N_STATES, i % N_STATES);
            let pa = if ia == CN2 { 0.5 } else { 0.5 / 3.0 };
            let pb = if ib == CN2 { 0.5 } else { 0.5 / 3.0 };
            *p = pa * pb;
            if ia != ib {
                *p *= 1.0 - same_prob;
            }
            norm += *p;
        }
        for p in &mut probs {
            *p /= norm;
        }
    }
    probs
}

fn norm_cdf(mean: f64, dev: f64) -> f64 {
    let top = 1.0 - 0.5 * erfc((1.0 - mean) / (dev * 2f64.sqrt()));
    let bot = 1.0 - 0.5 * erfc((0.0 - mean) / (dev * 2f64.sqrt()));
    top - bot
}

/// `norm_prob` with the C float/double promotions.
fn norm_prob(baf: f32, p: &GaussParam) -> f64 {
    let d = baf - p.mean;
    let d2 = -d * d;
    (d2 as f64 * 0.5 / p.dev2 as f64).exp()
        / p.norm as f64
        / (2.0 * std::f64::consts::PI * p.dev2 as f64).sqrt()
}

fn set_gauss_params(s: &mut Sample) {
    for g in &mut s.gauss {
        g.dev2 = s.baf_dev2;
    }
    let dev = (s.baf_dev2 as f64).sqrt();
    let means: [(usize, f64); 9] = [
        (CN1_PK_R, 0.0),
        (CN1_PK_A, 1.0),
        (CN2_PK_RR, 0.0),
        (CN2_PK_RA, 0.5),
        (CN2_PK_AA, 1.0),
        (CN3_PK_RRR, 0.0),
        // `2+smpl->cell_frac` is a float sum in C, `1.0+…` a double one.
        (CN3_PK_RRA, 1.0 / ((2.0f32 + s.cell_frac) as f64)),
        (
            CN3_PK_RAA,
            (1.0 + s.cell_frac as f64) / ((2.0f32 + s.cell_frac) as f64),
        ),
        (CN3_PK_AAA, 1.0),
    ];
    for (k, m) in means {
        s.gauss[k].mean = m as f32;
        s.gauss[k].norm = norm_cdf(s.gauss[k].mean as f64, dev) as f32;
    }
}

/// `smooth_data`: moving average of window `win` in single precision. The
/// C code primes the window with `win - win/2` values whatever the array
/// length, reading past short arrays; those reads count as zeros here.
fn smooth_data(dat: &mut [f32], win: i64) {
    if win <= 1 {
        return;
    }
    let win = win as usize;
    let (k1, k2) = (win / 2, win - win / 2);
    let ndat = dat.len();
    let mut buf: std::collections::VecDeque<f32> = std::collections::VecDeque::with_capacity(win);
    let mut sum: f32 = 0.0;
    for i in 0..k2 {
        let v = dat.get(i).copied().unwrap_or(0.0);
        sum += v;
        buf.push_back(v);
    }
    for i in 0..ndat {
        dat[i] = sum / buf.len() as f32;
        if i >= k1 {
            if let Some(v) = buf.pop_front() {
                sum -= v;
            }
        }
        if i + k2 < ndat {
            sum += dat[i + k2];
            buf.push_back(dat[i + k2]);
        }
    }
}

fn copy_number_state(two: bool, istate: usize, ismpl: usize) -> char {
    const CODE: [char; 5] = ['0', '1', '2', '3', '4'];
    if !two {
        return CODE[istate];
    }
    let idx = if ismpl != 0 {
        istate - (istate / N_STATES) * N_STATES
    } else {
        istate / N_STATES
    };
    CODE[idx]
}

fn avg_ii_prob(n: usize, mat: &[f64]) -> f64 {
    let mut avg = 0.0;
    for i in 0..n {
        avg += mat[n * i + i];
    }
    avg / n as f64
}

struct Cnv {
    cfg: CnvConfig,
    query: Sample,
    control: Sample,
    two: bool,
    nstates: usize,
    tprob: Vec<f64>,
    hmm: Hmm,
    eprob: Vec<f64>,
    sites: Vec<u32>,
    nonref_afs: Vec<f32>,
    f_rr: f32,
    f_ra: f32,
    f_aa: f32,
    n_rr: u64,
    n_ra: u64,
    n_aa: u64,
    af_file: Option<AfFile>,
    summary_fh: Option<BufWriter<File>>,
    prev_rid: Option<String>,
    ntot: usize,
    nused: usize,
}

fn open_out(dir: &Path, name: &str) -> Result<(BufWriter<File>, String)> {
    let p = dir.join(name);
    let f = File::create(&p).with_context(|| format!("Error opening file: {}", p.display()))?;
    Ok((BufWriter::new(f), p.display().to_string()))
}

fn init_sample_files(s: &mut Sample, dir: &Path) -> Result<()> {
    let name = s.name.clone().unwrap_or_default();
    let (mut dat, dat_fname) = open_out(dir, &format!("dat.{name}.tab"))?;
    let (mut cn, cn_fname) = open_out(dir, &format!("cn.{name}.tab"))?;
    let (mut sm, _) = open_out(dir, &format!("summary.{name}.tab"))?;
    writeln!(dat, "# [1]Chromosome\t[2]Position\t[3]BAF\t[4]LRR")?;
    writeln!(
        cn,
        "# [1]Chromosome\t[2]Position\t[3]CN\t[4]P(CN0)\t[5]P(CN1)\t[6]P(CN2)\t[7]P(CN3)"
    )?;
    writeln!(
        sm,
        "# RG, Regions [2]Chromosome\t[3]Start\t[4]End\t[5]Copy Number state\t[6]Quality\t[7]nSites\t[8]nHETs"
    )?;
    s.dat_fh = Some(dat);
    s.cn_fh = Some(cn);
    s.summary_fh = Some(sm);
    s.dat_fname = dat_fname;
    s.cn_fname = cn_fname;
    Ok(())
}

/// Per-sample BAF/LRR of one record: `(baf, lrr, present)`; missing BAF is
/// `-0.1`. A sample absent from the file reads as present, mirroring the
/// out-of-range read bcftools performs for a missing control.
fn parse_lrr_baf(
    rec: &VcfRecord,
    ibaf: usize,
    ilrr: Option<usize>,
    idx: Option<usize>,
) -> (f32, f32, bool) {
    let Some(idx) = idx else {
        return (-0.1, 0.0, true);
    };
    let field = |i: usize| -> Option<f32> {
        let s = rec.samples.get(idx)?.split(':').nth(i)?;
        if s == "." || s.is_empty() {
            return None;
        }
        s.parse::<f32>().ok().filter(|v| !v.is_nan())
    };
    let mut baf = field(ibaf).unwrap_or(-0.1);
    let lrr = match ilrr {
        Some(i) => match field(i) {
            Some(v) => v,
            None => {
                baf = -0.1;
                0.0
            }
        },
        None => 0.0,
    };
    (baf, lrr, baf >= 0.0)
}

impl Cnv {
    fn new(cfg: CnvConfig, headers: &[String]) -> Result<Self> {
        let samples = extract_samples(headers);
        let qname = match &cfg.query_sample {
            Some(n) => {
                if !samples.iter().any(|s| s == n) {
                    bail!("The sample \"{n}\" not found");
                }
                n.clone()
            }
            None => {
                if samples.len() > 1 {
                    bail!("Multi-sample VCF, missing the -s option");
                }
                samples.first().cloned().unwrap_or_default()
            }
        };
        let qidx = samples.iter().position(|s| *s == qname);
        let cidx = match &cfg.control_sample {
            Some(c) => match samples.iter().position(|s| s == c) {
                Some(i) => Some(i),
                None => bail!("The sample not found in the VCF: {c}"),
            },
            None => None,
        };
        let two = cfg.control_sample.is_some();
        let nstates = if two { N_STATES * N_STATES } else { N_STATES };
        let tprob = init_tprob_matrix(nstates, cfg.ij_prob, cfg.same_prob)?;
        let iprobs = init_iprobs(nstates, cfg.same_prob);
        let mut hmm = Hmm::new(nstates, &tprob, 10000, ());
        hmm.init_states(Some(&iprobs));
        std::fs::create_dir_all(&cfg.output_dir)
            .with_context(|| format!("Error opening file: {}", cfg.output_dir.display()))?;
        let mut query = Sample {
            name: Some(qname),
            idx: qidx,
            baf_dev2: cfg.baf_dev2.0,
            baf_dev2_dflt: cfg.baf_dev2.0,
            lrr_dev2: cfg.lrr_dev2.0,
            cell_frac: cfg.cell_frac.0,
            cell_frac_dflt: cfg.cell_frac.0,
            ..Default::default()
        };
        let mut control = Sample {
            name: cfg.control_sample.clone(),
            idx: cidx,
            baf_dev2: cfg.baf_dev2.1,
            baf_dev2_dflt: cfg.baf_dev2.1,
            lrr_dev2: cfg.lrr_dev2.1,
            cell_frac: cfg.cell_frac.1,
            cell_frac_dflt: cfg.cell_frac.1,
            ..Default::default()
        };
        init_sample_files(&mut query, &cfg.output_dir)?;
        let mut summary_fh = None;
        if two {
            init_sample_files(&mut control, &cfg.output_dir)?;
            summary_fh = Some(open_out(&cfg.output_dir, "summary.tab")?.0);
        }
        {
            let fh: &mut BufWriter<File> = match summary_fh.as_mut() {
                Some(f) => f,
                None => query.summary_fh.as_mut().unwrap(),
            };
            writeln!(
                fh,
                "# This file was produced by: kira-bt cnv({})",
                env!("CARGO_PKG_VERSION")
            )?;
            write!(fh, "# The command line was:\tkira-bt cnv {}", cfg.cmdline)?;
            if two {
                writeln!(
                    fh,
                    "\n#\n# RG, Regions\t[2]Chromosome\t[3]Start\t[4]End\t[5]Copy number:{}\t[6]Copy number:{}\t[7]Quality\t[8]nSites in (5)\t[9]nHETs in (5)\t[10]nSites in (6)\t[11]nHETs in(6)",
                    query.name.as_deref().unwrap_or(""),
                    control.name.as_deref().unwrap_or("")
                )?;
            } else {
                writeln!(
                    fh,
                    "\n#\n# RG, Regions\t[2]Chromosome\t[3]Start\t[4]End\t[5]Copy number:{}\t[6]Quality\t[7]nSites\t[8]nHETs",
                    query.name.as_deref().unwrap_or("")
                )?;
            }
        }
        if cfg.optimize_frac != 0.0 {
            writeln!(
                query.summary_fh.as_mut().unwrap(),
                "# CF, cell fraction estimate\t[2]Chromosome\t[3]Start\t[4]End\t[5]Cell fraction\t[6]BAF deviation"
            )?;
            if two {
                writeln!(
                    control.summary_fh.as_mut().unwrap(),
                    "# CF, cell fraction estimate\t[2]Chromosome\t[3]Start\t[4]End\t[5]Cell fraction\t[6]BAF deviation"
                )?;
                writeln!(
                    summary_fh.as_mut().unwrap(),
                    "# CF, cell fraction estimate\t[2]Chromosome\t[3]Start\t[4]End\t[5]Cell fraction:{}\t[6]Cell fraction:{}\t[7]BAF deviation:{}\t[8]BAF deviation:{}",
                    query.name.as_deref().unwrap_or(""),
                    control.name.as_deref().unwrap_or(""),
                    query.name.as_deref().unwrap_or(""),
                    control.name.as_deref().unwrap_or("")
                )?;
            }
        }
        let af_file = match &cfg.af_file {
            Some(p) => Some(AfFile::load(p)?),
            None => None,
        };
        Ok(Self {
            cfg,
            query,
            control,
            two,
            nstates,
            tprob,
            hmm,
            eprob: Vec::new(),
            sites: Vec::new(),
            nonref_afs: Vec::new(),
            f_rr: 0.0,
            f_ra: 0.0,
            f_aa: 0.0,
            n_rr: 0,
            n_ra: 0,
            n_aa: 0,
            af_file,
            summary_fh,
            prev_rid: None,
            ntot: 0,
            nused: 0,
        })
    }

    fn set_observed_prob(&self, s: &mut Sample, isite: usize) {
        let baf = s.baf[isite];
        let lrr: f32 = if self.cfg.lrr_bias > 0.0 {
            s.lrr[isite]
        } else {
            0.0
        };
        let (f_rr, f_ra, f_aa) = (self.f_rr, self.f_ra, self.f_aa);
        if baf < 0.0 {
            s.pobs[CN0] = 0.5;
            for i in 1..N_STATES {
                s.pobs[i] = (1.0 - s.pobs[CN0]) / (N_STATES as f64 - 1.0);
            }
            return;
        }
        let g = &s.gauss;
        let mut cn1_baf = norm_prob(baf, &g[CN1_PK_R]) * (f_rr as f64 + f_ra as f64 * 0.5)
            + norm_prob(baf, &g[CN1_PK_A]) * (f_aa as f64 + f_ra as f64 * 0.5);
        let mut cn2_baf = norm_prob(baf, &g[CN2_PK_RR]) * f_rr as f64
            + norm_prob(baf, &g[CN2_PK_RA]) * f_ra as f64
            + norm_prob(baf, &g[CN2_PK_AA]) * f_aa as f64;
        let mut cn3_baf = norm_prob(baf, &g[CN3_PK_RRR]) * f_rr as f64
            + norm_prob(baf, &g[CN3_PK_RRA]) * (f_ra as f64 * 0.5)
            + norm_prob(baf, &g[CN3_PK_RAA]) * (f_ra as f64 * 0.5)
            + norm_prob(baf, &g[CN3_PK_AAA]) * f_aa as f64;
        let norm = cn1_baf + cn2_baf + cn3_baf;
        cn1_baf /= norm;
        cn2_baf /= norm;
        cn3_baf /= norm;
        let lrr = lrr as f64;
        let ld = s.lrr_dev2 as f64;
        let cn1_lrr = (-(lrr + 0.45) * (lrr + 0.45) / ld).exp();
        let cn2_lrr = (-(lrr - 0.00) * (lrr - 0.00) / ld).exp();
        let cn3_lrr = (-(lrr - 0.30) * (lrr - 0.30) / ld).exp();
        let (bb, lb, e) = (self.cfg.baf_bias, self.cfg.lrr_bias, self.cfg.err_prob);
        s.pobs[CN0] = 0.0;
        s.pobs[CN1] = e + (1.0 - bb + bb * cn1_baf) * (1.0 - lb + lb * cn1_lrr);
        s.pobs[CN2] = e + (1.0 - bb + bb * cn2_baf) * (1.0 - lb + lb * cn2_lrr);
        s.pobs[CN3] = e + (1.0 - bb + bb * cn3_baf) * (1.0 - lb + lb * cn3_lrr);
    }

    fn set_emission_probs(&mut self) {
        if self.af_file.is_none() {
            self.f_rr = 0.76;
            self.f_ra = 0.14;
            self.f_aa = 0.098;
        }
        set_gauss_params(&mut self.query);
        if self.two {
            set_gauss_params(&mut self.control);
        }
        let n = self.sites.len();
        self.eprob.resize(n * self.nstates, 0.0);
        for i in 0..n {
            if self.af_file.is_some() {
                let af = self.nonref_afs[i];
                self.f_rr = (1.0 - af) * (1.0 - af);
                self.f_ra = 2.0 * af * (1.0 - af);
                self.f_aa = af * af;
            }
            let mut q = std::mem::take(&mut self.query);
            self.set_observed_prob(&mut q, i);
            self.query = q;
            if self.two {
                let mut c = std::mem::take(&mut self.control);
                self.set_observed_prob(&mut c, i);
                self.control = c;
                for a in 0..N_STATES {
                    for b in 0..N_STATES {
                        self.eprob[self.nstates * i + a * N_STATES + b] =
                            self.query.pobs[a] * self.control.pobs[b];
                    }
                }
            } else {
                for a in 0..N_STATES {
                    self.eprob[self.nstates * i + a] = self.query.pobs[a];
                }
            }
        }
    }

    /// `update_sample_args`: re-estimate the aberrant cell fraction from the
    /// CN3 posteriors; returns true when converged.
    fn update_sample_args(&mut self, is_control: bool) -> Result<bool> {
        let nstates = self.nstates;
        let fwd = self.hmm.fwd_bwd_prob().to_vec();
        let nsites = self.sites.len();
        let s = if is_control {
            &self.control
        } else {
            &self.query
        };
        let mut tmpf: Vec<f32> = Vec::with_capacity(nsites);
        for i in 0..nsites {
            let mut baf = s.baf[i];
            if baf as f64 > 4.0 / 5.0 {
                continue;
            }
            if baf > 0.5 {
                baf = 1.0 - baf;
            }
            if (baf as f64) < 1.0 / 5.0 {
                continue;
            }
            let probs = &fwd[i * nstates..(i + 1) * nstates];
            let mut prob_cn3 = 0.0;
            if !self.two {
                prob_cn3 = probs[CN3];
            } else if !is_control {
                for j in 0..N_STATES {
                    prob_cn3 += probs[CN3 * N_STATES + j];
                }
            } else {
                for j in 0..N_STATES {
                    prob_cn3 += probs[CN3 + j * N_STATES];
                }
            }
            tmpf.push(prob_cn3 as f32);
        }
        smooth_data(&mut tmpf, 50);
        let mut k = 0;
        let (mut mean_cn3, mut norm_cn3) = (0.0f64, 0.0f64);
        let (mut baf_aa_dev2, mut norm_baf_aa_dev2) = (0.0f64, 0.0f64);
        for i in 0..nsites {
            let mut baf = s.baf[i];
            if baf as f64 > 4.0 / 5.0 {
                baf_aa_dev2 += (1.0 - baf as f64) * (1.0 - baf as f64);
                norm_baf_aa_dev2 += 1.0;
                continue;
            }
            if baf > 0.5 {
                baf = 1.0 - baf;
            }
            if (baf as f64) < 1.0 / 5.0 {
                continue;
            }
            let prob_cn3 = tmpf[k] as f64;
            k += 1;
            mean_cn3 += prob_cn3 * baf as f64;
            norm_cn3 += prob_cn3;
        }
        let s = if is_control {
            &mut self.control
        } else {
            &mut self.query
        };
        if norm_cn3 == 0.0 {
            s.cell_frac = 1.0;
            return Ok(true);
        }
        mean_cn3 /= norm_cn3;
        let mut k = 0;
        let mut baf_dev2 = 0.0f64;
        for i in 0..nsites {
            let mut baf = s.baf[i];
            if baf > 0.5 {
                baf = 1.0 - baf;
            }
            if (baf as f64) < 1.0 / 5.0 {
                continue;
            }
            let prob_cn3 = tmpf[k] as f64;
            k += 1;
            baf_dev2 += prob_cn3 * (baf as f64 - mean_cn3) * (baf as f64 - mean_cn3);
        }
        baf_dev2 /= norm_cn3;
        baf_aa_dev2 /= norm_baf_aa_dev2;
        if baf_dev2 < baf_aa_dev2 {
            baf_dev2 = baf_aa_dev2;
        }
        let max_mean_cn3 = 0.5 - baf_dev2.sqrt() * 1.644854;
        if max_mean_cn3 <= 0.0 {
            bail!("cnv: BAF deviation too large to estimate the aberrant cell fraction");
        }
        let new_frac = 1.0 / mean_cn3 - 2.0;
        if mean_cn3 > max_mean_cn3 || new_frac < self.cfg.optimize_frac {
            s.cell_frac = 1.0;
            return Ok(true);
        }
        let new_frac = new_frac.min(1.0);
        let converged = (new_frac - s.cell_frac as f64).abs() < 1e-1;
        let mut bd = baf_dev2;
        if bd > 3.0 * s.baf_dev2_dflt as f64 {
            bd = 3.0 * s.baf_dev2_dflt as f64;
        } else if bd < 0.5 * s.baf_dev2_dflt as f64 {
            bd = 0.5 * s.baf_dev2_dflt as f64;
        }
        s.cell_frac = new_frac as f32;
        s.baf_dev2 = bd as f32;
        Ok(converged)
    }

    /// `update_args`: true when another iteration is needed.
    fn update_args(&mut self) -> Result<bool> {
        let mut converged = self.update_sample_args(false)? as i32;
        if self.two {
            converged += self.update_sample_args(true)? as i32;
            return Ok(converged != 2);
        }
        Ok(converged == 0)
    }

    fn write_rg(
        &mut self,
        chr: &str,
        start_pos: u32,
        end_pos: u32,
        start_cn: usize,
        qual: f64,
        counts: [usize; 4],
    ) -> Result<()> {
        let [smpl_ntot, smpl_nhet, ctrl_ntot, ctrl_nhet] = counts;
        let cq = copy_number_state(self.two, start_cn, 0);
        writeln!(
            self.query.summary_fh.as_mut().unwrap(),
            "RG\t{}\t{}\t{}\t{}\t{:.1}\t{}\t{}",
            chr,
            start_pos,
            end_pos,
            cq,
            qual,
            smpl_ntot,
            smpl_nhet
        )?;
        if self.two {
            let cc = copy_number_state(self.two, start_cn, 1);
            writeln!(
                self.control.summary_fh.as_mut().unwrap(),
                "RG\t{}\t{}\t{}\t{}\t{:.1}\t{}\t{}",
                chr,
                start_pos,
                end_pos,
                cc,
                qual,
                ctrl_ntot,
                ctrl_nhet
            )?;
            writeln!(
                self.summary_fh.as_mut().unwrap(),
                "RG\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{}\t{}\t{}\t{}",
                chr,
                start_pos,
                end_pos,
                cq,
                cc,
                qual,
                smpl_ntot,
                smpl_nhet,
                ctrl_ntot,
                ctrl_nhet
            )?;
        }
        Ok(())
    }

    /// `cnv_flush_viterbi`: run the HMM on the buffered chromosome.
    fn flush_viterbi(&mut self) -> Result<()> {
        let nsites = self.sites.len();
        if nsites == 0 {
            return Ok(());
        }
        let chr = self.prev_rid.clone().unwrap_or_default();
        let tprob = self.tprob.clone();
        self.hmm.set_tprob(&tprob, 10000);
        if self.cfg.lrr_bias > 0.0 {
            smooth_data(&mut self.query.lrr, self.cfg.lrr_smooth_win);
            if self.two {
                smooth_data(&mut self.control.lrr, self.cfg.lrr_smooth_win);
            }
        }
        self.query.cell_frac = self.query.cell_frac_dflt;
        self.control.cell_frac = self.control.cell_frac_dflt;
        self.query.baf_dev2 = self.query.baf_dev2_dflt;
        self.control.baf_dev2 = self.control.baf_dev2_dflt;
        set_gauss_params(&mut self.query);
        if self.two {
            set_gauss_params(&mut self.control);
        }
        if self.cfg.optimize_frac != 0.0 {
            let mut niter = 0;
            eprintln!("Attempting to estimate the fraction of aberrant cells (chr {chr}):");
            loop {
                eprint!(
                    "\t.. {:.6} {:.6}",
                    self.query.cell_frac, self.query.baf_dev2
                );
                if self.two {
                    eprint!(
                        "\t.. {:.6} {:.6}",
                        self.control.cell_frac, self.control.baf_dev2
                    );
                }
                eprintln!();
                self.set_emission_probs();
                let (e, s) = (
                    std::mem::take(&mut self.eprob),
                    std::mem::take(&mut self.sites),
                );
                self.hmm.run_fwd_bwd(&e, &s);
                self.eprob = e;
                self.sites = s;
                let again = self.update_args()?;
                niter += 1;
                if !(again && niter < 20) {
                    break;
                }
            }
            if niter >= 20 {
                self.query.cell_frac = self.query.cell_frac_dflt;
                self.control.cell_frac = self.control.cell_frac_dflt;
                self.query.baf_dev2 = self.query.baf_dev2_dflt;
                self.control.baf_dev2 = self.control.baf_dev2_dflt;
                set_gauss_params(&mut self.query);
                if self.two {
                    set_gauss_params(&mut self.control);
                }
            }
            eprint!(
                "\t.. {:.6} {:.6}",
                self.query.cell_frac, self.query.baf_dev2
            );
            if self.two {
                eprint!(
                    "\t.. {:.6} {:.6}",
                    self.control.cell_frac, self.control.baf_dev2
                );
            }
            eprintln!();
            let (first, last) = (self.sites[0] + 1, self.sites[nsites - 1] + 1);
            writeln!(
                self.query.summary_fh.as_mut().unwrap(),
                "CF\t{}\t{}\t{}\t{:.2}\t{:.6}",
                chr,
                first,
                last,
                self.query.cell_frac,
                (self.query.baf_dev2 as f64).sqrt()
            )?;
            if self.two {
                writeln!(
                    self.control.summary_fh.as_mut().unwrap(),
                    "CF\t{}\t{}\t{}\t{:.2}\t{:.6}",
                    chr,
                    first,
                    last,
                    self.control.cell_frac,
                    (self.control.baf_dev2 as f64).sqrt()
                )?;
                writeln!(
                    self.summary_fh.as_mut().unwrap(),
                    "CF\t{}\t{}\t{}\t{:.2}\t{:.2}\t{:.6}\t{:.6}",
                    chr,
                    first,
                    last,
                    self.query.cell_frac,
                    self.control.cell_frac,
                    (self.query.baf_dev2 as f64).sqrt(),
                    (self.control.baf_dev2 as f64).sqrt()
                )?;
            }
        }
        self.set_emission_probs();
        let (eprob, sites) = (
            std::mem::take(&mut self.eprob),
            std::mem::take(&mut self.sites),
        );
        while self.cfg.baum_welch_th != 0.0 {
            let ns = self.nstates;
            let ori_ii = avg_ii_prob(ns, self.hmm.tprob());
            self.hmm.run_baum_welch(&eprob, &sites);
            let new_ii = avg_ii_prob(ns, self.hmm.tprob());
            eprintln!(
                "{}\t{}\t{}",
                crate::util::fmt_e(ori_ii, 6),
                crate::util::fmt_e(new_ii, 6),
                crate::util::fmt_e(new_ii - ori_ii, 6)
            );
            let t = init_tprob_matrix(ns, 1.0 - new_ii, self.cfg.same_prob)?;
            self.hmm.set_tprob(&t, 10000);
            if (new_ii - ori_ii).abs() < self.cfg.baum_welch_th {
                let arr = self.hmm.tprob().to_vec();
                let mut out = std::io::stdout().lock();
                for i in 0..ns {
                    for j in 0..ns {
                        write!(out, " {:.15}", arr[ns * j + i])?;
                    }
                    writeln!(out)?;
                }
                break;
            }
        }
        self.hmm.run_viterbi(&eprob, &sites);
        self.hmm.run_fwd_bwd(&eprob, &sites);
        self.eprob = eprob;
        self.sites = sites;
        let ns = self.nstates;
        let fwd = self.hmm.fwd_bwd_prob().to_vec();
        let mut qual = 0.0f64;
        let mut start_cn = self.hmm.vstate(0) as usize;
        let mut start_pos = self.sites[0];
        let mut istart_pos = 0usize;
        let mut counts = [0usize; 4];
        for isite in 0..nsites {
            let state = self.hmm.vstate(isite) as usize;
            // bcftools reads row isite, the posterior of the previous site.
            let pval = &fwd[isite * ns..(isite + 1) * ns];
            qual += pval[start_cn];
            let pos1 = self.sites[isite] + 1;
            {
                let mut line = format!(
                    "{}\t{}\t{}",
                    chr,
                    pos1,
                    copy_number_state(self.two, state, 0)
                );
                if !self.two {
                    for v in pval.iter().take(ns) {
                        line.push_str(&format!("\t{:.6}", v));
                    }
                } else {
                    for i in 0..N_STATES {
                        let mut sum = 0.0;
                        for j in 0..N_STATES {
                            sum += pval[i * N_STATES + j];
                        }
                        line.push_str(&format!("\t{:.6}", sum));
                    }
                }
                writeln!(self.query.cn_fh.as_mut().unwrap(), "{line}")?;
                let b = self.query.baf[isite];
                if b >= 0.0 {
                    if b > 0.25 && b < 0.75 {
                        counts[1] += 1;
                    }
                    counts[0] += 1;
                }
            }
            if self.two {
                let mut line = format!(
                    "{}\t{}\t{}",
                    chr,
                    pos1,
                    copy_number_state(self.two, state, 1)
                );
                for i in 0..N_STATES {
                    let mut sum = 0.0;
                    for j in 0..N_STATES {
                        sum += pval[i + N_STATES * j];
                    }
                    line.push_str(&format!("\t{:.6}", sum));
                }
                writeln!(self.control.cn_fh.as_mut().unwrap(), "{line}")?;
                let b = self.control.baf[isite];
                if b >= 0.0 {
                    if b > 0.25 && b < 0.75 {
                        counts[3] += 1;
                    }
                    counts[2] += 1;
                }
            }
            if start_cn != state {
                let q = phred_score(1.0 - qual / (isite - istart_pos) as f64);
                self.write_rg(&chr, start_pos + 1, self.sites[isite], start_cn, q, counts)?;
                istart_pos = isite;
                start_pos = self.sites[isite];
                start_cn = state;
                qual = 0.0;
                counts = [0; 4];
            }
        }
        let q = phred_score(1.0 - qual / (nsites - istart_pos) as f64);
        self.write_rg(
            &chr,
            start_pos + 1,
            self.sites[nsites - 1] + 1,
            start_cn,
            q,
            counts,
        )?;
        Ok(())
    }

    /// `cnv_next_line`.
    fn next_line(&mut self, rec: &VcfRecord) -> Result<()> {
        if self.prev_rid.as_deref() != Some(rec.chrom.as_str()) {
            self.flush_viterbi()?;
            self.prev_rid = Some(rec.chrom.clone());
            self.sites.clear();
            self.eprob.clear();
            self.query.baf.clear();
            self.query.lrr.clear();
            self.control.baf.clear();
            self.control.lrr.clear();
            self.nonref_afs.clear();
            self.n_rr = 0;
            self.n_aa = 0;
            self.n_ra = 0;
        }
        self.ntot += 1;
        let Some(fmt) = rec.format.as_deref() else {
            return Ok(());
        };
        let keys: Vec<&str> = fmt.split(':').collect();
        let Some(ibaf) = keys.iter().position(|k| *k == "BAF") else {
            return Ok(());
        };
        let ilrr = keys.iter().position(|k| *k == "LRR");
        if self.cfg.lrr_bias > 0.0 && ilrr.is_none() {
            return Ok(());
        }
        let ilrr = if self.cfg.lrr_bias > 0.0 { ilrr } else { None };
        let (baf1, lrr1, ok1) = parse_lrr_baf(rec, ibaf, ilrr, self.query.idx);
        let (baf2, lrr2, ok2) = parse_lrr_baf(rec, ibaf, ilrr, self.control.idx);
        if !ok1 && !ok2 {
            return Ok(());
        }
        let pos = rec.pos.saturating_sub(1);
        self.sites.push(pos);
        self.query.lrr.push(lrr1);
        self.query.baf.push(baf1);
        if let Some(af) = &self.af_file {
            let mut alleles = vec![rec.ref_allele.clone()];
            if rec.alt != "." && !rec.alt.is_empty() {
                alleles.extend(rec.alt.split(',').map(str::to_string));
            }
            let v = af
                .read_af(&rec.chrom, rec.pos, &alleles)?
                .map(|x| x as f32)
                .unwrap_or(0.1);
            self.nonref_afs.push(v);
        }
        if self.two {
            self.control.lrr.push(lrr2);
            self.control.baf.push(baf2);
            if baf2 >= 0.0 {
                writeln!(
                    self.control.dat_fh.as_mut().unwrap(),
                    "{}\t{}\t{:.3}\t{:.3}",
                    rec.chrom,
                    rec.pos,
                    baf2,
                    lrr2
                )?;
            }
        }
        if baf1 >= 0.0 {
            writeln!(
                self.query.dat_fh.as_mut().unwrap(),
                "{}\t{}\t{:.3}\t{:.3}",
                rec.chrom,
                rec.pos,
                baf1,
                lrr1
            )?;
            if (baf1 as f64) < 1.0 / 5.0 {
                self.n_rr += 1;
            } else if baf1 as f64 > 4.0 / 5.0 {
                self.n_aa += 1;
            } else {
                self.n_ra += 1;
            }
        }
        self.nused += 1;
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        self.flush_viterbi()?;
        for fh in [
            &mut self.query.dat_fh,
            &mut self.query.cn_fh,
            &mut self.query.summary_fh,
            &mut self.control.dat_fh,
            &mut self.control.cn_fh,
            &mut self.control.summary_fh,
            &mut self.summary_fh,
        ]
        .into_iter()
        .flatten()
        {
            fh.flush()?;
        }
        self.write_plot_scripts()?;
        eprintln!(
            "Number of lines: total/processed: {}/{}",
            self.ntot, self.nused
        );
        Ok(())
    }

    /// `create_plots`: the matplotlib scripts bcftools leaves next to the
    /// tables; python is invoked only with `-p`.
    fn write_plot_scripts(&self) -> Result<()> {
        let dir = self.cfg.output_dir.display().to_string();
        let qn = self.query.name.clone().unwrap_or_default();
        let script = if !self.two {
            let p = self.cfg.output_dir.join(format!("plot.{qn}.py"));
            std::fs::write(
                &p,
                plot_single(&self.query.dat_fname, &self.query.cn_fname, &qn, &dir),
            )?;
            p
        } else {
            let cn = self.control.name.clone().unwrap_or_default();
            let p = self.cfg.output_dir.join(format!("plot.{cn}.{qn}.py"));
            std::fs::write(
                &p,
                plot_pair(
                    &cn,
                    &qn,
                    &dir,
                    &self.control.dat_fname,
                    &self.query.dat_fname,
                    &self.control.cn_fname,
                    &self.query.cn_fname,
                ),
            )?;
            p
        };
        if self.cfg.plot_th <= 100.0 {
            let status = std::process::Command::new("python")
                .arg(&script)
                .arg("-p")
                .arg(format!("{:.6}", self.cfg.plot_th))
                .status();
            match status {
                Ok(s) if s.success() => {}
                Ok(s) => eprintln!(
                    "The command returned non-zero status {}: python {} -p {:.6}",
                    s,
                    script.display(),
                    self.cfg.plot_th
                ),
                Err(e) => eprintln!("The command failed: python {} .. {}", script.display(), e),
            }
        }
        Ok(())
    }
}

fn plot_single(dat: &str, cn: &str, name: &str, dir: &str) -> String {
    format!(
        "import matplotlib as mpl\nmpl.use('Agg')\nimport matplotlib.pyplot as plt\nimport csv\nimport numpy as np\ncsv.register_dialect('tab', delimiter='\\t', quoting=csv.QUOTE_NONE)\n\ndat = {{}}\nwith open('{dat}', 'r') as f:\n    reader = csv.reader(f, 'tab')\n    for row in reader:\n        chr = row[0]\n        if chr[0]=='#': continue\n        if chr not in dat: dat[chr] = []\n        dat[chr].append([int(row[1]), float(row[2]), float(row[3])])\n\ncnv = {{}}\nwith open('{cn}', 'r') as f:\n    reader = csv.reader(f, 'tab')\n    for row in reader:\n        chr = row[0]\n        if chr[0]=='#': continue\n        if chr not in cnv: cnv[chr] = []\n        row[2] = int(row[2]) + 0.5\n        row[1] = int(row[1])\n        cnv[chr].append(row[1:])\n\nfor chr in dat:\n    fig,(ax1, ax2, ax3) = plt.subplots(3,1,figsize=(10,8),sharex=True)\n    ax1.plot([x[0] for x in dat[chr]],[x[2] for x in dat[chr]], '.', ms=3)\n    ax2.plot([x[0] for x in dat[chr]],[x[1] for x in dat[chr]], '.', ms=3)\n    cn_dat = cnv[chr]\n    xgrid = [float(x[0]) for x in cn_dat]\n    ygrid = np.linspace(0,5,6)\n    xgrid, ygrid = np.meshgrid(xgrid, ygrid)\n    heat = np.zeros_like(xgrid)\n    for x in range(len(heat[0])-1):\n       heat[0][x] = cn_dat[x][2]\n       heat[1][x] = cn_dat[x][3]\n       heat[2][x] = cn_dat[x][4]\n       heat[3][x] = cn_dat[x][5]\n    mesh = ax3.pcolormesh(xgrid, ygrid, heat, cmap='bwr_r', shading='auto', alpha=0)\n    mesh.set_clim(vmin=-1,vmax=1)\n    ax3.plot([x[0] for x in cn_dat],[x[1] for x in cn_dat],'.-',ms=3,color='black')\n    fig.suptitle('{name} (chr '+chr+')')\n    ax1.set_title('Log-R intensities Ratio',fontsize=10)\n    ax2.set_title('B-Allele Frequency',fontsize=10)\n    ax3.set_title('Copy Number Variation',fontsize=10)\n    ax1.set_ylabel('LRR')\n    ax2.set_ylabel('BAF')\n    ax3.set_ylabel('CN')\n    ax3.set_xlabel('Coordinate (chrom '+chr+')',fontsize=10)\n    ax3.set_ylim(-0.1,4.1)\n    ax3.set_yticks([0.5,1.5,2.5,3.5])\n    ax3.set_yticklabels(['CN0','CN1','CN2','CN3'])\n    plt.subplots_adjust(left=0.08,right=0.95,bottom=0.08,top=0.92)\n    plt.savefig('{dir}/plot.{name}.chr'+chr+'.png')\n    plt.close()\n\n"
    )
}

fn plot_pair(
    control: &str,
    query: &str,
    dir: &str,
    cdat: &str,
    qdat: &str,
    ccn: &str,
    qcn: &str,
) -> String {
    format!(
        "import matplotlib as mpl\nmpl.use('Agg')\nimport matplotlib.pyplot as plt\nimport csv,argparse\nimport numpy as np\ncsv.register_dialect('tab', delimiter='\\t', quoting=csv.QUOTE_NONE)\n\ncontrol_sample = '{control}'\nquery_sample   = '{query}'\n\nparser = argparse.ArgumentParser()\nparser.add_argument('-p', '--plot-threshold', type=float)\nparser.add_argument('-c', '--chromosome')\nargs = parser.parse_args()\nif args.plot_threshold==None: args.plot_threshold = 0\n\ndef chroms_to_plot(th):\n   dat = {{}}\n   with open('{dir}/summary.tab', 'r') as f:\n       reader = csv.reader(f, 'tab')\n       for row in reader:\n           if row[0]!='RG': continue\n           chr   = row[1]\n           start = row[2]\n           end   = row[3]\n           qual  = float(row[6])\n           if row[4]==row[5] and args.plot_threshold!=0: continue\n           if chr not in dat: dat[chr] = 0.0\n           if qual > dat[chr]: dat[chr] = qual\n   out = {{}}\n   for chr in dat:\n       if (chr not in dat) or dat[chr]<th: continue\n       out[chr] = 1\n   return out\nif args.chromosome!=None:\n   plot_chroms = {{ args.chromosome:1 }}\nelse:\n   plot_chroms = chroms_to_plot(args.plot_threshold)\n\ndef read_dat(file,dat,plot_chr):\n   with open(file, 'r') as f:\n       reader = csv.reader(f, 'tab')\n       for row in reader:\n           chr = row[0]\n           if chr != plot_chr: continue\n           dat.append([int(row[1]), float(row[2]), float(row[3])])\ndef read_cnv(file,cnv,plot_chr):\n   with open(file, 'r') as f:\n       reader = csv.reader(f, 'tab')\n       for row in reader:\n           chr = row[0]\n           if chr != plot_chr: continue\n           row[2] = int(row[2]) + 0.5\n           row[1] = int(row[1])\n           cnv.append(row[1:])\ndef find_diffs(a,b):\n    out = []\n    diff = []\n    for i in range(len(a)):\n        if a[i][1]!=b[i][1]:\n            if i>0: diff.append([b[i-1][0],b[i-1][1],a[i-1][1]])\n            diff.append([b[i][0],b[i][1],a[i][1]])\n        elif len(diff):\n            diff.append([b[i][0],b[i][1],a[i][1]])\n            out.append(diff)\n            diff = []\n    if len(diff): out.append(diff)\n    return out\n\nfor chr in sorted(plot_chroms.keys()):\n    control_dat = []\n    control_cnv = []\n    query_dat   = []\n    query_cnv   = []\n    read_dat('{cdat}',control_dat,chr)\n    read_dat('{qdat}',query_dat,chr)\n    read_cnv('{ccn}',control_cnv,chr)\n    read_cnv('{qcn}',query_cnv,chr)\n\n    fig,(ax1,ax2,ax3,ax4,ax5,ax6) = plt.subplots(6,1,figsize=(10,8),sharex=True)\n    ax1.plot([x[0] for x in control_dat],[x[2] for x in control_dat], '.', ms=3,color='red')\n    ax2.plot([x[0] for x in control_dat],[x[1] for x in control_dat], '.', ms=3,color='red')\n    cn_dat = control_cnv\n    xgrid = [float(x[0]) for x in cn_dat]\n    ygrid = np.linspace(0,5,6)\n    xgrid, ygrid = np.meshgrid(xgrid, ygrid)\n    heat = np.zeros_like(xgrid)\n    for x in range(len(heat[0])-1):\n       heat[0][x] = cn_dat[x][2]\n       heat[1][x] = cn_dat[x][3]\n       heat[2][x] = cn_dat[x][4]\n       heat[3][x] = cn_dat[x][5]\n    mesh = ax3.pcolormesh(xgrid, ygrid, heat, cmap='bwr', shading='auto', alpha=0)\n    mesh.set_clim(vmin=-1,vmax=1)\n    ax3.plot([x[0] for x in cn_dat],[x[1] for x in cn_dat],'-',ms=3,color='black',lw=1.7)\n\n    ax6.plot([x[0] for x in query_dat],[x[2] for x in query_dat], '.', ms=3)\n    ax5.plot([x[0] for x in query_dat],[x[1] for x in query_dat], '.', ms=3)\n    cn_dat = query_cnv\n    xgrid = [float(x[0]) for x in cn_dat]\n    ygrid = np.linspace(0,5,6)\n    xgrid, ygrid = np.meshgrid(xgrid, ygrid)\n    heat = np.zeros_like(xgrid)\n    for x in range(len(heat[0])-1):\n       heat[0][x] = cn_dat[x][2]\n       heat[1][x] = cn_dat[x][3]\n       heat[2][x] = cn_dat[x][4]\n       heat[3][x] = cn_dat[x][5]\n    mesh = ax4.pcolormesh(xgrid, ygrid, heat, cmap='bwr_r')\n    mesh.set_clim(vmin=-1,vmax=1)\n    ax4.plot([x[0] for x in cn_dat],[x[1] for x in cn_dat],'-',ms=3,color='black',lw=1.7)\n    ax3.annotate(control_sample, xy=(0.02,0.1), xycoords='axes fraction', color='red',fontsize=12, va='bottom',ha='left')\n    ax4.annotate(query_sample, xy=(0.02,0.9), xycoords='axes fraction', color='blue',fontsize=12, va='top',ha='left')\n\n    diffs = find_diffs(control_cnv,query_cnv)\n    for diff in diffs:\n        ax3.plot([x[0] for x in diff],[x[1] for x in diff],'-',ms=3,color='blue',lw=1.7)\n        ax4.plot([x[0] for x in diff],[x[2] for x in diff],'-',ms=3,color='red',lw=1.7)\n\n    fig.suptitle('chr '+chr+', '+control_sample+' vs '+query_sample)\n    ax1.tick_params(axis='both', labelsize=8)\n    ax2.tick_params(axis='both', labelsize=8)\n    ax3.tick_params(axis='both', labelsize=8)\n    ax4.tick_params(axis='both', labelsize=8)\n    ax5.tick_params(axis='both', labelsize=8)\n    ax6.tick_params(axis='both', labelsize=8)\n    ax6.set_xlabel('Coordinate (chrom '+chr+')',fontsize=8)\n    ax1.set_ylabel('LRR')\n    ax2.set_ylabel('BAF')\n    ax3.set_ylabel('CN')\n    ax6.set_ylabel('LRR')\n    ax5.set_ylabel('BAF')\n    ax4.set_ylabel('CN')\n    ax3.set_ylim(-0.1,4.1)\n    ax3.set_yticks([0.5,1.5,2.5,3.5])\n    ax3.set_yticklabels(['CN0','CN1','CN2','CN3'])\n    ax4.set_ylim(-0.1,4.1)\n    ax4.set_yticks([0.5,1.5,2.5,3.5])\n    ax4.set_yticklabels(['CN0','CN1','CN2','CN3'])\n    plt.subplots_adjust(left=0.08,right=0.95,bottom=0.08,top=0.92,hspace=0)\n    plt.savefig('{dir}/plot.{control}.{query}.chr'+chr+'.png')\n    plt.close()\n\n"
    )
}

/// Run `cnv` with a parsed configuration.
pub fn run(cfg: CnvConfig) -> Result<()> {
    let input = cfg.input.clone();
    let region_filter = match &cfg.regions {
        Some(s) if cfg.regions_is_file => Some(RegionSet::from_file(s)?),
        Some(s) => Some(RegionSet::from_cli(s)?),
        None => None,
    };
    let target_filter = match &cfg.targets {
        Some(s) if cfg.targets_is_file => Some(RegionSet::from_file(s)?),
        Some(s) => Some(RegionSet::from_cli(s.trim_start_matches('^'))?),
        None => None,
    };
    let targets_inverse = cfg
        .targets
        .as_deref()
        .is_some_and(|s| !cfg.targets_is_file && s.starts_with('^'));
    let (regions_overlap, targets_overlap) = (cfg.regions_overlap, cfg.targets_overlap);
    let mut reader = UnifiedVcfReader::open(&input).with_context(|| {
        format!(
            "Failed to read from {}",
            if input == Path::new("-") {
                "standard input".to_string()
            } else {
                input.display().to_string()
            }
        )
    })?;
    let headers = reader.header()?.to_vec();
    let mut cnv = Cnv::new(cfg, &headers)?;
    let mut contigs = ContigDict::new();
    let handle = |line: &str, cnv: &mut Cnv, contigs: &mut ContigDict| -> Result<()> {
        if line.is_empty() || line.as_bytes()[0] == b'#' {
            return Ok(());
        }
        if let Some(tf) = &target_filter {
            if tf.line_passes_mode(line, targets_overlap) == targets_inverse {
                return Ok(());
            }
        }
        let Some(rec) = parse_vcf_record(line, 0, contigs).map_err(|e| anyhow::anyhow!("{e}"))?
        else {
            return Ok(());
        };
        if let Some(af) = &cnv.af_file {
            if !af.has_site(&rec.chrom, rec.pos) {
                return Ok(());
            }
        }
        cnv.next_line(&rec)
    };
    let indexed = region_filter
        .as_ref()
        .filter(|_| input != Path::new("-") && crate::csi::find_index_for(&input).is_some());
    if let Some(rf) = indexed {
        rf.stream_with_index(&input, regions_overlap, |line| {
            handle(line, &mut cnv, &mut contigs)
        })?;
    } else {
        while let Some(line) = reader.read_line()? {
            if let Some(rf) = &region_filter {
                if line.is_empty()
                    || line.as_bytes()[0] == b'#'
                    || !rf.line_passes_mode(&line, regions_overlap)
                {
                    continue;
                }
            }
            handle(&line, &mut cnv, &mut contigs)?;
        }
    }
    cnv.finish()
}

#[cfg(test)]
#[path = "../../tests/unit/cnv.rs"]
mod tests;
