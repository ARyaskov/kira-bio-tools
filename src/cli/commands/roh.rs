use anyhow::Result;
use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use crate::VcfReader;
use crate::cli::args::RohArgs;
use crate::roh::{self, GtClass, RohOpts, RohSite, State};

/// `kira-bt roh`: runs-of-homozygosity detection via the two-state AZ/HW HMM
/// in [`crate::roh`]. Viterbi gives the per-site call, forward-backward gives the
/// fwd-bwd phred quality. Output mirrors `bcftools roh` (`ST` per-site, `RG`
/// per autozygous region).
pub fn cmd_roh(args: RohArgs) -> Result<()> {
    let mut run = RohRun::from_args(&args);
    run.apply_passthrough(&args.passthrough);

    let mut reader = VcfReader::open(&args.input)?;
    let headers = reader.header()?;
    let all_samples = sample_names(&headers);
    let sel: Vec<(usize, String)> = match &run.samples {
        Some(req) => all_samples
            .iter()
            .enumerate()
            .filter(|(_, n)| req.iter().any(|r| r == *n))
            .map(|(i, n)| (i, n.clone()))
            .collect(),
        None => all_samples
            .iter()
            .enumerate()
            .map(|(i, n)| (i, n.clone()))
            .collect(),
    };
    let sel = if sel.is_empty() {
        vec![(0usize, "sample".to_string())]
    } else {
        sel
    };

    let af_file = match &run.af_file {
        Some(p) => Some(load_af_file(p)?),
        None => None,
    };
    let gmap = match &run.genetic_map {
        Some(p) => Some(roh::parse_genetic_map(p)?),
        None => None,
    };

    let mut tracks: Vec<Vec<RohSite>> = vec![Vec::new(); sel.len()];
    while let Some(rec) = reader.next_record()? {
        if let Some(r) = &run.region {
            if !r.contains(&rec.chrom, rec.pos) {
                continue;
            }
        }
        if !run.include_noalt && (rec.alt == "." || rec.alt.is_empty() || rec.alt == rec.ref_allele) {
            continue;
        }
        if run.opts.skip_indels && is_indel(&rec.ref_allele, &rec.alt) {
            continue;
        }
        let af = resolve_af(&rec.chrom, rec.pos, &rec.info, &run, af_file.as_ref());
        let gpos = gmap
            .as_ref()
            .and_then(|m| m.get(&rec.chrom))
            .and_then(|v| interp_cm(v, rec.pos));
        let gts = parse_all_gts(&rec.format, &rec.samples);
        for (ti, (col, _)) in sel.iter().enumerate() {
            let gt = gts.get(*col).copied().unwrap_or(GtClass::Missing);
            if run.opts.ignore_homref && gt == GtClass::HomRef {
                continue;
            }
            tracks[ti].push(RohSite {
                chrom: rec.chrom.clone(),
                pos: rec.pos,
                genetic_pos: gpos,
                gt,
                af,
            });
        }
    }

    if run.estimate_af {
        for track in tracks.iter_mut() {
            let gts: Vec<GtClass> = track.iter().map(|s| s.gt).collect();
            let afs = roh::estimate_af(&gts, 100);
            for (s, a) in track.iter_mut().zip(afs) {
                s.af = a;
            }
        }
    }

    let mut out: BufWriter<Box<dyn Write>> = match &args.output {
        Some(p) => BufWriter::new(Box::new(std::fs::File::create(p)?)),
        None => BufWriter::new(Box::new(std::io::stdout())),
    };
    for (ti, (_, name)) in sel.iter().enumerate() {
        run_sample(name, &tracks[ti], &run, &mut out)?;
    }
    out.flush()?;
    Ok(())
}

/// Resolved run configuration: structured [`RohArgs`] overlaid with the
/// bcftools-style tokens passed after `--`.
struct RohRun {
    opts: RohOpts,
    af_tag: String,
    af_file: Option<PathBuf>,
    af_dflt: f64,
    estimate_af: bool,
    include_noalt: bool,
    region: Option<Region>,
    samples: Option<Vec<String>>,
    genetic_map: Option<PathBuf>,
    output_st: bool,
    output_rg: bool,
}

impl RohRun {
    fn from_args(args: &RohArgs) -> Self {
        Self {
            opts: RohOpts {
                hw_to_az: args.hw_to_az,
                az_to_hw: args.az_to_hw,
                rec_rate: args.rec_rate,
                af_dflt: args.af_dflt,
                ignore_homref: args.ignore_homref,
                skip_indels: args.skip_indels,
                viterbi_training: if args.viterbi_training > 0.0 { 10 } else { 0 },
            },
            af_tag: args.af_tag.clone(),
            af_file: args.af_file.clone(),
            af_dflt: args.af_dflt,
            estimate_af: args.estimate_af.is_some(),
            include_noalt: args.include_noalt,
            region: args.regions.as_deref().and_then(parse_region),
            samples: args
                .samples
                .as_ref()
                .map(|s| s.split(',').map(|x| x.trim().to_string()).collect()),
            genetic_map: args.genetic_map.clone(),
            output_st: args.output_type.contains('s'),
            output_rg: args.output_type.contains('r'),
        }
    }

    fn apply_passthrough(&mut self, toks: &[String]) {
        let mut i = 0usize;
        while i < toks.len() {
            let t = toks[i].as_str();
            macro_rules! val {
                () => {{
                    i += 1;
                    toks.get(i).cloned()
                }};
            }
            match t {
                "-O" | "--output-type" => {
                    if let Some(v) = val!() {
                        self.set_output(&v);
                    }
                }
                _ if t.starts_with("-O") => self.set_output(&t[2..]),
                "--AF-tag" => {
                    if let Some(v) = val!() {
                        self.af_tag = v;
                    }
                }
                "--AF-dflt" => {
                    if let Some(v) = val!() {
                        if let Ok(x) = v.parse() {
                            self.af_dflt = x;
                            self.opts.af_dflt = x;
                        }
                    }
                }
                "--AF-file" => {
                    if let Some(v) = val!() {
                        self.af_file = Some(v.into());
                    }
                }
                "-e" | "-E" | "--estimate-AF" => {
                    let _ = val!();
                    self.estimate_af = true;
                }
                "-a" | "--hw-to-az" => {
                    if let Some(v) = val!() {
                        if let Ok(x) = v.parse() {
                            self.opts.hw_to_az = x;
                        }
                    }
                }
                "-H" | "--az-to-hw" => {
                    if let Some(v) = val!() {
                        if let Ok(x) = v.parse() {
                            self.opts.az_to_hw = x;
                        }
                    }
                }
                "-M" | "--rec-rate" => {
                    if let Some(v) = val!() {
                        if let Ok(x) = v.parse() {
                            self.opts.rec_rate = x;
                        }
                    }
                }
                "-V" | "--viterbi-training" => {
                    let f = val!().and_then(|v| v.parse::<f64>().ok()).unwrap_or(0.0);
                    self.opts.viterbi_training = if f > 0.0 { 10 } else { 0 };
                }
                "-I" | "--skip-indels" => self.opts.skip_indels = true,
                "--ignore-homref" => self.opts.ignore_homref = true,
                "--include-noalt" => self.include_noalt = true,
                "-G" | "--GTs-only" => {
                    let _ = val!();
                }
                _ if t.starts_with("-G") => {}
                "-r" | "--regions" => {
                    if let Some(v) = val!() {
                        self.region = parse_region(&v);
                    }
                }
                "-s" | "--samples" => {
                    if let Some(v) = val!() {
                        self.samples =
                            Some(v.split(',').map(|x| x.trim().to_string()).collect());
                    }
                }
                "-m" | "--genetic-map" => {
                    if let Some(v) = val!() {
                        self.genetic_map = Some(v.into());
                    }
                }
                _ => {}
            }
            i += 1;
        }
    }

    /// Interpret a bcftools `-O` type string (`s`=ST, `r`=RG, `z`=bgzip).
    fn set_output(&mut self, s: &str) {
        self.output_st = s.contains('s');
        self.output_rg = s.contains('r');
    }
}

/// Run the HMM for one sample, processing each contiguous chromosome block
/// independently so transitions never cross a chromosome boundary.
fn run_sample(
    name: &str,
    sites: &[RohSite],
    run: &RohRun,
    out: &mut impl Write,
) -> Result<()> {
    let mut i = 0usize;
    while i < sites.len() {
        let mut j = i + 1;
        while j < sites.len() && sites[j].chrom == sites[i].chrom {
            j += 1;
        }
        let block = &sites[i..j];
        let mut opts = run.opts.clone();
        let iters = opts.viterbi_training;
        if iters > 0 {
            roh::baum_welch_train(block, &mut opts, iters);
        }
        let path = roh::viterbi(block, &opts);
        let post = roh::forward_backward(block, &opts);
        emit_block(name, block, &path, &post, run, out)?;
        i = j;
    }
    Ok(())
}

fn emit_block(
    name: &str,
    block: &[RohSite],
    path: &[State],
    post: &[[f64; 2]],
    run: &RohRun,
    out: &mut impl Write,
) -> Result<()> {
    let mut rg_beg: Option<u32> = None;
    let mut rg_end = 0u32;
    let mut rg_qsum = 0.0f64;
    let mut rg_n = 0u32;

    for k in 0..block.len() {
        let az = path[k] == State::AZ;
        // post[k] = [P(AZ), P(HW)]; quality is phred(1 - P(called state)).
        let p_called = if az { post[k][0] } else { post[k][1] };
        let qual = roh::phred_score(1.0 - p_called);

        if run.output_st {
            writeln!(
                out,
                "ST\t{}\t{}\t{}\t{}\t{:.1}",
                name,
                block[k].chrom,
                block[k].pos,
                if az { 1 } else { 0 },
                qual
            )?;
        }

        if az {
            if rg_beg.is_none() {
                rg_beg = Some(block[k].pos);
                rg_qsum = 0.0;
                rg_n = 0;
            }
            rg_end = block[k].pos;
            rg_qsum += qual;
            rg_n += 1;
        } else if let Some(beg) = rg_beg.take() {
            emit_rg(name, &block[k].chrom, beg, rg_end, rg_n, rg_qsum, run, out)?;
        }
    }
    if let Some(beg) = rg_beg.take() {
        let chrom = &block[block.len() - 1].chrom;
        emit_rg(name, chrom, beg, rg_end, rg_n, rg_qsum, run, out)?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn emit_rg(
    name: &str,
    chrom: &str,
    beg: u32,
    end: u32,
    n: u32,
    qsum: f64,
    run: &RohRun,
    out: &mut impl Write,
) -> Result<()> {
    if run.output_rg && n > 0 {
        writeln!(
            out,
            "RG\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}",
            name,
            chrom,
            beg,
            end,
            end - beg + 1,
            n,
            qsum / n as f64
        )?;
    }
    Ok(())
}

#[derive(Clone)]
struct Region {
    chrom: String,
    start: Option<u32>,
    end: Option<u32>,
}

impl Region {
    fn contains(&self, chrom: &str, pos: u32) -> bool {
        if self.chrom != chrom {
            return false;
        }
        if let Some(s) = self.start {
            if pos < s {
                return false;
            }
        }
        if let Some(e) = self.end {
            if pos > e {
                return false;
            }
        }
        true
    }
}

fn parse_region(s: &str) -> Option<Region> {
    if s.contains(':') {
        let (chrom, beg, end) = crate::regions::parse_region_spec(s).ok()?;
        return Some(Region { chrom, start: Some(beg), end: (end != u32::MAX).then_some(end) });
    }
    Some(Region {
        chrom: s.to_string(),
        start: None,
        end: None,
    })
}

fn sample_names(headers: &[String]) -> Vec<String> {
    headers
        .iter()
        .find(|h| h.starts_with("#CHROM\t"))
        .map(|h| h.split('\t').skip(9).map(|s| s.trim().to_string()).collect())
        .unwrap_or_default()
}

fn parse_all_gts(format: &Option<String>, samples: &[String]) -> Vec<GtClass> {
    let Some(fmt) = format else {
        return vec![GtClass::Missing; samples.len()];
    };
    let Some(gt_idx) = fmt.split(':').position(|k| k == "GT") else {
        return vec![GtClass::Missing; samples.len()];
    };
    samples
        .iter()
        .map(|s| classify_gt(s.split(':').nth(gt_idx)))
        .collect()
}

fn classify_gt(gt: Option<&str>) -> GtClass {
    let Some(gt) = gt else {
        return GtClass::Missing;
    };
    if gt.contains('.') {
        return GtClass::Missing;
    }
    let sep = if gt.contains('|') { '|' } else { '/' };
    let mut vals = Vec::new();
    for p in gt.split(sep) {
        match p.parse::<u32>() {
            Ok(v) => vals.push(v),
            Err(_) => return GtClass::Missing,
        }
    }
    if vals.is_empty() {
        GtClass::Missing
    } else if vals.iter().all(|v| *v == 0) {
        GtClass::HomRef
    } else if vals.iter().all(|v| *v == vals[0]) {
        GtClass::HomAlt
    } else {
        GtClass::Het
    }
}

fn is_indel(ref_allele: &str, alt: &str) -> bool {
    if ref_allele.len() != 1 {
        return true;
    }
    alt.split(',')
        .any(|a| a.len() != 1 || a.starts_with('<') || a == "*")
}

/// Alt-allele frequency for a site: `--AF-file` first, then the INFO `--AF-tag`,
/// then `--AF-dflt`.
fn resolve_af(
    chrom: &str,
    pos: u32,
    info: &str,
    run: &RohRun,
    af_file: Option<&HashMap<(String, u32), f64>>,
) -> f64 {
    if let Some(m) = af_file {
        if let Some(a) = m.get(&(chrom.to_string(), pos)) {
            return *a;
        }
    }
    if let Some(a) = info_af(info, &run.af_tag) {
        return a;
    }
    run.af_dflt
}

fn info_af(info: &str, tag: &str) -> Option<f64> {
    for field in info.split(';') {
        if let Some((k, v)) = field.split_once('=') {
            if k == tag {
                return v.split(',').next().unwrap_or(v).parse::<f64>().ok();
            }
        }
    }
    None
}

fn load_af_file(path: &PathBuf) -> Result<HashMap<(String, u32), f64>> {
    use std::io::{BufRead, BufReader, Read};
    let file = std::fs::File::open(path)?;
    let reader: Box<dyn Read> = if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        Box::new(flate2::read::MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let mut map = HashMap::new();
    for line in BufReader::new(reader).lines() {
        let l = line?;
        let t = l.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = t.split_whitespace().collect();
        if cols.len() < 3 {
            continue;
        }
        let chrom = cols[0].to_string();
        let Ok(pos) = cols[1].parse::<u32>() else {
            continue;
        };
        // AF is the last column (CHROM POS [REF ALT] AF).
        if let Ok(af) = cols[cols.len() - 1].parse::<f64>() {
            map.insert((chrom, pos), af);
        }
    }
    Ok(map)
}

/// Linear interpolation of cumulative cM at `pos` from a sorted `(pos, cM)` map.
fn interp_cm(map: &[(u32, f64)], pos: u32) -> Option<f64> {
    if map.is_empty() {
        return None;
    }
    match map.binary_search_by_key(&pos, |x| x.0) {
        Ok(i) => Some(map[i].1),
        Err(i) => {
            if i == 0 {
                Some(map[0].1)
            } else if i >= map.len() {
                Some(map[map.len() - 1].1)
            } else {
                let (p0, c0) = map[i - 1];
                let (p1, c1) = map[i];
                let frac = (pos - p0) as f64 / (p1 - p0).max(1) as f64;
                Some(c0 + frac * (c1 - c0))
            }
        }
    }
}
