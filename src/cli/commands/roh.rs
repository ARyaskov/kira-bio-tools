use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use flate2::Compression;

use crate::bgzf::{BgzfWriter, FILE_BUFFER_SIZE, STREAM_BUFFER_SIZE};
use crate::cli::args::RohArgs;
use crate::cli::getopt::{OptSpec, opt, parse_f64};
use crate::filter::FilterEngine;
use crate::regions::RegionSet;
use crate::roh::{RohCaller, RohOpts, sample_list};
use crate::vcf::header::extract_samples;
use crate::vcf::{ContigDict, UnifiedVcfReader, parse_vcf_record};

const SPECS: &[OptSpec] = &[
    opt(None, "AF-tag", true),
    opt(None, "AF-file", true),
    opt(None, "AF-dflt", true),
    opt(None, "include", true),
    opt(None, "exclude", true),
    opt(None, "include-noalt", false),
    opt(Some('b'), "buffer-size", true),
    opt(Some('i'), "ignore-homref", false),
    opt(Some('e'), "estimate-AF", true),
    opt(Some('o'), "output", true),
    opt(Some('O'), "output-type", true),
    opt(Some('G'), "GTs-only", true),
    opt(Some('s'), "samples", true),
    opt(Some('S'), "samples-file", true),
    opt(Some('a'), "hw-to-az", true),
    opt(Some('H'), "az-to-hw", true),
    opt(Some('V'), "viterbi-training", true),
    opt(Some('t'), "targets", true),
    opt(Some('T'), "targets-file", true),
    opt(None, "targets-overlap", true),
    opt(Some('r'), "regions", true),
    opt(Some('R'), "regions-file", true),
    opt(None, "regions-overlap", true),
    opt(Some('m'), "genetic-map", true),
    opt(Some('M'), "rec-rate", true),
    opt(Some('I'), "skip-indels", false),
    opt(None, "threads", true),
];

fn overlap_opt(v: &str, what: &str) -> Result<u8> {
    match v {
        "0" | "pos" => Ok(0),
        "1" | "record" => Ok(1),
        "2" | "variant" => Ok(2),
        _ => bail!("Could not parse: {what} {v}"),
    }
}

/// Fold the tokens after `--` into the clap-parsed arguments.
fn apply_passthrough(args: &mut RohArgs) -> Result<()> {
    let toks = std::mem::take(&mut args.passthrough);
    let p = crate::cli::getopt::parse(&toks, SPECS)?;
    for (name, val) in p.opts {
        let v = val.unwrap_or_default();
        match name {
            "AF-tag" => args.af_tag = Some(v),
            "AF-file" => args.af_file = Some(PathBuf::from(v)),
            "AF-dflt" => args.af_dflt = Some(parse_f64(&v, "--AF-dflt")?),
            "include" => args.include = Some(v),
            "exclude" => args.exclude = Some(v),
            "include-noalt" => args.include_noalt = true,
            "buffer-size" => args.buffer_size = Some(v),
            "ignore-homref" => args.ignore_homref = true,
            "estimate-AF" => args.estimate_af = Some(v),
            "output" => args.output = Some(PathBuf::from(v)),
            "output-type" => args.output_type = Some(v),
            "GTs-only" => args.gts_only = Some(parse_f64(&v, "-G")?),
            "samples" => args.samples = Some(v),
            "samples-file" => args.samples_file = Some(PathBuf::from(v)),
            "hw-to-az" => args.hw_to_az = Some(parse_f64(&v, "-a")?),
            "az-to-hw" => args.az_to_hw = Some(parse_f64(&v, "-H")?),
            "viterbi-training" => {
                args.viterbi_training = Some(parse_f64(&v, "--viterbi-training")?)
            }
            "targets" => args.targets = Some(v),
            "targets-file" => args.targets_file = Some(PathBuf::from(v)),
            "targets-overlap" => args.targets_overlap = overlap_opt(&v, "--targets-overlap")?,
            "regions" => args.regions = Some(v),
            "regions-file" => args.regions_file = Some(PathBuf::from(v)),
            "regions-overlap" => args.regions_overlap = overlap_opt(&v, "--regions-overlap")?,
            "genetic-map" => args.genetic_map = Some(v),
            "rec-rate" => args.rec_rate = Some(parse_f64(&v, "-M")?),
            "skip-indels" => args.skip_indels = true,
            "threads" => args.threads = v.parse().unwrap_or(0),
            _ => {}
        }
    }
    if args.input.is_none()
        && let Some(f) = p.positional.first() {
            args.input = Some(PathBuf::from(f));
        }
    Ok(())
}

pub fn cmd_roh(mut args: RohArgs) -> Result<()> {
    apply_passthrough(&mut args)?;
    let cmdline = std::env::args().skip(1).collect::<Vec<_>>().join(" ");
    let Some(input) = args.input.clone() else {
        bail!("Usage: kira-bt roh [options] <in.vcf.gz>")
    };
    let mut naf_opts = 0;
    if args.af_tag.is_some() {
        naf_opts += 1;
    }
    if args.af_file.is_some() {
        naf_opts += 1;
    }
    if args.estimate_af.is_some() {
        naf_opts += 1;
    }
    if naf_opts > 1 {
        bail!("Error: The options --AF-tag, --AF-file and -e are mutually exclusive");
    }
    if args.af_file.is_some() && (args.targets.is_some() || args.targets_file.is_some()) {
        bail!("Error: The options --AF-file and -t are mutually exclusive");
    }
    if args.include.is_some() && args.exclude.is_some() {
        bail!(
            "Error: only one --include or --exclude expression can be given, and they cannot be combined"
        );
    }
    let t2az = args.hw_to_az.unwrap_or(6.7e-8);
    let t2hw = args.az_to_hw.unwrap_or(5e-9);
    if !(0.0..=1.0).contains(&t2az) {
        bail!(
            "Error: The parameter --hw-to-az is not in [0,1] .. {}",
            crate::util::fmt_e(t2az, 6)
        );
    }
    if !(0.0..=1.0).contains(&t2hw) {
        bail!(
            "Error: The parameter --az-to-hw is not in [0,1] .. {}",
            crate::util::fmt_e(t2hw, 6)
        );
    }
    if args.viterbi_training.is_some() && args.buffer_size.is_some() {
        bail!("Error: cannot use -b with -V");
    }
    let ot = args.output_type.clone().unwrap_or_default();
    let mut output_st = ot.contains(['s', 'S']);
    let mut output_rg = ot.contains(['r', 'R']);
    let output_gz = ot.contains('z');
    if !output_st && !output_rg {
        output_st = true;
        output_rg = true;
    }

    let region_filter = if let Some(s) = &args.regions {
        Some(RegionSet::from_cli(s)?)
    } else if let Some(p) = &args.regions_file {
        Some(RegionSet::from_file(p)?)
    } else {
        None
    };
    let target_filter = if let Some(s) = &args.targets {
        Some(RegionSet::from_cli(s.trim_start_matches('^'))?)
    } else if let Some(p) = &args.targets_file {
        Some(RegionSet::from_file(p)?)
    } else {
        None
    };
    let targets_inverse = args.targets.as_deref().is_some_and(|s| s.starts_with('^'));

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
    let samples = extract_samples(&headers);
    let roh_smpl: Vec<usize> = if let Some(s) = &args.samples {
        sample_list(&samples, s, false)?
    } else if let Some(p) = &args.samples_file {
        sample_list(&samples, &p.display().to_string(), true)?
    } else {
        (0..samples.len()).collect()
    };
    let filter = match (&args.include, &args.exclude) {
        (Some(e), _) => Some((
            FilterEngine::new(&headers, Some(e.as_str()), false).context("--include expression")?,
            false,
        )),
        (_, Some(e)) => Some((
            FilterEngine::new(&headers, Some(e.as_str()), false).context("--exclude expression")?,
            true,
        )),
        _ => None,
    };

    let out_path = args
        .output
        .clone()
        .filter(|p| p != Path::new("-") && p != Path::new("stdout"));
    let raw: Box<dyn Write + Send> = match &out_path {
        Some(p) => Box::new(
            std::fs::File::create(p).with_context(|| format!("Failed to open {}", p.display()))?,
        ),
        None => Box::new(std::io::stdout()),
    };
    let out: Box<dyn Write> = if output_gz {
        let buf = if out_path.is_none() {
            STREAM_BUFFER_SIZE
        } else {
            FILE_BUFFER_SIZE
        };
        Box::new(BgzfWriter::from_writer_buffered(
            raw,
            Compression::new(6),
            buf,
        )?)
    } else {
        Box::new(std::io::BufWriter::with_capacity(1 << 20, raw))
    };

    let opts = RohOpts {
        t2az,
        t2hw,
        rec_rate: args.rec_rate.unwrap_or(0.0),
        genmap_fname: args.genetic_map.clone(),
        af_tag: args.af_tag.clone(),
        af_file: args.af_file.as_ref().map(|p| p.display().to_string()),
        dflt_af: args.af_dflt.unwrap_or(0.0),
        estimate_af: args.estimate_af.clone(),
        gts_only: args.gts_only,
        snps_only: args.skip_indels,
        skip_homref: args.ignore_homref,
        include_noalt: args.include_noalt,
        vi_training: args.viterbi_training,
        buffer_size: args.buffer_size.clone(),
        output_st,
        output_rg,
        cmdline,
    };
    let mut caller = RohCaller::new(opts, &headers, samples, roh_smpl, out)?;
    let mut contigs = ContigDict::new();
    let handle = |line: &str, caller: &mut RohCaller, contigs: &mut ContigDict| -> Result<()> {
        if line.is_empty() || line.as_bytes()[0] == b'#' {
            return Ok(());
        }
        if let Some(tf) = &target_filter
            && tf.line_passes_mode(line, args.targets_overlap) == targets_inverse {
                return Ok(());
            }
        let Some(rec) = parse_vcf_record(line, 0, contigs).map_err(|e| anyhow::anyhow!("{e}"))?
        else {
            return Ok(());
        };
        if !caller.af_file_has(&rec.chrom, rec.pos) {
            return Ok(());
        }
        let pass = match &filter {
            Some((engine, exclude)) => {
                let p = engine.eval(&rec)?.pass_site;
                if *exclude { !p } else { p }
            }
            None => true,
        };
        caller.process(&rec, pass)
    };
    let indexed = region_filter
        .as_ref()
        .filter(|_| input != Path::new("-") && crate::csi::find_index_for(&input).is_some());
    if let Some(rf) = indexed {
        rf.stream_with_index(&input, args.regions_overlap, |line| {
            handle(line, &mut caller, &mut contigs)
        })?;
    } else {
        while let Some(line) = reader.read_line()? {
            if let Some(rf) = &region_filter
                && (line.is_empty()
                    || line.as_bytes()[0] == b'#'
                    || !rf.line_passes_mode(&line, args.regions_overlap))
                {
                    continue;
                }
            handle(&line, &mut caller, &mut contigs)?;
        }
    }
    caller.finish()
}
