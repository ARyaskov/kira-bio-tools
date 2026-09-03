use anyhow::{Context, Result, bail};
use std::fs::File;
use std::io::{BufWriter, Write};

use crate::VcfReader;
use crate::call::mcall::ConstrainMode;
use crate::call::pedigree::{parse_ploidy_file, parse_prior_freqs, parse_sex_file};
use crate::call::stream::{CallConfig, CallMode, CallStream};
use crate::cli::args::CallArgs;

pub fn cmd_call(args: CallArgs) -> Result<()> {
    let cfg = config_from_args(&args)?;
    let out_path = args.output.clone().unwrap_or_else(|| {
        let mut p = args.input.clone();
        p.set_extension("call.vcf");
        p
    });

    let mut reader = VcfReader::open(&args.input)?;
    let headers = reader.header()?;

    let out = File::create(&out_path)?;
    let mut w = BufWriter::new(out);
    let mut stream = CallStream::new(&mut w, cfg)?;
    stream.write_header(&headers)?;

    while let Some(rec) = reader.next_record()? {
        stream.push(&rec)?;
    }
    stream.finish()?;

    w.flush()?;
    Ok(())
}

fn config_from_args(args: &CallArgs) -> Result<CallConfig> {
    let mut emit_gq = false;
    let mut emit_gp = false;
    if let Some(a) = &args.annotate {
        for t in a.split(',').map(str::trim) {
            if t.eq_ignore_ascii_case("GQ") {
                emit_gq = true;
            }
            if t.eq_ignore_ascii_case("GP") {
                emit_gp = true;
            }
        }
    }

    // `-S` is both the sample list and the sex table a ploidy file is keyed on.
    let sex_map = match &args.samples_file {
        Some(p) if p.as_os_str() != "-" && p.exists() => parse_sex_file(p)?,
        _ => Default::default(),
    };
    let ploidy_regions = match &args.ploidy_file {
        Some(p) => parse_ploidy_file(p).with_context(|| format!("--ploidy-file {}", p.display()))?,
        None => Vec::new(),
    };

    let prior_freqs = match (&args.prior_freqs, &args.prior_af) {
        (Some(spec), _) => Some(parse_prior_freqs(spec)?),
        (None, Some(tag)) => Some(parse_prior_freqs(tag.trim_start_matches("INFO/"))?),
        (None, None) => None,
    };
    let constrain = match args.constrain.as_deref().map(str::trim) {
        None => ConstrainMode::None,
        Some("trio") => ConstrainMode::Trio,
        Some("alleles") => ConstrainMode::Alleles,
        Some(other) => bail!("call: -C expects 'trio' or 'alleles' (got {other:?})"),
    };
    let mut novel_rate = [1e-8, 1e-9, 1e-9];
    if let Some(spec) = &args.novel_rate {
        let vals: Vec<f64> = spec
            .split(',')
            .map(|s| s.trim().parse::<f64>().with_context(|| format!("-X/--novel-rate: bad value {s:?}")))
            .collect::<Result<_>>()?;
        match vals.len() {
            1 => novel_rate = [vals[0]; 3],
            2 => novel_rate = [vals[0], vals[1], vals[1]],
            3 => novel_rate = [vals[0], vals[1], vals[2]],
            n => bail!("-X/--novel-rate: expected 1 to 3 rates, got {n}"),
        }
    }
    let ped_path = if constrain == ConstrainMode::Trio { args.samples_file.clone() } else { None };

    let cfg = CallConfig {
        mode: if args.consensus_caller {
            CallMode::Consensus
        } else {
            CallMode::Multiallelic
        },
        variants_only: args.variants_only,
        emit_gq,
        emit_gp,
        theta: args.prior,
        ploidy: parse_ploidy(args.ploidy.as_deref())?,
        ploidy_regions,
        sex_map,
        samples: args.samples.clone(),
        samples_file: args.samples_file.clone(),
        prior_freqs,
        constrain,
        ped_path,
        novel_rate,
        groups_path: args.group_samples.clone(),
        group_tag: args.group_samples_tag.clone().unwrap_or_else(|| "AD".to_string()),
        gvcf: match &args.gvcf {
            Some(spec) => Some(
                spec.split(',')
                    .map(|s| s.trim().parse::<u32>().with_context(|| format!("-g/--gvcf: bad depth {s:?}")))
                    .collect::<Result<Vec<u32>>>()?,
            ),
            None => None,
        },
    };
    cfg.validate()?;
    Ok(cfg)
}

/// `--ploidy` takes a plain number; named assemblies are region-dependent and
/// belong in `--ploidy-file`, so they are rejected instead of defaulting to 2.
fn parse_ploidy(spec: Option<&str>) -> Result<u8> {
    let Some(s) = spec else { return Ok(2) };
    match s.trim() {
        "0" => Ok(0),
        "1" => Ok(1),
        "2" | "-" => Ok(2),
        other => bail!(
            "--ploidy: expected 0, 1, 2 or '-' (got {other:?}); \
             use --ploidy-file for region-dependent ploidy"
        ),
    }
}
