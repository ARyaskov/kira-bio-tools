use anyhow::{Context, Result, bail};
use std::fs::File;
use std::io::{BufWriter, Write};

use crate::VcfReader;
use crate::call::pedigree::{parse_ploidy_file, parse_sex_file};
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
