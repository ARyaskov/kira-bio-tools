use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result, bail};

use crate::annotate::postproc::{
    HeaderOptions, OutputKind, apply_to_header, parse_output_type, parse_samples_cli,
    read_samples_file, resolve_samples_keep,
};
use crate::cli::args::CsqArgs;
use crate::csi::{IndexKind, build_index};
use crate::csq::gff::GffOpts;
use crate::csq::{Caller, CsqOpts, CsqOutput, Phase, ncsq2_to_nfmt};
use crate::fasta::IndexedFasta;
use crate::filter::FilterEngine;
use crate::regions::RegionSet;
use crate::vcf::header::extract_samples;
use crate::vcf::{ContigDict, UnifiedVcfReader, VcfSink, parse_vcf_record};

pub fn cmd_csq(args: CsqArgs) -> Result<()> {
    let Some(fasta_path) = &args.fasta_ref else {
        bail!("Error: missing the -f, --fasta-ref option")
    };
    let Some(gff_path) = &args.gff else {
        bail!("Error: missing the -g, --gff-annot option")
    };
    let input = args
        .input
        .clone()
        .unwrap_or_else(|| Path::new("-").to_path_buf());
    let mut phase = match args.phase.as_str() {
        "a" => Phase::AsIs,
        "m" => Phase::Merge,
        "r" => Phase::Require,
        "R" => Phase::NonRef,
        "s" => Phase::Skip,
        s => bail!("The -p argument not recognised, expected one of a,m,r,R,s: {s}"),
    };
    if args.ncsq == 0 {
        bail!("Expected positive integer with -n, got {}", args.ncsq);
    }
    let ncsq2_max = 2 * args.ncsq;
    let brief = match args.trim_protein_seq {
        Some(0) => bail!("Could not parse argument: --trim-protein-seq 0"),
        Some(n) => n,
        None if args.brief_predictions => 1,
        None => 0,
    };
    let unify = match args.unify_chr_names {
        0 => false,
        1 => true,
        n => bail!("Could not parse: --unify-chr-names {n}"),
    };
    if args.include.is_some() && args.exclude.is_some() {
        bail!("Only one of -i or -e can be given.");
    }
    let (text_out, kind) = match args.output_type.as_deref() {
        Some(s) if s.starts_with('t') => (true, OutputKind::Vcf),
        Some(s) => (false, parse_output_type(s)?),
        None => (false, OutputKind::Vcf),
    };

    let gff = crate::csq::gff::parse(
        gff_path,
        &GffOpts {
            strip_chr_names: unify,
            force: args.force,
            verbosity: args.verbosity,
        },
    )?;
    let fasta = IndexedFasta::open(fasta_path)
        .with_context(|| format!("Failed to load the fai index: {}", fasta_path.display()))?;

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

    let mut reader = UnifiedVcfReader::open(&input).context("open input")?;
    let headers = reader.header()?.to_vec();
    let samples = extract_samples(&headers);
    let smpl: Vec<usize> = if samples.is_empty() || args.samples.as_deref() == Some("-") {
        phase = Phase::DropGt;
        Vec::new()
    } else if let Some(s) = &args.samples {
        let (names, inv) = parse_samples_cli(s);
        for n in &names {
            if !samples.iter().any(|x| x == n) {
                bail!("No such sample: \"{n}\"");
            }
        }
        resolve_samples_keep(&samples, &names, inv)
    } else if let Some(p) = &args.samples_file {
        let (names, inv) = read_samples_file(p)?;
        for n in &names {
            if !samples.iter().any(|x| x == n) {
                bail!("No such sample: \"{n}\"");
            }
        }
        resolve_samples_keep(&samples, &names, inv)
    } else {
        (0..samples.len()).collect()
    };
    let hdr_nsmpl = if phase == Phase::DropGt {
        0
    } else {
        samples.len()
    };

    let filter = match (&args.include, &args.exclude) {
        (Some(e), _) => Some((
            FilterEngine::new(&headers, Some(e.as_str()), false)
                .context("-i/--include expression")?,
            false,
        )),
        (_, Some(e)) => Some((
            FilterEngine::new(&headers, Some(e.as_str()), false)
                .context("-e/--exclude expression")?,
            true,
        )),
        _ => None,
    };

    let tag = args.custom_tag.clone();
    let cmdline = std::env::args().skip(1).collect::<Vec<_>>().join(" ");
    let out: CsqOutput = if text_out {
        let mut w: Box<dyn Write> = match args.output.as_deref() {
            Some(p) if p != Path::new("-") => Box::new(std::io::BufWriter::with_capacity(
                1 << 20,
                std::fs::File::create(p)
                    .with_context(|| format!("Failed to write to {}", p.display()))?,
            )),
            _ => Box::new(std::io::BufWriter::with_capacity(
                1 << 20,
                std::io::stdout(),
            )),
        };
        writeln!(
            w,
            "# This file was produced by: kira-bt csq({})",
            env!("CARGO_PKG_VERSION")
        )?;
        writeln!(w, "# The command line was:\tkira-bt {cmdline}")?;
        writeln!(w, "# LOG\t[2]Message")?;
        writeln!(
            w,
            "# CSQ\t[2]Sample\t[3]Haplotype\t[4]Chromosome\t[5]Position\t[6]Consequence"
        )?;
        CsqOutput::Text(w)
    } else {
        let mut extra = Vec::new();
        if !headers
            .iter()
            .any(|h| h.starts_with(&format!("##INFO=<ID={tag},")))
        {
            extra.push(format!(
                "##INFO=<ID={tag},Number=.,Type=String,Description=\"{} consequence annotation from BCFtools/csq, see http://samtools.github.io/bcftools/howtos/csq-calling.html for details. Format: Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change\">",
                if args.local_csq { "Local" } else { "Haplotype-aware" }
            ));
        }
        if hdr_nsmpl > 0
            && !headers
                .iter()
                .any(|h| h.starts_with(&format!("##FORMAT=<ID={tag},")))
        {
            extra.push(format!(
                "##FORMAT=<ID={tag},Number=.,Type=Integer,Description=\"Bitmask of indexes to INFO/BCSQ, with interleaved first/second haplotype. Use \\\"bcftools query -f'[%CHROM\\t%POS\\t%SAMPLE\\t%TBCSQ\\n]'\\\" to translate.\">"
            ));
        }
        if !args.no_version {
            extra.push(format!(
                "##kira_bt_csqVersion={}+htslib-compat",
                env!("CARGO_PKG_VERSION")
            ));
            extra.push(format!(
                "##kira_bt_csqCommand={cmdline}; Date={}",
                crate::util::now_string()
            ));
        }
        let opts = HeaderOptions {
            no_version: true,
            extra_header_lines: &extra,
            remove: None,
            rename_chrs: None,
            rename_annots: None,
            mark_sites: None,
            set_id: false,
            samples_keep: None,
            version_line: None,
        };
        let out_headers = apply_to_header(headers.clone(), &opts);
        let mut sink = VcfSink::open(args.output.as_deref(), kind, &out_headers)?;
        sink.write_header(&out_headers)?;
        CsqOutput::Vcf(sink)
    };

    let opts = CsqOpts {
        phase,
        ncsq2_max,
        nfmt_bcsq: ncsq2_to_nfmt(ncsq2_max),
        local_csq: args.local_csq,
        brief_predictions: brief,
        unify_chr_names: unify,
        verbosity: args.verbosity,
        bcsq_tag: tag,
        hdr_nsmpl,
        fa_fname: fasta_path.display().to_string(),
    };
    if args.verbosity > 0 {
        eprintln!("Calling...");
    }
    let mut caller = Caller::new(gff, fasta, opts, smpl, samples, out);
    let mut contigs = ContigDict::new();
    let handle = |line: &str, caller: &mut Caller, contigs: &mut ContigDict| -> Result<()> {
        if line.is_empty() || line.as_bytes()[0] == b'#' {
            return Ok(());
        }
        if let Some(tf) = &target_filter
            && tf.line_passes_mode(line, args.targets_overlap) == targets_inverse
        {
            return Ok(());
        }
        let Some(rec) = parse_vcf_record(line, 0, contigs).map_err(|e| anyhow::anyhow!("{e}"))?
        else {
            return Ok(());
        };
        let pass = match &filter {
            Some((engine, exclude)) => {
                let p = engine.eval(&rec)?.pass_site;
                Some(if *exclude { !p } else { p })
            }
            None => None,
        };
        caller.process(rec, pass)
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
    caller.finish()?;
    match caller.into_output() {
        CsqOutput::Vcf(sink) => sink.finish()?,
        CsqOutput::Text(mut w) => w.flush()?,
    }
    if args.write_index
        && let Some(out) = args.output.as_deref()
    {
        if matches!(kind, OutputKind::VcfGz(_) | OutputKind::Bcf(_)) && out != Path::new("-") {
            let idx = std::path::PathBuf::from(format!("{}.csi", out.display()));
            build_index(out, &idx, IndexKind::Csi, None)
                .with_context(|| format!("--write-index: write index {}", idx.display()))?;
        } else {
            eprintln!(
                "Warning: --write-index requires a compressed VCF or BCF written to a file; skipped"
            );
        }
    }
    Ok(())
}
