//! `consensus`: bcftools-compatible driver over [`crate::consensus`].

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};

use crate::cli::args::ConsensusArgs;
use crate::consensus::{Consensus, Mark, Mask, MaskWith, Opts, PICK_ALT, PICK_IUPAC, PICK_LONG, PICK_REF, PICK_SHORT, Rec};
use crate::filter::FilterEngine;
use crate::regions::RegionSet;
use crate::vcf::{HeaderInfo, UnifiedVcfReader, parse_vcf_line};

pub fn cmd_consensus(args: ConsensusArgs) -> Result<()> {
    let Some(ref_path) = args.fasta_ref.as_ref() else { bail!("missing -f <ref.fa>") };
    let Some(input) = args.input.as_ref() else { bail!("missing input VCF/BCF for consensus") };

    let (haplotype, allele) = parse_haplotype(args.haplotype.as_deref())?;
    let mut masks = Vec::new();
    for (i, path) in args.mask.iter().enumerate() {
        let with = match args.mask_with.get(i) {
            Some(w) => MaskWith::parse(w)?,
            None => MaskWith::Char(b'N'),
        };
        masks.push(Mask::load(path, with)?);
    }
    let one_char = |v: Option<&String>, opt: &str| -> Result<Option<u8>> {
        match v {
            None => Ok(None),
            Some(s) if s.len() == 1 => Ok(Some(s.as_bytes()[0])),
            Some(s) => bail!("Expected single character with -{opt}, got \"{s}\""),
        }
    };

    // Records per chromosome, filtered like next_vcf_line().
    let mut reader = UnifiedVcfReader::open(input).with_context(|| format!("open {}", input.display()))?;
    let headers: Vec<String> = reader.header()?.to_vec();
    let hdr = HeaderInfo::parse(&headers);
    let (filter, exclude) = match (&args.include, &args.exclude) {
        (Some(_), Some(_)) => bail!("Error: only one -i or -e expression can be given, and they cannot be combined"),
        (Some(i), None) => (Some(i.as_str()), false),
        (None, Some(e)) => (Some(e.as_str()), true),
        (None, None) => (None, false),
    };
    let engine = filter.map(|f| FilterEngine::new(&headers, Some(f), exclude)).transpose()?;
    let regions = RegionSet::from_args(args.regions.as_deref(), args.regions_file.as_deref())?;

    let (isample, iupac_gts) = resolve_samples(&hdr.samples, args.samples.as_deref(), args.samples_file.as_deref(), haplotype != 0 || allele != 0)?;

    let mut vcf: HashMap<String, Vec<Rec>> = HashMap::new();
    while let Some(line) = reader.read_line()? {
        if line.is_empty() || line.as_bytes()[0] == b'#' {
            continue;
        }
        let Some(rec) = parse_vcf_line(&line) else { bail!("cannot parse VCF line: {line}") };
        if let Some(rf) = &regions {
            if !rf.record_passes(&rec.chrom, rec.pos, &rec.ref_allele, &rec.alt, args.regions_overlap) {
                continue;
            }
        }
        if let Some(e) = &engine {
            if !e.eval(&rec)?.pass_site {
                continue;
            }
        }
        vcf.entry(rec.chrom.clone()).or_default().push(Rec::from_record(&rec));
    }
    let contigs: HashMap<String, ()> = hdr.contigs.names().iter().map(|n| (n.clone(), ())).collect();

    let opts = Opts {
        output_iupac: args.iupac_codes,
        haplotype,
        allele,
        isample,
        iupac_gts,
        missing_allele: one_char(args.missing.as_ref(), "M")?,
        absent_allele: one_char(args.absent.as_ref(), "a")?,
        mark_del: args.mark_del.as_deref().and_then(|s| s.bytes().next()),
        mark_ins: args.mark_ins.as_deref().map(|s| Mark::parse(s, "mark-ins")).transpose()?.unwrap_or(Mark::None),
        mark_snv: args.mark_snv.as_deref().map(|s| Mark::parse(s, "mark-snv")).transpose()?.unwrap_or(Mark::None),
        chr_prefix: args.prefix.clone(),
        masks,
    };

    let out: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(BufWriter::with_capacity(1 << 20, File::create(p).with_context(|| format!("create {}", p.display()))?)),
        None => Box::new(BufWriter::with_capacity(1 << 20, std::io::stdout())),
    };
    let chain_out: Option<Box<dyn Write>> = match &args.chain {
        Some(p) => Some(Box::new(BufWriter::new(File::create(p).with_context(|| format!("create {}", p.display()))?))),
        None => None,
    };
    let fasta = crate::util::open_text_maybe_gz(ref_path)?;
    let mut c = Consensus::new(&opts, &vcf, &contigs, out, chain_out);
    c.run(fasta)?;
    if args.verbosity > 0 {
        eprintln!("Applied {} variants", c.applied());
    }
    Ok(())
}

/// `-H`: `R`, `A`, `L`, `S`, `LR`, `LA`, `SR`, `SA`, `I`, `N` or `NpIu`.
fn parse_haplotype(h: Option<&str>) -> Result<(i64, u32)> {
    let Some(h) = h else { return Ok((0, 0)) };
    let up = h.to_ascii_uppercase();
    let allele = match up.as_str() {
        "R" => PICK_REF,
        "A" => PICK_ALT,
        "L" | "LR" => PICK_LONG | PICK_REF,
        "S" | "SR" => PICK_SHORT | PICK_REF,
        "LA" => PICK_LONG | PICK_ALT,
        "SA" => PICK_SHORT | PICK_ALT,
        "I" => PICK_IUPAC,
        _ => {
            let digits = h.bytes().take_while(u8::is_ascii_digit).count();
            let (num, rest) = h.split_at(digits);
            if digits == 0 || !(rest.is_empty() || rest.eq_ignore_ascii_case("pIu")) {
                bail!("Error: Could not parse \"--haplotype {h}\", expected number of number followed with \"pIu\"");
            }
            let n: i64 = num.parse()?;
            if n <= 0 {
                bail!("Error: Expected positive integer with --haplotype");
            }
            return Ok((n, if rest.is_empty() { 0 } else { PICK_IUPAC }));
        }
    };
    Ok((0, allele))
}

/// `-s LIST` / `-S FILE` (`^` negates); without either every sample is used,
/// `-s -` uses REF/ALT without samples. With `-H` one sample is used,
/// otherwise the selected samples' genotypes become IUPAC codes (as
/// bcftools 1.19 does, its usage text notwithstanding).
fn resolve_samples(
    names: &[String],
    list: Option<&str>,
    file: Option<&Path>,
    with_haplotype: bool,
) -> Result<(Option<usize>, Option<Vec<usize>>)> {
    let mut wanted: Vec<String> = Vec::new();
    let mut negate = false;
    let mut idx: Vec<usize> = Vec::new();
    match (list, file) {
        (Some("-"), _) => return Ok((None, None)),
        (Some(l), _) => {
            let l = l.strip_prefix('^').inspect(|_| negate = true).unwrap_or(l);
            wanted.extend(l.split(',').filter(|s| !s.is_empty()).map(str::to_string));
        }
        (None, Some(f)) => {
            let text = std::fs::read_to_string(f).with_context(|| format!("read {}", f.display()))?;
            for line in text.lines() {
                let s = line.trim().split_whitespace().next().unwrap_or("");
                if !s.is_empty() {
                    wanted.push(s.strip_prefix('^').inspect(|_| negate = true).unwrap_or(s).to_string());
                }
            }
        }
        (None, None) => idx.extend(0..names.len()),
    }
    let listed = list.is_some() || file.is_some();
    if listed {
        if negate {
            idx.extend((0..names.len()).filter(|&i| !wanted.contains(&names[i])));
        } else {
            for w in &wanted {
                match names.iter().position(|n| n == w) {
                    Some(i) => idx.push(i),
                    None => bail!("No such sample: {w}"),
                }
            }
        }
        if idx.is_empty() {
            bail!("No matching sample found");
        }
    }
    if with_haplotype {
        if idx.len() > 1 {
            bail!("Too many samples, only one can be used with -H");
        }
        Ok((idx.first().copied(), None))
    } else {
        Ok((None, Some(idx)))
    }
}
