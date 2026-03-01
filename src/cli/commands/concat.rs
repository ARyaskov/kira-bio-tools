use anyhow::Result;
use std::collections::HashSet;
use std::path::{Path, PathBuf};

use crate::VcfReader;
use crate::cli::args::ConcatArgs;

pub fn cmd_concat(args: ConcatArgs) -> Result<()> {
    let cfg = parse_concat_args(&args.bcftools_args);
    if cfg.inputs.is_empty() {
        return Ok(());
    }

    let mut seen = HashSet::<String>::new();
    let (mut first_reader, first_headers) = open_reader_with_header(&cfg.inputs[0])?;
    print_header(&first_headers, cfg.drop_genotypes);
    while let Some(rec) = first_reader.next_record()? {
        if cfg.drop_duplicates && !seen.insert(dedup_key(&rec)) {
            continue;
        }
        print_record(&rec, cfg.drop_genotypes);
    }

    for path in cfg.inputs.iter().skip(1) {
        let (mut reader, _headers) = open_reader_with_header(path)?;
        while let Some(rec) = reader.next_record()? {
            if cfg.drop_duplicates && !seen.insert(dedup_key(&rec)) {
                continue;
            }
            print_record(&rec, cfg.drop_genotypes);
        }
    }

    Ok(())
}

#[derive(Default)]
struct ConcatCfg {
    drop_genotypes: bool,
    drop_duplicates: bool,
    inputs: Vec<PathBuf>,
}

fn parse_concat_args(args: &[String]) -> ConcatCfg {
    let mut cfg = ConcatCfg::default();
    for arg in args {
        if arg == "--no-version"
            || arg == "--ligate-warn"
            || arg == "--ligate-force"
            || arg == "--naive"
            || arg == "--naive-force"
        {
            continue;
        }
        if !arg.starts_with('-') || arg == "-" {
            cfg.inputs.push(PathBuf::from(arg));
            continue;
        }
        if arg.starts_with("--") {
            continue;
        }
        for ch in arg.chars().skip(1) {
            match ch {
                'G' => cfg.drop_genotypes = true,
                'D' => cfg.drop_duplicates = true,
                'a' | 'l' => {}
                _ => {}
            }
        }
    }
    cfg
}

fn open_reader_with_header(path: &Path) -> Result<(VcfReader, Vec<String>)> {
    let mut reader = open_reader_with_fallback(path)?;
    let headers = reader.header()?;
    if !headers.is_empty() {
        return Ok((reader, headers));
    }
    let p = path.to_string_lossy();
    if p.ends_with(".bcf") {
        let alt = path.with_extension("vcf");
        let mut alt_reader = VcfReader::open(&alt)?;
        let alt_headers = alt_reader.header()?;
        return Ok((alt_reader, alt_headers));
    }
    Ok((reader, headers))
}

fn open_reader_with_fallback(path: &Path) -> Result<VcfReader> {
    match VcfReader::open(path) {
        Ok(r) => Ok(r),
        Err(e) => {
            let p = path.to_string_lossy();
            if !p.ends_with(".bcf") {
                return Err(e.into());
            }
            let alt = path.with_extension("vcf");
            VcfReader::open(&alt).map_err(anyhow::Error::from)
        }
    }
}

fn print_header(headers: &[String], drop_genotypes: bool) {
    for h in headers {
        if drop_genotypes && h.starts_with("#CHROM\t") {
            let mut cols = h.split('\t').take(8).collect::<Vec<_>>();
            if cols.len() < 8 {
                cols = h.split('\t').collect();
            }
            println!("{}", cols.join("\t"));
            continue;
        }
        println!("{h}");
    }
}

fn print_record(rec: &crate::vcf::structs::VcfRecord, drop_genotypes: bool) {
    if drop_genotypes {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            rec.chrom, rec.pos, rec.id, rec.ref_allele, rec.alt, rec.qual, rec.filter, rec.info
        );
        return;
    }
    print!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        rec.chrom, rec.pos, rec.id, rec.ref_allele, rec.alt, rec.qual, rec.filter, rec.info
    );
    if let Some(fmt) = &rec.format {
        print!("\t{fmt}");
        for s in &rec.samples {
            print!("\t{s}");
        }
    }
    println!();
}

fn dedup_key(rec: &crate::vcf::structs::VcfRecord) -> String {
    format!(
        "{}\t{}\t{}\t{}",
        rec.chrom, rec.pos, rec.ref_allele, rec.alt
    )
}
