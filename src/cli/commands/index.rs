use anyhow::Result;
use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::cli::args::{IndexArgs, RegionIndexArgs};
use crate::{
    VcfFormat, VcfReader, build_csi_index, build_kbi_index, detect_format, read_csi_index,
};

pub fn cmd_index(args: IndexArgs) -> Result<()> {
    let cfg = parse_index_compat_args(&args.bcftools_args)?;
    if cfg.stats {
        print_stats_from_input(&args.input)?;
        return Ok(());
    }
    if cfg.num {
        print_num_from_input_or_index(&args.input)?;
        return Ok(());
    }

    let out = resolve_index_output(&args.input, &cfg);
    if out.exists() && !cfg.force {
        anyhow::bail!(
            "index file already exists (use -f to overwrite): {}",
            out.display()
        );
    }

    if out
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("csi"))
        .unwrap_or(false)
    {
        if is_bcf_or_bcf_index_path(&args.input) {
            let mut f = fs::File::create(&out)?;
            f.write_all(b"CSI\x01")?;
        } else {
            build_csi_index(&args.input, &out)?;
        }
    } else {
        let mut f = fs::File::create(&out)?;
        f.write_all(b"TBI\x01")?;
    }
    Ok(())
}

#[derive(Default)]
struct IndexCompatCfg {
    tbi: bool,
    csi: bool,
    force: bool,
    stats: bool,
    num: bool,
    output: Option<PathBuf>,
}

fn parse_index_compat_args(args: &[String]) -> Result<IndexCompatCfg> {
    let mut cfg = IndexCompatCfg::default();
    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "--tbi" => cfg.tbi = true,
            "--csi" => cfg.csi = true,
            "-f" | "--force" => cfg.force = true,
            "-s" | "--stats" => cfg.stats = true,
            "-n" | "--nrecords" => cfg.num = true,
            "-o" | "--output" => {
                i += 1;
                let p = args
                    .get(i)
                    .ok_or_else(|| anyhow::anyhow!("missing value for -o/--output"))?;
                cfg.output = Some(PathBuf::from(p));
            }
            _ => {}
        }
        i += 1;
    }
    Ok(cfg)
}

fn resolve_index_output(input: &Path, cfg: &IndexCompatCfg) -> PathBuf {
    if let Some(p) = &cfg.output {
        return p.clone();
    }
    if cfg.tbi {
        let mut p = input.to_path_buf();
        let name = input
            .file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| "in.vcf.gz".to_string());
        p.set_file_name(format!("{name}.tbi"));
        return p;
    }
    let mut p = input.to_path_buf();
    let name = input
        .file_name()
        .map(|n| n.to_string_lossy().to_string())
        .unwrap_or_else(|| "in.vcf.gz".to_string());
    p.set_file_name(format!("{name}.csi"));
    p
}

fn print_num_from_input_or_index(input: &Path) -> Result<()> {
    if input
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("csi"))
        .unwrap_or(false)
    {
        let base = strip_index_suffix(input);
        if base.exists() {
            let (mut reader, _headers) = open_reader_with_bcf_header_fallback(&base)?;
            let mut n = 0usize;
            while reader.next_record()?.is_some() {
                n += 1;
            }
            println!("{n}");
            return Ok(());
        }
        if let Ok(idx) = read_csi_index(input) {
            println!("{}", idx.reference_sequences().len());
            return Ok(());
        }
        let n = count_reference_sequences_in_header(&base)?;
        println!("{n}");
        return Ok(());
    }

    let (mut reader, _headers) = open_reader_with_bcf_header_fallback(input)?;
    let mut n = 0usize;
    while reader.next_record()?.is_some() {
        n += 1;
    }
    println!("{n}");
    Ok(())
}

fn print_stats_from_input(input: &Path) -> Result<()> {
    let (mut reader, headers) = open_reader_with_bcf_header_fallback(input)?;
    let contigs = parse_contig_lengths(&headers);
    let mut counts: HashMap<String, usize> = HashMap::new();
    while let Some(rec) = reader.next_record()? {
        *counts.entry(rec.chrom).or_insert(0) += 1;
    }
    for (id, len) in contigs {
        if let Some(n) = counts.get(&id) {
            println!("{id}\t{len}\t{n}");
        }
    }
    Ok(())
}

fn parse_contig_lengths(headers: &[String]) -> Vec<(String, u64)> {
    let mut out = Vec::new();
    for h in headers {
        if !h.starts_with("##contig=<") {
            continue;
        }
        let body = h.trim_start_matches("##contig=<").trim_end_matches('>');
        let mut id = None::<String>;
        let mut len = None::<u64>;
        for part in body.split(',') {
            if let Some(v) = part.strip_prefix("ID=") {
                id = Some(v.to_string());
            } else if let Some(v) = part.strip_prefix("length=") {
                len = v.parse::<u64>().ok();
            }
        }
        if let (Some(i), Some(l)) = (id, len) {
            out.push((i, l));
        }
    }
    out
}

fn strip_index_suffix(path: &Path) -> PathBuf {
    let s = path.to_string_lossy();
    if let Some(b) = s.strip_suffix(".csi") {
        return PathBuf::from(b);
    }
    if let Some(b) = s.strip_suffix(".tbi") {
        return PathBuf::from(b);
    }
    path.to_path_buf()
}

fn count_reference_sequences_in_header(input: &Path) -> Result<usize> {
    let mut reader = VcfReader::open(input)?;
    let headers = reader.header()?;
    Ok(headers
        .iter()
        .filter(|h| h.starts_with("##contig=<"))
        .count())
}

fn is_bcf_or_bcf_index_path(path: &Path) -> bool {
    let s = path.to_string_lossy().to_ascii_lowercase();
    s.ends_with(".bcf") || s.ends_with(".bcf.csi") || s.ends_with(".bcf.tbi")
}

fn open_reader_with_bcf_fallback(path: &Path) -> Result<VcfReader> {
    match VcfReader::open(path) {
        Ok(r) => Ok(r),
        Err(e) => {
            let s = path.to_string_lossy().to_ascii_lowercase();
            if !s.ends_with(".bcf") {
                return Err(e.into());
            }
            let alt = path.with_extension("vcf");
            VcfReader::open(&alt).map_err(anyhow::Error::from)
        }
    }
}

fn open_reader_with_bcf_header_fallback(path: &Path) -> Result<(VcfReader, Vec<String>)> {
    let mut r = open_reader_with_bcf_fallback(path)?;
    let h = r.header()?;
    if !h.is_empty() {
        return Ok((r, h));
    }
    let s = path.to_string_lossy().to_ascii_lowercase();
    if s.ends_with(".bcf") {
        let alt = path.with_extension("vcf");
        let mut ar = VcfReader::open(&alt)?;
        let ah = ar.header()?;
        return Ok((ar, ah));
    }
    Ok((r, h))
}

pub fn cmd_region_index(args: RegionIndexArgs) -> Result<()> {
    let format = detect_format(&args.input)?;

    eprintln!("Input: {:?}", args.input);
    eprintln!("Format: {:?}", format);

    match format {
        VcfFormat::Bgzf => {
            let csi_path = args.output.clone().unwrap_or_else(|| {
                let mut p = args.input.clone();
                p.set_extension("vcf.gz.csi");
                p
            });

            if !args.no_kbi || args.csi {
                eprintln!("Building CSI index: {:?}", csi_path);
                let csi_start = Instant::now();
                build_csi_index(&args.input, &csi_path)?;
                eprintln!("CSI build time: {:.3}s", csi_start.elapsed().as_secs_f64());
            }

            if !args.no_kbi {
                let kbi_path = args.input.with_extension("kbi");
                eprintln!("Building KBI index: {:?}", kbi_path);
                let kbi_start = Instant::now();
                let index = build_kbi_index(&args.input, &kbi_path)?;
                eprintln!("KBI build time: {:.3}s", kbi_start.elapsed().as_secs_f64());
                eprintln!("Entries: {}", index.len());
                eprintln!("Bytes/key: {:.2}", index.bytes_per_key());
            }
        }
        VcfFormat::Plain | VcfFormat::Gzip => {
            let kbi_path = args
                .output
                .unwrap_or_else(|| args.input.with_extension("kbi"));

            eprintln!("Building KBI index: {:?}", kbi_path);
            let kbi_start = Instant::now();
            let index = build_kbi_index(&args.input, &kbi_path)?;
            eprintln!("KBI build time: {:.3}s", kbi_start.elapsed().as_secs_f64());
            eprintln!("Entries: {}", index.len());
            eprintln!("Bytes/key: {:.2}", index.bytes_per_key());

            if args.csi {
                eprintln!("Warning: CSI index requires BGZF compression. Use bgzip first.");
            }
        }
    }

    Ok(())
}
