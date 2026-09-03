use anyhow::{Context, Result, bail};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::cli::args::{IndexArgs, RegionIndexArgs};
use crate::csi::{BinIndex, IndexKind, build_index, build_index_in_memory, find_index_for};
use crate::vcf::header::HeaderInfo;
use crate::{VcfFormat, VcfReader, build_kbi_index, detect_format};

pub fn cmd_index(args: IndexArgs) -> Result<()> {
    let mut cfg = IndexCompatCfg {
        tbi: args.tbi,
        force: args.force,
        stats: args.stats,
        num: args.nrecords,
        all: args.all,
        min_shift: args.min_shift as u8,
        output: args.output.clone(),
    };
    apply_passthrough(&mut cfg, &args.passthrough)?;

    if cfg.stats {
        return print_stats(&args.input, cfg.all);
    }
    if cfg.num {
        return print_num(&args.input);
    }

    let kind = if cfg.tbi { IndexKind::Tbi } else { IndexKind::Csi };
    let out = resolve_index_output(&args.input, &cfg, kind);
    if out.exists() && !cfg.force {
        bail!("index file already exists (use -f to overwrite): {}", out.display());
    }
    if !(1..=30).contains(&cfg.min_shift) {
        bail!("--min-shift must be between 1 and 30");
    }
    build_index(&args.input, &out, kind, Some(cfg.min_shift))?;
    Ok(())
}

#[derive(Default)]
struct IndexCompatCfg {
    tbi: bool,
    force: bool,
    stats: bool,
    num: bool,
    all: bool,
    min_shift: u8,
    output: Option<PathBuf>,
}

fn apply_passthrough(cfg: &mut IndexCompatCfg, args: &[String]) -> Result<()> {
    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "-t" | "--tbi" => cfg.tbi = true,
            "-c" | "--csi" => cfg.tbi = false,
            "-f" | "--force" => cfg.force = true,
            "-s" | "--stats" => cfg.stats = true,
            "-n" | "--nrecords" => cfg.num = true,
            "-a" | "--all" => cfg.all = true,
            "-m" | "--min-shift" => {
                i += 1;
                let v = args.get(i).ok_or_else(|| anyhow::anyhow!("missing value for -m/--min-shift"))?;
                cfg.min_shift = v.parse().with_context(|| format!("bad --min-shift {v:?}"))?;
            }
            "-o" | "--output" => {
                i += 1;
                let p = args.get(i).ok_or_else(|| anyhow::anyhow!("missing value for -o/--output"))?;
                cfg.output = Some(PathBuf::from(p));
            }
            "--threads" => {
                i += 1;
            }
            other => {
                if let Some(v) = other.strip_prefix("-m") {
                    cfg.min_shift = v.parse().with_context(|| format!("bad --min-shift {v:?}"))?;
                } else {
                    bail!("index: unknown option {other:?}");
                }
            }
        }
        i += 1;
    }
    Ok(())
}

fn resolve_index_output(input: &Path, cfg: &IndexCompatCfg, kind: IndexKind) -> PathBuf {
    if let Some(p) = &cfg.output {
        return p.clone();
    }
    let mut s = input.as_os_str().to_os_string();
    s.push(match kind {
        IndexKind::Csi => ".csi",
        IndexKind::Tbi => ".tbi",
    });
    PathBuf::from(s)
}

fn strip_index_suffix(path: &Path) -> Option<PathBuf> {
    let s = path.to_string_lossy();
    for suf in [".csi", ".tbi"] {
        if let Some(b) = s.strip_suffix(suf) {
            return Some(PathBuf::from(b));
        }
    }
    None
}

/// The data file and its index for `-n`/`-s`: `input` may be the data file
/// or the index itself. Builds the index in memory when none is on disk.
fn load_data_and_index(input: &Path) -> Result<(PathBuf, BinIndex)> {
    if let Some(base) = strip_index_suffix(input) {
        let idx = BinIndex::load(input).with_context(|| format!("read index {}", input.display()))?;
        return Ok((base, idx));
    }
    if let Some(p) = find_index_for(input) {
        let idx = BinIndex::load(&p).with_context(|| format!("read index {}", p.display()))?;
        return Ok((input.to_path_buf(), idx));
    }
    let idx = build_index_in_memory(input, IndexKind::Csi, None)
        .with_context(|| format!("no index for {} and it could not be built", input.display()))?;
    Ok((input.to_path_buf(), idx))
}

/// Contig names in index order with their header lengths.
fn contig_table(data: &Path, idx: &BinIndex) -> Result<Vec<(String, Option<u64>)>> {
    let mut names: Vec<String> = idx.names().to_vec();
    let info = if data.exists() {
        let mut r = VcfReader::open(data)?;
        let h = r.header()?;
        Some(HeaderInfo::parse(&h))
    } else {
        None
    };
    if names.is_empty() {
        // BCF index: names come from the header in rid order.
        if let Some(i) = &info {
            names = i.contigs.names().to_vec();
        }
    }
    Ok(names
        .into_iter()
        .map(|n| {
            let len = info.as_ref().and_then(|i| i.contigs.id(&n).and_then(|id| i.contigs.length(id)));
            (n, len)
        })
        .collect())
}

fn print_num(input: &Path) -> Result<()> {
    let (_, idx) = load_data_and_index(input)?;
    println!("{}", idx.total_records());
    Ok(())
}

fn print_stats(input: &Path, show_all: bool) -> Result<()> {
    let (data, idx) = load_data_and_index(input)?;
    let table = contig_table(&data, &idx)?;
    let mut printed: HashMap<usize, bool> = HashMap::new();
    for (i, (name, len)) in table.iter().enumerate() {
        let n = idx.n_records(i).unwrap_or(0);
        printed.insert(i, true);
        if show_all || n > 0 {
            let l = len.map(|v| v.to_string()).unwrap_or_else(|| ".".into());
            println!("{name}\t{l}\t{n}");
        }
    }
    // References present in the index but not named (should not happen).
    for i in table.len()..idx.n_refs() {
        let n = idx.n_records(i).unwrap_or(0);
        if show_all || n > 0 {
            println!("{i}\t.\t{n}");
        }
    }
    Ok(())
}

pub fn cmd_region_index(args: RegionIndexArgs) -> Result<()> {
    let format = detect_format(&args.input)?;

    eprintln!("Input: {:?}", args.input);
    eprintln!("Format: {:?}", format);

    match format {
        VcfFormat::Bgzf => {
            let csi_path = args.output.clone().unwrap_or_else(|| {
                let mut s = args.input.as_os_str().to_os_string();
                s.push(".csi");
                PathBuf::from(s)
            });

            if !args.no_kbi || args.csi {
                eprintln!("Building CSI index: {:?}", csi_path);
                let csi_start = Instant::now();
                build_index(&args.input, &csi_path, IndexKind::Csi, Some(args.min_shift))?;
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
