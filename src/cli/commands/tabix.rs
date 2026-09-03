use anyhow::{Context, Result, bail};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::cli::args::TabixArgs;
use crate::csi::{BinIndex, IndexKind, IndexedVcfReader, build_index, find_index_for};
use crate::util::Region;
use crate::vcf::VcfReader;

pub fn cmd_tabix(args: TabixArgs) -> Result<()> {
    if args.list_chroms {
        return cmd_list_tabix(&args);
    }
    if args.only_header {
        return print_vcf_header(&args.input);
    }
    if args.regions.is_empty() && args.regions_file.is_none() && args.targets_file.is_none() {
        return cmd_index_tabix(&args);
    }
    cmd_query_tabix(&args)
}

fn index_path_for(input: &Path, kind: IndexKind) -> PathBuf {
    let mut s = input.as_os_str().to_os_string();
    s.push(match kind {
        IndexKind::Csi => ".csi",
        IndexKind::Tbi => ".tbi",
    });
    PathBuf::from(s)
}

fn cmd_index_tabix(args: &TabixArgs) -> Result<()> {
    if !args.input.exists() {
        bail!("Input file does not exist: {}", args.input.display());
    }
    if let Some(p) = &args.preset {
        if !p.eq_ignore_ascii_case("vcf") {
            bail!("tabix preset {p:?} is not supported; only VCF/BCF inputs can be indexed");
        }
    }
    let custom_cols = args.sequence_col.is_some_and(|c| c != 1)
        || args.begin_col.is_some_and(|c| c != 2)
        || args.end_col.is_some_and(|c| c != 0 && c != 2);
    if custom_cols {
        bail!("custom -s/-b/-e columns are not supported; only the VCF layout can be indexed");
    }
    let kind = if args.csi { IndexKind::Csi } else { IndexKind::Tbi };
    let out = index_path_for(&args.input, kind);
    if out.exists() && !args.force {
        bail!("the index file exists; use -f to overwrite: {}", out.display());
    }
    build_index(&args.input, &out, kind, args.min_shift)?;
    Ok(())
}

/// Regions from a `-R`/`-T` file: BED (0-based half-open, by suffix), two
/// columns (`CHROM POS`), three columns (`CHROM BEG END`, 1-based inclusive) or
/// region strings.
pub fn read_region_file(path: &Path) -> Result<Vec<(String, u32, u32)>> {
    let is_bed = path
        .to_string_lossy()
        .to_ascii_lowercase()
        .trim_end_matches(".gz")
        .ends_with(".bed");
    let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let reader: Box<dyn BufRead> = if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') || t.starts_with("track") || t.starts_with("browser") {
            continue;
        }
        let cols: Vec<&str> = t.split('\t').collect();
        if cols.len() >= 3 {
            let b: u32 = cols[1].trim().parse().with_context(|| format!("region file: bad start {:?}", cols[1]))?;
            let e: u32 = cols[2].trim().parse().with_context(|| format!("region file: bad end {:?}", cols[2]))?;
            let beg = if is_bed { b + 1 } else { b };
            out.push((cols[0].to_string(), beg.max(1), e));
        } else if cols.len() == 2 {
            let p: u32 = cols[1].trim().parse().with_context(|| format!("region file: bad position {:?}", cols[1]))?;
            let beg = if is_bed { p + 1 } else { p };
            out.push((cols[0].to_string(), beg.max(1), if is_bed { u32::MAX } else { p }));
        } else {
            let r = Region::parse(t).ok_or_else(|| anyhow::anyhow!("invalid region {t:?}"))?;
            let (b, e) = r.bounds();
            out.push((r.chr, b, e));
        }
    }
    Ok(out)
}

fn cmd_query_tabix(args: &TabixArgs) -> Result<()> {
    let mut reader = IndexedVcfReader::open(&args.input)
        .with_context(|| format!("[tabix] could not load index for {}", args.input.display()))?;

    let mut regions: Vec<(String, u32, u32)> = Vec::new();
    for s in &args.regions {
        let r = Region::parse_with(s, false).ok_or_else(|| anyhow::anyhow!("invalid region {s:?}"))?;
        let (mut b, mut e) = r.bounds();
        if args.zero_based && r.start.is_some() {
            b += 1;
            if e != u32::MAX {
                e = e.max(b);
            }
        }
        regions.push((r.chr, b, e));
    }
    if let Some(p) = &args.regions_file {
        regions.extend(read_region_file(p)?);
    }
    if let Some(p) = &args.targets_file {
        regions.extend(read_region_file(p)?);
    }
    if regions.is_empty() {
        bail!("No regions specified");
    }

    let stdout = std::io::stdout();
    let mut out = BufWriter::with_capacity(1 << 20, stdout.lock());
    if let Some(reheader_file) = &args.reheader {
        let content = std::fs::read_to_string(reheader_file)?;
        out.write_all(content.as_bytes())?;
        if !content.ends_with('\n') {
            out.write_all(b"\n")?;
        }
    } else if args.print_header {
        for h in reader.headers() {
            out.write_all(h.as_bytes())?;
            out.write_all(b"\n")?;
        }
    }

    for (chrom, beg, end) in &regions {
        let label = format!("{chrom}:{beg}-{end}");
        match args.regions_overlap {
            Some(1) => writeln!(out, "#{label}")?,
            _ => {}
        }
        let prefix = args.regions_overlap == Some(2);
        reader.query(chrom, *beg, *end, |line| {
            if prefix {
                out.write_all(label.as_bytes())?;
                out.write_all(b"\t")?;
            }
            out.write_all(line.as_bytes())?;
            out.write_all(b"\n")?;
            Ok(true)
        })?;
    }
    out.flush()?;
    Ok(())
}

fn cmd_list_tabix(args: &TabixArgs) -> Result<()> {
    if let Some(idx_path) = find_index_for(&args.input) {
        let idx = BinIndex::load(&idx_path)?;
        if !idx.names().is_empty() {
            for name in idx.names() {
                println!("{name}");
            }
            return Ok(());
        }
        // BCF index carries no names: fall back to the header in rid order.
    }
    let mut reader = VcfReader::open(&args.input)?;
    let _ = reader.header()?;
    for name in reader.reference_sequences()? {
        println!("{name}");
    }
    Ok(())
}

fn print_vcf_header(path: &Path) -> Result<()> {
    let mut reader = VcfReader::open(path)?;
    let headers = reader.header()?;
    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    for line in headers {
        out.write_all(line.as_bytes())?;
        out.write_all(b"\n")?;
    }
    out.flush()?;
    Ok(())
}
