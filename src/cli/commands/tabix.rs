use anyhow::{Context, Result, bail};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::cli::args::TabixArgs;
use crate::csi::{
    BinIndex, IndexKind, IndexedVcfReader, TabixConf, build_index, build_tabix_index,
    find_index_for,
};
use crate::util::Region;
use crate::vcf::VcfReader;

pub fn cmd_tabix(args: TabixArgs) -> Result<()> {
    if args.list_chroms {
        return cmd_list_tabix(&args);
    }
    if args.only_header {
        return print_header_lines(&args.input, header_conf(&args)?);
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

/// The column layout for indexing: `-p`, else custom `-s/-b/-e` (generic,
/// gff-like defaults), else guessed from the file name; `-c`, `-S` and `-0`
/// adjust whichever applies.
fn index_conf(args: &TabixArgs) -> Result<TabixConf> {
    let custom = args.sequence_col.is_some() || args.begin_col.is_some() || args.end_col.is_some();
    let mut conf = if let Some(p) = &args.preset {
        TabixConf::preset(p)
            .ok_or_else(|| anyhow::anyhow!("unknown preset {p:?}; expected gff, bed, sam or vcf"))?
    } else if custom {
        TabixConf::GFF
    } else {
        TabixConf::from_extension(&args.input).ok_or_else(|| {
            anyhow::anyhow!(
                "cannot guess the format of {}; use -p gff|bed|sam|vcf or -s/-b/-e",
                args.input.display()
            )
        })?
    };
    if let Some(c) = args.sequence_col {
        conf.sc = c as i32;
    }
    if let Some(c) = args.begin_col {
        conf.bc = c as i32;
    }
    if let Some(c) = args.end_col {
        conf.ec = c as i32;
    }
    if let Some(c) = args.comment_char {
        conf.meta = c as u8;
    }
    if let Some(n) = args.skip_lines {
        conf.skip = n as i32;
    }
    if args.zero_based {
        conf.preset |= crate::csi::TBX_UCSC;
    }
    Ok(conf)
}

/// The layout for reading headers: the existing index's, else as for indexing.
fn header_conf(args: &TabixArgs) -> Result<TabixConf> {
    if let Some(idx) = find_index_for(&args.input)
        && let Ok(index) = BinIndex::load(&idx)
        && let Some(h) = index.tabix_header()
    {
        return Ok(TabixConf::from_header(h));
    }
    Ok(index_conf(args).unwrap_or(TabixConf::VCF))
}

fn cmd_index_tabix(args: &TabixArgs) -> Result<()> {
    if !args.input.exists() {
        bail!("Input file does not exist: {}", args.input.display());
    }
    let conf = index_conf(args)?;
    let kind = if args.csi {
        IndexKind::Csi
    } else {
        IndexKind::Tbi
    };
    let out = index_path_for(&args.input, kind);
    if out.exists() && !args.force {
        bail!(
            "the index file exists; use -f to overwrite: {}",
            out.display()
        );
    }
    if conf == TabixConf::VCF {
        // VCF and BCF (CSI only) through the VCF/BCF driver.
        build_index(&args.input, &out, kind, args.min_shift)?;
    } else {
        let idx = build_tabix_index(&args.input, &conf, kind, args.min_shift)?;
        idx.save(&out)
            .with_context(|| format!("write {}", out.display()))?;
    }
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
        if t.is_empty() || t.starts_with('#') || t.starts_with("track") || t.starts_with("browser")
        {
            continue;
        }
        let cols: Vec<&str> = t.split('\t').collect();
        if cols.len() >= 3 {
            let b: u32 = cols[1]
                .trim()
                .parse()
                .with_context(|| format!("region file: bad start {:?}", cols[1]))?;
            let e: u32 = cols[2]
                .trim()
                .parse()
                .with_context(|| format!("region file: bad end {:?}", cols[2]))?;
            let beg = if is_bed { b + 1 } else { b };
            out.push((cols[0].to_string(), beg.max(1), e));
        } else if cols.len() == 2 {
            let p: u32 = cols[1]
                .trim()
                .parse()
                .with_context(|| format!("region file: bad position {:?}", cols[1]))?;
            let beg = if is_bed { p + 1 } else { p };
            out.push((
                cols[0].to_string(),
                beg.max(1),
                if is_bed { u32::MAX } else { p },
            ));
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
        let r =
            Region::parse_with(s, false).ok_or_else(|| anyhow::anyhow!("invalid region {s:?}"))?;
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
        if let Some(1) = args.regions_overlap {
            writeln!(out, "#{label}")?
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

/// `-H`: the leading lines a layout treats as header (skipped lines and
/// lines starting with the comment character), for any tabix-indexed text.
fn print_header_lines(path: &Path, conf: TabixConf) -> Result<()> {
    let reader = crate::util::open_text_maybe_gz(path)?;
    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim_end_matches('\r');
        if !((i as i64) < conf.skip as i64 || line.as_bytes().first() == Some(&conf.meta)) {
            break;
        }
        out.write_all(line.as_bytes())?;
        out.write_all(b"\n")?;
    }
    out.flush()?;
    Ok(())
}
