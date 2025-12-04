use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use std::time::Instant;

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};

use kira_bio_tools::{
    build_csi_index, build_kbi_index, chr_id_to_name, chr_name_to_id, detect_format, fetch_line,
    CsiQuery, KbiIndex, Region, VcfFormat, VcfReader,
};

use kira_bio_tools::norm::normalize;

use memmap2::Mmap;
use rayon::prelude::*;

use kira_bio_tools::annotate::annotate_vcf_ani;
#[cfg(feature = "gpu")]
use kira_bio_tools::annotate_gpu::{annotate_vcf_ani_gpu, GpuAni};

use kira_bio_tools::filter_args::FilterArgs;
use kira_bio_tools::norm::turbo_norm_vcf;

use kira_bio_tools::annotate_index::AniIndex;

#[cfg(feature = "opencl")]
use kira_bio_tools::annotate_opencl::annotate_vcf_ani_opencl;
#[cfg(feature = "opencl")]
use kira_bio_tools::annotate_opencl::OpenCLAni;

#[derive(Parser)]
#[command(name = "kira-bt")]
#[command(about = "High-performance bioinformatics tools with full tabix compatibility")]
#[command(version, author)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    #[command(about = "Build annotation index (.ani) from VCF")]
    AnnotateIndex(AnnotateIndexArgs),

    #[command(about = "Annotate VCF using ANI index (bcftools annotate -a style)")]
    Annotate(AnnotateArgs),

    #[command(about = "Tabix-compatible indexer and query tool")]
    Tabix(TabixArgs),

    #[command(
        about = "Index a VCF file (extended functionality)",
        visible_alias = "idx"
    )]
    Index(IndexArgs),

    #[command(about = "Query regions from indexed VCF (extended functionality)")]
    Query(QueryArgs),

    #[command(about = "Display index statistics")]
    Stat(StatArgs),

    #[command(about = "List chromosome names from index")]
    List(ListArgs),

    #[command(about = "Print VCF header", visible_alias = "H")]
    Header(HeaderArgs),

    #[command(about = "Normalization", visible_alias = "N")]
    Norm(NormArgs),

    #[command(about = "Build ANI annotation index")]
    DbBuild(DbBuildArgs),

    #[command(about = "Filter")]
    Filter(FilterArgs),
}

#[derive(Parser)]
struct TabixArgs {
    #[arg(help = "Input file (for indexing) or indexed file (for querying)")]
    input: PathBuf,

    #[arg(help = "Regions to query (chr:start-end or chr:pos)")]
    regions: Vec<String>,

    #[arg(
        short = '0',
        long = "zero-based",
        help = "Position is 0-based half-open"
    )]
    zero_based: bool,

    #[arg(
        short = 'b',
        long = "begin",
        help = "Column of start chromosomal position [4]"
    )]
    begin_col: Option<usize>,

    #[arg(
        short = 'c',
        long = "comment",
        help = "Skip lines started with character [#]"
    )]
    comment_char: Option<char>,

    #[arg(
        short = 'C',
        long = "csi",
        help = "Produce CSI format index instead of TBI"
    )]
    csi: bool,

    #[arg(
        short = 'e',
        long = "end",
        help = "Column of end chromosomal position [5]"
    )]
    end_col: Option<usize>,

    #[arg(short = 'f', long = "force", help = "Force overwrite of index file")]
    force: bool,

    #[arg(
        short = 'm',
        long = "min-shift",
        help = "Set minimal interval size for CSI to 2^INT [14]"
    )]
    min_shift: Option<u8>,

    #[arg(
        short = 'p',
        long = "preset",
        help = "Input format: gff, bed, sam, vcf"
    )]
    preset: Option<String>,

    #[arg(short = 's', long = "sequence", help = "Column of sequence name [1]")]
    sequence_col: Option<usize>,

    #[arg(short = 'S', long = "skip-lines", help = "Skip first INT lines [0]")]
    skip_lines: Option<usize>,

    #[arg(
        short = 'h',
        long = "print-header",
        help = "Print also the header/meta lines"
    )]
    print_header: bool,

    #[arg(short = 'H', long = "only-header", help = "Print only the header")]
    only_header: bool,

    #[arg(short = 'l', long = "list-chroms", help = "List chromosome names")]
    list_chroms: bool,

    #[arg(
        short = 'r',
        long = "reheader",
        help = "Replace header with content of FILE"
    )]
    reheader: Option<PathBuf>,

    #[arg(
        short = 'R',
        long = "regions",
        help = "Restrict to regions listed in FILE"
    )]
    regions_file: Option<PathBuf>,

    #[arg(short = 't', long = "targets", help = "Similar to -R but sequential")]
    targets_file: Option<PathBuf>,

    #[arg(
        short = 'D',
        long = "no-download",
        help = "Do not download index for remote files"
    )]
    no_download: bool,

    #[arg(
        long = "cache",
        help = "Set BGZF block cache size in MB [10]",
        default_value = "10"
    )]
    cache_size: usize,

    #[arg(
        long = "regions-overlap",
        help = "Output region names before records (0=off, 1=on, 2=tab-separated)"
    )]
    regions_overlap: Option<u8>,

    #[arg(
        long = "verbosity",
        help = "Set log verbosity (0=silent, 3=default, 4+=debug)"
    )]
    verbosity: Option<u8>,

    #[arg(short = '@', long = "threads", help = "Number of threads [0]")]
    threads: Option<usize>,
}

#[derive(Parser)]
struct IndexArgs {
    #[arg(help = "Input VCF file (.vcf, .vcf.gz with BGZF)")]
    input: PathBuf,

    #[arg(short, long, help = "Output index file")]
    output: Option<PathBuf>,

    #[arg(
        short,
        long,
        default_value = "vcf",
        help = "Format preset (vcf, bed, gff, sam)"
    )]
    preset: String,

    #[arg(short, long, help = "Force overwrite existing index")]
    force: bool,

    #[arg(
        short = 's',
        long,
        default_value = "14",
        help = "Minimum bin shift (CSI)"
    )]
    min_shift: u8,

    #[arg(short = 'd', long, default_value = "5", help = "Bin depth (CSI)")]
    depth: u8,

    #[arg(short = 'C', long, help = "Generate CSI index (for BGZF files)")]
    csi: bool,

    #[arg(long, help = "Skip KBI index generation")]
    no_kbi: bool,
}

#[derive(Parser)]
struct QueryArgs {
    #[arg(help = "Indexed VCF file")]
    file: PathBuf,

    #[arg(help = "Regions to query (chr:start-end or chr:pos)")]
    regions: Vec<String>,

    #[arg(short = 'R', long, help = "File with regions (one per line)")]
    regions_file: Option<PathBuf>,

    #[arg(short, long, help = "Print only matching record count")]
    count: bool,

    #[arg(short = 'h', long, help = "Include header in output")]
    print_header: bool,

    #[arg(short = 'H', long, help = "Print only header")]
    only_header: bool,
}

#[derive(Parser)]
struct StatArgs {
    #[arg(help = "Index file (.kbi or .csi)")]
    index: PathBuf,
}

#[derive(Parser)]
struct ListArgs {
    #[arg(help = "Indexed VCF file")]
    file: PathBuf,
}

#[derive(Parser)]
struct HeaderArgs {
    #[arg(help = "VCF file")]
    file: PathBuf,
}

#[derive(Parser)]
struct NormArgs {
    #[arg(help = "Input VCF file (.vcf or .vcf.gz)")]
    input: PathBuf,

    #[arg(short, long, help = "Output VCF file (optional)")]
    output: Option<PathBuf>,
}

#[derive(Parser)]
struct AnnotateIndexArgs {
    #[arg(help = "Annotation VCF to index")]
    input: PathBuf,

    #[arg(short, long, help = "Output .ani file")]
    output: Option<PathBuf>,
}

#[derive(Parser)]
struct AnnotateArgs {
    #[arg(help = "Input VCF file (.vcf)")]
    input: PathBuf,

    #[arg(
        short = 'a',
        long = "annotations",
        help = "Annotation DB (.vcf or .ani)"
    )]
    annotations: PathBuf,

    #[arg(short, long, help = "Output file (optional)")]
    output: Option<PathBuf>,

    #[arg(long, help = "Use GPU (CUDA) for annotation")]
    gpu: bool,

    #[arg(long, help = "Use OpenCL backend (AMD/Intel/Apple/GPU/CPU)")]
    opencl: bool,
}

#[derive(Parser)]
struct DbBuildArgs {
    #[arg(help = "Input annotation database VCF")]
    input: PathBuf,

    #[arg(short, long, help = "Output ANI file")]
    output: Option<PathBuf>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let start = Instant::now();

    let result = match cli.command {
        Commands::Tabix(args) => cmd_tabix(args),
        Commands::Index(args) => cmd_index(args),
        Commands::Query(args) => cmd_query(args),
        Commands::Stat(args) => cmd_stat(args),
        Commands::List(args) => cmd_list(args),
        Commands::Header(args) => cmd_header(args),
        Commands::Norm(args) => cmd_norm(args),
        Commands::AnnotateIndex(args) => cmd_annotate_index(args),
        Commands::Annotate(args) => cmd_annotate(args),
        Commands::DbBuild(args) => cmd_db_build(args),
        Commands::Filter(args) => kira_bio_tools::filter::run_filter(&args),
    };

    if std::env::var("KIRA_BT_TIMING").is_ok() {
        let elapsed = start.elapsed();
        eprintln!("Total time: {:.3}s", elapsed.as_secs_f64());
    }

    result
}

fn cmd_tabix(args: TabixArgs) -> Result<()> {
    if args.list_chroms {
        return cmd_list_tabix(&args);
    }

    if args.only_header {
        return cmd_header_tabix(&args);
    }

    if args.regions.is_empty() && args.regions_file.is_none() && args.targets_file.is_none() {
        return cmd_index_tabix(args);
    }

    cmd_query_tabix(args)
}

fn cmd_index_tabix(args: TabixArgs) -> Result<()> {
    let format = detect_format(&args.input)?;

    if !args.input.exists() {
        anyhow::bail!("Input file does not exist: {:?}", args.input);
    }

    match format {
        VcfFormat::Bgzf => {
            let index_path = if args.csi {
                let mut p = args.input.clone();
                let name = p.file_name().unwrap().to_string_lossy().to_string();
                p.set_file_name(format!("{}.csi", name));
                p
            } else {
                let mut p = args.input.clone();
                let name = p.file_name().unwrap().to_string_lossy().to_string();
                p.set_file_name(format!("{}.tbi", name));
                p
            };

            if index_path.exists() && !args.force {
                anyhow::bail!(
                    "Index file already exists: {:?}. Use -f to overwrite.",
                    index_path
                );
            }

            eprintln!("Building index: {:?}", index_path);
            build_csi_index(&args.input, &index_path)?;

            let kbi_path = args.input.with_extension("kbi");
            if !kbi_path.exists() {
                eprintln!("Building KBI index: {:?}", kbi_path);
                build_kbi_index(&args.input, &kbi_path)?;
            }
        }
        VcfFormat::Plain | VcfFormat::Gzip => {
            anyhow::bail!(
                "File must be BGZF-compressed. Use bgzip to compress: bgzip {:?}",
                args.input
            );
        }
    }

    Ok(())
}

fn cmd_query_tabix(args: TabixArgs) -> Result<()> {
    let format = detect_format(&args.input)?;

    if format != VcfFormat::Bgzf {
        anyhow::bail!("Query requires BGZF-compressed file: {:?}", args.input);
    }

    let csi_path = {
        let mut p = args.input.clone();
        let name = p.file_name().unwrap().to_string_lossy().to_string();
        p.set_file_name(format!("{}.csi", name));
        p
    };

    let tbi_path = {
        let mut p = args.input.clone();
        let name = p.file_name().unwrap().to_string_lossy().to_string();
        p.set_file_name(format!("{}.tbi", name));
        p
    };

    let kbi_path = args.input.with_extension("kbi");

    let index_path = if csi_path.exists() {
        csi_path
    } else if tbi_path.exists() {
        tbi_path
    } else if kbi_path.exists() {
        kbi_path.clone()
    } else {
        anyhow::bail!(
            "No index found for {:?}. Run 'kira-bt tabix {:?}' first.",
            args.input,
            args.input
        );
    };

    let use_kbi = index_path.extension().map(|e| e == "kbi").unwrap_or(false);

    let mut regions = args.regions.clone();

    if let Some(ref regions_file) = args.regions_file {
        let file = File::open(regions_file)?;
        for line in BufReader::new(file).lines() {
            let line = line?;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            regions.push(line);
        }
    }

    if let Some(ref targets_file) = args.targets_file {
        let file = File::open(targets_file)?;
        for line in BufReader::new(file).lines() {
            let line = line?;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            regions.push(line);
        }
    }

    if regions.is_empty() {
        anyhow::bail!("No regions specified");
    }

    if let Some(ref reheader_file) = args.reheader {
        let reheader_content = std::fs::read_to_string(reheader_file)?;
        print!("{}", reheader_content);
    } else if args.print_header {
        print_vcf_header(&args.input)?;
    }

    if use_kbi {
        query_with_kbi(&args, &regions)?;
    } else {
        query_with_csi(&args, &regions, &index_path)?;
    }

    Ok(())
}

fn query_with_kbi(args: &TabixArgs, regions: &[String]) -> Result<()> {
    let kbi_path = args.input.with_extension("kbi");
    let index = KbiIndex::load(&kbi_path)?;

    for (_region_idx, region_str) in regions.iter().enumerate() {
        if let Some(overlap_mode) = args.regions_overlap {
            match overlap_mode {
                1 => println!("#{}", region_str),
                2 => {
                    print!("{}\t", region_str);
                }
                _ => {}
            }
        }

        let region = Region::parse(region_str)
            .ok_or_else(|| anyhow::anyhow!("Invalid region: {}", region_str))?;

        let chr_id = chr_name_to_id(&region.chr)
            .ok_or_else(|| anyhow::anyhow!("Unknown chromosome: {}", region.chr))?;

        let (start, end) = if args.zero_based {
            let s = region.start.unwrap_or(0);
            let e = region.end.unwrap_or(u32::MAX);
            (s, e.saturating_sub(1))
        } else {
            let s = region.start.unwrap_or(1);
            let e = region.end.unwrap_or(u32::MAX);
            (s, e)
        };

        let results = index.range(chr_id, start, end);

        for (_pos, offset) in results {
            let line = fetch_line(&args.input, offset)?;
            println!("{}", line);
        }
    }

    Ok(())
}

fn query_with_csi(args: &TabixArgs, regions: &[String], index_path: &PathBuf) -> Result<()> {
    let csi = CsiQuery::open(index_path)?;

    for region_str in regions {
        if let Some(overlap_mode) = args.regions_overlap {
            match overlap_mode {
                1 => println!("#{}", region_str),
                2 => {
                    print!("{}\t", region_str);
                }
                _ => {}
            }
        }

        let region = Region::parse(region_str)
            .ok_or_else(|| anyhow::anyhow!("Invalid region: {}", region_str))?;

        let chr_id = chr_name_to_id(&region.chr)
            .ok_or_else(|| anyhow::anyhow!("Unknown chromosome: {}", region.chr))?;

        let (start, end) = if args.zero_based {
            let s = region.start.unwrap_or(0);
            let e = region.end.unwrap_or(u32::MAX);
            (s, e.saturating_sub(1))
        } else {
            let s = region.start.unwrap_or(1);
            let e = region.end.unwrap_or(u32::MAX);
            (s, e)
        };

        let chunks = csi.query((chr_id - 1) as usize, start, end);

        for (chunk_start, _chunk_end) in chunks {
            let line = fetch_line(&args.input, chunk_start)?;
            println!("{}", line);
        }
    }

    Ok(())
}

fn cmd_list_tabix(args: &TabixArgs) -> Result<()> {
    let kbi_path = args.input.with_extension("kbi");

    if kbi_path.exists() {
        let index = KbiIndex::load(&kbi_path)?;

        for chr_id in 1..=25u8 {
            if let Some(name) = chr_id_to_name(chr_id) {
                let results = index.range(chr_id, 0, u32::MAX);
                if !results.is_empty() {
                    println!("{}", name);
                }
            }
        }
    } else {
        let mut reader = VcfReader::open(&args.input)?;
        let _ = reader.header()?;

        for name in reader.reference_sequences() {
            println!("{}", name);
        }
    }

    Ok(())
}

fn cmd_header_tabix(args: &TabixArgs) -> Result<()> {
    print_vcf_header(&args.input)
}

fn cmd_index(args: IndexArgs) -> Result<()> {
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

fn cmd_query(args: QueryArgs) -> Result<()> {
    if args.only_header {
        return cmd_header(HeaderArgs { file: args.file });
    }

    let kbi_path = args.file.with_extension("kbi");
    let csi_path = {
        let mut p = args.file.clone();
        let name = p.file_name().unwrap().to_string_lossy().to_string();
        p.set_file_name(format!("{}.csi", name));
        p
    };

    let use_kbi = kbi_path.exists();
    let use_csi = csi_path.exists() && !use_kbi;

    if !use_kbi && !use_csi {
        anyhow::bail!("No index found. Run 'kira-bt index {:?}' first.", args.file);
    }

    let mut regions = args.regions.clone();
    if let Some(ref regions_file) = args.regions_file {
        let file = File::open(regions_file)?;
        for line in BufReader::new(file).lines() {
            regions.push(line?);
        }
    }

    if regions.is_empty() {
        anyhow::bail!("No regions specified");
    }

    if args.print_header {
        print_vcf_header(&args.file)?;
    }

    let mut total_count = 0usize;

    if use_kbi {
        let index = KbiIndex::load(&kbi_path)?;

        for region_str in &regions {
            let region = Region::parse(region_str)
                .ok_or_else(|| anyhow::anyhow!("Invalid region: {}", region_str))?;

            let chr_id = chr_name_to_id(&region.chr)
                .ok_or_else(|| anyhow::anyhow!("Unknown chromosome: {}", region.chr))?;

            let start = region.start.unwrap_or(0);
            let end = region.end.unwrap_or(u32::MAX);

            let results = index.range(chr_id, start, end);
            total_count += results.len();

            if !args.count {
                for (_pos, offset) in results {
                    let line = fetch_line(&args.file, offset)?;
                    println!("{}", line);
                }
            }
        }
    } else if use_csi {
        let csi = CsiQuery::open(&csi_path)?;

        for region_str in &regions {
            let region = Region::parse(region_str)
                .ok_or_else(|| anyhow::anyhow!("Invalid region: {}", region_str))?;

            let chr_id = chr_name_to_id(&region.chr)
                .ok_or_else(|| anyhow::anyhow!("Unknown chromosome: {}", region.chr))?;

            let start = region.start.unwrap_or(0);
            let end = region.end.unwrap_or(u32::MAX);

            let chunks = csi.query((chr_id - 1) as usize, start, end);

            for (chunk_start, _chunk_end) in chunks {
                let line = fetch_line(&args.file, chunk_start)?;
                if !args.count {
                    println!("{}", line);
                }
                total_count += 1;
            }
        }
    }

    if args.count {
        println!("{}", total_count);
    }

    Ok(())
}

fn cmd_stat(args: StatArgs) -> Result<()> {
    let file_size = fs::metadata(&args.index)?.len();

    if args.index.extension().map(|e| e == "kbi").unwrap_or(false) {
        let index = KbiIndex::load(&args.index)?;

        println!("Index Statistics (KBI)");
        println!("======================");
        println!("File:          {:?}", args.index);
        println!(
            "File size:     {} bytes ({:.2} MB)",
            file_size,
            file_size as f64 / 1024.0 / 1024.0
        );
        println!("Entries:       {}", index.len());
        println!(
            "Memory usage:  {} bytes ({:.2} MB)",
            index.memory_usage(),
            index.memory_usage() as f64 / 1024.0 / 1024.0
        );
        println!("Bytes/key:     {:.2}", index.bytes_per_key());
    } else if args
        .index
        .extension()
        .map(|e| e == "csi" || e == "tbi")
        .unwrap_or(false)
    {
        let _csi = CsiQuery::open(&args.index)?;

        println!("Index Statistics (CSI/TBI)");
        println!("==========================");
        println!("File:          {:?}", args.index);
        println!(
            "File size:     {} bytes ({:.2} MB)",
            file_size,
            file_size as f64 / 1024.0 / 1024.0
        );
        println!("Format:        CSI/TBI (tabix-compatible)");
    } else {
        anyhow::bail!("Unknown index format");
    }

    Ok(())
}

fn cmd_list(args: ListArgs) -> Result<()> {
    let kbi_path = args.file.with_extension("kbi");

    if kbi_path.exists() {
        let index = KbiIndex::load(&kbi_path)?;

        for chr_id in 1..=25u8 {
            if let Some(name) = chr_id_to_name(chr_id) {
                let results = index.range(chr_id, 0, u32::MAX);
                if !results.is_empty() {
                    println!("{}", name);
                }
            }
        }
    } else {
        let mut reader = VcfReader::open(&args.file)?;
        let _ = reader.header()?;

        for name in reader.reference_sequences() {
            println!("{}", name);
        }
    }

    Ok(())
}

fn cmd_header(args: HeaderArgs) -> Result<()> {
    print_vcf_header(&args.file)
}

macro_rules! stage {
    ($name:expr, $block:block) => {{
        let __s = Instant::now();
        let __r = { $block };
        eprintln!("[norm] {}: {:.6}s", $name, __s.elapsed().as_secs_f64());
        __r
    }};
}

fn cmd_annotate(args: AnnotateArgs) -> Result<()> {
    let out = args.output.clone().unwrap_or_else(|| {
        let mut p = args.input.clone();
        p.set_extension("annot.vcf");
        p
    });

    let ani_path = if args.annotations.extension().unwrap_or_default() == "ani" {
        args.annotations.clone()
    } else {
        let mut p = args.annotations.clone();
        p.set_extension("ani");
        p
    };

    if !ani_path.exists() {
        anyhow::bail!("Annotation index not found: {:?}", ani_path);
    }

    eprintln!("[annotate] ANI = {:?}", ani_path);

    // ------------------------------
    // GPU backend (CUDA)
    // ------------------------------
    #[cfg(feature = "gpu")]
    if args.gpu {
        eprintln!("[annotate] Using CUDA GPU backend…");

        let ani = AniIndex::open(&ani_path)?;
        let gpu = GpuAni::load(&ani)?;

        annotate_vcf_ani_gpu(&gpu, &ani, &args.input, &out)?;
        return Ok(());
    }

    // ------------------------------
    // OpenCL backend
    // ------------------------------
    #[cfg(feature = "opencl")]
    if args.opencl {
        eprintln!("[annotate] Using OpenCL backend…");

        let ani = AniIndex::open(&ani_path)?;
        let gpu = OpenCLAni::new(&ani)?;

        annotate_vcf_ani_opencl(&gpu, &ani, &args.input, &out)?;
        return Ok(());
    }

    // ------------------------------
    // CPU fallback
    // ------------------------------
    annotate_vcf_ani(&ani_path, &args.input, &out)?;
    Ok(())
}

fn cmd_annotate_index(args: AnnotateIndexArgs) -> Result<()> {
    let out = args.output.clone().unwrap_or_else(|| {
        let mut p = args.input.clone();
        p.set_extension("ani");
        p
    });

    eprintln!("[annotate-index] Input  = {:?}", args.input);
    eprintln!("[annotate-index] Output = {:?}", out);

    kira_bio_tools::annotate_index::build_ani_index(&args.input, &out)?;

    Ok(())
}

fn cmd_db_build(args: DbBuildArgs) -> Result<()> {
    let out = args.output.clone().unwrap_or_else(|| {
        let mut p = args.input.clone();
        p.set_extension("ani");
        p
    });

    eprintln!("[db-build] Input: {:?}", args.input);
    eprintln!("[db-build] Output: {:?}", out);

    kira_bio_tools::annotate_index::build_ani_index(&args.input, &out)?;

    eprintln!("[db-build] Done");
    Ok(())
}

fn cmd_norm(args: NormArgs) -> Result<()> {
    let fmt = detect_format(&args.input)?;
    if fmt != VcfFormat::Plain {
        anyhow::bail!("Turbo mode currently supports only plain VCF");
    }

    let out_path = args.output.clone().unwrap_or_else(|| {
        let mut p = args.input.clone();
        p.set_extension("norm.vcf");
        p
    });

    turbo_norm_vcf(&args.input, &out_path)?;

    Ok(())
}

fn print_vcf_header(path: &PathBuf) -> Result<()> {
    let mut reader = VcfReader::open(path)?;
    let headers = reader.header()?;

    for line in headers {
        println!("{}", line);
    }

    Ok(())
}
