use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use std::time::Instant;

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};

use kira_bio_tools::{
    build_csi_index, build_kbi_index, chr_id_to_name, chr_name_to_id, detect_format,
    fetch_line, CsiQuery, KbiIndex, Region, VcfFormat, VcfReader,
};

#[derive(Parser)]
#[command(name = "kira-vcf")]
#[command(about = "High-performance VCF indexer with tabix compatibility")]
#[command(version, author)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    #[command(
        about = "Index a VCF file (tabix-compatible for BGZF files)",
        visible_alias = "idx"
    )]
    Index(IndexArgs),

    #[command(about = "Query regions from indexed VCF")]
    Query(QueryArgs),

    #[command(about = "Display index statistics")]
    Stat(StatArgs),

    #[command(about = "List chromosome names from index")]
    List(ListArgs),

    #[command(about = "Print VCF header", visible_alias = "H")]
    Header(HeaderArgs),

    #[command(about = "Print indexed regions (like tabix -R)")]
    Regions(RegionsArgs),
}

#[derive(Parser)]
struct IndexArgs {
    #[arg(help = "Input VCF file (.vcf, .vcf.gz with BGZF)")]
    input: PathBuf,

    #[arg(short, long, help = "Output index file")]
    output: Option<PathBuf>,

    #[arg(short, long, default_value = "vcf", help = "Format preset (vcf, bed, gff, sam)")]
    preset: String,

    #[arg(short, long, help = "Force overwrite existing index")]
    force: bool,

    #[arg(short = 's', long, default_value = "14", help = "Minimum bin shift (CSI)")]
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
struct RegionsArgs {
    #[arg(help = "Indexed VCF file")]
    file: PathBuf,

    #[arg(short = 'R', long, help = "Regions file")]
    regions_file: PathBuf,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let start = Instant::now();

    let result = match cli.command {
        Commands::Index(args) => cmd_index(args),
        Commands::Query(args) => cmd_query(args),
        Commands::Stat(args) => cmd_stat(args),
        Commands::List(args) => cmd_list(args),
        Commands::Header(args) => cmd_header(args),
        Commands::Regions(args) => cmd_regions(args),
    };

    let elapsed = start.elapsed();
    eprintln!("Total time: {:.3}s", elapsed.as_secs_f64());

    result
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
            let kbi_path = args.output.unwrap_or_else(|| args.input.with_extension("kbi"));

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
        anyhow::bail!("No index found. Run 'kira-vcf index {:?}' first.", args.file);
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
                for (pos, offset) in results {
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
        println!("File size:     {} bytes ({:.2} MB)", file_size, file_size as f64 / 1024.0 / 1024.0);
        println!("Entries:       {}", index.len());
        println!("Memory usage:  {} bytes ({:.2} MB)", index.memory_usage(), index.memory_usage() as f64 / 1024.0 / 1024.0);
        println!("Bytes/key:     {:.2}", index.bytes_per_key());
    } else if args.index.extension().map(|e| e == "csi").unwrap_or(false) {
        let _csi = CsiQuery::open(&args.index)?;

        println!("Index Statistics (CSI)");
        println!("======================");
        println!("File:          {:?}", args.index);
        println!("File size:     {} bytes ({:.2} MB)", file_size, file_size as f64 / 1024.0 / 1024.0);
        println!("Format:        CSI v1 (tabix-compatible)");
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

fn cmd_regions(args: RegionsArgs) -> Result<()> {
    let query_args = QueryArgs {
        file: args.file,
        regions: Vec::new(),
        regions_file: Some(args.regions_file),
        count: false,
        print_header: false,
        only_header: false,
    };
    cmd_query(query_args)
}

fn print_vcf_header(path: &PathBuf) -> Result<()> {
    let mut reader = VcfReader::open(path)?;
    let headers = reader.header()?;

    for line in headers {
        println!("{}", line);
    }

    Ok(())
}