//! CLI tool for VCF indexing and querying.
//!
//! Commands:
//! - index: Build index from VCF file
//! - query: Look up positions in index
//! - stat: Show index statistics
//! - range: Query by genomic interval

use std::fs::File;
use std::io::{BufRead, BufReader, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::time::Instant;

use clap::{Parser, Subcommand};
use memmap2::MmapOptions;

use kira_bio_tools::{GenomicKey, VcfIndex, VcfIndexBuilder, chr_name_to_id};

/// CLI definition
#[derive(Parser)]
#[command(name = "kira-vcf")]
#[command(about = "High-performance VCF indexer using learned indexes")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Build index from VCF file
    Index {
        /// Input VCF file
        #[arg(short, long)]
        input: PathBuf,

        /// Output index file (.kbi)
        #[arg(short, long)]
        output: PathBuf,

        /// Show progress every N records
        #[arg(long, default_value = "100000")]
        progress: usize,
    },

    /// Query positions in index
    Query {
        /// Index file (.kbi)
        #[arg(short, long)]
        index: PathBuf,

        /// VCF file (optional, for retrieving lines)
        #[arg(short, long)]
        vcf: Option<PathBuf>,

        /// Positions to query (chr[:pos] format, comma-separated)
        #[arg(short, long)]
        positions: String,
    },

    /// Show index statistics
    Stat {
        /// Index file (.kbi)
        #[arg(short, long)]
        index: PathBuf,
    },

    /// Query range of positions
    Range {
        /// Index file (.kbi)
        #[arg(short, long)]
        index: PathBuf,

        /// Chromosome
        #[arg(short, long)]
        chr: String,

        /// Start position
        #[arg(short, long)]
        start: u32,

        /// End position
        #[arg(short, long)]
        end: u32,
    },
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Index {
            input,
            output,
            progress,
        } => {
            cmd_index(&input, &output, progress)?;
        }
        Commands::Query {
            index,
            vcf,
            positions,
        } => {
            cmd_query(&index, vcf.as_deref(), &positions)?;
        }
        Commands::Stat { index } => {
            cmd_stat(&index)?;
        }
        Commands::Range {
            index,
            chr,
            start,
            end,
        } => {
            cmd_range(&index, &chr, start, end)?;
        }
    }

    Ok(())
}

/// Parse chromosome name from raw bytes without allocations.
///
/// Supports:
/// - "1".."22"
/// - "chr1".."chr22" (any case, e.g. "CHR1")
/// - "X", "Y", "M", "MT" with optional "chr" prefix
fn parse_chr_id_bytes(bytes: &[u8]) -> Option<u8> {
    if bytes.is_empty() {
        return None;
    }

    // Strip optional "chr" / "CHR" prefix
    let mut i = 0usize;
    if bytes.len() >= 3 {
        let c0 = bytes[0];
        let c1 = bytes[1];
        let c2 = bytes[2];
        if (c0 == b'c' || c0 == b'C') && (c1 == b'h' || c1 == b'H') && (c2 == b'r' || c2 == b'R') {
            i = 3;
        }
    }

    let rem = &bytes[i..];
    if rem.is_empty() {
        return None;
    }

    // Special chromosomes
    if rem.len() == 1 {
        match rem[0] {
            b'X' | b'x' => return Some(23),
            b'Y' | b'y' => return Some(24),
            b'M' | b'm' => return Some(25),
            _ => {}
        }
    } else if rem.len() == 2 {
        // MT / mt
        if (rem[0] == b'M' || rem[0] == b'm') && (rem[1] == b'T' || rem[1] == b't') {
            return Some(25);
        }
    }

    // Numeric chromosomes 1..22
    let mut n: u32 = 0;
    for &b in rem {
        if b < b'0' || b > b'9' {
            return None;
        }
        n = n.saturating_mul(10).saturating_add((b - b'0') as u32);
        if n > 22 {
            return None;
        }
    }

    if (1..=22).contains(&n) {
        Some(n as u8)
    } else {
        None
    }
}

/// Parse POS field from raw bytes as u32.
fn parse_pos_bytes(bytes: &[u8]) -> Option<u32> {
    if bytes.is_empty() {
        return None;
    }

    let mut n: u64 = 0;
    for &b in bytes {
        if b < b'0' || b > b'9' {
            return None;
        }
        n = n * 10 + (b - b'0') as u64;
        if n > u32::MAX as u64 {
            return None;
        }
    }

    Some(n as u32)
}

/// Build index from VCF file using mmap + zero-allocation ASCII parser.
fn cmd_index(input: &PathBuf, output: &PathBuf, progress: usize) -> anyhow::Result<()> {
    eprintln!("[index] Building index from: {}", input.display());
    let start_time = Instant::now();

    eprintln!("[index] Opening VCF...");
    let file = File::open(input)?;
    let file_size = file.metadata()?.len();
    eprintln!(
        "[index] VCF opened, file size: {:.2} MB",
        file_size as f64 / 1024.0 / 1024.0
    );

    eprintln!("[index] Memory-mapping VCF...");
    let mmap = unsafe { MmapOptions::new().map(&file)? };
    let data: &[u8] = &mmap;
    let len = data.len();

    eprintln!("[index] Initializing builder...");
    let mut builder = VcfIndexBuilder::with_capacity(10_000_000);

    let mut byte_offset: u64 = 0;
    let mut count: usize = 0;
    let mut line_count: usize = 0;

    eprintln!("[index] Starting VCF scan (mmap + byte parser)...");
    let scan_start = Instant::now();

    let mut pos_in_file: usize = 0;
    while pos_in_file < len {
        let line_start = pos_in_file;

        // Find end of line ('\n') or EOF
        let rel = match data[line_start..].iter().position(|&b| b == b'\n') {
            Some(p) => p,
            None => len - line_start,
        };
        let line_end = line_start + rel;
        let next_pos = if line_end < len && data[line_end] == b'\n' {
            line_end + 1
        } else {
            line_end
        };

        let line = &data[line_start..line_end];
        line_count += 1;
        byte_offset = line_start as u64;

        if line_count % 10_000 == 0 {
            let elapsed = scan_start.elapsed().as_secs_f64();
            eprintln!(
                "[index] Scan progress: {:>8} lines, {:>8} variants, offset={}, {:.2}s elapsed",
                line_count, count, byte_offset, elapsed
            );
        }

        // Skip empty lines
        if line.is_empty() {
            pos_in_file = next_pos;
            continue;
        }

        // Header lines start with '#'
        if line[0] == b'#' {
            pos_in_file = next_pos;
            continue;
        }

        // Parse VCF line: CHROM\tPOS\t...
        // We only need first two fields, so split at first two tabs
        let mut field_iter = line.splitn(3, |&b| b == b'\t');
        let chrom_bytes = field_iter.next().unwrap_or(&[]);
        let pos_bytes = field_iter.next().unwrap_or(&[]);

        if !chrom_bytes.is_empty() && !pos_bytes.is_empty() {
            if let (Some(chr_id), Some(pos)) =
                (parse_chr_id_bytes(chrom_bytes), parse_pos_bytes(pos_bytes))
            {
                let key = GenomicKey::new(chr_id, pos);
                builder.add(key, byte_offset)?;
                count += 1;

                if count % progress == 0 {
                    let chrom_str = std::str::from_utf8(chrom_bytes).unwrap_or("<?>");
                    eprintln!(
                        "[index] Records: {:>8} (line {}), last key: {}:{} offset {}",
                        count, line_count, chrom_str, pos, byte_offset
                    );
                }
            }
        }

        pos_in_file = next_pos;
    }

    let scan_time = scan_start.elapsed();
    eprintln!(
        "[index] Finished VCF scan: {} lines, {} records, {:.2}s",
        line_count,
        count,
        scan_time.as_secs_f64()
    );

    eprintln!(
        "[index] Building MPH index for {} entries (last_offset={}, capacity_hint={})...",
        count, byte_offset, 10_000_000
    );
    let build_start = Instant::now();
    let index = builder.build()?;
    let build_time = build_start.elapsed();
    eprintln!(
        "[index] MPH build finished in {:.2}s",
        build_time.as_secs_f64()
    );

    eprintln!("[index] Saving index to: {}", output.display());
    let save_start = Instant::now();
    index.save(output)?;
    let save_time = save_start.elapsed();
    eprintln!("[index] Save finished in {:.2}s", save_time.as_secs_f64());

    let total_time = start_time.elapsed();
    eprintln!("\n[index] Statistics:");
    eprintln!("  Total entries:    {}", index.len());
    eprintln!(
        "  Memory usage:     {:.2} MB",
        index.memory_usage() as f64 / 1024.0 / 1024.0
    );
    eprintln!("  Bytes per key:    {:.2}", index.bytes_per_key());
    eprintln!("  Scan time:        {:.2}s", scan_time.as_secs_f64());
    eprintln!("  Build time:       {:.2}s", build_time.as_secs_f64());
    eprintln!("  Save time:        {:.2}s", save_time.as_secs_f64());
    eprintln!("  Total time:       {:.2}s", total_time.as_secs_f64());
    if build_time.as_secs_f64() > 0.0 {
        eprintln!(
            "  Build rate:       {:.0} entries/s",
            count as f64 / build_time.as_secs_f64()
        );
    }

    Ok(())
}

fn cmd_query(index_path: &Path, vcf_path: Option<&Path>, positions: &str) -> anyhow::Result<()> {
    eprintln!("Loading index: {}", index_path.display());
    let index = VcfIndex::load_mmap(index_path)?;

    let mut vcf_file = match vcf_path {
        Some(path) => Some(File::open(path)?),
        None => None,
    };

    for pos_str in positions.split(',') {
        let pos_str = pos_str.trim();

        // allow "chr:1234" or just "1234"
        let (chr_id, pos) = if let Some((chr, p)) = pos_str.split_once(':') {
            let chr_id =
                chr_name_to_id(chr).ok_or_else(|| anyhow::anyhow!("Unknown chromosome {}", chr))?;
            let p = p.parse::<u32>()?;
            (chr_id, p)
        } else {
            (0u8, pos_str.parse::<u32>()?)
        };

        let key = GenomicKey::new(chr_id, pos);

        if let Some(offset) = index.get(key) {
            println!("{} → offset {}", pos_str, offset);

            if let Some(file) = vcf_file.as_mut() {
                file.seek(SeekFrom::Start(offset))?;
                let mut reader = BufReader::new(file);
                let mut line = String::new();
                reader.read_line(&mut line)?;
                print!("{}", line);
            }
        } else {
            println!("{} not found", pos_str);
        }
    }

    Ok(())
}

fn cmd_stat(index_path: &PathBuf) -> anyhow::Result<()> {
    let load_start = Instant::now();
    let index = VcfIndex::load_mmap(index_path)?;
    let load_time = load_start.elapsed();

    let file_size = std::fs::metadata(index_path)?.len();

    println!("Index Statistics");
    println!("================");
    println!("File:             {}", index_path.display());
    println!(
        "File size:        {:.2} MB",
        file_size as f64 / 1024.0 / 1024.0
    );
    println!("Entries:          {}", index.len());
    println!(
        "Memory usage:     {:.2} MB",
        index.memory_usage() as f64 / 1024.0 / 1024.0
    );
    println!("Bytes per key:    {:.2}", index.bytes_per_key());
    println!(
        "Load time:        {:.2}ms",
        load_time.as_secs_f64() * 1000.0
    );

    Ok(())
}

fn cmd_range(index_path: &PathBuf, chr: &str, start: u32, end: u32) -> anyhow::Result<()> {
    let index = VcfIndex::load_mmap(index_path)?;

    let chr_id =
        chr_name_to_id(chr).ok_or_else(|| anyhow::anyhow!("Unknown chromosome: {}", chr))?;

    let results = index.range(chr_id, start, end);

    println!(
        "Found {} positions in {}:{}-{}",
        results.len(),
        chr,
        start,
        end
    );
    for (pos, offset) in results {
        println!("{}\t{}\t{}", chr, pos, offset);
    }

    Ok(())
}
