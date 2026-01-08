use anyhow::Result;
use clap::{Parser, Subcommand};
use std::time::Instant;

use kira_bio_tools::cli::args::*;
use kira_bio_tools::cli::commands::*;

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

    #[command(about = "OpenCL annotate server (stdin commands)")]
    AnnotateServe(AnnotateServeArgs),

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
        Commands::AnnotateServe(args) => cmd_annotate_serve(args),
        Commands::DbBuild(args) => cmd_db_build(args),
        Commands::Filter(args) => cmd_filter(&args),
    };

    if std::env::var("KIRA_BT_TIMING").is_ok() {
        let elapsed = start.elapsed();
        eprintln!("Total time: {:.3}s", elapsed.as_secs_f64());
    }

    result
}
