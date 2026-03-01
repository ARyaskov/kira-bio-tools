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

    #[command(about = "Index VCF/BCF files", visible_alias = "idx")]
    Index(IndexArgs),

    #[command(
        about = "Index a VCF file (extended functionality)",
        visible_alias = "kindex"
    )]
    RegionIndex(RegionIndexArgs),

    #[command(about = "Transform VCF/BCF into user-defined formats")]
    Query(QueryCompatArgs),

    #[command(
        about = "Query regions from indexed VCF (extended functionality)",
        visible_alias = "rquery"
    )]
    RegionQuery(RegionQueryArgs),

    #[command(about = "Display index statistics")]
    Stat(StatArgs),

    #[command(about = "List chromosome names from index")]
    List(ListArgs),

    #[command(about = "Print VCF header", visible_alias = "H")]
    Header(HeaderArgs),

    #[command(about = "View VCF/BCF file headers")]
    Head(HeadArgs),

    #[command(about = "Normalization", visible_alias = "N")]
    Norm(NormArgs),

    #[command(about = "Build ANI annotation index")]
    DbBuild(DbBuildArgs),

    #[command(about = "Filter")]
    Filter(FilterArgs),

    #[command(about = "Check sample concordance, detect sample swaps and contamination")]
    Gtcheck(GtcheckArgs),

    #[command(about = "Haplotype aware consequence caller")]
    Csq(CsqArgs),

    #[command(about = "Convert VCF/BCF to other formats and back")]
    Convert(ConvertArgs),

    #[command(about = "Create consensus sequence by applying VCF variants")]
    Consensus(ConsensusArgs),

    #[command(about = "Concatenate VCF/BCF files from the same set of samples")]
    Concat(ConcatArgs),

    #[command(about = "Copy Number Variation caller")]
    Cnv(CnvArgs),

    #[command(about = "SNP/indel calling")]
    Call(CallArgs),

    #[command(about = "Multi-way pileup producing genotype likelihoods")]
    Mpileup(MpileupArgs),

    #[command(about = "Sort VCF/BCF files")]
    Sort(SortArgs),

    #[command(about = "Produce VCF/BCF stats (former vcfcheck)")]
    Stats(StatsArgs),

    #[command(about = "Identify runs of homo/auto-zygosity")]
    Roh(RohArgs),

    #[command(about = "Modify VCF/BCF header, change sample names")]
    Reheader(ReheaderArgs),

    #[command(about = "Detect contaminations and whole-chromosome aberrations")]
    Polysomy(PolysomyArgs),

    #[command(about = "Merge VCF/BCF files from non-overlapping sample sets")]
    Merge(MergeArgs),

    #[command(about = "Intersections of VCF/BCF files")]
    Isec(IsecArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let start = Instant::now();

    let result = match cli.command {
        Commands::Tabix(args) => cmd_tabix(args),
        Commands::Index(args) => cmd_index(args),
        Commands::RegionIndex(args) => cmd_region_index(args),
        Commands::Query(args) => cmd_query(args),
        Commands::RegionQuery(args) => cmd_region_query(args),
        Commands::Stat(args) => cmd_stat(args),
        Commands::List(args) => cmd_list(args),
        Commands::Header(args) => cmd_header(args),
        Commands::Head(args) => cmd_head(args),
        Commands::Norm(args) => cmd_norm(args),
        Commands::AnnotateIndex(args) => cmd_annotate_index(args),
        Commands::Annotate(args) => cmd_annotate(args),
        Commands::AnnotateServe(args) => cmd_annotate_serve(args),
        Commands::DbBuild(args) => cmd_db_build(args),
        Commands::Filter(args) => cmd_filter(&args),
        Commands::Gtcheck(args) => cmd_gtcheck(args),
        Commands::Csq(args) => cmd_csq(args),
        Commands::Convert(args) => cmd_convert(args),
        Commands::Consensus(args) => cmd_consensus(args),
        Commands::Concat(args) => cmd_concat(args),
        Commands::Cnv(args) => cmd_cnv(args),
        Commands::Call(args) => cmd_call(args),
        Commands::Mpileup(args) => cmd_mpileup(args),
        Commands::Sort(args) => cmd_sort(args),
        Commands::Stats(args) => cmd_stats(args),
        Commands::Roh(args) => cmd_roh(args),
        Commands::Reheader(args) => cmd_reheader(args),
        Commands::Polysomy(args) => cmd_polysomy(args),
        Commands::Merge(args) => cmd_merge(args),
        Commands::Isec(args) => cmd_isec(args),
    };

    if std::env::var("KIRA_BT_TIMING").is_ok() {
        let elapsed = start.elapsed();
        eprintln!("Total time: {:.3}s", elapsed.as_secs_f64());
    }

    result
}
