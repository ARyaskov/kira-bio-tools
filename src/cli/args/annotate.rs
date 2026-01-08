use clap::{ArgGroup, Parser};
use std::path::PathBuf;

#[derive(Parser)]
pub struct AnnotateIndexArgs {
    #[arg(help = "Annotation VCF to index")]
    pub input: PathBuf,

    #[arg(short, long, help = "Output .ani file")]
    pub output: Option<PathBuf>,
}

#[derive(Parser)]
#[command(group(
    ArgGroup::new("anno_source")
        .required(true)
        .args(["annotations", "ani"])
))]
pub struct AnnotateArgs {
    #[arg(help = "Input VCF file (.vcf)")]
    pub input: PathBuf,

    #[arg(
        short = 'a',
        long = "annotations",
        help = "Annotation DB (.vcf, .tab, or .ani)"
    )]
    pub annotations: Option<PathBuf>,

    #[arg(long = "ani", help = "Annotation index (.ani), no auto-build")]
    pub ani: Option<PathBuf>,

    #[arg(short = 'o', long = "output", help = "Output file (optional)")]
    pub output: Option<PathBuf>,

    #[arg(
        short = 'c',
        long = "columns",
        help = "Column specification for TAB files (e.g., CHROM,POS,REF,ALT,ID,+INFO/AC)"
    )]
    pub columns: Option<String>,

    #[arg(
        short = 'h',
        long = "header-lines",
        help = "Header file (.hdr) with VCF header lines"
    )]
    pub header_lines: Option<PathBuf>,

    #[arg(long, help = "Use GPU (CUDA) for annotation")]
    pub gpu: bool,

    #[arg(long, help = "Use OpenCL backend (AMD/Intel/Apple/GPU/CPU)")]
    pub opencl: bool,

    #[arg(long, help = "BGZF compression level (0-9, default 1)", value_parser = clap::value_parser!(u32))]
    pub bgzf_level: Option<u32>,

    #[arg(long, help = "Use cached plain VCF in CWD or create it from .gz input")]
    pub cache_plain: bool,

    #[arg(long, help = "Write plain VCF first, then BGZF-compress")]
    pub bgzf_after: bool,

    #[arg(long, help = "Use memory-mapped output for plain VCF")]
    pub mmap_output: bool,

    #[arg(long, help = "Do not flush mmap output explicitly")]
    pub mmap_no_flush: bool,

    #[arg(long, help = "Write output to RAM only (no disk)")]
    pub ram_output: bool,

    #[arg(
        long,
        help = "Max RAM output size in MB (unused)",
        default_value_t = 1024
    )]
    pub ram_max_mb: u32,
}

#[derive(Parser)]
#[command(group(
    ArgGroup::new("anno_source")
        .required(true)
        .args(["annotations", "ani"])
))]
pub struct AnnotateServeArgs {
    #[arg(
        short = 'a',
        long = "annotations",
        help = "Annotation DB (.vcf, .tab, or .ani)"
    )]
    pub annotations: Option<PathBuf>,

    #[arg(long = "ani", help = "Annotation index (.ani), no auto-build")]
    pub ani: Option<PathBuf>,

    #[arg(short = 'c', long = "columns", help = "Default column specification")]
    pub columns: Option<String>,

    #[arg(long, help = "Use GPU (CUDA) for annotation")]
    pub gpu: bool,

    #[arg(long, help = "Use OpenCL backend")]
    pub opencl: bool,

    #[arg(long, help = "BGZF compression level (0-9, default 1)", value_parser = clap::value_parser!(u32))]
    pub bgzf_level: Option<u32>,

    #[arg(long, help = "Use cached plain VCF in CWD or create it from .gz input")]
    pub cache_plain: bool,

    #[arg(long, help = "Write plain VCF first, then BGZF-compress")]
    pub bgzf_after: bool,

    #[arg(long, help = "Use memory-mapped output for plain VCF")]
    pub mmap_output: bool,

    #[arg(long, help = "Do not flush mmap output explicitly")]
    pub mmap_no_flush: bool,

    #[arg(long, help = "Write output to RAM only (no disk)")]
    pub ram_output: bool,

    #[arg(
        long,
        help = "Max RAM output size in MB (unused)",
        default_value_t = 1024
    )]
    pub ram_max_mb: u32,
}

#[derive(Parser)]
pub struct DbBuildArgs {
    #[arg(help = "Input annotation database VCF")]
    pub input: PathBuf,

    #[arg(short, long, help = "Output ANI file")]
    pub output: Option<PathBuf>,
}
