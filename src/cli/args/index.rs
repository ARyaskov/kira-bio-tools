use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct IndexArgs {
    pub input: PathBuf,

    #[arg(short = 't', long = "tbi")]
    pub tbi: bool,

    #[arg(short = 'c', long = "csi", default_value_t = true)]
    pub csi: bool,

    #[arg(short = 'm', long = "min-shift", default_value_t = 14)]
    pub min_shift: u32,

    #[arg(short = 'f', long = "force")]
    pub force: bool,

    #[arg(short = 'n', long = "nrecords")]
    pub nrecords: bool,

    #[arg(short = 's', long = "stats")]
    pub stats: bool,

    #[arg(long = "all")]
    pub all: bool,

    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    #[arg(long = "threads", default_value_t = 0)]
    pub threads: usize,

    #[arg(short = 'v', long = "verbosity", default_value_t = 1)]
    pub verbosity: u8,

    #[arg(last = true, allow_hyphen_values = true)]
    pub passthrough: Vec<String>,
}

#[derive(Parser)]
pub struct RegionIndexArgs {
    #[arg(help = "Input VCF file (.vcf, .vcf.gz with BGZF)")]
    pub input: PathBuf,

    #[arg(short, long, help = "Output index file")]
    pub output: Option<PathBuf>,

    #[arg(short, long, default_value = "vcf", help = "Format preset (vcf, bed, gff, sam)")]
    pub preset: String,

    #[arg(short, long, help = "Force overwrite existing index")]
    pub force: bool,

    #[arg(short = 's', long, default_value = "14", help = "Minimum bin shift (CSI)")]
    pub min_shift: u8,

    #[arg(short = 'd', long, default_value = "5", help = "Bin depth (CSI)")]
    pub depth: u8,

    #[arg(short = 'C', long, help = "Generate CSI index (for BGZF files)")]
    pub csi: bool,

    #[arg(long, help = "Skip KBI index generation")]
    pub no_kbi: bool,
}
