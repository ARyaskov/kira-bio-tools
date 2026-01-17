// Command-line arguments for VCF filtering

use clap::Args;
use std::path::PathBuf;

#[derive(Args, Debug)]
pub struct FilterArgs {
    #[arg(short = 'i', long = "include", value_name = "EXPR", conflicts_with = "exclude")]
    pub include: Option<String>,

    #[arg(short = 'e', long = "exclude", value_name = "EXPR", conflicts_with = "include")]
    pub exclude: Option<String>,

    #[arg(long = "expr", value_name = "EXPR", conflicts_with_all = ["include", "exclude"])]
    pub expr: Option<String>,

    #[arg(short = 's', long = "soft-filter")]
    pub soft_filter: Option<String>,

    #[arg(short = 'm', long = "mode")]
    pub mode: Option<String>,

    #[arg(short = 'S', long = "set-GTs")]
    pub set_gts: Option<String>,

    #[arg(short = 'g', long = "SnpGap")]
    pub snp_gap: Option<String>,

    #[arg(short = 'G', long = "IndelGap")]
    pub indel_gap: Option<u32>,

    #[arg(long = "mask")]
    pub mask: Option<String>,

    #[arg(short = 'M', long = "mask-file")]
    pub mask_file: Option<PathBuf>,

    #[arg(long = "mask-overlap")]
    pub mask_overlap: Option<u8>,

    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    #[arg(short = 'O', long = "output-type")]
    pub output_type: Option<String>,

    #[arg(long = "threads")]
    pub threads: Option<usize>,

    #[arg(long = "gpu", default_value_t = false)]
    pub gpu: bool,

    #[arg(long = "opencl", default_value_t = false)]
    pub opencl: bool,

    #[arg(long = "no-version", default_value_t = false)]
    pub no_version: bool,

    pub input: PathBuf,
}
