// Command-line arguments for VCF filtering

use clap::Args;
use std::path::PathBuf;

#[derive(Args, Debug)]
pub struct FilterArgs {
    /// Filter expression (bcftools-compatible)
    #[arg(long = "expr")]
    pub expr: String,

    /// Name to use when marking filtered-out records (like bcftools --soft-filter)
    #[arg(long = "soft-filter")]
    pub soft_filter: Option<String>,

    /// Preload specific INFO fields for faster evaluation
    #[arg(long = "preload-info")]
    pub preload_info: Option<String>,

    /// Output only records passing the filter (like bcftools -i)
    #[arg(long = "pass-only", default_value_t = false)]
    pub pass_only: bool,

    /// Input VCF file
    pub input: PathBuf,

    /// Output VCF file
    #[arg(short = 'o', long = "output")]
    pub output: PathBuf,
}
