use clap::Parser;
use std::path::PathBuf;

#[derive(Parser)]
pub struct MpileupArgs {
    #[arg(help = "Input alignment files (BAM/CRAM/SAM)", required = true)]
    pub inputs: Vec<PathBuf>,

    #[arg(short, long, help = "Output VCF file (optional)")]
    pub output: Option<PathBuf>,

    #[arg(last = true, trailing_var_arg = true)]
    pub bcftools_args: Vec<String>,
}
