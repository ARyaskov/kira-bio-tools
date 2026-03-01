use clap::Parser;
use std::path::PathBuf;

#[derive(Parser)]
pub struct StatsArgs {
    #[arg(help = "Input VCF/BCF file(s)", required = true)]
    pub inputs: Vec<PathBuf>,

    #[arg(last = true, trailing_var_arg = true)]
    pub bcftools_args: Vec<String>,
}
