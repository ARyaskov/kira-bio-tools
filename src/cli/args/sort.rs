use clap::Parser;
use std::path::PathBuf;

#[derive(Parser)]
pub struct SortArgs {
    #[arg(help = "Input VCF/BCF file")]
    pub input: PathBuf,

    #[arg(short, long, help = "Output VCF file (optional)")]
    pub output: Option<PathBuf>,

    #[arg(last = true, trailing_var_arg = true)]
    pub bcftools_args: Vec<String>,
}
