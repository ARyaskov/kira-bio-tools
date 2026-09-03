use clap::Parser;
use std::path::PathBuf;

#[derive(Parser)]
pub struct PolysomyArgs {
    #[arg(help = "Input VCF/BCF file")]
    pub input: PathBuf,

    #[arg(last = true)]
    pub bcftools_args: Vec<String>,
}
