use clap::Parser;
use std::path::PathBuf;

#[derive(Parser)]
pub struct RohArgs {
    #[arg(help = "Input VCF/BCF file")]
    pub input: PathBuf,

    #[arg(last = true, trailing_var_arg = true)]
    pub bcftools_args: Vec<String>,
}
