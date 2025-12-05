use clap::Parser;
use std::path::PathBuf;

#[derive(Parser)]
pub struct NormArgs {
    #[arg(help = "Input VCF file (.vcf or .vcf.gz)")]
    pub input: PathBuf,

    #[arg(short, long, help = "Output VCF file (optional)")]
    pub output: Option<PathBuf>,
}
