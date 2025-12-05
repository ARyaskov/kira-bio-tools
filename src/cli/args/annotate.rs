use clap::Parser;
use std::path::PathBuf;

#[derive(Parser)]
pub struct AnnotateIndexArgs {
    #[arg(help = "Annotation VCF to index")]
    pub input: PathBuf,

    #[arg(short, long, help = "Output .ani file")]
    pub output: Option<PathBuf>,
}

#[derive(Parser)]
pub struct AnnotateArgs {
    #[arg(help = "Input VCF file (.vcf)")]
    pub input: PathBuf,

    #[arg(
        short = 'a',
        long = "annotations",
        help = "Annotation DB (.vcf or .ani)"
    )]
    pub annotations: PathBuf,

    #[arg(short, long, help = "Output file (optional)")]
    pub output: Option<PathBuf>,

    #[arg(long, help = "Use GPU (CUDA) for annotation")]
    pub gpu: bool,

    #[arg(long, help = "Use OpenCL backend (AMD/Intel/Apple/GPU/CPU)")]
    pub opencl: bool,
}

#[derive(Parser)]
pub struct DbBuildArgs {
    #[arg(help = "Input annotation database VCF")]
    pub input: PathBuf,

    #[arg(short, long, help = "Output ANI file")]
    pub output: Option<PathBuf>,
}
