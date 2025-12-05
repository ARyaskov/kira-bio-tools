use clap::Parser;
use std::path::PathBuf;

#[derive(Parser)]
pub struct IndexArgs {
    #[arg(help = "Input VCF file (.vcf, .vcf.gz with BGZF)")]
    pub input: PathBuf,

    #[arg(short, long, help = "Output index file")]
    pub output: Option<PathBuf>,

    #[arg(
        short,
        long,
        default_value = "vcf",
        help = "Format preset (vcf, bed, gff, sam)"
    )]
    pub preset: String,

    #[arg(short, long, help = "Force overwrite existing index")]
    pub force: bool,

    #[arg(
        short = 's',
        long,
        default_value = "14",
        help = "Minimum bin shift (CSI)"
    )]
    pub min_shift: u8,

    #[arg(short = 'd', long, default_value = "5", help = "Bin depth (CSI)")]
    pub depth: u8,

    #[arg(short = 'C', long, help = "Generate CSI index (for BGZF files)")]
    pub csi: bool,

    #[arg(long, help = "Skip KBI index generation")]
    pub no_kbi: bool,
}
