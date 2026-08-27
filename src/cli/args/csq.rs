use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct CsqArgs {
    pub input: Option<PathBuf>,

    #[arg(short = 'f', long = "fasta-ref")]
    pub fasta_ref: Option<PathBuf>,

    #[arg(short = 'g', long = "gff-annot")]
    pub gff: Option<PathBuf>,

    #[arg(short = 'l', long = "local-csq")]
    pub local_csq: bool,

    #[arg(short = 'p', long = "phase", default_value = "a")]
    pub phase: String,

    #[arg(short = 'c', long = "custom-tag", default_value = "BCSQ")]
    pub custom_tag: String,

    #[arg(short = 's', long = "samples")]
    pub samples: Option<String>,

    #[arg(short = 'S', long = "samples-file")]
    pub samples_file: Option<PathBuf>,

    #[arg(short = 'B', long = "brief-predictions")]
    pub brief_predictions: bool,

    #[arg(short = 'H', long = "haplotypes")]
    pub haplotypes: bool,

    #[arg(short = 'q', long = "quiet")]
    pub quiet: bool,

    #[arg(long = "ncsq", default_value_t = 5)]
    pub ncsq: u32,

    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    #[arg(short = 'O', long = "output-type")]
    pub output_type: Option<String>,

    #[arg(short = 'r', long = "regions")]
    pub regions: Option<String>,

    #[arg(short = 'R', long = "regions-file")]
    pub regions_file: Option<PathBuf>,

    #[arg(long = "regions-overlap", default_value = "1")]
    pub regions_overlap: u8,

    #[arg(short = 't', long = "targets")]
    pub targets: Option<String>,

    #[arg(short = 'T', long = "targets-file")]
    pub targets_file: Option<PathBuf>,

    #[arg(short = 'i', long = "include")]
    pub include: Option<String>,

    #[arg(short = 'e', long = "exclude")]
    pub exclude: Option<String>,

    #[arg(long = "no-version")]
    pub no_version: bool,

    #[arg(long = "threads", default_value_t = 0)]
    pub threads: usize,

    #[arg(short = 'v', long = "verbosity", default_value_t = 1)]
    pub verbosity: u8,
}
