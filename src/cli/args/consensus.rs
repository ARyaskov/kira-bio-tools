use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct ConsensusArgs {
    pub input: Option<PathBuf>,

    #[arg(short = 'f', long = "fasta-ref")]
    pub fasta_ref: Option<PathBuf>,

    #[arg(short = 'c', long = "chain")]
    pub chain: Option<PathBuf>,

    #[arg(short = 'H', long = "haplotype", default_value = "1")]
    pub haplotype: String,

    #[arg(short = 'I', long = "iupac-codes")]
    pub iupac_codes: bool,

    #[arg(short = 'a', long = "absent")]
    pub absent: Option<String>,

    #[arg(long = "mark-del")]
    pub mark_del: Option<String>,

    #[arg(long = "mark-ins")]
    pub mark_ins: Option<String>,

    #[arg(long = "mark-snv")]
    pub mark_snv: Option<String>,

    #[arg(short = 'm', long = "mask", action = clap::ArgAction::Append)]
    pub mask: Vec<PathBuf>,

    #[arg(long = "mask-with", action = clap::ArgAction::Append)]
    pub mask_with: Vec<String>,

    #[arg(short = 'M', long = "missing")]
    pub missing: Option<String>,

    #[arg(short = 's', long = "samples")]
    pub samples: Option<String>,

    #[arg(short = 'S', long = "samples-file")]
    pub samples_file: Option<PathBuf>,

    #[arg(short = 'p', long = "prefix")]
    pub prefix: Option<String>,

    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    #[arg(short = 'i', long = "include")]
    pub include: Option<String>,

    #[arg(short = 'e', long = "exclude")]
    pub exclude: Option<String>,

    #[arg(short = 'r', long = "regions")]
    pub regions: Option<String>,

    #[arg(short = 'R', long = "regions-file")]
    pub regions_file: Option<PathBuf>,

    #[arg(long = "regions-overlap", default_value = "1")]
    pub regions_overlap: u8,

    #[arg(short = 'v', long = "verbosity", default_value_t = 1)]
    pub verbosity: u8,
}
