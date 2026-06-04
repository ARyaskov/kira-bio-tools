use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct GtcheckArgs {
    pub input: Option<PathBuf>,

    #[arg(short = 'g', long = "genotypes")]
    pub genotypes: Option<PathBuf>,

    #[arg(short = 'p', long = "pairs")]
    pub pairs: Option<String>,

    #[arg(short = 'P', long = "pairs-file")]
    pub pairs_file: Option<PathBuf>,

    #[arg(short = 's', long = "samples")]
    pub samples: Option<String>,

    #[arg(short = 'S', long = "samples-file")]
    pub samples_file: Option<PathBuf>,

    #[arg(short = 'e', long = "error-probability", default_value_t = 40)]
    pub error_probability: u32,

    #[arg(short = 'n', long = "n-matches", default_value_t = 0)]
    pub n_matches: usize,

    #[arg(short = 'H', long = "homs-only")]
    pub homs_only: bool,

    #[arg(short = 'u', long = "use", default_value = "GT")]
    pub use_tag: String,

    #[arg(long = "no-HWE-prob")]
    pub no_hwe_prob: bool,

    #[arg(long = "keep-refs")]
    pub keep_refs: bool,

    #[arg(long = "distinctive-sites")]
    pub distinctive_sites: Option<String>,

    #[arg(long = "dry-run")]
    pub dry_run: bool,

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

    #[arg(short = 'v', long = "verbosity", default_value_t = 1)]
    pub verbosity: u8,
}
