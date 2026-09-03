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

    /// `[qry|gt:]LIST`, may be given twice (once per file).
    #[arg(short = 's', long = "samples", action = clap::ArgAction::Append)]
    pub samples: Vec<String>,

    /// `[qry|gt:]FILE`, may be given twice (once per file).
    #[arg(short = 'S', long = "samples-file", action = clap::ArgAction::Append)]
    pub samples_file: Vec<String>,

    /// Phred-scaled genotyping error probability; 0 counts mismatches instead [40].
    #[arg(short = 'E', long = "error-probability")]
    pub error_probability: Option<u32>,

    /// Print only the top INT matches per sample; negative sorts by HWE probability.
    #[arg(long = "n-matches", default_value_t = 0, allow_hyphen_values = true)]
    pub n_matches: i64,

    #[arg(short = 'H', long = "homs-only")]
    pub homs_only: bool,

    /// `TAG1[,TAG2]`: tag used in the query file and in the -g file [PL,GT].
    #[arg(short = 'u', long = "use", default_value = "PL,GT")]
    pub use_tag: String,

    #[arg(long = "no-HWE-prob")]
    pub no_hwe_prob: bool,

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

    #[arg(long = "targets-overlap", default_value = "0")]
    pub targets_overlap: u8,

    /// `[qry|gt:]EXPR`
    #[arg(short = 'i', long = "include", action = clap::ArgAction::Append)]
    pub include: Vec<String>,

    /// `[qry|gt:]EXPR`; a bare integer is the old `--error-probability`.
    #[arg(short = 'e', long = "exclude", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    pub exclude: Vec<String>,

    #[arg(short = 'v', long = "verbosity", default_value_t = 1)]
    pub verbosity: u8,
}
