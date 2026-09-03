use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct StatsArgs {
    #[arg(required = true)]
    pub inputs: Vec<PathBuf>,

    #[arg(short = 's', long = "samples")]
    pub samples: Option<String>,

    #[arg(short = 'S', long = "samples-file")]
    pub samples_file: Option<PathBuf>,

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

    #[arg(short = 'f', long = "apply-filters")]
    pub apply_filters: Option<String>,

    #[arg(short = 'F', long = "fasta-ref")]
    pub fasta_ref: Option<PathBuf>,

    #[arg(short = 'd', long = "depth", default_value = "0,500,1")]
    pub depth: String,

    #[arg(long = "af-bins", default_value = "0.01,0.05,0.1,1")]
    pub af_bins: String,

    #[arg(long = "af-tag")]
    pub af_tag: Option<String>,

    #[arg(short = '1', long = "1st-allele-only")]
    pub first_allele_only: bool,

    #[arg(short = 'I', long = "split-by-ID")]
    pub split_by_id: bool,

    #[arg(short = 'c', long = "collapse")]
    pub collapse: Option<String>,

    #[arg(short = 'E', long = "exons")]
    pub exons: Option<PathBuf>,

    #[arg(short = 'u', long = "user-tstv")]
    pub user_tstv: Option<String>,

    #[arg(long = "threads", default_value_t = 0)]
    pub threads: usize,

    #[arg(short = 'v', long = "verbose")]
    pub verbose: bool,
}
