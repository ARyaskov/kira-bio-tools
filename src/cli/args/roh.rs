use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct RohArgs {
    pub input: PathBuf,

    #[arg(long = "AF-tag", default_value = "AF")]
    pub af_tag: String,

    #[arg(long = "AF-file")]
    pub af_file: Option<PathBuf>,

    #[arg(long = "AF-dflt", default_value_t = 0.0)]
    pub af_dflt: f64,

    #[arg(short = 'G', long = "GTs-only", default_value_t = 30.0)]
    pub gts_only: f64,

    #[arg(short = 'I', long = "ignore-homref")]
    pub ignore_homref: bool,

    #[arg(long = "include-noalt")]
    pub include_noalt: bool,

    #[arg(short = 'X', long = "skip-indels")]
    pub skip_indels: bool,

    #[arg(short = 'V', long = "viterbi-training", default_value_t = 0.0)]
    pub viterbi_training: f64,

    #[arg(short = 'M', long = "rec-rate", default_value_t = 1e-8)]
    pub rec_rate: f64,

    #[arg(short = 'a', long = "hw-to-az", default_value_t = 6.7e-8)]
    pub hw_to_az: f64,

    #[arg(short = 'H', long = "az-to-hw", default_value_t = 5e-9)]
    pub az_to_hw: f64,

    #[arg(short = 'b', long = "buffer-size", default_value_t = 0)]
    pub buffer_size: i64,

    #[arg(short = 'm', long = "genetic-map")]
    pub genetic_map: Option<PathBuf>,

    #[arg(short = 'E', long = "estimate-AF")]
    pub estimate_af: Option<String>,

    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    #[arg(short = 'O', long = "output-type", default_value = "rsz")]
    pub output_type: String,

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

    #[arg(long = "threads", default_value_t = 0)]
    pub threads: usize,

    #[arg(short = 'v', long = "verbosity", default_value_t = 1)]
    pub verbosity: u8,

    #[arg(last = true, allow_hyphen_values = true)]
    pub passthrough: Vec<String>,
}
