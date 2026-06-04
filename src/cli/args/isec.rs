use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct IsecArgs {
    #[arg(required = true)]
    pub inputs: Vec<PathBuf>,

    #[arg(short = 'p', long = "prefix")]
    pub prefix: Option<PathBuf>,

    #[arg(short = 'n', long = "nfiles")]
    pub nfiles: Option<String>,

    #[arg(short = 'w', long = "write")]
    pub write: Option<String>,

    #[arg(short = 'c', long = "collapse")]
    pub collapse: Option<String>,

    #[arg(short = 'C', long = "complement")]
    pub complement: bool,

    #[arg(short = 'f', long = "apply-filters")]
    pub apply_filters: Option<String>,

    #[arg(short = 'l', long = "file-list")]
    pub file_list: Option<PathBuf>,

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

    #[arg(short = 'W', long = "write-index", num_args = 0..=1, default_missing_value = "csi")]
    pub write_index: Option<String>,
}
