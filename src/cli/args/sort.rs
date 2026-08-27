use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct SortArgs {
    pub input: PathBuf,

    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    #[arg(short = 'O', long = "output-type")]
    pub output_type: Option<String>,

    #[arg(short = 'm', long = "max-mem", default_value = "768M")]
    pub max_mem: String,

    #[arg(short = 'T', long = "temp-dir")]
    pub temp_dir: Option<PathBuf>,

    #[arg(long = "no-version")]
    pub no_version: bool,

    #[arg(long = "threads", default_value_t = 0)]
    pub threads: usize,

    #[arg(short = 'W', long = "write-index", num_args = 0..=1, default_missing_value = "csi")]
    pub write_index: Option<String>,

    #[arg(short = 'v', long = "verbosity", default_value_t = 1)]
    pub verbosity: u8,

    #[arg(last = true, trailing_var_arg = true, allow_hyphen_values = true)]
    pub passthrough: Vec<String>,
}
