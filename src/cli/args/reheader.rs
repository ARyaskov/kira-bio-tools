use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct ReheaderArgs {
    pub input: PathBuf,

    #[arg(short = 'h', long = "header")]
    pub header: Option<PathBuf>,

    #[arg(short = 's', long = "samples")]
    pub samples_file: Option<PathBuf>,

    #[arg(long = "samples-list")]
    pub samples_list: Option<String>,

    #[arg(short = 'f', long = "fai")]
    pub fai: Option<PathBuf>,

    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    #[arg(long = "temp-prefix")]
    pub temp_prefix: Option<String>,

    #[arg(long = "threads", default_value_t = 0)]
    pub threads: usize,

    #[arg(short = 'v', long = "verbosity", default_value_t = 1)]
    pub verbosity: u8,

    #[arg(last = true, trailing_var_arg = true, allow_hyphen_values = true)]
    pub passthrough: Vec<String>,
}
