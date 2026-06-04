use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct HeadArgs {
    pub input: PathBuf,

    #[arg(short = 'h', long = "headers", default_value_t = -1)]
    pub headers: i64,

    #[arg(short = 'n', long = "records", default_value_t = -1)]
    pub records: i64,

    #[arg(short = 's', long = "samples", default_value_t = -1)]
    pub samples: i64,

    #[arg(short = 'v', long = "verbosity", default_value_t = 1)]
    pub verbosity: u8,

    #[arg(last = true, trailing_var_arg = true, allow_hyphen_values = true)]
    pub passthrough: Vec<String>,
}
