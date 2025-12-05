use clap::Args;
use std::path::PathBuf;

#[derive(Args, Debug)]
pub struct FilterArgs {
    #[arg(long = "expr")]
    pub expr: String,

    #[arg(long = "soft-filter")]
    pub soft_filter: Option<String>,

    #[arg(long = "preload-info")]
    pub preload_info: Option<String>,

    #[arg(long = "pass-only", default_value_t = false)]
    pub pass_only: bool,

    pub input: PathBuf,

    #[arg(short = 'o', long = "output")]
    pub output: PathBuf,
}
