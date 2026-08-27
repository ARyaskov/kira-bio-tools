use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct SamViewArgs {
    pub input: PathBuf,

    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    #[arg(short = 'b', long = "bam")]
    pub bam: bool,

    #[arg(short = 'C', long = "cram")]
    pub cram: bool,

    #[arg(short = 'h', long = "with-header")]
    pub with_header: bool,

    #[arg(short = 'H', long = "header-only")]
    pub header_only: bool,

    #[arg(short = 'c', long = "count")]
    pub count: bool,

    #[arg(short = 'f', long = "require-flags", default_value_t = 0)]
    pub require_flags: u16,

    #[arg(short = 'F', long = "exclude-flags", default_value_t = 0)]
    pub exclude_flags: u16,

    #[arg(short = 'q', long = "min-MQ", default_value_t = 0)]
    pub min_mq: u32,

    #[arg(short = 'r', long = "region")]
    pub region: Option<String>,

    #[arg(short = 'L', long = "regions-file")]
    pub regions_file: Option<PathBuf>,

    #[arg(short = 's', long = "subsample", default_value_t = 1.0)]
    pub subsample: f64,

    #[arg(long = "threads", default_value_t = 0)]
    pub threads: usize,

    #[arg(short = 'v', long = "verbosity", default_value_t = 1)]
    pub verbosity: u8,
}
