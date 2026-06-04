use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct ConcatArgs {
    pub inputs: Vec<PathBuf>,

    #[arg(short = 'f', long = "file-list")]
    pub file_list: Option<PathBuf>,

    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    #[arg(short = 'O', long = "output-type")]
    pub output_type: Option<String>,

    #[arg(short = 'a', long = "allow-overlaps")]
    pub allow_overlaps: bool,

    #[arg(short = 'l', long = "ligate")]
    pub ligate: bool,

    #[arg(long = "ligate-warn")]
    pub ligate_warn: bool,

    #[arg(long = "ligate-force")]
    pub ligate_force: bool,

    #[arg(short = 'D', long = "remove-duplicates")]
    pub remove_duplicates: bool,

    #[arg(short = 'd', long = "rm-dups")]
    pub rm_dups: Option<String>,

    #[arg(short = 'n', long = "naive")]
    pub naive: bool,

    #[arg(long = "naive-force")]
    pub naive_force: bool,

    #[arg(short = 'G', long = "drop-genotypes")]
    pub drop_genotypes: bool,

    #[arg(short = 'c', long = "compact-PS")]
    pub compact_ps: bool,

    #[arg(short = 'q', long = "min-PQ", default_value_t = 30)]
    pub min_pq: u32,

    #[arg(short = 'r', long = "regions")]
    pub regions: Option<String>,

    #[arg(short = 'R', long = "regions-file")]
    pub regions_file: Option<PathBuf>,

    #[arg(long = "regions-overlap", default_value = "1")]
    pub regions_overlap: u8,

    #[arg(long = "no-version")]
    pub no_version: bool,

    #[arg(long = "threads", default_value_t = 0)]
    pub threads: usize,

    #[arg(short = 'W', long = "write-index", num_args = 0..=1, default_missing_value = "csi")]
    pub write_index: Option<String>,
}
