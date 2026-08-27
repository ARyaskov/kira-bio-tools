use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct MergeArgs {
    #[arg(required = true)]
    pub inputs: Vec<PathBuf>,

    #[arg(short = 'l', long = "file-list")]
    pub file_list: Option<PathBuf>,

    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    #[arg(short = 'O', long = "output-type")]
    pub output_type: Option<String>,

    #[arg(short = 'm', long = "merge", default_value = "both")]
    pub merge: String,

    #[arg(short = 'f', long = "apply-filters")]
    pub apply_filters: Option<String>,

    #[arg(short = 'F', long = "filter-logic", default_value = "+")]
    pub filter_logic: String,

    #[arg(long = "force-samples")]
    pub force_samples: bool,

    #[arg(long = "force-single")]
    pub force_single: bool,

    #[arg(long = "force-no-index")]
    pub force_no_index: bool,

    #[arg(long = "no-index")]
    pub no_index: bool,

    #[arg(short = 'i', long = "info-rules")]
    pub info_rules: Option<String>,

    #[arg(long = "missing-rules")]
    pub missing_rules: Option<String>,

    #[arg(long = "missing-to-ref", short = '0')]
    pub missing_to_ref: bool,

    #[arg(short = 'g', long = "gvcf")]
    pub gvcf: Option<PathBuf>,

    #[arg(short = 'L', long = "local-alleles", default_value_t = 0)]
    pub local_alleles: u32,

    #[arg(long = "print-header")]
    pub print_header: bool,

    #[arg(long = "use-header")]
    pub use_header: Option<PathBuf>,

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
