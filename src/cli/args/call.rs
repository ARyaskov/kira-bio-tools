use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct CallArgs {
    pub input: PathBuf,

    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    #[arg(short = 'O', long = "output-type")]
    pub output_type: Option<String>,

    #[arg(short = 'c', long = "consensus-caller")]
    pub consensus_caller: bool,

    #[arg(short = 'm', long = "multiallelic-caller")]
    pub multiallelic_caller: bool,

    #[arg(short = 'v', long = "variants-only")]
    pub variants_only: bool,

    #[arg(short = 'A', long = "keep-alts")]
    pub keep_alts: bool,

    #[arg(short = 'M', long = "keep-masked-ref")]
    pub keep_masked_ref: bool,

    #[arg(long = "keep-unseen-allele")]
    pub keep_unseen_allele: bool,

    #[arg(long = "insert-missed")]
    pub insert_missed: Option<u32>,

    #[arg(short = 'p', long = "pval-threshold", default_value_t = 0.5)]
    pub pval_threshold: f64,

    #[arg(short = 'P', long = "prior", default_value_t = 1.1e-3)]
    pub prior: f64,

    #[arg(short = 'F', long = "prior-freqs")]
    pub prior_freqs: Option<String>,

    /// Site-specific allele-frequency prior from an INFO tag (e.g. `INFO/AF`
    /// after annotating from gnomAD); shorthand for `-F` with a single tag.
    #[arg(long = "prior-af")]
    pub prior_af: Option<String>,

    #[arg(short = 'X', long = "novel-rate")]
    pub novel_rate: Option<String>,

    #[arg(short = 'C', long = "constrain")]
    pub constrain: Option<String>,

    /// gVCF output: `INT[,INT...]` minimum per-sample depths of the blocks.
    #[arg(short = 'g', long = "gvcf")]
    pub gvcf: Option<String>,

    #[arg(short = 'a', long = "annotate")]
    pub annotate: Option<String>,

    #[arg(long = "ploidy")]
    pub ploidy: Option<String>,

    #[arg(long = "ploidy-file")]
    pub ploidy_file: Option<PathBuf>,

    #[arg(short = 'G', long = "group-samples")]
    pub group_samples: Option<PathBuf>,

    #[arg(long = "group-samples-tag")]
    pub group_samples_tag: Option<String>,

    #[arg(short = 's', long = "samples")]
    pub samples: Option<String>,

    #[arg(short = 'S', long = "samples-file")]
    pub samples_file: Option<PathBuf>,

    #[arg(short = 'V', long = "skip-variants")]
    pub skip_variants: Option<String>,

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

    #[arg(long = "no-version")]
    pub no_version: bool,

    #[arg(long = "threads", default_value_t = 0)]
    pub threads: usize,

    #[arg(long = "verbosity", default_value_t = 1)]
    pub verbosity: u8,

    #[arg(short = 'W', long = "write-index", num_args = 0..=1, default_missing_value = "csi")]
    pub write_index: Option<String>,
}
