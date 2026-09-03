use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(disable_help_flag = true)]
pub struct ViewArgs {
    #[arg(long = "help", action = clap::ArgAction::Help, help = "Print help")]
    pub help: Option<bool>,

    pub input: PathBuf,

    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    #[arg(short = 'O', long = "output-type")]
    pub output_type: Option<String>,

    #[arg(short = 'l', long = "compression-level", default_value_t = -1)]
    pub compression_level: i32,

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

    #[arg(long = "targets-overlap", default_value = "0")]
    pub targets_overlap: u8,

    #[arg(short = 's', long = "samples")]
    pub samples: Option<String>,

    #[arg(short = 'S', long = "samples-file")]
    pub samples_file: Option<PathBuf>,

    #[arg(long = "force-samples")]
    pub force_samples: bool,

    #[arg(short = 'i', long = "include")]
    pub include: Option<String>,

    #[arg(short = 'e', long = "exclude")]
    pub exclude: Option<String>,

    #[arg(short = 'a', long = "trim-alt-alleles")]
    pub trim_alt_alleles: bool,

    #[arg(long = "trim-unseen-allele")]
    pub trim_unseen_allele: bool,

    #[arg(short = 'I', long = "no-update")]
    pub no_update: bool,

    #[arg(short = 'G', long = "drop-genotypes")]
    pub drop_genotypes: bool,

    #[arg(short = 'h', long = "header-only")]
    pub header_only: bool,

    #[arg(short = 'H', long = "no-header")]
    pub no_header: bool,

    #[arg(long = "with-header", default_value_t = true)]
    pub with_header: bool,

    #[arg(short = 'f', long = "apply-filters")]
    pub apply_filters: Option<String>,

    #[arg(short = 'v', long = "types")]
    pub types: Option<String>,

    #[arg(short = 'V', long = "exclude-types")]
    pub exclude_types: Option<String>,

    #[arg(long = "known")]
    pub known: bool,

    #[arg(long = "novel")]
    pub novel: bool,

    #[arg(short = 'm', long = "min-alleles")]
    pub min_alleles: Option<u32>,

    #[arg(short = 'M', long = "max-alleles")]
    pub max_alleles: Option<u32>,

    #[arg(short = 'g', long = "genotype")]
    pub genotype: Option<String>,

    #[arg(short = 'p', long = "phased")]
    pub phased: bool,

    #[arg(short = 'P', long = "exclude-phased")]
    pub exclude_phased: bool,

    #[arg(short = 'u', long = "uncalled")]
    pub uncalled: bool,

    #[arg(short = 'U', long = "exclude-uncalled")]
    pub exclude_uncalled: bool,

    #[arg(short = 'x', long = "private")]
    pub private: bool,

    #[arg(short = 'X', long = "exclude-private")]
    pub exclude_private: bool,

    #[arg(short = 'c', long = "min-ac")]
    pub min_ac: Option<String>,

    #[arg(short = 'C', long = "max-ac")]
    pub max_ac: Option<String>,

    #[arg(short = 'q', long = "min-af")]
    pub min_af: Option<String>,

    #[arg(short = 'Q', long = "max-af")]
    pub max_af: Option<String>,

    #[arg(long = "no-version")]
    pub no_version: bool,

    #[arg(long = "threads", default_value_t = 0)]
    pub threads: usize,

    #[arg(short = 'W', long = "write-index", num_args = 0..=1, default_missing_value = "csi")]
    pub write_index: Option<String>,
}
