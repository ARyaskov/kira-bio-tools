use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct NormArgs {
    pub input: PathBuf,

    #[arg(short = 'f', long = "fasta-ref")]
    pub fasta_ref: Option<PathBuf>,

    #[arg(long = "fasta-ref-fai")]
    pub fasta_ref_fai: Option<PathBuf>,

    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    #[arg(short = 'O', long = "output-type")]
    pub output_type: Option<String>,

    #[arg(short = 'm', long = "multiallelics")]
    pub multiallelics: Option<String>,

    #[arg(long = "atomize")]
    pub atomize: bool,

    #[arg(long = "atom-overlaps", default_value = "*")]
    pub atom_overlaps: String,

    #[arg(short = 'c', long = "check-ref", default_value = "e")]
    pub check_ref: String,

    #[arg(short = 'N', long = "do-not-normalize")]
    pub do_not_normalize: bool,

    #[arg(short = 'd', long = "rm-dup")]
    pub rm_dup: Option<String>,

    #[arg(long = "site-win", default_value_t = 1000)]
    pub site_win: u32,

    #[arg(long = "strict-filter")]
    pub strict_filter: bool,

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

    #[arg(long = "keep-sum")]
    pub keep_sum: Option<String>,

    #[arg(long = "old-rec-tag")]
    pub old_rec_tag: Option<String>,

    #[arg(long = "multi-overlaps", default_value = "0")]
    pub multi_overlaps: String,

    #[arg(short = 'S', long = "sort", default_value = "lex")]
    pub sort: String,

    #[arg(long = "force")]
    pub force: bool,

    #[arg(long = "no-version")]
    pub no_version: bool,

    #[arg(long = "threads", default_value_t = 0)]
    pub threads: usize,

    #[arg(short = 'W', long = "write-index", num_args = 0..=1, default_missing_value = "csi")]
    pub write_index: Option<String>,
}
