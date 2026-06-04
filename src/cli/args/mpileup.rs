use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct MpileupArgs {
    #[arg(required = true)]
    pub inputs: Vec<PathBuf>,

    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    #[arg(short = 'O', long = "output-type")]
    pub output_type: Option<String>,

    #[arg(short = 'f', long = "fasta-ref")]
    pub fasta_ref: Option<PathBuf>,

    #[arg(short = 'g', long = "gvcf")]
    pub gvcf: Option<String>,

    #[arg(long = "no-reference")]
    pub no_reference: bool,

    #[arg(short = 'a', long = "annotate")]
    pub annotate: Option<String>,

    #[arg(short = 'b', long = "bam-list")]
    pub bam_list: Option<PathBuf>,

    #[arg(short = 'd', long = "max-depth", default_value_t = 250)]
    pub max_depth: u32,

    #[arg(short = 'D', long = "max-idepth")]
    pub max_idepth: Option<u32>,

    #[arg(short = 'q', long = "min-MQ", default_value_t = 0)]
    pub min_mq: u32,

    #[arg(short = 'Q', long = "min-BQ", default_value_t = 13)]
    pub min_bq: u32,

    #[arg(short = 'm', long = "min-ireads", default_value_t = 1)]
    pub min_ireads: u32,

    #[arg(short = 'I', long = "skip-indels")]
    pub skip_indels: bool,

    /// Fused with `call -m -v`: emit only sites whose posterior best GT is non-ref.
    #[arg(long = "variants-only")]
    pub variants_only: bool,

    /// Per-allele variant prior for Bayesian GT call (bcftools default 1.1e-3).
    #[arg(long = "prior", default_value_t = 1.1e-3)]
    pub prior: f64,

    /// Minimum ALT reads to consider an SNV/indel candidate (variants-only mode).
    #[arg(long = "min-alt-reads", default_value_t = 2)]
    pub min_alt_reads: u32,

    /// Minimum allele fraction for an ALT candidate (variants-only mode).
    #[arg(long = "min-af", default_value_t = 0.20)]
    pub min_af: f64,

    /// Minimum site QUAL to emit in variants-only mode (matches bcftools-call default).
    #[arg(long = "min-qual", default_value_t = 10)]
    pub min_qual: u8,

    #[arg(short = 'x', long = "no-BAQ")]
    pub no_baq: bool,

    #[arg(short = 'B', long = "no-BAQ-old")]
    pub no_baq_old: bool,

    #[arg(short = 'E', long = "redo-BAQ")]
    pub redo_baq: bool,

    #[arg(short = 'X', long = "config")]
    pub config: Option<String>,

    #[arg(short = 'F', long = "gap-frac", default_value_t = 0.002)]
    pub gap_frac: f64,

    #[arg(short = 'h', long = "tandem-qual", default_value_t = 500)]
    pub tandem_qual: u32,

    #[arg(short = 'L', long = "max-idepth-orphan")]
    pub max_idepth_orphan: Option<u32>,

    #[arg(short = 'o', long = "open-prob")]
    pub open_prob_dup: Option<u32>,

    #[arg(short = 's', long = "samples")]
    pub samples: Option<String>,

    #[arg(short = 'S', long = "samples-file")]
    pub samples_file: Option<PathBuf>,

    #[arg(short = 'r', long = "regions")]
    pub regions: Option<String>,

    #[arg(short = 'R', long = "regions-file")]
    pub regions_file: Option<PathBuf>,

    #[arg(short = 't', long = "targets")]
    pub targets: Option<String>,

    #[arg(short = 'T', long = "targets-file")]
    pub targets_file: Option<PathBuf>,

    #[arg(long = "ef")]
    pub ef: Option<u32>,

    #[arg(long = "df")]
    pub df: Option<u32>,

    #[arg(long = "if")]
    pub if_: Option<u32>,

    #[arg(long = "nf")]
    pub nf: Option<u32>,

    #[arg(long = "rf", aliases = ["incl-flags"])]
    pub rf: Option<String>,

    #[arg(long = "ff", aliases = ["excl-flags"], default_value = "UNMAP,SECONDARY,QCFAIL,DUP")]
    pub ff: Option<String>,

    #[arg(short = 'M', long = "min-ireads-frac")]
    pub min_ireads_frac: Option<f64>,

    #[arg(long = "indel-size", default_value_t = 110)]
    pub indel_size: u32,

    #[arg(long = "threads", default_value_t = 0)]
    pub threads: usize,

    #[arg(long = "stream", help = "Streaming mode: spawn reader thread per BAM (no full Vec<Record> in RAM)")]
    pub stream: bool,

    #[arg(short = 'v', long = "verbosity", default_value_t = 1)]
    pub verbosity: u8,

    #[arg(last = true, trailing_var_arg = true, allow_hyphen_values = true)]
    pub passthrough: Vec<String>,
}
