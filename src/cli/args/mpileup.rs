use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(disable_help_flag = true)]
pub struct MpileupArgs {
    #[arg(long = "help", action = clap::ArgAction::Help, help = "Print help")]
    pub help: Option<bool>,

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

    /// Comma list of tags; a leading `-` removes a tag (`-a -AD`); may repeat.
    #[arg(short = 'a', long = "annotate", allow_hyphen_values = true, action = clap::ArgAction::Append)]
    pub annotate: Vec<String>,

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
    #[arg(short = 'v', long = "variants-only")]
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

    #[arg(short = 'B', long = "no-BAQ")]
    pub no_baq: bool,

    #[arg(short = 'E', long = "redo-BAQ")]
    pub redo_baq: bool,

    /// Disable read-pair overlap detection (mates sharing a base then count twice).
    #[arg(short = 'x', long = "ignore-overlaps")]
    pub ignore_overlaps: bool,

    /// Down-weight qualities of reads with a high mismatch rate (paralog
    /// suspects): `off`, `auto`, or `FULL,SLOPE`. A kira extension, off by default.
    #[arg(long = "nm-weight", default_value = "off")]
    pub nm_weight: String,

    /// Require more support for short indels inside homopolymer/STR tracts and
    /// lower their quality by run length. A kira extension, off by default.
    #[arg(long = "hp-indel")]
    pub hp_indel: bool,

    /// Confirm indel candidates by pair-HMM realignment of the covering reads:
    /// `off`, `ins` (insertions) or `all`. A kira extension, off by default.
    #[arg(long = "indel-realign", default_value = "off")]
    pub indel_realign: String,

    /// Recover indels the aligner hid as mismatches by local assembly in active
    /// regions; records may need sorting afterwards. A kira extension.
    #[arg(long = "assemble")]
    pub assemble: bool,

    /// Recalibrate base qualities from the data (reported Q × reference
    /// trinucleotide × strand × cycle, learned at non-variant sites) before
    /// BAQ. A kira extension; needs `-f` and is not available with `--stream`.
    #[arg(long = "recal")]
    pub recal: bool,

    #[arg(short = 'X', long = "config")]
    pub config: Option<String>,

    #[arg(short = 'F', long = "gap-frac", default_value_t = 0.002)]
    pub gap_frac: f64,

    #[arg(short = 'h', long = "tandem-qual", default_value_t = 500)]
    pub tandem_qual: u32,

    #[arg(short = 'L', long = "max-idepth-orphan")]
    pub max_idepth_orphan: Option<u32>,

    #[arg(long = "open-prob")]
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

    #[arg(long = "rf", aliases = ["incl-flags", "skip-any-unset"])]
    pub rf: Option<String>,

    #[arg(long = "ff", aliases = ["excl-flags", "skip-any-set"], default_value = "UNMAP,SECONDARY,QCFAIL,DUP")]
    pub ff: Option<String>,

    #[arg(long = "skip-all-set")]
    pub skip_all_set: Option<String>,

    #[arg(long = "skip-all-unset")]
    pub skip_all_unset: Option<String>,

    // Accepted for bcftools compatibility; they do not change the pileup.
    #[arg(short = 'A', long = "count-orphans")]
    pub count_orphans: bool,

    #[arg(long = "ambig-reads")]
    pub ambig_reads: Option<String>,

    #[arg(short = 'G', long = "read-groups")]
    pub read_groups: Option<PathBuf>,

    #[arg(short = 'p', long = "per-sample-mF")]
    pub per_sample_mf: bool,

    #[arg(short = 'P', long = "platforms")]
    pub platforms: Option<String>,

    #[arg(short = 'C', long = "adjust-MQ", default_value_t = 0)]
    pub adjust_mq: u32,

    #[arg(short = 'e', long = "ext-prob", default_value_t = 20)]
    pub ext_prob: u32,

    #[arg(long = "seed")]
    pub seed: Option<u64>,

    #[arg(long = "indels-2.0")]
    pub indels_20: bool,

    #[arg(long = "indels-cns")]
    pub indels_cns: bool,

    #[arg(long = "indel-bias")]
    pub indel_bias: Option<f64>,

    #[arg(long = "max-BQ")]
    pub max_bq: Option<u32>,

    #[arg(long = "delta-BQ")]
    pub delta_bq: Option<u32>,

    #[arg(short = 'M', long = "min-ireads-frac")]
    pub min_ireads_frac: Option<f64>,

    #[arg(long = "indel-size", default_value_t = 110)]
    pub indel_size: u32,

    #[arg(long = "threads", default_value_t = 0)]
    pub threads: usize,

    #[arg(long = "stream", help = "Streaming mode: spawn reader thread per BAM (no full Vec<Record> in RAM)")]
    pub stream: bool,

    #[arg(long = "verbosity", default_value_t = 1)]
    pub verbosity: u8,

    #[arg(last = true, allow_hyphen_values = true)]
    pub passthrough: Vec<String>,
}
