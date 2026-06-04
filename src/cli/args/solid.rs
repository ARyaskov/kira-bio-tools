use std::path::PathBuf;

use clap::Parser;

/// Fused single-binary pipeline: align → coordinate-sort → markdup → mpileup,
/// entirely in memory. The ONLY file written to disk is the output VCF.
///
/// Aligner-stage options use the `--aligner-*` prefix (mirroring kira-ls-aligner
/// `mem`), the in-memory BAM stage uses `--bam-*`, and mpileup uses `--mpileup-*`.
#[derive(Parser, Debug)]
pub struct SolidArgs {
    /// Reference FASTA (shared by the aligner and mpileup/BAQ).
    #[arg(long = "aligner-ref", value_name = "REF")]
    pub aligner_ref: PathBuf,

    /// Prebuilt kira-ls-aligner index (`.kiraidx`). Built in memory if omitted.
    #[arg(long = "aligner-index", value_name = "IDX")]
    pub aligner_index: Option<PathBuf>,

    /// R1 FASTQ(.gz).
    #[arg(long = "aligner-r1", value_name = "R1")]
    pub aligner_r1: PathBuf,

    /// R2 FASTQ(.gz) for paired-end. Omit for single-end.
    #[arg(long = "aligner-r2", value_name = "R2")]
    pub aligner_r2: Option<PathBuf>,

    /// Read group fields, e.g. `ID:hg002\tSM:HG002\tLB:lib1\tPL:ILLUMINA`.
    #[arg(long = "aligner-rg", value_name = "RG")]
    pub aligner_rg: Option<String>,

    /// Insert-size spec `MIN,MAX,MEAN,SD`.
    #[arg(long = "aligner-insert-size", default_value = "0,1000,200,50")]
    pub aligner_insert_size: String,

    /// Aligner batch size in bases.
    #[arg(long = "aligner-batch", default_value_t = 1_000_000)]
    pub aligner_batch: usize,

    /// Mark PCR/optical duplicates during the in-memory coordinate sort.
    #[arg(long = "bam-markdup", default_value_t = false)]
    pub bam_markdup: bool,

    /// mpileup FORMAT/INFO annotations.
    #[arg(long = "mpileup-annotate", default_value = "AD,DP,SP")]
    pub mpileup_annotate: String,

    /// mpileup max depth (0 = unlimited).
    #[arg(long = "mpileup-max-depth", default_value_t = 0)]
    pub mpileup_max_depth: u32,

    /// Emit variant sites only (skip ref-only positions).
    #[arg(long = "mpileup-variants-only", default_value_t = true)]
    pub mpileup_variants_only: bool,

    /// Minimum read mapping quality for the mpileup stage. Applied IN-STREAM during the in-memory
    /// pileup (no intermediate BAM is written). 20 matches the validated `samtools view -q 20` step.
    #[arg(long = "mpileup-min-mq", default_value_t = 0)]
    pub mpileup_min_mq: u32,

    /// Minimum base quality for the mpileup stage. 6 is the validated chr20 optimum
    /// (with NM-aware weighting) for SNP recall; 13 is the samtools default.
    #[arg(long = "mpileup-min-bq", default_value_t = 13)]
    pub mpileup_min_bq: u32,

    /// Number of threads.
    #[arg(short = 't', long = "threads", default_value_t = 8)]
    pub threads: usize,

    /// Output VCF path — the only artifact written to disk.
    #[arg(short = 'o', long = "output", value_name = "VCF")]
    pub output: PathBuf,
}
