use std::path::PathBuf;

use clap::Parser;

/// Fused single-binary pipeline: align → coordinate-sort → markdup → mpileup →
/// call, entirely in memory. The ONLY file written to disk is the output VCF.
///
/// Aligner-stage options use the `--aligner-*` prefix (mirroring kira-ls-aligner
/// `mem`), the in-memory BAM stage uses `--bam-*`, mpileup uses `--mpileup-*`,
/// and the optional caller uses `--call-*`.
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

    /// NM-aware quality weighting for the mpileup stage: `off`, `auto`, or
    /// `FULL,SLOPE`. Suppresses paralog misplacements the aligner still lets
    /// through. `KIRA_NM_WEIGHT` sets it when the flag is absent.
    #[arg(long = "mpileup-nm-weight", env = "KIRA_NM_WEIGHT", default_value = "off")]
    pub mpileup_nm_weight: String,

    /// Pair-HMM realignment of indel candidates in the mpileup stage: `off`,
    /// `ins` or `all`. `KIRA_INDEL_REALIGN` sets it when the flag is absent.
    #[arg(long = "mpileup-indel-realign", env = "KIRA_INDEL_REALIGN", default_value = "off")]
    pub mpileup_indel_realign: String,

    /// Run the variant caller on the pileup, in-stream, so the pipeline emits a
    /// called VCF instead of the pileup's per-site maximum-likelihood genotypes.
    /// Off by default: it changes what the output VCF means.
    #[arg(long = "call", default_value_t = false)]
    pub call: bool,

    /// Ploidy for the call stage: 1 (haploid), 2 (diploid) or 0 (no call).
    /// Ignored unless `--call` is set.
    #[arg(long = "call-ploidy", value_name = "N", default_value_t = 2)]
    pub call_ploidy: u8,

    /// Region-dependent ploidy for the call stage: `chrom beg end sex ploidy`,
    /// as in `bcftools call --ploidy-file`. The fused pipeline has one sample
    /// and no sex table, so only `*` sex rows apply here.
    #[arg(long = "call-ploidy-file", value_name = "FILE")]
    pub call_ploidy_file: Option<PathBuf>,

    /// Use the consensus caller (`call -c`) instead of the multiallelic one.
    #[arg(long = "call-consensus", default_value_t = false)]
    pub call_consensus: bool,

    /// Extra FORMAT tags from the call stage, e.g. `GQ` or `GQ,GP`.
    #[arg(long = "call-annotate", value_name = "TAGS")]
    pub call_annotate: Option<String>,

    /// Mutation rate prior for the call stage (`bcftools call -P`).
    #[arg(long = "call-prior", value_name = "FLOAT", default_value_t = 1.1e-3)]
    pub call_prior: f64,

    /// Number of threads.
    #[arg(short = 't', long = "threads", default_value_t = 8)]
    pub threads: usize,

    /// Process the genome in reference windows of this many megabases instead of
    /// holding every alignment in memory at once.
    ///
    /// The default (0) keeps the whole run resident — fastest, but peak memory
    /// scales with the input (~12 GB for a 30x human chr20). With a window size
    /// set, alignments are spilled to temporary BAMs bucketed by reference window
    /// and each window is sorted, deduplicated and called on its own, bounding
    /// peak memory by one window's depth. Smaller windows mean lower peak memory
    /// and more temporary I/O; 32 is a reasonable starting point.
    #[arg(long = "window-mb", value_name = "MB", default_value_t = 0)]
    pub window_mb: u32,

    /// Directory for the windowed mode's temporary BAMs (default: alongside the
    /// output VCF).
    #[arg(long = "window-tmpdir", value_name = "DIR")]
    pub window_tmpdir: Option<PathBuf>,

    /// Output VCF path — the only artifact written to disk.
    #[arg(short = 'o', long = "output", value_name = "VCF")]
    pub output: PathBuf,
}
