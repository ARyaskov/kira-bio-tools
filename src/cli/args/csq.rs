use clap::Parser;
use std::path::PathBuf;

/// Haplotype-aware consequence caller (`bcftools csq`).
#[derive(Parser, Debug)]
pub struct CsqArgs {
    pub input: Option<PathBuf>,

    /// Reference file in fasta format
    #[arg(short = 'f', long = "fasta-ref")]
    pub fasta_ref: Option<PathBuf>,

    /// GFF3 annotation file
    #[arg(short = 'g', long = "gff-annot")]
    pub gff: Option<PathBuf>,

    /// Abbreviate protein-changing predictions to max INT aminoacids
    #[arg(short = 'B', long = "trim-protein-seq")]
    pub trim_protein_seq: Option<usize>,

    #[arg(short = 'b', long = "brief-predictions", hide = true)]
    pub brief_predictions: bool,

    /// Use this tag instead of the default BCSQ
    #[arg(short = 'c', long = "custom-tag", default_value = "BCSQ")]
    pub custom_tag: String,

    /// Localized predictions, consider only one VCF record at a time
    #[arg(short = 'l', long = "local-csq")]
    pub local_csq: bool,

    /// Maximum number of per-haplotype consequences to consider for each site
    #[arg(short = 'n', long = "ncsq", default_value_t = 15)]
    pub ncsq: usize,

    /// How to handle unphased heterozygous genotypes: a|m|r|R|s
    #[arg(short = 'p', long = "phase", default_value = "r")]
    pub phase: String,

    /// Run even if some sanity checks fail
    #[arg(long = "force")]
    pub force: bool,

    /// Automatically unify chromosome naming (e.g. chrX vs X) in GFF, fasta, and VCF
    #[arg(long = "unify-chr-names", default_value_t = 1)]
    pub unify_chr_names: u8,

    /// Exclude sites for which the expression is true
    #[arg(short = 'e', long = "exclude")]
    pub exclude: Option<String>,

    /// Select sites for which the expression is true
    #[arg(short = 'i', long = "include")]
    pub include: Option<String>,

    /// Do not append version and command line to the header
    #[arg(long = "no-version")]
    pub no_version: bool,

    /// Write output to a file [standard output]
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    /// b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF, t: plain tab-delimited text
    #[arg(short = 'O', long = "output-type")]
    pub output_type: Option<String>,

    /// Restrict to comma-separated list of regions
    #[arg(short = 'r', long = "regions")]
    pub regions: Option<String>,

    /// Restrict to regions listed in a file
    #[arg(short = 'R', long = "regions-file")]
    pub regions_file: Option<PathBuf>,

    /// Include if POS in the region (0), record overlaps (1), variant overlaps (2)
    #[arg(long = "regions-overlap", default_value_t = 1)]
    pub regions_overlap: u8,

    /// Samples to include or "-" to apply all variants and ignore samples
    #[arg(short = 's', long = "samples")]
    pub samples: Option<String>,

    /// Samples to include
    #[arg(short = 'S', long = "samples-file")]
    pub samples_file: Option<PathBuf>,

    /// Similar to -r but streams rather than index-jumps
    #[arg(short = 't', long = "targets")]
    pub targets: Option<String>,

    /// Similar to -R but streams rather than index-jumps
    #[arg(short = 'T', long = "targets-file")]
    pub targets_file: Option<PathBuf>,

    /// Include if POS in the region (0), record overlaps (1), variant overlaps (2)
    #[arg(long = "targets-overlap", default_value_t = 0)]
    pub targets_overlap: u8,

    /// Use multithreading with INT worker threads
    #[arg(long = "threads", default_value_t = 0)]
    pub threads: usize,

    /// Verbosity level 0-2
    #[arg(
        short = 'v',
        long = "verbose",
        alias = "verbosity",
        default_value_t = 1
    )]
    pub verbosity: u8,

    /// Automatically index the output files
    #[arg(long = "write-index")]
    pub write_index: bool,
}
