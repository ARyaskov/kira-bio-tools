use clap::Parser;
use std::path::PathBuf;

/// HMM model for detecting runs of autozygosity (`bcftools roh`).
#[derive(Parser, Debug)]
pub struct RohArgs {
    pub input: Option<PathBuf>,

    /// If AF is not known, use this allele frequency [skip]
    #[arg(long = "AF-dflt")]
    pub af_dflt: Option<f64>,

    /// Use TAG for allele frequency
    #[arg(long = "AF-tag")]
    pub af_tag: Option<String>,

    /// Read allele frequencies from file (CHR\tPOS\tREF,ALT\tAF)
    #[arg(long = "AF-file")]
    pub af_file: Option<PathBuf>,

    /// Buffer size and the number of overlapping sites, 0 for unlimited
    #[arg(short = 'b', long = "buffer-size")]
    pub buffer_size: Option<String>,

    /// Estimate AF from FORMAT/TAG (GT or PL) of all samples ("-") or samples listed in FILE
    #[arg(short = 'e', long = "estimate-AF")]
    pub estimate_af: Option<String>,

    /// Exclude sites for which the expression is true
    #[arg(long = "exclude")]
    pub exclude: Option<String>,

    /// Use GTs and ignore PLs, instead using FLOAT for PL of the two least likely genotypes
    #[arg(short = 'G', long = "GTs-only")]
    pub gts_only: Option<f64>,

    /// Select sites for which the expression is true
    #[arg(long = "include")]
    pub include: Option<String>,

    /// Skip hom-ref genotypes (0/0)
    #[arg(short = 'i', long = "ignore-homref")]
    pub ignore_homref: bool,

    /// Include sites with no ALT allele (ignored by default)
    #[arg(long = "include-noalt")]
    pub include_noalt: bool,

    /// Skip indels as their genotypes are enriched for errors
    #[arg(short = 'I', long = "skip-indels")]
    pub skip_indels: bool,

    /// Genetic map in IMPUTE2 format, single file or mask with "{CHROM}"
    #[arg(short = 'm', long = "genetic-map")]
    pub genetic_map: Option<String>,

    /// Constant recombination rate per bp
    #[arg(short = 'M', long = "rec-rate")]
    pub rec_rate: Option<f64>,

    /// Write output to a file [standard output]
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    /// Output s:per-site, r:regions, z:compressed [sr]
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

    /// List of samples to analyze [all samples]
    #[arg(short = 's', long = "samples")]
    pub samples: Option<String>,

    /// File of samples to analyze [all samples]
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

    /// P(AZ|HW) transition probability from HW (Hardy-Weinberg) to AZ (autozygous) state
    #[arg(short = 'a', long = "hw-to-az")]
    pub hw_to_az: Option<f64>,

    /// P(HW|AZ) transition probability from AZ to HW state
    #[arg(short = 'H', long = "az-to-hw")]
    pub az_to_hw: Option<f64>,

    /// Estimate HMM parameters, FLOAT is the convergence threshold (experimental)
    #[arg(short = 'V', long = "viterbi-training")]
    pub viterbi_training: Option<f64>,

    /// bcftools-style options given after `--`
    #[arg(last = true, allow_hyphen_values = true)]
    pub passthrough: Vec<String>,
}
