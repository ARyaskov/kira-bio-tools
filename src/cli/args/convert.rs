use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct ConvertArgs {
    pub input: Option<PathBuf>,

    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    #[arg(short = 'O', long = "output-type")]
    pub output_type: Option<String>,

    #[arg(long = "gvcf2vcf")]
    pub gvcf2vcf: bool,

    #[arg(long = "tsv2vcf")]
    pub tsv2vcf: Option<PathBuf>,

    #[arg(short = 'G', long = "gensample")]
    pub gensample: Option<String>,

    #[arg(short = 'g', long = "gensample2vcf")]
    pub gen2vcf: Option<String>,

    #[arg(long = "haplegendsample")]
    pub haplegendsample: Option<String>,

    #[arg(long = "haplegendsample2vcf")]
    pub haplegend2vcf: Option<String>,

    #[arg(short = 'H', long = "hapsample")]
    pub hapsample: Option<String>,

    #[arg(long = "hapsample2vcf")]
    pub hap2vcf: Option<String>,

    #[arg(long = "haploid")]
    pub haploid: bool,

    #[arg(short = 'c', long = "columns")]
    pub columns: Option<String>,

    #[arg(short = 'f', long = "fasta-ref")]
    pub fasta_ref: Option<PathBuf>,

    #[arg(long = "chrom")]
    pub chrom: bool,

    #[arg(long = "vcf-ids")]
    pub vcf_ids: bool,

    #[arg(long = "sex")]
    pub sex: Option<PathBuf>,

    #[arg(long = "keep-duplicates")]
    pub keep_duplicates: bool,

    #[arg(long = "tag")]
    pub tag: Option<String>,

    #[arg(short = 's', long = "samples")]
    pub samples: Option<String>,

    #[arg(short = 'S', long = "samples-file")]
    pub samples_file: Option<PathBuf>,

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

    #[arg(short = 'i', long = "include")]
    pub include: Option<String>,

    #[arg(short = 'e', long = "exclude")]
    pub exclude: Option<String>,

    #[arg(long = "no-version")]
    pub no_version: bool,

    #[arg(long = "threads", default_value_t = 0)]
    pub threads: usize,

    #[arg(short = 'v', long = "verbosity", default_value_t = 1)]
    pub verbosity: u8,

    #[arg(short = 'W', long = "write-index", num_args = 0..=1, default_missing_value = "csi")]
    pub write_index: Option<String>,

    #[arg(last = true, allow_hyphen_values = true)]
    pub passthrough: Vec<String>,
}
