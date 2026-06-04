use clap::{ArgGroup, Parser};
use std::path::PathBuf;

#[derive(Parser)]
pub struct AnnotateIndexArgs {
    #[arg(help = "Annotation VCF to index")]
    pub input: PathBuf,

    #[arg(short, long, help = "Output .ani file")]
    pub output: Option<PathBuf>,
}

#[derive(Parser)]
#[command(group(
    ArgGroup::new("anno_source")
        .required(true)
        .args(["annotations", "ani"])
))]
pub struct AnnotateArgs {
    #[arg(help = "Input VCF file (.vcf)")]
    pub input: PathBuf,

    #[arg(short = 'a', long = "annotations", help = "Annotation DB (.vcf, .tab, or .ani)")]
    pub annotations: Option<PathBuf>,

    #[arg(long = "ani", help = "Annotation index (.ani), no auto-build")]
    pub ani: Option<PathBuf>,

    #[arg(short = 'o', long = "output", help = "Output file (optional)")]
    pub output: Option<PathBuf>,

    #[arg(
        short = 'O',
        long = "output-type",
        help = "Output type: v=VCF, z=BGZF-VCF, u=uncompressed BCF, b=compressed BCF; suffix 0-9 for level (e.g. z6)"
    )]
    pub output_type: Option<String>,

    #[arg(short = 'c', long = "columns", help = "Column spec: CHROM,POS,REF,ALT,ID,+INFO/AC")]
    pub columns: Option<String>,

    #[arg(short = 'C', long = "columns-file", help = "Read -c columns from file (NAME[ TYPE] per line)")]
    pub columns_file: Option<PathBuf>,

    #[arg(short = 'h', long = "header-lines", help = "File with header lines to append")]
    pub header_lines: Option<PathBuf>,

    #[arg(short = 'H', long = "header-line", help = "Single header line to append (repeatable)")]
    pub header_line: Vec<String>,

    #[arg(short = 'i', long = "include", help = "Include sites where EXPR is true")]
    pub include: Option<String>,

    #[arg(short = 'e', long = "exclude", help = "Exclude sites where EXPR is true")]
    pub exclude: Option<String>,

    #[arg(short = 'k', long = "keep-sites", help = "Keep -i/-e excluded sites unchanged instead of dropping")]
    pub keep_sites: bool,

    #[arg(short = 'I', long = "set-id", help = "Set ID column using %CHROM/%POS/%REF/%ALT/%INFO_X template; prefix + to fill only missing")]
    pub set_id: Option<String>,

    #[arg(short = 'm', long = "mark-sites", help = "Add INFO/TAG (+TAG) for matched or (-TAG) for unmatched sites")]
    pub mark_sites: Option<String>,

    #[arg(short = 'x', long = "remove", help = "Comma list: ID, QUAL, FILTER, INFO[/TAG], FORMAT[/TAG]; ^prefix keeps listed")]
    pub remove: Option<String>,

    #[arg(short = 'r', long = "regions", help = "Comma list of regions: CHR or CHR:beg or CHR:beg-end")]
    pub regions: Option<String>,

    #[arg(short = 'R', long = "regions-file", help = "Regions file (BED 0-based or TSV 1-based: CHR[\\tBEG[\\tEND]])")]
    pub regions_file: Option<PathBuf>,

    #[arg(long = "regions-overlap", help = "Overlap mode: 0=pos in region, 1=record overlaps, 2=variant overlaps", default_value = "1")]
    pub regions_overlap: u8,

    #[arg(short = 's', long = "samples", help = "Sample list: name1,name2; prefix ^ inverts")]
    pub samples: Option<String>,

    #[arg(short = 'S', long = "samples-file", help = "Samples file (one per line); prefix ^ inverts")]
    pub samples_file: Option<PathBuf>,

    #[arg(long = "pair-logic", help = "REF/ALT match: snps|indels|both|all|some|exact|id", default_value = "some")]
    pub pair_logic: String,

    #[arg(long = "rename-chrs", help = "File with `old\\tnew` per line")]
    pub rename_chrs: Option<PathBuf>,

    #[arg(long = "rename-annots", help = "File with `TYPE/old\\tnew` (TYPE=INFO|FORMAT|FILTER)")]
    pub rename_annots: Option<PathBuf>,

    #[arg(long = "no-version", help = "Do not append the kira_bt version line to the output header")]
    pub no_version: bool,

    #[arg(long = "force", help = "Continue on parse errors")]
    pub force: bool,

    #[arg(short = 'v', long = "verbosity", help = "Verbosity 0-3", default_value_t = 1)]
    pub verbosity: u8,

    #[arg(long = "threads", help = "Output compression worker threads", default_value_t = 0)]
    pub threads: usize,

    #[arg(short = 'W', long = "write-index", num_args = 0..=1, default_missing_value = "csi", help = "Auto-index output [csi|tbi]")]
    pub write_index: Option<String>,

    #[arg(long, help = "Use GPU (CUDA) for annotation")]
    pub gpu: bool,

    #[arg(long, help = "Use OpenCL backend")]
    pub opencl: bool,

    #[arg(long, help = "BGZF compression level (0-9, default 1)", value_parser = clap::value_parser!(u32))]
    pub bgzf_level: Option<u32>,

    #[arg(long, help = "Use cached plain VCF in CWD")]
    pub cache_plain: bool,

    #[arg(long, help = "Write plain VCF first, then BGZF-compress")]
    pub bgzf_after: bool,

    #[arg(long, help = "Use memory-mapped output for plain VCF")]
    pub mmap_output: bool,

    #[arg(long, help = "Do not flush mmap output explicitly")]
    pub mmap_no_flush: bool,

    #[arg(long, help = "Write output to RAM only")]
    pub ram_output: bool,

    #[arg(long, help = "Max RAM output size in MB", default_value_t = 1024)]
    pub ram_max_mb: u32,

    #[arg(long = "no-ktile", help = "Disable .ktile sidecar entirely")]
    pub no_ktile: bool,

    #[arg(long = "no-build-ktile", help = "Use .ktile if present but don't auto-build")]
    pub no_build_ktile: bool,

    #[arg(long = "force-build-ktile", help = "Override auto-skip heuristic")]
    pub force_build_ktile: bool,
}

#[derive(Parser)]
#[command(group(
    ArgGroup::new("anno_source")
        .required(true)
        .args(["annotations", "ani"])
))]
pub struct AnnotateServeArgs {
    #[arg(short = 'a', long = "annotations", help = "Annotation DB (.vcf, .tab, or .ani)")]
    pub annotations: Option<PathBuf>,

    #[arg(long = "ani", help = "Annotation index (.ani), no auto-build")]
    pub ani: Option<PathBuf>,

    #[arg(short = 'c', long = "columns", help = "Default column specification")]
    pub columns: Option<String>,

    #[arg(long, help = "Use GPU (CUDA) for annotation")]
    pub gpu: bool,

    #[arg(long, help = "Use OpenCL backend")]
    pub opencl: bool,

    #[arg(long, help = "BGZF compression level (0-9, default 1)", value_parser = clap::value_parser!(u32))]
    pub bgzf_level: Option<u32>,

    #[arg(long, help = "Use cached plain VCF in CWD or create it from .gz input")]
    pub cache_plain: bool,

    #[arg(long, help = "Write plain VCF first, then BGZF-compress")]
    pub bgzf_after: bool,

    #[arg(long, help = "Use memory-mapped output for plain VCF")]
    pub mmap_output: bool,

    #[arg(long, help = "Do not flush mmap output explicitly")]
    pub mmap_no_flush: bool,

    #[arg(long, help = "Write output to RAM only (no disk)")]
    pub ram_output: bool,

    #[arg(long, help = "Max RAM output size in MB (unused)", default_value_t = 1024)]
    pub ram_max_mb: u32,
}

#[derive(Parser)]
pub struct DbBuildArgs {
    #[arg(help = "Input annotation database VCF")]
    pub input: PathBuf,

    #[arg(short, long, help = "Output ANI file")]
    pub output: Option<PathBuf>,
}
