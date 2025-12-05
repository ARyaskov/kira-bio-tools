use clap::Parser;
use std::path::PathBuf;

#[derive(Parser)]
pub struct TabixArgs {
    #[arg(help = "Input file (for indexing) or indexed file (for querying)")]
    pub input: PathBuf,

    #[arg(help = "Regions to query (chr:start-end or chr:pos)")]
    pub regions: Vec<String>,

    #[arg(
        short = '0',
        long = "zero-based",
        help = "Position is 0-based half-open"
    )]
    pub zero_based: bool,

    #[arg(
        short = 'b',
        long = "begin",
        help = "Column of start chromosomal position [4]"
    )]
    pub begin_col: Option<usize>,

    #[arg(
        short = 'c',
        long = "comment",
        help = "Skip lines started with character [#]"
    )]
    pub comment_char: Option<char>,

    #[arg(
        short = 'C',
        long = "csi",
        help = "Produce CSI format index instead of TBI"
    )]
    pub csi: bool,

    #[arg(
        short = 'e',
        long = "end",
        help = "Column of end chromosomal position [5]"
    )]
    pub end_col: Option<usize>,

    #[arg(short = 'f', long = "force", help = "Force overwrite of index file")]
    pub force: bool,

    #[arg(
        short = 'm',
        long = "min-shift",
        help = "Set minimal interval size for CSI to 2^INT [14]"
    )]
    pub min_shift: Option<u8>,

    #[arg(
        short = 'p',
        long = "preset",
        help = "Input format: gff, bed, sam, vcf"
    )]
    pub preset: Option<String>,

    #[arg(short = 's', long = "sequence", help = "Column of sequence name [1]")]
    pub sequence_col: Option<usize>,

    #[arg(short = 'S', long = "skip-lines", help = "Skip first INT lines [0]")]
    pub skip_lines: Option<usize>,

    #[arg(
        short = 'h',
        long = "print-header",
        help = "Print also the header/meta lines"
    )]
    pub print_header: bool,

    #[arg(short = 'H', long = "only-header", help = "Print only the header")]
    pub only_header: bool,

    #[arg(short = 'l', long = "list-chroms", help = "List chromosome names")]
    pub list_chroms: bool,

    #[arg(
        short = 'r',
        long = "reheader",
        help = "Replace header with content of FILE"
    )]
    pub reheader: Option<PathBuf>,

    #[arg(
        short = 'R',
        long = "regions",
        help = "Restrict to regions listed in FILE"
    )]
    pub regions_file: Option<PathBuf>,

    #[arg(short = 't', long = "targets", help = "Similar to -R but sequential")]
    pub targets_file: Option<PathBuf>,

    #[arg(
        short = 'D',
        long = "no-download",
        help = "Do not download index for remote files"
    )]
    pub no_download: bool,

    #[arg(
        long = "cache",
        help = "Set BGZF block cache size in MB [10]",
        default_value = "10"
    )]
    pub cache_size: usize,

    #[arg(
        long = "regions-overlap",
        help = "Output region names before records (0=off, 1=on, 2=tab-separated)"
    )]
    pub regions_overlap: Option<u8>,

    #[arg(
        long = "verbosity",
        help = "Set log verbosity (0=silent, 3=default, 4+=debug)"
    )]
    pub verbosity: Option<u8>,

    #[arg(short = '@', long = "threads", help = "Number of threads [0]")]
    pub threads: Option<usize>,
}
