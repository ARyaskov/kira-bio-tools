use clap::Parser;
use std::path::PathBuf;

#[derive(Parser)]
pub struct QueryCompatArgs {
    #[arg(allow_hyphen_values = true, trailing_var_arg = true)]
    pub bcftools_args: Vec<String>,
}

#[derive(Parser)]
#[command(disable_help_flag = true)]
pub struct RegionQueryArgs {
    #[arg(long = "help", action = clap::ArgAction::Help, help = "Print help")]
    pub help: Option<bool>,

    #[arg(help = "Indexed VCF file")]
    pub file: PathBuf,

    #[arg(help = "Regions to query (chr:start-end or chr:pos)")]
    pub regions: Vec<String>,

    #[arg(short = 'R', long, help = "File with regions (one per line)")]
    pub regions_file: Option<PathBuf>,

    #[arg(short, long, help = "Print only matching record count")]
    pub count: bool,

    #[arg(short = 'h', long, help = "Include header in output")]
    pub print_header: bool,

    #[arg(short = 'H', long, help = "Print only header")]
    pub only_header: bool,
}

#[derive(Parser)]
pub struct StatArgs {
    #[arg(help = "Index file (.kbi or .csi)")]
    pub index: PathBuf,
}

#[derive(Parser)]
pub struct ListArgs {
    #[arg(help = "Indexed VCF file")]
    pub file: PathBuf,
}

#[derive(Parser)]
pub struct HeaderArgs {
    #[arg(help = "VCF file")]
    pub file: PathBuf,
}
