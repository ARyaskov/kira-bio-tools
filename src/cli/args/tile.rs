use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
pub struct TileArgs {
    #[command(subcommand)]
    pub action: TileAction,
}

#[derive(Subcommand)]
pub enum TileAction {
    /// Build a `.ktile` sidecar from a VCF/BGZF input. The sidecar speeds
    /// up subsequent annotate runs by skipping the BGZF decode + per-line
    /// parse pass.
    Build(TileBuildArgs),
}

#[derive(Parser)]
pub struct TileBuildArgs {
    /// Input VCF or VCF.gz file.
    #[arg(short = 'i', long)]
    pub input: PathBuf,

    /// Output `.ktile` path. Defaults to `<input>.ktile`.
    #[arg(short = 'o', long)]
    pub output: Option<PathBuf>,
}
