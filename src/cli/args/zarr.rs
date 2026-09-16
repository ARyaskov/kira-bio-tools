use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct ZarrArgs {
    /// Input VCF/BCF.
    pub input: PathBuf,
    /// Output Zarr store directory.
    #[arg(short = 'o', long = "output")]
    pub output: PathBuf,
    /// Variants per chunk.
    #[arg(long = "chunk-variants", default_value_t = 10_000)]
    pub chunk_variants: usize,
    /// Samples per chunk.
    #[arg(long = "chunk-samples", default_value_t = 1_000)]
    pub chunk_samples: usize,
    /// zlib compression level (0-9).
    #[arg(short = 'l', long = "level", default_value_t = 6)]
    pub level: u32,
    /// `INFO/TAG,FORMAT/TAG,...` to export besides GT (default: every
    /// Number=1 Integer/Float field and every INFO Flag).
    #[arg(long = "fields")]
    pub fields: Option<String>,
    /// Overwrite an existing store.
    #[arg(short = 'f', long = "force")]
    pub force: bool,
}
