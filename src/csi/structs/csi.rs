use std::io;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum CsiError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
    #[error("VCF error: {0}")]
    Vcf(#[from] crate::vcf::VcfError),
    #[error("CSI format error: {0}")]
    Format(String),
}

pub type Result<T> = std::result::Result<T, CsiError>;
