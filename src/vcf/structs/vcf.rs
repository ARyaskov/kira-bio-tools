use std::collections::HashMap;
use std::io;
use thiserror::Error;

use crate::util::{ChrId, GenomicKey};

#[derive(Debug, Error)]
pub enum VcfError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
    #[error("BGZF error: {0}")]
    Bgzf(#[from] crate::bgzf::BgzfError),
    #[error("Invalid VCF format: {0}")]
    InvalidFormat(String),
    #[error("Parse error at line {line}: {message}")]
    ParseError { line: usize, message: String },
}

pub type Result<T> = std::result::Result<T, VcfError>;

#[derive(Debug, Clone)]
pub struct VcfRecord {
    pub chr_id: ChrId,
    pub position: u32,
    pub offset: u64,
}

impl VcfRecord {
    pub fn key(&self) -> GenomicKey {
        GenomicKey::new(self.chr_id, self.position)
    }
}

#[derive(Debug, Clone)]
pub struct VcfParsedRecord {
    pub chrom: String,
    pub pos: u32,
    pub filter: String,
    pub info: HashMap<String, String>,
    pub raw_line: String,
}

impl VcfParsedRecord {
    pub fn to_line(&self) -> &str {
        &self.raw_line
    }
}
