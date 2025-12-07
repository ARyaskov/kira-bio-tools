use std::collections::HashMap;
use std::fmt;
use std::io;

use crate::bgzf::BgzfError;
use crate::util::{ChrId, GenomicKey};

#[derive(Debug)]
pub enum VcfError {
    Io(io::Error),
    ParseError(String),
    InvalidFormat,
    MissingField(String),
}

impl fmt::Display for VcfError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            VcfError::Io(e) => write!(f, "IO error: {}", e),
            VcfError::ParseError(s) => write!(f, "Parse error: {}", s),
            VcfError::InvalidFormat => write!(f, "Invalid VCF format"),
            VcfError::MissingField(s) => write!(f, "Missing field: {}", s),
        }
    }
}

impl std::error::Error for VcfError {}

impl From<io::Error> for VcfError {
    fn from(error: io::Error) -> Self {
        VcfError::Io(error)
    }
}

impl From<BgzfError> for VcfError {
    fn from(error: BgzfError) -> Self {
        match error {
            BgzfError::Io(e) => VcfError::Io(e),
            _ => VcfError::InvalidFormat,
        }
    }
}

pub type Result<T> = std::result::Result<T, VcfError>;

#[derive(Debug, Clone)]
pub struct VcfRecord {
    pub chrom: String,
    pub pos: u32,
    pub id: String,
    pub ref_allele: String,
    pub alt: String,
    pub qual: String,
    pub filter: String,
    pub info: String,
    pub format: Option<String>,
    pub samples: Vec<String>,
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

#[derive(Clone, Copy)]
pub struct VcfFields<'a> {
    pub chrom: &'a str,
    pub pos: &'a str,
    pub id: &'a str,
    pub ref_allele: &'a str,
    pub alt: &'a str,
    pub qual: &'a str,
    pub filter: &'a str,
    pub info: &'a str,
}

impl<'a> VcfFields<'a> {
    #[inline]
    pub fn position(&self) -> Option<u32> {
        fast_parse_u32(self.pos.as_bytes())
    }
}

#[derive(Clone)]
pub struct VcfFieldsFull<'a> {
    pub chrom: &'a str,
    pub pos: &'a str,
    pub id: &'a str,
    pub ref_allele: &'a str,
    pub alt: &'a str,
    pub qual: &'a str,
    pub filter: &'a str,
    pub info: &'a str,
    pub format: Option<&'a str>,
    pub samples: Vec<&'a str>,
}

impl<'a> VcfFieldsFull<'a> {
    #[inline]
    pub fn position(&self) -> Option<u32> {
        fast_parse_u32(self.pos.as_bytes())
    }

    #[inline]
    pub fn as_basic_fields(&self) -> VcfFields<'a> {
        VcfFields {
            chrom: self.chrom,
            pos: self.pos,
            id: self.id,
            ref_allele: self.ref_allele,
            alt: self.alt,
            qual: self.qual,
            filter: self.filter,
            info: self.info,
        }
    }
}

#[inline]
pub fn fast_parse_u32(bytes: &[u8]) -> Option<u32> {
    if bytes.is_empty() || bytes.len() > 10 {
        return None;
    }

    let mut result = 0u32;
    for &byte in bytes {
        if !byte.is_ascii_digit() {
            return None;
        }
        result = result.wrapping_mul(10).wrapping_add((byte - b'0') as u32);
    }

    Some(result)
}
