use std::fs::File;
use std::io::Read;
use std::path::Path;

use crate::bgzf::structs::{Result, is_bgzf_header};

pub fn is_bgzf<P: AsRef<Path>>(path: P) -> Result<bool> {
    let mut file = File::open(path)?;
    let mut header = [0u8; 18];
    let n = file.read(&mut header)?;
    Ok(is_bgzf_header(&header[..n]))
}
