use std::fs::File;
use std::io::Read;
use std::path::Path;

use crate::bgzf::structs::Result;

pub fn is_bgzf<P: AsRef<Path>>(path: P) -> Result<bool> {
    let mut file = File::open(path)?;
    let mut header = [0u8; 18];
    let n = file.read(&mut header)?;

    if n < 18 {
        return Ok(false);
    }

    Ok(header[0] == 0x1f
        && header[1] == 0x8b
        && header[2] == 0x08
        && header[3] == 0x04
        && header[12] == b'B'
        && header[13] == b'C')
}
