use std::fs::File;
use std::io::{BufRead, BufReader, Seek, SeekFrom};
use std::path::Path;

use crate::bgzf::{BgzfReader, VirtualPosition};
use crate::util::{detect_format, VcfFormat};
use crate::vcf::structs::{Result, VcfError};

pub fn fetch_line<P: AsRef<Path>>(path: P, offset: u64) -> Result<String> {
    let path = path.as_ref();
    let format = detect_format(path)?;

    match format {
        VcfFormat::Plain => {
            let mut file = File::open(path)?;
            file.seek(SeekFrom::Start(offset))?;
            let mut reader = BufReader::new(file);
            let mut line = String::new();
            reader.read_line(&mut line)?;
            Ok(line.trim_end().to_string())
        }
        VcfFormat::Gzip => Err(VcfError::InvalidFormat(
            "Cannot seek in gzip file. Use BGZF compression.".into(),
        )),
        VcfFormat::Bgzf => {
            let file = File::open(path)?;
            let mut bgzf_reader = BgzfReader::new(file);
            bgzf_reader.seek(VirtualPosition::from_u64(offset))?;

            let mut line = String::new();
            bgzf_reader.read_line(&mut line)?;

            if line.ends_with('\n') {
                line.pop();
            }
            if line.ends_with('\r') {
                line.pop();
            }
            Ok(line)
        }
    }
}
