use std::io::{BufRead, Seek};
use std::path::Path;

use crate::bgzf::{BgzfReader, VirtualPosition};
use crate::util::{detect_format, VcfFormat};
use crate::vcf::structs::Result;

pub fn fetch_line<P: AsRef<Path>>(path: P, offset: u64) -> Result<String> {
    let path = path.as_ref();
    let format = detect_format(path)?;

    match format {
        VcfFormat::Plain => {
            let file = std::fs::File::open(path)?;
            let mut reader = std::io::BufReader::new(file);
            reader.seek(std::io::SeekFrom::Start(offset))?;

            let mut line = String::new();
            reader.read_line(&mut line)?;
            Ok(line)
        }
        VcfFormat::Gzip => Err(crate::vcf::structs::VcfError::InvalidFormat),
        VcfFormat::Bgzf => {
            let mut bgzf_reader = BgzfReader::open(path)?;
            bgzf_reader.seek(VirtualPosition::from_raw(offset))?;

            let mut line = String::new();
            bgzf_reader.read_line(&mut line)?;
            Ok(line)
        }
    }
}
