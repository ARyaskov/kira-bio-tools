//! Random access to single VCF lines by byte offset (plain text) or BGZF
//! virtual position. The file stays open between fetches and the last
//! decompressed block is cached, so consecutive hits in one block cost no I/O.

use std::fs::File;
use std::io::{BufRead, BufReader, Seek, SeekFrom};
use std::path::Path;

use crate::bgzf::{BgzfReader, VirtualPosition};
use crate::util::{VcfFormat, detect_format};
use crate::vcf::structs::{Result, VcfError};

enum Inner {
    Plain(BufReader<File>),
    Bgzf { reader: BgzfReader<File>, block_cpos: Option<u64>, block: Vec<u8> },
}

pub struct LineFetcher {
    inner: Inner,
    buf: String,
}

impl LineFetcher {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let inner = match detect_format(path)? {
            VcfFormat::Plain => Inner::Plain(BufReader::new(File::open(path)?)),
            VcfFormat::Gzip => return Err(VcfError::InvalidFormat),
            VcfFormat::Bgzf => Inner::Bgzf { reader: BgzfReader::open(path)?, block_cpos: None, block: Vec::new() },
        };
        Ok(Self { inner, buf: String::new() })
    }

    /// The line starting at `offset`, without its line terminator.
    pub fn fetch(&mut self, offset: u64) -> Result<&str> {
        self.buf.clear();
        match &mut self.inner {
            Inner::Plain(r) => {
                r.seek(SeekFrom::Start(offset))?;
                r.read_line(&mut self.buf)?;
            }
            Inner::Bgzf { reader, block_cpos, block } => {
                let cpos = offset >> 16;
                let upos = (offset & 0xffff) as usize;
                if *block_cpos != Some(cpos) {
                    reader.seek(VirtualPosition::from_raw(cpos << 16))?;
                    block.clear();
                    block.extend_from_slice(reader.fill_buf()?);
                    *block_cpos = Some(cpos);
                }
                let tail = block.get(upos..).unwrap_or(&[]);
                match memchr::memchr(b'\n', tail) {
                    Some(nl) => self.buf.push_str(&String::from_utf8_lossy(&tail[..=nl])),
                    None => {
                        // The line continues in the next block: read it through the reader.
                        reader.seek(VirtualPosition::from_raw(offset))?;
                        reader.read_line(&mut self.buf)?;
                        *block_cpos = None;
                    }
                }
            }
        }
        Ok(self.buf.trim_end_matches(['\r', '\n']))
    }
}

pub fn fetch_line<P: AsRef<Path>>(path: P, offset: u64) -> Result<String> {
    let mut f = LineFetcher::open(path)?;
    Ok(f.fetch(offset)?.to_string())
}
