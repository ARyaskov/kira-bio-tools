use std::fs::File;
use std::io::{self, BufRead, Read, Write};
use std::path::Path;

use noodles_bgzf as bgzf;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum BgzfError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
    #[error("Invalid BGZF format")]
    InvalidFormat,
}

pub type Result<T> = std::result::Result<T, BgzfError>;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct VirtualPosition(u64);

impl VirtualPosition {
    #[inline]
    pub fn new(compressed_offset: u64, uncompressed_offset: u16) -> Self {
        Self((compressed_offset << 16) | (uncompressed_offset as u64))
    }

    #[inline]
    pub fn compressed(&self) -> u64 {
        self.0 >> 16
    }

    #[inline]
    pub fn uncompressed(&self) -> u16 {
        (self.0 & 0xFFFF) as u16
    }

    #[inline]
    pub fn as_u64(&self) -> u64 {
        self.0
    }

    #[inline]
    pub fn from_u64(v: u64) -> Self {
        Self(v)
    }
}

impl From<bgzf::VirtualPosition> for VirtualPosition {
    fn from(vp: bgzf::VirtualPosition) -> Self {
        Self(u64::from(vp))
    }
}

impl From<VirtualPosition> for bgzf::VirtualPosition {
    fn from(vp: VirtualPosition) -> Self {
        bgzf::VirtualPosition::try_from(vp.0).unwrap_or_default()
    }
}

pub struct BgzfReader<R: Read> {
    inner: bgzf::Reader<R>,
}

impl BgzfReader<File> {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        Ok(Self::new(file))
    }
}

impl<R: Read> BgzfReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            inner: bgzf::Reader::new(reader),
        }
    }

    pub fn virtual_position(&self) -> VirtualPosition {
        self.inner.virtual_position().into()
    }
}

impl<R: Read> Read for BgzfReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.inner.read(buf)
    }
}

impl<R: Read> BufRead for BgzfReader<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        self.inner.fill_buf()
    }

    fn consume(&mut self, amt: usize) {
        self.inner.consume(amt)
    }
}

impl<R: Read + io::Seek> BgzfReader<R> {
    pub fn seek(&mut self, pos: VirtualPosition) -> Result<VirtualPosition> {
        let noodles_pos: bgzf::VirtualPosition = pos.into();
        self.inner.seek(noodles_pos).map_err(BgzfError::Io)?;
        Ok(self.virtual_position())
    }
}

pub struct BgzfWriter<W: Write> {
    inner: bgzf::Writer<W>,
}

impl BgzfWriter<File> {
    pub fn create<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::create(path)?;
        Ok(Self::new(file))
    }
}

impl<W: Write> BgzfWriter<W> {
    pub fn new(writer: W) -> Self {
        Self {
            inner: bgzf::Writer::new(writer),
        }
    }

    pub fn finish(self) -> Result<W> {
        self.inner.finish().map_err(BgzfError::Io)
    }
}

impl<W: Write> Write for BgzfWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.inner.write(buf)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

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

pub struct BgzfLineReader<R: Read> {
    reader: BgzfReader<R>,
    buf: String,
}

impl<R: Read> BgzfLineReader<R> {
    pub fn new(reader: BgzfReader<R>) -> Self {
        Self {
            reader,
            buf: String::with_capacity(4096),
        }
    }

    pub fn read_line(&mut self) -> Result<Option<(&str, VirtualPosition)>> {
        self.buf.clear();
        let vpos = self.reader.virtual_position();

        match self.reader.read_line(&mut self.buf) {
            Ok(0) => Ok(None),
            Ok(_) => {
                if self.buf.ends_with('\n') {
                    self.buf.pop();
                }
                if self.buf.ends_with('\r') {
                    self.buf.pop();
                }
                Ok(Some((&self.buf, vpos)))
            }
            Err(e) => Err(BgzfError::Io(e)),
        }
    }

    pub fn virtual_position(&self) -> VirtualPosition {
        self.reader.virtual_position()
    }
}

impl<R: Read + io::Seek> BgzfLineReader<R> {
    pub fn seek(&mut self, pos: VirtualPosition) -> Result<()> {
        self.reader.seek(pos)?;
        Ok(())
    }
}
