use std::fs::File;
use std::io::{self, BufRead, Read};
use std::path::Path;

use noodles_bgzf as bgzf;

use crate::bgzf::structs::{Result, VirtualPosition};

pub struct BgzfReader<R: Read> {
    inner: bgzf::io::Reader<R>,
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
            inner: bgzf::io::Reader::new(reader),
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
        self.inner
            .seek(noodles_pos)
            .map_err(crate::bgzf::structs::BgzfError::Io)?;
        Ok(self.virtual_position())
    }
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
            Err(e) => Err(crate::bgzf::structs::BgzfError::Io(e)),
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
