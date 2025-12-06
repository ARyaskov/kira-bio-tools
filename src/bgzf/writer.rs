use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use noodles_bgzf as bgzf;

use crate::bgzf::structs::Result;

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
        self.inner
            .finish()
            .map_err(crate::bgzf::structs::BgzfError::Io)
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
