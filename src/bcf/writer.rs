use super::BCF_MAGIC;
use super::header::{BcfHeaderDict, parse_header_to_dict, serialize_header};
use super::record::encode_record;
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use crate::bgzf::{BgzfWriter, FILE_BUFFER_SIZE, STREAM_BUFFER_SIZE};

enum BcfOut {
    Bgzf(BgzfWriter),
    Plain(BufWriter<Box<dyn Write + Send>>),
}

impl Write for BcfOut {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            BcfOut::Bgzf(w) => w.write(buf),
            BcfOut::Plain(w) => w.write(buf),
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        match self {
            BcfOut::Bgzf(w) => w.flush(),
            BcfOut::Plain(w) => w.flush(),
        }
    }
}

pub struct BcfWriter {
    inner: BcfOut,
    pub dict: BcfHeaderDict,
}

impl BcfWriter {
    pub fn create<P: AsRef<Path>>(p: P, compressed: bool, level: u32, headers: &[String]) -> Result<Self> {
        let path = p.as_ref();
        let f = File::create(path).with_context(|| format!("create {:?}", path))?;
        Self::from_writer(Box::new(f), compressed, level, headers, false)
    }

    /// `-O u` (`compressed == false`) writes raw BCF; `-O b` wraps it in BGZF.
    /// `streaming` selects a small output buffer so pipe consumers see data early.
    pub fn from_writer(
        w: Box<dyn Write + Send>,
        compressed: bool,
        level: u32,
        headers: &[String],
        streaming: bool,
    ) -> Result<Self> {
        let inner = if compressed {
            let buf = if streaming { STREAM_BUFFER_SIZE } else { FILE_BUFFER_SIZE };
            BcfOut::Bgzf(BgzfWriter::from_writer_buffered(w, flate2::Compression::new(level.min(9)), buf)?)
        } else {
            BcfOut::Plain(BufWriter::with_capacity(1 << 20, w))
        };
        let mut out = Self { inner, dict: parse_header_to_dict(headers) };
        out.write_header()?;
        Ok(out)
    }

    fn write_header(&mut self) -> Result<()> {
        self.inner.write_all(BCF_MAGIC)?;
        let text = serialize_header(&self.dict);
        let mut text_bytes = text.into_bytes();
        text_bytes.push(0);
        let l_text = text_bytes.len() as u32;
        self.inner.write_all(&l_text.to_le_bytes())?;
        self.inner.write_all(&text_bytes)?;
        Ok(())
    }

    pub fn write_vcf_line(&mut self, line: &str) -> Result<()> {
        if line.is_empty() || line.as_bytes()[0] == b'#' { return Ok(()); }
        encode_record(&mut self.inner, line, &self.dict)
    }

    pub fn finish(self) -> Result<()> {
        match self.inner {
            BcfOut::Bgzf(w) => w.finish().context("finalize BGZF-BCF"),
            BcfOut::Plain(mut w) => w.flush().context("flush BCF"),
        }
    }
}

#[cfg(test)]
#[path = "../../tests/unit/bcf_writer.rs"]
mod tests;
