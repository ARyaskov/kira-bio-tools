use super::header::{BcfHeaderDict, parse_header_to_dict, serialize_header};
use super::record::encode_record;
use super::BCF_MAGIC;
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub struct BcfWriter {
    inner: Box<dyn Write>,
    pub dict: BcfHeaderDict,
}

impl BcfWriter {
    pub fn create<P: AsRef<Path>>(p: P, compressed: bool, level: u32, headers: &[String]) -> Result<Self> {
        let path = p.as_ref();
        let inner: Box<dyn Write> = if compressed {
            let w = crate::bgzf::BgzfWriter::with_compression(path, flate2::Compression::new(level))
                .with_context(|| format!("open bgzf BCF {:?}", path))?;
            Box::new(w)
        } else {
            Box::new(BufWriter::with_capacity(1 << 20, File::create(path)?))
        };
        let mut w = Self { inner, dict: parse_header_to_dict(headers) };
        w.write_header()?;
        Ok(w)
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
        let mut inner = self.inner;
        inner.flush()?;
        Ok(())
    }
}

#[cfg(test)]
#[path = "../../tests/unit/bcf_writer.rs"]
mod tests;
