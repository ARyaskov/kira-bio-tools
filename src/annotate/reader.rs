use anyhow::Result;
use std::path::Path;

use crate::vcf::UnifiedVcfReader;

pub struct VcfAnnotationReader {
    inner: UnifiedVcfReader,
}

impl VcfAnnotationReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        Ok(Self {
            inner: UnifiedVcfReader::open(path)?,
        })
    }

    pub fn read_line(&mut self) -> Result<Option<String>> {
        Ok(self.inner.read_line()?)
    }

    pub fn read_header(&mut self) -> Result<Vec<String>> {
        Ok(self.inner.header()?)
    }
}

pub struct BatchVcfReader {
    reader: UnifiedVcfReader,
    batch_size: usize,
}

impl BatchVcfReader {
    pub fn new(reader: VcfAnnotationReader, batch_size: usize) -> Self {
        Self {
            reader: reader.inner,
            batch_size,
        }
    }

    pub fn read_batch(&mut self) -> Result<Vec<String>> {
        Ok(self.reader.read_batch(self.batch_size)?)
    }

    pub fn into_headers_and_self(mut self) -> Result<(Vec<String>, Self)> {
        let headers = self.reader.header()?;
        Ok((headers, self))
    }
}
