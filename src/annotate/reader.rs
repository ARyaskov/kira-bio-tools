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

pub struct StreamingVcfReader {
    reader: UnifiedVcfReader,
}

impl StreamingVcfReader {
    pub fn new(reader: VcfAnnotationReader) -> Self {
        Self {
            reader: reader.inner,
        }
    }

    pub fn read_line(&mut self) -> Result<Option<String>> {
        Ok(self.reader.read_line()?)
    }

    /// Phase 3 fast-path: read a line with the source's pre-parsed
    /// (chr_id, pos) metadata if available. Returns `None` metadata for
    /// sources without it (BGZF, Plain). Consumer can skip its own
    /// fast_parse_min when meta is `Some`.
    pub fn read_line_with_meta(
        &mut self,
    ) -> Result<Option<(String, Option<(u32, u32)>)>> {
        Ok(self.reader.read_line_with_meta()?)
    }

    /// Zero-copy fast path: appends the next line into `batch`
    /// without `String` intermediate. Returns `true` on success,
    /// `false` at EOF.
    pub fn read_line_into_batch(
        &mut self,
        batch: &mut crate::annotate::cpu_v2::ReadBatch,
    ) -> Result<bool> {
        Ok(self.reader.read_line_into_batch(batch)?)
    }

    pub fn into_headers_and_self(self) -> Result<(Vec<String>, Self)> {
        let headers = self.reader.header()?;
        Ok((headers, self))
    }
}
