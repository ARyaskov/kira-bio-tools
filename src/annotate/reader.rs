use anyhow::Result;
use std::path::Path;

use crate::vcf::UnifiedVcfReader;

pub struct VcfAnnotationReader {
    inner: UnifiedVcfReader,
    indexed_names: Vec<String>,
}

impl VcfAnnotationReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let inner = UnifiedVcfReader::open(path.as_ref())?;
        let source = if path
            .as_ref()
            .extension()
            .is_some_and(|extension| extension == "ktile")
        {
            path.as_ref().with_extension("")
        } else {
            path.as_ref().to_path_buf()
        };
        let indexed_names = crate::csi::find_index_for(&source)
            .map(crate::csi::BinIndex::load)
            .transpose()?
            .map(|index| index.names().to_vec())
            .unwrap_or(inner.reference_sequences()?);
        Ok(Self {
            inner,
            indexed_names,
        })
    }

    pub fn read_line(&mut self) -> Result<Option<String>> {
        Ok(self.inner.read_line()?)
    }

    pub fn read_header(&mut self) -> Result<Vec<String>> {
        Ok(complete_contig_header(
            self.inner.header()?,
            &self.indexed_names,
        ))
    }
}

pub struct StreamingVcfReader {
    reader: UnifiedVcfReader,
    indexed_names: Vec<String>,
}

impl StreamingVcfReader {
    pub fn new(reader: VcfAnnotationReader) -> Self {
        Self {
            reader: reader.inner,
            indexed_names: reader.indexed_names,
        }
    }

    pub fn read_line(&mut self) -> Result<Option<String>> {
        Ok(self.reader.read_line()?)
    }

    /// Phase 3 fast-path: read a line with the source's pre-parsed
    /// (chr_id, pos) metadata if available. Returns `None` metadata for
    /// sources without it (BGZF, Plain). Consumer can skip its own
    /// fast_parse_min when meta is `Some`.
    pub fn read_line_with_meta(&mut self) -> Result<Option<(String, Option<(u32, u32)>)>> {
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
        let headers = complete_contig_header(self.reader.header()?, &self.indexed_names);
        Ok((headers, self))
    }
}

pub(crate) fn complete_contig_header(mut headers: Vec<String>, names: &[String]) -> Vec<String> {
    let known: std::collections::HashSet<String> = headers
        .iter()
        .filter_map(|header| crate::vcf::header::parse_contig_line(header).map(|(name, _)| name))
        .collect();
    let offset = headers
        .iter()
        .position(|header| header.starts_with("#CHROM"))
        .unwrap_or(headers.len());
    headers.splice(
        offset..offset,
        names
            .iter()
            .filter(|name| !known.contains(*name))
            .map(|name| format!("##contig=<ID={name}>")),
    );
    headers
}
