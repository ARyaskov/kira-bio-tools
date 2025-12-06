use std::path::Path;

use crate::bgzf::VirtualPosition;
use crate::util::parse_vcf_line_fast;
use crate::vcf::structs::{Result, VcfRecord};
use crate::vcf::unified_reader::UnifiedVcfReader;

pub struct VcfReader {
    inner: UnifiedVcfReader,
    offset: u64,
}

impl VcfReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        Ok(Self {
            inner: UnifiedVcfReader::open(path)?,
            offset: 0,
        })
    }

    pub fn open_for_indexing<P: AsRef<Path>>(path: P) -> Result<Self> {
        Ok(Self {
            inner: UnifiedVcfReader::open_for_indexing(path)?,
            offset: 0,
        })
    }

    pub fn header(&mut self) -> Result<Vec<String>> {
        let headers = self.inner.header()?;
        for h in &headers {
            self.offset += (h.len() + 1) as u64;
        }
        Ok(headers)
    }

    pub fn next_record(&mut self) -> Result<Option<VcfRecord>> {
        let start_offset = self.offset;

        match self.inner.read_line()? {
            Some(line) => {
                self.offset += (line.len() + 1) as u64;

                if let Some((chr_id, position)) = parse_vcf_line_fast(line.as_bytes()) {
                    Ok(Some(VcfRecord {
                        chr_id,
                        position,
                        offset: start_offset,
                    }))
                } else {
                    self.next_record()
                }
            }
            None => Ok(None),
        }
    }

    pub fn next_raw_line(&mut self) -> Result<Option<(String, u64)>> {
        let start_offset = self.offset;

        match self.inner.read_line()? {
            Some(line) => {
                self.offset += (line.len() + 1) as u64;
                Ok(Some((line, start_offset)))
            }
            None => Ok(None),
        }
    }

    pub fn records(&mut self) -> RecordIterator<'_> {
        RecordIterator { reader: self }
    }

    pub fn reference_sequences(&self) -> &[String] {
        self.inner.reference_sequences()
    }

    pub fn virtual_position(&self) -> Option<VirtualPosition> {
        self.inner.virtual_position()
    }

    pub fn next_record_with_vpos(&mut self) -> Result<Option<(VcfRecord, VirtualPosition)>> {
        self.inner.next_record_with_vpos()
    }
}

pub struct RecordIterator<'a> {
    reader: &'a mut VcfReader,
}

impl<'a> Iterator for RecordIterator<'a> {
    type Item = Result<VcfRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.next_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

pub use VcfReader as PlainVcfReader;
pub use VcfReader as GzipVcfReader;
pub use VcfReader as BgzfVcfReader;
