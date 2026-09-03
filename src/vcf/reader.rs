use std::path::Path;

use crate::bgzf::VirtualPosition;
use crate::vcf::header::ContigDict;
use crate::vcf::structs::{Result, VcfRecord};
use crate::vcf::unified_reader::UnifiedVcfReader;

pub struct VcfReader {
    inner: UnifiedVcfReader,
}

impl VcfReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        Ok(Self {
            inner: UnifiedVcfReader::open(path)?,
        })
    }

    pub fn open_for_indexing<P: AsRef<Path>>(path: P) -> Result<Self> {
        Ok(Self {
            inner: UnifiedVcfReader::open_for_indexing(path)?,
        })
    }

    pub fn header(&mut self) -> Result<Vec<String>> {
        self.inner.header()
    }

    pub fn next_record(&mut self) -> Result<Option<VcfRecord>> {
        self.inner.read_record()
    }

    pub fn next_raw_line(&mut self) -> Result<Option<(String, u64)>> {
        match self.inner.read_record()? {
            Some(rec) => Ok(Some((format_vcf_record(&rec), rec.offset))),
            None => Ok(None),
        }
    }

    pub fn records(&mut self) -> RecordIterator<'_> {
        RecordIterator { reader: self }
    }

    pub fn reference_sequences(&self) -> Result<Vec<String>> {
        self.inner.reference_sequences()
    }

    pub fn contigs(&self) -> &ContigDict {
        self.inner.contigs()
    }

    pub fn virtual_position(&self) -> Option<VirtualPosition> {
        self.inner.virtual_position()
    }

    pub fn next_record_with_vpos(&mut self) -> Result<Option<(VcfRecord, VirtualPosition)>> {
        self.inner.next_record_with_vpos()
    }

    pub fn next_line_with_vpos(&mut self) -> Result<Option<(String, VirtualPosition, VirtualPosition)>> {
        self.inner.next_line_with_vpos()
    }

    pub fn into_inner(self) -> UnifiedVcfReader {
        self.inner
    }
}

pub struct RecordIterator<'a> {
    reader: &'a mut VcfReader,
}

impl<'a> Iterator for RecordIterator<'a> {
    type Item = Result<VcfRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.next_record() {
            Ok(Some(rec)) => Some(Ok(rec)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

pub fn format_vcf_record(rec: &VcfRecord) -> String {
    let mut line = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        rec.chrom, rec.pos, rec.id, rec.ref_allele, rec.alt, rec.qual, rec.filter, rec.info
    );

    if let Some(fmt) = &rec.format {
        line.push('\t');
        line.push_str(fmt);

        for sample in &rec.samples {
            line.push('\t');
            line.push_str(sample);
        }
    }

    line
}
