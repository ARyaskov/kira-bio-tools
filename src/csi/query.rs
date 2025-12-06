use std::fs::File;
use std::path::Path;

use noodles_csi as csi;

use crate::csi::structs::Result;

pub fn read_csi_index<P: AsRef<Path>>(path: P) -> Result<csi::Index> {
    let file = File::open(path)?;
    let mut reader = csi::Reader::new(file);
    reader
        .read_index()
        .map_err(|e| crate::csi::structs::CsiError::Format(e.to_string()))
}

pub struct CsiQuery {
    index: csi::Index,
}

impl CsiQuery {
    pub fn new(index: csi::Index) -> Self {
        Self { index }
    }

    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let index = read_csi_index(path)?;
        Ok(Self::new(index))
    }

    pub fn query(&self, ref_id: usize, _start: u32, _end: u32) -> Vec<(u64, u64)> {
        let ref_seqs = self.index.reference_sequences();

        if ref_id >= ref_seqs.len() {
            return Vec::new();
        }

        let ref_seq = &ref_seqs[ref_id];
        let mut chunks = Vec::new();

        for (_bin_id, bin) in ref_seq.bins() {
            for chunk in bin.chunks() {
                let chunk_start = u64::from(chunk.start());
                let chunk_end = u64::from(chunk.end());
                chunks.push((chunk_start, chunk_end));
            }
        }

        chunks.sort_by_key(|c| c.0);
        chunks
    }

    pub fn num_reference_sequences(&self) -> usize {
        self.index.reference_sequences().len()
    }
}
