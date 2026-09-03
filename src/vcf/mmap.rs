use memmap2::Mmap;
use rayon::prelude::*;

use crate::vcf::header::ContigDict;
use crate::vcf::structs::VcfRecord;

pub struct MmapVcfParser<'a> {
    data: &'a [u8],
}

impl<'a> MmapVcfParser<'a> {
    pub fn new(mmap: &'a Mmap) -> Self {
        Self { data: mmap }
    }

    pub fn parse_parallel(&self, num_threads: usize) -> Vec<VcfRecord> {
        let chunk_size = self.data.len() / num_threads;
        let chunks: Vec<_> = (0..num_threads)
            .map(|i| {
                let start = i * chunk_size;
                let end = if i == num_threads - 1 {
                    self.data.len()
                } else {
                    (i + 1) * chunk_size
                };
                (start, end)
            })
            .collect();

        chunks
            .into_par_iter()
            .flat_map(|(start, end)| {
                let adjusted_start = if start == 0 {
                    0
                } else {
                    self.data[start..]
                        .iter()
                        .position(|&b| b == b'\n')
                        .map(|p| start + p + 1)
                        .unwrap_or(end)
                };

                let adjusted_end = if end >= self.data.len() {
                    self.data.len()
                } else {
                    self.data[..end]
                        .iter()
                        .rposition(|&b| b == b'\n')
                        .map(|p| p + 1)
                        .unwrap_or(end)
                };

                self.parse_chunk(adjusted_start, adjusted_end)
            })
            .collect()
    }

    fn parse_chunk(&self, start: usize, end: usize) -> Vec<VcfRecord> {
        let mut records = Vec::new();
        let mut pos = start;
        // Contig ids are local to the chunk; callers needing global ids resolve by name.
        let mut contigs = ContigDict::new();

        while pos < end {
            let line_end = self.data[pos..end]
                .iter()
                .position(|&b| b == b'\n')
                .map(|p| pos + p)
                .unwrap_or(end);

            let line = &self.data[pos..line_end];

            if !line.is_empty() && line[0] != b'#' {
                if let Ok(line_str) = std::str::from_utf8(line) {
                    if let Ok(Some(rec)) =
                        crate::vcf::unified_reader::parse_vcf_record(line_str, pos as u64, &mut contigs)
                    {
                        records.push(rec);
                    }
                }
            }

            pos = line_end + 1;
        }

        records
    }
}
