use memmap2::Mmap;
use rayon::prelude::*;

use crate::util::chr_name_to_id;
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

        while pos < end {
            let line_end = self.data[pos..end]
                .iter()
                .position(|&b| b == b'\n')
                .map(|p| pos + p)
                .unwrap_or(end);

            let line = &self.data[pos..line_end];

            if !line.is_empty() && line[0] != b'#' {
                if let Ok(line_str) = std::str::from_utf8(line) {
                    if let Some(rec) = parse_full_vcf_record(line_str, pos as u64) {
                        records.push(rec);
                    }
                }
            }

            pos = line_end + 1;
        }

        records
    }
}

fn parse_full_vcf_record(line: &str, offset: u64) -> Option<VcfRecord> {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 8 {
        return None;
    }

    let chrom = cols[0];
    let pos = cols[1].parse::<u32>().ok()?;
    let chr_id = chr_name_to_id(chrom).unwrap_or(0);

    let format = if cols.len() > 8 {
        Some(cols[8].to_string())
    } else {
        None
    };

    let samples = if cols.len() > 9 {
        cols[9..].iter().map(|s| s.to_string()).collect()
    } else {
        Vec::new()
    };

    Some(VcfRecord {
        chrom: chrom.to_string(),
        pos,
        id: cols[2].to_string(),
        ref_allele: cols[3].to_string(),
        alt: cols[4].to_string(),
        qual: cols[5].to_string(),
        filter: cols[6].to_string(),
        info: cols[7].to_string(),
        format,
        samples,
        chr_id,
        position: pos,
        offset,
    })
}
