use std::fs::File;
use std::io::{self, BufWriter};
use std::path::Path;

use indexmap::IndexMap;
use noodles_bgzf::VirtualPosition as NoodlesVirtualPosition;
use noodles_csi as csi;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles_csi::binning_index::index::reference_sequence::Bin;
use thiserror::Error;

use crate::bgzf::VirtualPosition as OurVirtualPosition;
use crate::vcf::{BgzfVcfReader, VcfError};

#[derive(Debug, Error)]
pub enum CsiError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
    #[error("VCF error: {0}")]
    Vcf(#[from] VcfError),
    #[error("CSI format error: {0}")]
    Format(String),
}

pub type Result<T> = std::result::Result<T, CsiError>;

const MIN_SHIFT: u8 = 14;
const DEPTH: u8 = 5;

pub fn build_csi_index<P: AsRef<Path>>(vcf_path: P, output_path: P) -> Result<()> {
    use std::collections::BTreeMap;

    let mut reader = BgzfVcfReader::open(vcf_path.as_ref())?;
    let _header = reader.header()?;

    let mut ref_seqs: BTreeMap<usize, BTreeMap<usize, Vec<Chunk>>> = BTreeMap::new();
    let mut records_indexed = 0u64;

    loop {
        let current_vpos: OurVirtualPosition = reader.virtual_position();

        match reader.next_record()? {
            Some(record) => {
                let ref_id = (record.chr_id.saturating_sub(1)) as usize;
                let start = record.position;
                let end = record.position + 1;

                let bin_id = reg2bin(start as usize, end as usize, MIN_SHIFT, DEPTH);

                let start_vp = NoodlesVirtualPosition::try_from(record.offset).unwrap_or_default();
                let end_vp = NoodlesVirtualPosition::try_from(current_vpos.as_u64()).unwrap_or_default();

                let chunk = Chunk::new(start_vp, end_vp);

                ref_seqs
                    .entry(ref_id)
                    .or_default()
                    .entry(bin_id)
                    .or_default()
                    .push(chunk);

                records_indexed += 1;
            }
            None => break,
        }
    }

    let max_ref_id = ref_seqs.keys().max().copied().unwrap_or(0);
    let mut reference_sequences = Vec::with_capacity(max_ref_id + 1);

    for ref_id in 0..=max_ref_id {
        let bins: IndexMap<usize, Bin> = if let Some(bins_map) = ref_seqs.get(&ref_id) {
            bins_map
                .iter()
                .map(|(&id, chunks)| (id, Bin::new(chunks.clone())))
                .collect()
        } else {
            IndexMap::new()
        };

        let linear_index: IndexMap<usize, NoodlesVirtualPosition> = IndexMap::new();
        reference_sequences.push(csi::binning_index::index::ReferenceSequence::new(
            bins,
            linear_index,
            None,
        ));
    }

    let index = csi::Index::builder()
        .set_min_shift(MIN_SHIFT)
        .set_depth(DEPTH)
        .set_reference_sequences(reference_sequences)
        .build();

    let file = File::create(output_path)?;
    let mut writer = csi::Writer::new(BufWriter::new(file));
    writer
        .write_index(&index)
        .map_err(|e| CsiError::Format(e.to_string()))?;

    eprintln!("CSI index created: {} records indexed", records_indexed);
    Ok(())
}

fn reg2bin(start: usize, end: usize, min_shift: u8, depth: u8) -> usize {
    let end = end.saturating_sub(1);
    let mut s = min_shift as usize;
    let mut t = ((1 << (depth as usize * 3)) - 1) / 7;

    for _ in 0..depth {
        if start >> s == end >> s {
            return t + (start >> s);
        }
        s += 3;
        t = (t - 1) / 8;
    }

    0
}

pub fn read_csi_index<P: AsRef<Path>>(path: P) -> Result<csi::Index> {
    let file = File::open(path)?;
    let mut reader = csi::Reader::new(file);
    reader
        .read_index()
        .map_err(|e| CsiError::Format(e.to_string()))
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