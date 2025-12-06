use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use indexmap::IndexMap;
use noodles_bgzf::VirtualPosition as NoodlesVirtualPosition;
use noodles_csi as csi;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles_csi::binning_index::index::reference_sequence::Bin;

use crate::csi::structs::Result;
use crate::csi::utils::reg2bin;
use crate::vcf::VcfReader;

const MIN_SHIFT: u8 = 14;
const DEPTH: u8 = 5;

pub fn build_csi_index<P: AsRef<Path>>(vcf_path: P, output_path: P) -> Result<()> {
    let mut reader = VcfReader::open_for_indexing(vcf_path.as_ref())?;
    let _header = reader.header()?;

    let mut ref_seqs: BTreeMap<usize, BTreeMap<usize, Vec<Chunk>>> = BTreeMap::new();
    let mut records_indexed = 0u64;

    loop {
        match reader.next_record_with_vpos()? {
            Some((record, current_vpos)) => {
                let ref_id = (record.chr_id.saturating_sub(1)) as usize;
                let start = record.position;
                let end = record.position + 1;

                let bin_id = reg2bin(start as usize, end as usize, MIN_SHIFT, DEPTH);

                let start_vp = NoodlesVirtualPosition::try_from(record.offset).unwrap_or_default();
                let end_vp =
                    NoodlesVirtualPosition::try_from(current_vpos.as_u64()).unwrap_or_default();

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
        .map_err(|e| crate::csi::structs::CsiError::Format(e.to_string()))?;

    eprintln!("CSI index created: {} records indexed", records_indexed);
    Ok(())
}
