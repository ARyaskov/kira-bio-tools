use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use indexmap::IndexMap;
use noodles_bgzf::VirtualPosition as NoodlesVirtualPosition;
use noodles_csi as csi;
use noodles_csi::binning_index::index::reference_sequence::Bin;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;

use crate::csi::structs::Result;
use crate::csi::utils::reg2bin;
use crate::vcf::VcfReader;

const MIN_SHIFT: u8 = 14;
const DEPTH: u8 = 5;

pub fn build_csi_index<P: AsRef<Path>>(vcf_path: P, output_path: P) -> Result<()> {
    use noodles_csi::binning_index::index::header::{Format, Header};

    let mut reader = VcfReader::open_for_indexing(vcf_path.as_ref())?;
    let _header = reader.header()?;
    // Reference ids MUST follow the file's own @contig order (what htslib/bcftools expect), not a
    // global chromosome table.
    let contigs = reader.reference_sequences()?;
    let name_to_id: std::collections::HashMap<String, usize> =
        contigs.iter().enumerate().map(|(i, n)| (n.clone(), i)).collect();

    let mut ref_seqs: BTreeMap<usize, BTreeMap<usize, Vec<Chunk>>> = BTreeMap::new();
    // Per-reference (count, min start vpos, max end vpos) for the index metadata block (so that
    // `bcftools index -n` and stats report the right counts).
    let mut ref_meta: BTreeMap<usize, (u64, NoodlesVirtualPosition, NoodlesVirtualPosition)> =
        BTreeMap::new();
    // Per-reference linear index: 16kb window -> smallest record vpos in it. htslib uses this
    // `loffset` to pick the first BGZF block for a region query; without it, region queries skip
    // blocks and miss records.
    let mut ref_linear: BTreeMap<usize, BTreeMap<usize, u64>> = BTreeMap::new();
    let mut records_indexed = 0u64;

    loop {
        match reader.next_record_with_vpos()? {
            Some((record, current_vpos)) => {
                let ref_id = name_to_id
                    .get(&record.chrom)
                    .copied()
                    .unwrap_or_else(|| (record.chr_id.saturating_sub(1)) as usize);
                // CSI/tabix binning is 0-based half-open: a variant at 1-based POS occupies [POS-1, POS).
                let start = record.position.saturating_sub(1);
                let end = record.position;

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

                let m = ref_meta.entry(ref_id).or_insert((0, start_vp, end_vp));
                m.0 += 1;
                if start_vp < m.1 {
                    m.1 = start_vp;
                }
                if end_vp > m.2 {
                    m.2 = end_vp;
                }

                // First (lowest-vpos) record in each 16kb window sets its linear offset.
                let w = (start >> MIN_SHIFT) as usize;
                ref_linear
                    .entry(ref_id)
                    .or_default()
                    .entry(w)
                    .or_insert(record.offset);

                records_indexed += 1;
            }
            None => break,
        }
    }

    let max_ref_id = ref_seqs.keys().max().copied().unwrap_or(0);
    let n_ref = contigs.len().max(max_ref_id + 1);
    let mut reference_sequences = Vec::with_capacity(n_ref);

    for ref_id in 0..n_ref {
        let bins: IndexMap<usize, Bin> = if let Some(bins_map) = ref_seqs.get(&ref_id) {
            bins_map
                .iter()
                .map(|(&id, chunks)| (id, Bin::new(chunks.clone())))
                .collect()
        } else {
            IndexMap::new()
        };

        let linear_index: IndexMap<usize, NoodlesVirtualPosition> = {
            let mut out = IndexMap::new();
            if let Some(lin) = ref_linear.get(&ref_id) {
                if let Some(&max_w) = lin.keys().max() {
                    let mut carry = 0u64;
                    for w in 0..=max_w {
                        carry = lin.get(&w).copied().unwrap_or(carry);
                        out.insert(w, NoodlesVirtualPosition::try_from(carry).unwrap_or_default());
                    }
                }
            }
            out
        };
        let metadata = ref_meta.get(&ref_id).map(|&(count, start, end)| {
            csi::binning_index::index::reference_sequence::Metadata::new(start, end, count, 0)
        });
        reference_sequences.push(csi::binning_index::index::ReferenceSequence::new(
            bins,
            linear_index,
            metadata,
        ));
    }

    // Tabix/VCF header embedded in the CSI aux block — without this, htslib/bcftools report
    // "Invalid index header" and cannot use the index.
    let names: noodles_csi::binning_index::index::header::ReferenceSequenceNames =
        contigs.iter().map(|c| c.as_bytes().into()).collect();
    let header = Header::builder()
        .set_format(Format::Vcf)
        .set_reference_sequence_name_index(0)
        .set_start_position_index(1)
        .set_reference_sequence_names(names)
        .build();

    let index = csi::Index::builder()
        .set_min_shift(MIN_SHIFT)
        .set_depth(DEPTH)
        .set_header(header)
        .set_reference_sequences(reference_sequences)
        .build();

    let file = File::create(output_path)?;
    let mut writer = csi::io::Writer::new(BufWriter::new(file));
    writer
        .write_index(&index)
        .map_err(|e| crate::csi::structs::CsiError::Format(e.to_string()))?;

    eprintln!("CSI index created: {} records indexed", records_indexed);
    Ok(())
}
