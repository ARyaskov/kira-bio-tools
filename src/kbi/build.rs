use std::path::Path;

use crate::kbi::bitmap::PosBitmap;
use crate::kbi::builder::KbiBuilder;
use crate::kbi::index::KbiIndex;
use crate::kbi::structs::Result;
use crate::vcf::VcfReader;

pub fn build_kbi_index<P: AsRef<Path>>(vcf_path: P, output_path: P) -> Result<KbiIndex> {
    build(vcf_path, output_path, None)
}

/// The `.kbi` plus the position bitmap of the same records, in one pass.
pub fn build_kbi_index_with_bitmap<P: AsRef<Path>>(
    vcf_path: P,
    output_path: P,
) -> Result<(KbiIndex, PosBitmap)> {
    let mut bitmap = PosBitmap::default();
    let index = build(vcf_path, output_path, Some(&mut bitmap))?;
    Ok((index, bitmap))
}

fn build<P: AsRef<Path>>(
    vcf_path: P,
    output_path: P,
    mut bitmap: Option<&mut PosBitmap>,
) -> Result<KbiIndex> {
    let mut reader = VcfReader::open(vcf_path)?;
    reader.header()?;

    let estimated_capacity = 10_000_000;
    let mut builder = KbiBuilder::with_capacity(estimated_capacity);

    let mut count = 0usize;
    for record in reader.records() {
        let record = record?;
        builder.add_record(&record);
        if let Some(bm) = bitmap.as_deref_mut() {
            bm.insert_span(&record.chrom, record.pos, record.ref_allele.len() as u32);
        }
        count += 1;

        if count.is_multiple_of(1_000_000) {
            eprintln!("Processed {} records...", count);
        }
    }
    builder.set_contigs(reader.contigs().clone());

    eprintln!("Building index for {} entries...", builder.len());
    let index = builder.build()?;

    eprintln!("Saving index to {:?}...", output_path.as_ref());
    index.save(&output_path)?;

    Ok(index)
}
