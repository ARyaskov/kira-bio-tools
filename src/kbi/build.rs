use std::path::Path;

use crate::kbi::builder::KbiBuilder;
use crate::kbi::index::KbiIndex;
use crate::kbi::structs::Result;
use crate::vcf::VcfReader;

pub fn build_kbi_index<P: AsRef<Path>>(vcf_path: P, output_path: P) -> Result<KbiIndex> {
    let mut reader = VcfReader::open(vcf_path)?;
    reader.header()?;

    let estimated_capacity = 10_000_000;
    let mut builder = KbiBuilder::with_capacity(estimated_capacity);

    let mut count = 0usize;
    for record in reader.records() {
        let record = record?;
        builder.add_record(&record);
        count += 1;

        if count % 1_000_000 == 0 {
            eprintln!("Processed {} records...", count);
        }
    }

    eprintln!("Building index for {} entries...", builder.len());
    let index = builder.build()?;

    eprintln!("Saving index to {:?}...", output_path.as_ref());
    index.save(&output_path)?;

    Ok(index)
}
