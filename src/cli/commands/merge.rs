use anyhow::Result;

use crate::VcfReader;
use crate::cli::args::MergeArgs;

pub fn cmd_merge(args: MergeArgs) -> Result<()> {
    let _ = &args.bcftools_args;
    if args.inputs.is_empty() {
        return Ok(());
    }

    let mut first = VcfReader::open(&args.inputs[0])?;
    let headers = first.header()?;
    for h in &headers {
        println!("{h}");
    }
    while let Some(rec) = first.next_record()? {
        print_record(&rec);
    }

    for path in args.inputs.iter().skip(1) {
        let mut r = VcfReader::open(path)?;
        let _ = r.header()?;
        while let Some(rec) = r.next_record()? {
            print_record(&rec);
        }
    }

    Ok(())
}

fn print_record(rec: &crate::vcf::structs::VcfRecord) {
    print!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        rec.chrom, rec.pos, rec.id, rec.ref_allele, rec.alt, rec.qual, rec.filter, rec.info
    );
    if let Some(fmt) = &rec.format {
        print!("\t{fmt}");
        for s in &rec.samples {
            print!("\t{s}");
        }
    }
    println!();
}
