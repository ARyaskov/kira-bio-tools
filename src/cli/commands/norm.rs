use anyhow::Result;

use crate::cli::args::NormArgs;
use crate::{VcfFormat, detect_format};

pub fn cmd_norm(args: NormArgs) -> Result<()> {
    let _ = &args.fasta_ref;

    let fmt = detect_format(&args.input)?;
    if fmt != VcfFormat::Plain {
        anyhow::bail!("Turbo mode currently supports only plain VCF");
    }

    let out_path = args.output.clone().unwrap_or_else(|| {
        let mut p = args.input.clone();
        p.set_extension("norm.vcf");
        p
    });

    crate::norm::turbo_norm_vcf(&args.input, &out_path)?;

    Ok(())
}
