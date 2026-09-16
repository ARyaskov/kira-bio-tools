use crate::cli::args::ZarrArgs;
use crate::zarr::{ZarrOpts, export_vcf_zarr};
use anyhow::Result;

pub fn cmd_zarr(args: ZarrArgs) -> Result<()> {
    let opts = ZarrOpts {
        chunk_variants: args.chunk_variants,
        chunk_samples: args.chunk_samples,
        level: args.level,
        fields: args.fields.as_deref().map(|f| {
            f.split(',')
                .map(|t| t.trim().to_string())
                .filter(|t| !t.is_empty())
                .collect()
        }),
    };
    let start = std::time::Instant::now();
    let s = export_vcf_zarr(&args.input, &args.output, &opts, args.force)?;
    eprintln!(
        "[KIRA_BT] zarr: {} variants, {} samples, {} arrays -> {} ({:.2}s)",
        s.variants,
        s.samples,
        s.arrays.len(),
        args.output.display(),
        start.elapsed().as_secs_f64()
    );
    Ok(())
}
