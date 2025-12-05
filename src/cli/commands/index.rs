use anyhow::Result;
use std::time::Instant;

use crate::cli::args::IndexArgs;
use crate::{build_csi_index, build_kbi_index, detect_format, VcfFormat};

pub fn cmd_index(args: IndexArgs) -> Result<()> {
    let format = detect_format(&args.input)?;

    eprintln!("Input: {:?}", args.input);
    eprintln!("Format: {:?}", format);

    match format {
        VcfFormat::Bgzf => {
            let csi_path = args.output.clone().unwrap_or_else(|| {
                let mut p = args.input.clone();
                p.set_extension("vcf.gz.csi");
                p
            });

            if !args.no_kbi || args.csi {
                eprintln!("Building CSI index: {:?}", csi_path);
                let csi_start = Instant::now();
                build_csi_index(&args.input, &csi_path)?;
                eprintln!("CSI build time: {:.3}s", csi_start.elapsed().as_secs_f64());
            }

            if !args.no_kbi {
                let kbi_path = args.input.with_extension("kbi");
                eprintln!("Building KBI index: {:?}", kbi_path);
                let kbi_start = Instant::now();
                let index = build_kbi_index(&args.input, &kbi_path)?;
                eprintln!("KBI build time: {:.3}s", kbi_start.elapsed().as_secs_f64());
                eprintln!("Entries: {}", index.len());
                eprintln!("Bytes/key: {:.2}", index.bytes_per_key());
            }
        }
        VcfFormat::Plain | VcfFormat::Gzip => {
            let kbi_path = args
                .output
                .unwrap_or_else(|| args.input.with_extension("kbi"));

            eprintln!("Building KBI index: {:?}", kbi_path);
            let kbi_start = Instant::now();
            let index = build_kbi_index(&args.input, &kbi_path)?;
            eprintln!("KBI build time: {:.3}s", kbi_start.elapsed().as_secs_f64());
            eprintln!("Entries: {}", index.len());
            eprintln!("Bytes/key: {:.2}", index.bytes_per_key());

            if args.csi {
                eprintln!("Warning: CSI index requires BGZF compression. Use bgzip first.");
            }
        }
    }

    Ok(())
}
