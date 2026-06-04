use std::time::Instant;

use anyhow::Result;

use crate::annotate::ktile::write_ktile_from_vcf;
use crate::cli::args::{TileAction, TileArgs};

pub fn cmd_tile(args: TileArgs) -> Result<()> {
    match args.action {
        TileAction::Build(b) => {
            let output = b.output.clone().unwrap_or_else(|| {
                let mut p = b.input.clone();
                let mut name = p
                    .file_name()
                    .and_then(|n| n.to_str())
                    .unwrap_or("input")
                    .to_string();
                // Append `.ktile` regardless of existing extension so
                // `foo.vcf.gz` → `foo.vcf.gz.ktile`. This makes sidecar
                // detection a simple `<input>.ktile` exists check later.
                name.push_str(".ktile");
                p.set_file_name(name);
                p
            });

            let start = Instant::now();
            let stats = write_ktile_from_vcf(&b.input, &output)?;
            let elapsed = start.elapsed().as_secs_f64();
            let mb = stats.bytes_written as f64 / (1024.0 * 1024.0);

            eprintln!(
                "[tile] Built {} from {} records in {:.2}s ({:.1} MB, {:.1} MB/s)",
                output.display(),
                stats.n_records,
                elapsed,
                mb,
                mb / elapsed.max(0.001),
            );
            Ok(())
        }
    }
}
