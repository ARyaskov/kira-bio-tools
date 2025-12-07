use anyhow::Result;
use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::cli::args::{HeaderArgs, ListArgs, QueryArgs, StatArgs};
use crate::{chr_id_to_name, chr_name_to_id, fetch_line, CsiQuery, KbiIndex, Region, VcfReader};

pub fn cmd_query(args: QueryArgs) -> Result<()> {
    if args.only_header {
        return cmd_header(HeaderArgs { file: args.file });
    }

    let kbi_path = args.file.with_extension("kbi");
    let csi_path = {
        let mut p = args.file.clone();
        let name = p.file_name().unwrap().to_string_lossy().to_string();
        p.set_file_name(format!("{}.csi", name));
        p
    };

    let use_kbi = kbi_path.exists();
    let use_csi = csi_path.exists() && !use_kbi;

    if !use_kbi && !use_csi {
        anyhow::bail!("No index found. Run 'kira-bt index {:?}' first.", args.file);
    }

    let mut regions = args.regions.clone();
    if let Some(ref regions_file) = args.regions_file {
        let file = File::open(regions_file)?;
        for line in BufReader::new(file).lines() {
            regions.push(line?);
        }
    }

    if regions.is_empty() {
        anyhow::bail!("No regions specified");
    }

    if args.print_header {
        print_vcf_header(&args.file)?;
    }

    let mut total_count = 0usize;

    if use_kbi {
        let index = KbiIndex::load(&kbi_path)?;

        for region_str in &regions {
            let region = Region::parse(region_str)
                .ok_or_else(|| anyhow::anyhow!("Invalid region: {}", region_str))?;

            let chr_id = chr_name_to_id(&region.chr)
                .ok_or_else(|| anyhow::anyhow!("Unknown chromosome: {}", region.chr))?;

            let start = region.start.unwrap_or(0);
            let end = region.end.unwrap_or(u32::MAX);

            let results = index.range(chr_id, start, end);
            total_count += results.len();

            if !args.count {
                for (_pos, offset) in results {
                    let line = fetch_line(&args.file, offset)?;
                    println!("{}", line);
                }
            }
        }
    } else if use_csi {
        let csi = CsiQuery::open(&csi_path)?;

        for region_str in &regions {
            let region = Region::parse(region_str)
                .ok_or_else(|| anyhow::anyhow!("Invalid region: {}", region_str))?;

            let chr_id = chr_name_to_id(&region.chr)
                .ok_or_else(|| anyhow::anyhow!("Unknown chromosome: {}", region.chr))?;

            let start = region.start.unwrap_or(0);
            let end = region.end.unwrap_or(u32::MAX);

            let chunks = csi.query((chr_id - 1) as usize, start, end);

            for (chunk_start, _chunk_end) in chunks {
                let line = fetch_line(&args.file, chunk_start)?;
                if !args.count {
                    println!("{}", line);
                }
                total_count += 1;
            }
        }
    }

    if args.count {
        println!("{}", total_count);
    }

    Ok(())
}

pub fn cmd_stat(args: StatArgs) -> Result<()> {
    use std::fs;

    let file_size = fs::metadata(&args.index)?.len();

    if args.index.extension().map(|e| e == "kbi").unwrap_or(false) {
        let index = KbiIndex::load(&args.index)?;

        println!("Index Statistics (KBI)");
        println!("======================");
        println!("File:          {:?}", args.index);
        println!(
            "File size:     {} bytes ({:.2} MB)",
            file_size,
            file_size as f64 / 1024.0 / 1024.0
        );
        println!("Entries:       {}", index.len());
        println!(
            "Memory usage:  {} bytes ({:.2} MB)",
            index.memory_usage(),
            index.memory_usage() as f64 / 1024.0 / 1024.0
        );
        println!("Bytes/key:     {:.2}", index.bytes_per_key());
    } else if args
        .index
        .extension()
        .map(|e| e == "csi" || e == "tbi")
        .unwrap_or(false)
    {
        let _csi = CsiQuery::open(&args.index)?;

        println!("Index Statistics (CSI/TBI)");
        println!("==========================");
        println!("File:          {:?}", args.index);
        println!(
            "File size:     {} bytes ({:.2} MB)",
            file_size,
            file_size as f64 / 1024.0 / 1024.0
        );
        println!("Format:        CSI/TBI (tabix-compatible)");
    } else {
        anyhow::bail!("Unknown index format");
    }

    Ok(())
}

pub fn cmd_list(args: ListArgs) -> Result<()> {
    let kbi_path = args.file.with_extension("kbi");

    if kbi_path.exists() {
        let index = KbiIndex::load(&kbi_path)?;

        for chr_id in 1..=25u8 {
            if let Some(name) = chr_id_to_name(chr_id) {
                let results = index.range(chr_id, 0, u32::MAX);
                if !results.is_empty() {
                    println!("{}", name);
                }
            }
        }
    } else {
        let mut reader = VcfReader::open(&args.file)?;
        let _ = reader.header()?;

        for name in reader.reference_sequences()? {
            println!("{}", name);
        }
    }

    Ok(())
}

pub fn cmd_header(args: HeaderArgs) -> Result<()> {
    print_vcf_header(&args.file)
}

fn print_vcf_header(path: &std::path::Path) -> Result<()> {
    let mut reader = VcfReader::open(path)?;
    let headers = reader.header()?;

    for line in headers {
        println!("{}", line);
    }

    Ok(())
}
