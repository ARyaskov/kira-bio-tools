use anyhow::Result;
use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::cli::args::TabixArgs;
use crate::{
    build_csi_index, build_kbi_index, chr_id_to_name, chr_name_to_id, detect_format, fetch_line,
    CsiQuery, KbiIndex, Region, VcfFormat, VcfReader,
};

pub fn cmd_tabix(args: TabixArgs) -> Result<()> {
    if args.list_chroms {
        return cmd_list_tabix(&args);
    }

    if args.only_header {
        return cmd_header_tabix(&args);
    }

    if args.regions.is_empty() && args.regions_file.is_none() && args.targets_file.is_none() {
        return cmd_index_tabix(args);
    }

    cmd_query_tabix(args)
}

fn cmd_index_tabix(args: TabixArgs) -> Result<()> {
    let format = detect_format(&args.input)?;

    if !args.input.exists() {
        anyhow::bail!("Input file does not exist: {:?}", args.input);
    }

    match format {
        VcfFormat::Bgzf => {
            let index_path = if args.csi {
                let mut p = args.input.clone();
                let name = p.file_name().unwrap().to_string_lossy().to_string();
                p.set_file_name(format!("{}.csi", name));
                p
            } else {
                let mut p = args.input.clone();
                let name = p.file_name().unwrap().to_string_lossy().to_string();
                p.set_file_name(format!("{}.tbi", name));
                p
            };

            if index_path.exists() && !args.force {
                anyhow::bail!(
                    "Index file already exists: {:?}. Use -f to overwrite.",
                    index_path
                );
            }

            eprintln!("Building index: {:?}", index_path);
            build_csi_index(&args.input, &index_path)?;

            let kbi_path = args.input.with_extension("kbi");
            if !kbi_path.exists() {
                eprintln!("Building KBI index: {:?}", kbi_path);
                build_kbi_index(&args.input, &kbi_path)?;
            }
        }
        VcfFormat::Plain | VcfFormat::Gzip => {
            anyhow::bail!(
                "File must be BGZF-compressed. Use bgzip to compress: bgzip {:?}",
                args.input
            );
        }
    }

    Ok(())
}

fn cmd_query_tabix(args: TabixArgs) -> Result<()> {
    let format = detect_format(&args.input)?;

    if format != VcfFormat::Bgzf {
        anyhow::bail!("Query requires BGZF-compressed file: {:?}", args.input);
    }

    let csi_path = {
        let mut p = args.input.clone();
        let name = p.file_name().unwrap().to_string_lossy().to_string();
        p.set_file_name(format!("{}.csi", name));
        p
    };

    let tbi_path = {
        let mut p = args.input.clone();
        let name = p.file_name().unwrap().to_string_lossy().to_string();
        p.set_file_name(format!("{}.tbi", name));
        p
    };

    let kbi_path = args.input.with_extension("kbi");

    let index_path = if csi_path.exists() {
        csi_path
    } else if tbi_path.exists() {
        tbi_path
    } else if kbi_path.exists() {
        kbi_path.clone()
    } else {
        anyhow::bail!(
            "No index found for {:?}. Run 'kira-bt tabix {:?}' first.",
            args.input,
            args.input
        );
    };

    let use_kbi = index_path.extension().map(|e| e == "kbi").unwrap_or(false);

    let mut regions = args.regions.clone();

    if let Some(ref regions_file) = args.regions_file {
        let file = File::open(regions_file)?;
        for line in BufReader::new(file).lines() {
            let line = line?;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            regions.push(line);
        }
    }

    if let Some(ref targets_file) = args.targets_file {
        let file = File::open(targets_file)?;
        for line in BufReader::new(file).lines() {
            let line = line?;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            regions.push(line);
        }
    }

    if regions.is_empty() {
        anyhow::bail!("No regions specified");
    }

    if let Some(ref reheader_file) = args.reheader {
        let reheader_content = std::fs::read_to_string(reheader_file)?;
        print!("{}", reheader_content);
    } else if args.print_header {
        print_vcf_header(&args.input)?;
    }

    if use_kbi {
        query_with_kbi(&args, &regions)?;
    } else {
        query_with_csi(&args, &regions, &index_path)?;
    }

    Ok(())
}

fn query_with_kbi(args: &TabixArgs, regions: &[String]) -> Result<()> {
    let kbi_path = args.input.with_extension("kbi");
    let index = KbiIndex::load(&kbi_path)?;

    for (_region_idx, region_str) in regions.iter().enumerate() {
        if let Some(overlap_mode) = args.regions_overlap {
            match overlap_mode {
                1 => println!("#{}", region_str),
                2 => {
                    print!("{}\t", region_str);
                }
                _ => {}
            }
        }

        let region = Region::parse(region_str)
            .ok_or_else(|| anyhow::anyhow!("Invalid region: {}", region_str))?;

        let chr_id = chr_name_to_id(&region.chr)
            .ok_or_else(|| anyhow::anyhow!("Unknown chromosome: {}", region.chr))?;

        let (start, end) = if args.zero_based {
            let s = region.start.unwrap_or(0);
            let e = region.end.unwrap_or(u32::MAX);
            (s, e.saturating_sub(1))
        } else {
            let s = region.start.unwrap_or(1);
            let e = region.end.unwrap_or(u32::MAX);
            (s, e)
        };

        let results = index.range(chr_id, start, end);

        for (_pos, offset) in results {
            let line = fetch_line(&args.input, offset)?;
            println!("{}", line);
        }
    }

    Ok(())
}

fn query_with_csi(
    args: &TabixArgs,
    regions: &[String],
    index_path: &std::path::PathBuf,
) -> Result<()> {
    let csi = CsiQuery::open(index_path)?;

    for region_str in regions {
        if let Some(overlap_mode) = args.regions_overlap {
            match overlap_mode {
                1 => println!("#{}", region_str),
                2 => {
                    print!("{}\t", region_str);
                }
                _ => {}
            }
        }

        let region = Region::parse(region_str)
            .ok_or_else(|| anyhow::anyhow!("Invalid region: {}", region_str))?;

        let chr_id = chr_name_to_id(&region.chr)
            .ok_or_else(|| anyhow::anyhow!("Unknown chromosome: {}", region.chr))?;

        let (start, end) = if args.zero_based {
            let s = region.start.unwrap_or(0);
            let e = region.end.unwrap_or(u32::MAX);
            (s, e.saturating_sub(1))
        } else {
            let s = region.start.unwrap_or(1);
            let e = region.end.unwrap_or(u32::MAX);
            (s, e)
        };

        let chunks = csi.query((chr_id - 1) as usize, start, end);

        for (chunk_start, _chunk_end) in chunks {
            let line = fetch_line(&args.input, chunk_start)?;
            println!("{}", line);
        }
    }

    Ok(())
}

fn cmd_list_tabix(args: &TabixArgs) -> Result<()> {
    let kbi_path = args.input.with_extension("kbi");

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
        let mut reader = VcfReader::open(&args.input)?;
        let _ = reader.header()?;

        for name in reader.reference_sequences()? {
            println!("{}", name);
        }
    }

    Ok(())
}

fn cmd_header_tabix(args: &TabixArgs) -> Result<()> {
    print_vcf_header(&args.input)
}

fn print_vcf_header(path: &std::path::Path) -> Result<()> {
    let mut reader = VcfReader::open(path)?;
    let headers = reader.header()?;

    for line in headers {
        println!("{}", line);
    }

    Ok(())
}
