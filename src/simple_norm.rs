use anyhow::Result;
use memmap2::{Mmap, MmapMut, MmapOptions};
use rayon::prelude::*;
use std::fs::{File, OpenOptions};
use std::io::Write;
use std::time::Instant;

use crate::norm::NormContext;

/// Fast VCF normalization with SIMD and mmap write
pub fn turbo_norm_vcf(input: &std::path::Path, output: &std::path::Path) -> Result<()> {
    let total_start = Instant::now();

    // Detect CPU features once
    let norm_ctx = NormContext::detect();

    eprintln!("[CPU] Features: {}", norm_ctx.features());

    // Load input file
    let load_start = Instant::now();
    let file = File::open(input)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let data = &mmap[..];

    eprintln!(
        "[Timing] File load: {:.3}s",
        load_start.elapsed().as_secs_f64()
    );

    // Parse line positions
    let parse_start = Instant::now();
    let line_positions: Vec<usize> = data
        .par_iter()
        .enumerate()
        .filter_map(|(i, &byte)| if byte == b'\n' { Some(i) } else { None })
        .collect();

    eprintln!(
        "[Timing] Parse: {:.3}s | Lines: {}",
        parse_start.elapsed().as_secs_f64(),
        line_positions.len()
    );

    // Process lines in parallel with SIMD normalization
    let process_start = Instant::now();

    let normalized: Vec<Vec<u8>> = line_positions
        .par_iter()
        .enumerate()
        .map(|(i, &line_end)| {
            let line_start = if i == 0 { 0 } else { line_positions[i - 1] + 1 };
            let line = &data[line_start..line_end];

            // Skip empty or header lines
            if line.is_empty() || line[0] == b'#' {
                return line.to_vec();
            }

            // Find REF and ALT columns (tab-separated)
            let mut tabs = 0;
            let mut ref_start = 0;
            let mut ref_end = 0;
            let mut alt_start = 0;
            let mut alt_end = 0;

            for (idx, &byte) in line.iter().enumerate() {
                if byte == b'\t' {
                    tabs += 1;
                    match tabs {
                        3 => ref_start = idx + 1,
                        4 => {
                            ref_end = idx;
                            alt_start = idx + 1;
                        }
                        5 => {
                            alt_end = idx;
                            break;
                        }
                        _ => {}
                    }
                }
            }

            if tabs < 4 {
                return line.to_vec();
            }

            if alt_end == 0 {
                alt_end = line.len();
            }

            // Normalize with SIMD
            let ref_allele = unsafe { std::str::from_utf8_unchecked(&line[ref_start..ref_end]) };
            let alt_allele = unsafe { std::str::from_utf8_unchecked(&line[alt_start..alt_end]) };

            let (norm_ref, norm_alt, _, _) = norm_ctx.normalize(ref_allele, alt_allele);

            // Reconstruct line
            let mut result = Vec::with_capacity(line.len() + 8);
            result.extend_from_slice(&line[..ref_start]);
            result.extend_from_slice(norm_ref.as_bytes());
            result.push(b'\t');
            result.extend_from_slice(norm_alt.as_bytes());
            result.extend_from_slice(&line[alt_end..]);

            result
        })
        .collect();

    eprintln!(
        "[Timing] Normalize: {:.3}s",
        process_start.elapsed().as_secs_f64()
    );

    // Write output using mmap for better performance
    let write_start = Instant::now();

    // Calculate exact output size
    let total_output_size: usize = normalized.iter().map(|v| v.len() + 1).sum(); // +1 for newline

    // Create output file with exact size
    let out_file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(true)
        .open(output)?;

    out_file.set_len(total_output_size as u64)?;

    // Memory-map output file
    let mut mmap_out = unsafe { MmapOptions::new().map_mut(&out_file)? };

    // Write in parallel chunks
    let chunk_size = 10_000;
    let chunks: Vec<_> = normalized.chunks(chunk_size).collect();

    let mut offset = 0;
    for chunk in chunks {
        for line in chunk {
            let len = line.len();
            mmap_out[offset..offset + len].copy_from_slice(line);
            offset += len;
            mmap_out[offset] = b'\n';
            offset += 1;
        }
    }

    // Ensure data is written to disk
    mmap_out.flush()?;

    eprintln!(
        "[Timing] Write: {:.3}s",
        write_start.elapsed().as_secs_f64()
    );
    eprintln!(
        "[Timing] Total: {:.3}s",
        total_start.elapsed().as_secs_f64()
    );

    Ok(())
}
