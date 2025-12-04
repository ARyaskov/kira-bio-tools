#![cfg(feature = "opencl")]

use crate::annotate_index::{AniEntryCL, AniIndex};
use crate::chr_name_to_id;
use anyhow::{Context, Result};
use ocl::{Buffer, ProQue};
use std::path::Path;

use std::fs::File;
use std::io::{BufWriter, Write};

const BATCH: usize = 20000;

pub struct OpenCLAni {
    pq: ProQue,
    buf_g: Buffer<u32>,
    buf_entries: Buffer<AniEntryCL>,
    m: u32,
    n: usize,
    string_block: Vec<u8>,
}

impl OpenCLAni {
    pub fn new(ani: &AniIndex) -> Result<Self> {
        println!("=== DEBUG ANI METADATA ===");
        println!("entries.len() = {}", ani.entries.len());
        println!("string_block.len() = {}", ani.string_block.len());
        println!("mph.n = {}", ani.mph.n);
        println!("mph.m = {}", ani.mph.m);
        println!("g.len() = {}", ani.mph.g.len());
        println!("sizeof(AniEntryCL) = {}", std::mem::size_of::<AniEntryCL>());
        println!(
            "sizeof(AniEntry)   = {}",
            std::mem::size_of::<crate::annotate_index::AniEntry>()
        );

        let src = include_str!("ani_opencl.cl");

        // Build OpenCL program
        let pq = ProQue::builder()
            .src(src)
            .dims(1)
            .build()
            .context("Failed to build OpenCL program")?;

        // Upload MPH g[]
        let buf_g = pq
            .buffer_builder::<u32>()
            .len(ani.mph.g.len())
            .copy_host_slice(&ani.mph.g)
            .build()
            .context("Failed to create g buffer")?;

        // Convert entries to AniEntryCL
        let entries_cl: Vec<AniEntryCL> = ani.entries.iter().map(|e| e.to_cl()).collect();

        let buf_entries = pq
            .buffer_builder::<AniEntryCL>()
            .len(entries_cl.len())
            .copy_host_slice(&entries_cl)
            .build()
            .context("Failed to create entry buffer")?;

        Ok(Self {
            pq,
            buf_g,
            buf_entries,
            m: ani.mph.m,
            n: ani.entries.len(),
            string_block: ani.string_block.clone(),
        })
    }

    pub fn lookup_batch(&self, keys: &[u64]) -> Result<Vec<i64>> {
        let n = keys.len();

        // Allocate buffers
        let buf_keys = self
            .pq
            .buffer_builder::<u64>()
            .len(n)
            .copy_host_slice(keys)
            .build()
            .context("Failed to upload batch keys")?;

        let buf_out = self
            .pq
            .buffer_builder::<i64>()
            .len(n)
            .build()
            .context("Failed to allocate batch output")?;

        // Build kernel
        let kernel = self
            .pq
            .kernel_builder("ani_lookup_kernel")
            .arg(&buf_keys)
            .arg(&self.buf_g)
            .arg(self.m)
            .arg(&self.buf_entries)
            .arg(&buf_out)
            .arg(n as i32)
            .global_work_size(n)
            .build()?;

        unsafe {
            kernel.enq()?;
        }

        let mut out = vec![0i64; n];
        buf_out
            .read(&mut out)
            .enq()
            .context("Failed to read GPU results")?;

        Ok(out)
    }
}

/// -----------
/// High-throughput annotate
/// -----------
pub fn annotate_vcf_ani_opencl(
    gpu: &OpenCLAni,
    ani: &AniIndex,
    input_vcf: &Path,
    output_vcf: &Path,
) -> Result<()> {
    use std::fs::File;
    use std::io::{BufRead, BufReader, BufWriter, Write};

    let fin = File::open(input_vcf)?;
    let rdr = BufReader::new(fin);

    let fout = File::create(output_vcf)?;
    let mut bw = BufWriter::new(fout);

    let mut batch_lines: Vec<String> = Vec::with_capacity(BATCH);
    let mut batch_keys: Vec<u64> = Vec::with_capacity(BATCH);

    for line in rdr.lines() {
        let line = line?;

        if line.starts_with('#') {
            bw.write_all(line.as_bytes())?;
            bw.write_all(b"\n")?;
            continue;
        }

        if let Some((chr, pos, ref_, alt)) = crate::annotate::parse_fields(&line) {
            // Build key (same as CPU version)
            let chr_id = chr_name_to_id(chr).unwrap_or(0);

            let mut h = fxhash::hash64(&[chr_id]);
            h ^= fxhash::hash64(pos.to_le_bytes().as_ref());
            h ^= fxhash::hash64(ref_.as_bytes());
            h ^= fxhash::hash64(alt.as_bytes());

            batch_keys.push(h);
            batch_lines.push(line);
        } else {
            // write as-is
            bw.write_all(line.as_bytes())?;
            bw.write_all(b"\n")?;
        }

        // If batch full → flush to GPU
        if batch_keys.len() >= BATCH {
            process_batch(&batch_lines, &batch_keys, gpu, ani, &mut bw)?;
            batch_keys.clear();
            batch_lines.clear();
        }
    }

    // Flush last batch
    if !batch_keys.is_empty() {
        process_batch(&batch_lines, &batch_keys, gpu, ani, &mut bw)?;
    }

    Ok(())
}

/// Process one batch using GPU lookup
fn process_batch(
    lines: &[String],
    keys: &[u64],
    gpu: &OpenCLAni,
    ani: &AniIndex,
    bw: &mut BufWriter<File>,
) -> Result<()> {
    let idxs = gpu.lookup_batch(keys)?;

    for (line, idx) in lines.iter().zip(idxs.iter()) {
        if *idx >= 0 && (*idx as usize) < ani.entries.len() {
            let entry = ani.entries[*idx as usize];
            let info =
                crate::annotate_index::read_cstring(&ani.string_block, entry.info_ofs as usize);

            let base = crate::annotate::extract_info(line);
            let merged = crate::annotate::merge_info(base, info);

            let mut cols: Vec<&str> = line.split('\t').collect();
            cols[7] = &merged;

            let newline = cols.join("\t");
            bw.write_all(newline.as_bytes())?;
            bw.write_all(b"\n")?;
        } else {
            // fallback
            bw.write_all(line.as_bytes())?;
            bw.write_all(b"\n")?;
        }
    }

    Ok(())
}
