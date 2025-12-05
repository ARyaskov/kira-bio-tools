#![cfg(feature = "opencl")]

use crate::annotate_index::{AniEntry, AniEntryCL, AniIndex};
use crate::chr_name_to_id;
use anyhow::{Context, Result};
use crossbeam::channel;
use ocl::{Buffer, ProQue};
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

pub struct OpenCLv2 {
    pq: ProQue,
    buf_g: Buffer<u32>,
    buf_chrpos: Buffer<u32>,
    buf_refhash: Buffer<u64>,
    buf_althash: Buffer<u64>,
    buf_out: Buffer<i64>,
    batch_cap: usize,
    m: u32,
}

impl OpenCLv2 {
    pub fn new(ani: &AniIndex, batch_cap: usize) -> Result<Self> {
        let src = include_str!("ani_opencl_v2.cl");

        let pq = ProQue::builder().src(src).dims(1).build()?;

        // --- Upload MPH g[] ---
        let buf_g = pq
            .buffer_builder::<u32>()
            .len(ani.mph.g.len())
            .copy_host_slice(&ani.mph.g)
            .build()
            .context("Failed to create buf_g")?;

        let buf_chrpos = pq.buffer_builder::<u32>().len(batch_cap).build()?;
        let buf_refhash = pq.buffer_builder::<u64>().len(batch_cap).build()?;
        let buf_althash = pq.buffer_builder::<u64>().len(batch_cap).build()?;
        let buf_out = pq.buffer_builder::<i64>().len(batch_cap).build()?;

        Ok(Self {
            pq,
            buf_g,
            buf_chrpos,
            buf_refhash,
            buf_althash,
            buf_out,
            batch_cap,
            m: ani.mph.m,
        })
    }

    #[inline]
    pub fn run_batch(
        &self,
        chr_pos: &[u32],
        ref_hash: &[u64],
        alt_hash: &[u64],
    ) -> Result<Vec<i64>> {
        let batch = chr_pos.len();

        self.buf_chrpos.write(chr_pos).enq()?;
        self.buf_refhash.write(ref_hash).enq()?;
        self.buf_althash.write(alt_hash).enq()?;

        let kernel = self
            .pq
            .kernel_builder("ani_lookup_kernel_v2")
            .arg(&self.buf_chrpos)
            .arg(&self.buf_refhash)
            .arg(&self.buf_althash)
            .arg(&self.buf_g)
            .arg(self.m)
            .arg(&self.buf_out)
            .arg(batch as i32)
            .global_work_size(batch)
            .build()?;

        unsafe {
            kernel.enq()?;
        }

        let mut out = vec![0i64; batch];
        self.buf_out.read(&mut out).enq()?;

        Ok(out)
    }
}

/// ==== NEW PIPELINE —_GPU CALLED IN ONE THREAD_====
///
/// All parsing + hashing = Rayon threads
/// GPU = only one thread.
/// Communication via crossbeam bounded channel.

pub fn annotate_vcf_opencl_v2(
    gpu: &OpenCLv2,
    ani: &AniIndex,
    input: &Path,
    output: &Path,
) -> Result<()> {
    const BATCH: usize = 200_000;

    let fin = File::open(input)?;
    let rdr = BufReader::new(fin);

    let fout = File::create(output)?;
    let mut bw = BufWriter::new(fout);

    let mut first_data_line: Option<String> = None;

    {
        let fin2 = File::open(input)?;
        let mut rdr2 = BufReader::new(fin2);

        let mut buf = String::new();
        loop {
            buf.clear();
            let n = rdr2.read_line(&mut buf)?;
            if n == 0 {
                break;
            }
            if buf.starts_with('#') {
                bw.write_all(buf.as_bytes())?;
            } else {
                first_data_line = Some(buf.trim_end().to_string());
                break;
            }
        }
    }

    let (tx, rx) = channel::bounded::<(String, u32, u64, u64)>(BATCH * 2);

    let parse_thread = std::thread::spawn({
        let tx = tx.clone();
        move || {
            if let Some(line) = first_data_line {
                if let Some((chr, pos, ref_, alt)) = crate::annotate::parse_fields(&line) {
                    let chrid = chr_name_to_id(chr).unwrap() as u32;
                    let chrpos_u32 = (chrid << 24) | pos;
                    let h_ref = fxhash::hash64(ref_.as_bytes());
                    let h_alt = fxhash::hash64(alt.as_bytes());
                    tx.send((line, chrpos_u32, h_ref, h_alt)).unwrap();
                } else {
                    tx.send((line, 0, 0, 0)).unwrap();
                }
            }

            rdr.lines().par_bridge().for_each(|raw| {
                let line = raw.unwrap();

                if line.starts_with('#') {
                    return;
                }

                if let Some((chr, pos, ref_, alt)) = crate::annotate::parse_fields(&line) {
                    let chrid = chr_name_to_id(chr).unwrap() as u32;
                    let chrpos_u32 = (chrid << 24) | pos;
                    let h_ref = fxhash::hash64(ref_.as_bytes());
                    let h_alt = fxhash::hash64(alt.as_bytes());
                    tx.send((line, chrpos_u32, h_ref, h_alt)).unwrap();
                } else {
                    tx.send((line, 0, 0, 0)).unwrap();
                }
            });
        }
    });

    drop(tx);

    let mut chrpos = Vec::<u32>::with_capacity(BATCH);
    let mut refhash = Vec::<u64>::with_capacity(BATCH);
    let mut althash = Vec::<u64>::with_capacity(BATCH);
    let mut lines = Vec::<String>::with_capacity(BATCH);

    for (line, cp, rh, ah) in rx {
        chrpos.push(cp);
        refhash.push(rh);
        althash.push(ah);
        lines.push(line);

        if chrpos.len() >= BATCH {
            flush_batch(
                gpu,
                ani,
                &mut bw,
                &mut chrpos,
                &mut refhash,
                &mut althash,
                &mut lines,
            )?;
        }
    }

    flush_batch(
        gpu,
        ani,
        &mut bw,
        &mut chrpos,
        &mut refhash,
        &mut althash,
        &mut lines,
    )?;

    parse_thread.join().unwrap();

    Ok(())
}

fn flush_batch(
    gpu: &OpenCLv2,
    ani: &AniIndex,
    bw: &mut BufWriter<File>,
    chrpos: &mut Vec<u32>,
    refhash: &mut Vec<u64>,
    althash: &mut Vec<u64>,
    lines: &mut Vec<String>,
) -> Result<()> {
    let n = chrpos.len();
    if n == 0 {
        return Ok(());
    }

    let out = gpu.run_batch(&chrpos[..], &refhash[..], &althash[..])?;

    for i in 0..n {
        let idx = out[i];
        let line = &lines[i];

        if idx >= 0 {
            let e = ani.entries[idx as usize];
            let info = crate::annotate_index::read_cstring(&ani.string_block, e.info_ofs as usize);

            let base = crate::annotate::extract_info(line);
            let merged = crate::annotate::merge_info(base, info);

            let mut cols: Vec<&str> = line.split('\t').collect();
            cols[7] = &merged;

            let newline = cols.join("\t");
            bw.write_all(newline.as_bytes())?;
            bw.write_all(b"\n")?;
        } else {
            bw.write_all(line.as_bytes())?;
            bw.write_all(b"\n")?;
        }
    }

    chrpos.clear();
    refhash.clear();
    althash.clear();
    lines.clear();

    Ok(())
}
