#![cfg(feature = "opencl")]

use crate::util::fast_hash64;
use anyhow::{Context, Result};
use ocl::{Buffer, ProQue};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::annotate::constants::OUTPUT_BUFFER_SIZE;
use crate::annotate::cpu_v2::field_metadata::{iter_ani_header_lines, load_and_infer_metadata};
use crate::annotate::cpu_v2::vcf_parsing::{parse_vcf_record, ParsedVcfRecord};
use crate::annotate::cpu_v2::{
    annotate_record_with_bundles, build_sample_map, expand_column_specs,
    extract_samples_from_headers, merge_annotation_headers, ColumnSpec,
};
use crate::annotate::reader::{StreamingVcfReader, VcfAnnotationReader};
use crate::annotate::structs::ani::AniIndex;
use crate::annotate::structs::annotate_mode::AnnotateMode;
use crate::annotate::structs::bundle::{AnnotationBundle, FieldNumber};
use crate::bgzf::BgzfWriter;
use crate::util::{chr_name_to_id, detect_format, VcfFormat};

const LINE_BATCH: usize = 200_000;

pub struct OpenCLv2 {
    pq: ProQue,
    buf_g: Buffer<u32>,
    buf_keys: Buffer<u64>,
    buf_entry_keys: Buffer<u64>,
    buf_out: Buffer<u32>,
    batch_cap: usize,
    m: u32,
    n: u32,
    salt: u64,
}

impl OpenCLv2 {
    pub fn new(ani: &AniIndex, batch_cap: usize) -> Result<Self> {
        let src = include_str!("ani_opencl_v2.cl");

        let pq = ProQue::builder().src(src).dims(1).build()?;

        let buf_g = pq
            .buffer_builder::<u32>()
            .len(ani.mph.g.len())
            .copy_host_slice(&ani.mph.g)
            .build()
            .context("Failed to create buf_g")?;

        let buf_keys = pq.buffer_builder::<u64>().len(batch_cap).build()?;
        let entry_keys = build_entry_keys(ani);
        let buf_entry_keys = pq
            .buffer_builder::<u64>()
            .len(entry_keys.len())
            .copy_host_slice(&entry_keys)
            .build()
            .context("Failed to create buf_entry_keys")?;
        let buf_out = pq.buffer_builder::<u32>().len(batch_cap).build()?;

        Ok(Self {
            pq,
            buf_g,
            buf_keys,
            buf_entry_keys,
            buf_out,
            batch_cap,
            m: ani.mph.m,
            n: ani.mph.n as u32,
            salt: ani.mph.salt,
        })
    }

    #[inline]
    pub fn run_batch(&self, keys: &[u64]) -> Result<Vec<u32>> {
        let batch = keys.len();
        if batch > self.batch_cap {
            anyhow::bail!("OpenCL batch size {} exceeds cap {}", batch, self.batch_cap);
        }

        self.buf_keys.write(keys).enq()?;

        let kernel = self
            .pq
            .kernel_builder("ani_lookup_kernel_v2")
            .arg(&self.buf_keys)
            .arg(&self.buf_g)
            .arg(self.m)
            .arg(self.n)
            .arg(self.salt)
            .arg(&self.buf_entry_keys)
            .arg(&self.buf_out)
            .arg(batch as i32)
            .global_work_size(batch)
            .build()?;

        unsafe {
            kernel.enq()?;
        }

        let mut out = vec![0u32; batch];
        self.buf_out.read(&mut out).enq()?;
        Ok(out)
    }
}

pub fn annotate_vcf_opencl_v2(
    gpu: &OpenCLv2,
    ani: &AniIndex,
    input: &Path,
    output: &Path,
    columns: &[String],
) -> Result<()> {
    let mut column_specs = ColumnSpec::parse_all(columns);
    let info_overwrite_all = column_specs
        .iter()
        .any(|c| c.key.eq_ignore_ascii_case("INFO"));
    let format_overwrite_all = column_specs
        .iter()
        .any(|c| c.key.eq_ignore_ascii_case("FMT") || c.key.eq_ignore_ascii_case("FORMAT"));

    let field_meta = load_and_infer_metadata(ani, false)?;
    let ani_headers = iter_ani_header_lines(ani);
    column_specs = expand_column_specs(&column_specs, &ani_headers, &field_meta);
    let column_modes: Vec<(String, AnnotateMode)> = column_specs
        .iter()
        .map(|c| (c.key.clone(), c.mode))
        .collect();

    let input_format = detect_format(input)?;
    let output_ext = output.extension().and_then(|s| s.to_str()).unwrap_or("");
    let output_wants_bgzf = matches!(output_ext, "gz" | "bgz" | "bgzf");
    let use_bgzf = matches!(input_format, VcfFormat::Bgzf) || output_wants_bgzf;

    let input_reader = VcfAnnotationReader::open(input)?;
    let streaming_reader = StreamingVcfReader::new(input_reader);
    let (headers, mut reader) = streaming_reader.into_headers_and_self()?;

    let merged_headers = merge_annotation_headers(&headers, &ani_headers, &column_specs)?;
    let input_samples = extract_samples_from_headers(&headers);
    let db_samples = extract_samples_from_headers(&ani_headers);
    let sample_map = build_sample_map(&input_samples, &db_samples);

    if use_bgzf {
        let mut writer = BgzfWriter::create(output)?;
        for h in &merged_headers {
            writeln!(writer, "{}", h)?;
        }
        process_records(
            &mut writer,
            &mut reader,
            gpu,
            ani,
            &field_meta,
            &column_modes,
            &sample_map,
            info_overwrite_all,
            format_overwrite_all,
        )?;
        writer.finish()?;
    } else {
        let file = File::create(output)?;
        let mut writer = BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, file);
        for h in &merged_headers {
            writeln!(writer, "{}", h)?;
        }
        process_records(
            &mut writer,
            &mut reader,
            gpu,
            ani,
            &field_meta,
            &column_modes,
            &sample_map,
            info_overwrite_all,
            format_overwrite_all,
        )?;
        writer.flush()?;
    }

    Ok(())
}

struct ParsedLine {
    raw: String,
    parsed: Option<ParsedVcfRecord>,
    alt_alleles: Vec<String>,
    chr_id: Option<u8>,
}

fn process_records<W: Write>(
    writer: &mut W,
    reader: &mut StreamingVcfReader,
    gpu: &OpenCLv2,
    ani: &AniIndex,
    field_meta: &HashMap<String, FieldNumber>,
    column_modes: &[(String, AnnotateMode)],
    sample_map: &[Option<usize>],
    info_overwrite_all: bool,
    format_overwrite_all: bool,
) -> Result<()> {
    let mut lines: Vec<String> = Vec::with_capacity(LINE_BATCH);
    loop {
        lines.clear();
        while lines.len() < LINE_BATCH {
            let Some(line) = reader.read_line()? else {
                break;
            };
            if line.starts_with('#') {
                continue;
            }
            lines.push(line);
        }

        if lines.is_empty() {
            break;
        }

        let parsed_lines: Vec<ParsedLine> = lines
            .par_iter()
            .map(|line| {
                let parsed = parse_vcf_record(line);
                let (alt_alleles, chr_id) = if let Some(ref p) = parsed {
                    let chr_id = chr_name_to_id(&p.chrom);
                    let alts = p.alt.split(',').map(|s| s.to_string()).collect::<Vec<_>>();
                    (alts, chr_id)
                } else {
                    (Vec::new(), None)
                };

                ParsedLine {
                    raw: line.clone(),
                    parsed,
                    alt_alleles,
                    chr_id,
                }
            })
            .collect();

        flush_batch(
            writer,
            gpu,
            ani,
            field_meta,
            column_modes,
            sample_map,
            info_overwrite_all,
            format_overwrite_all,
            &parsed_lines,
        )?;
    }

    Ok(())
}

fn flush_batch<W: Write>(
    writer: &mut W,
    gpu: &OpenCLv2,
    ani: &AniIndex,
    field_meta: &HashMap<String, FieldNumber>,
    column_modes: &[(String, AnnotateMode)],
    sample_map: &[Option<usize>],
    info_overwrite_all: bool,
    format_overwrite_all: bool,
    lines: &[ParsedLine],
) -> Result<()> {
    if lines.is_empty() {
        return Ok(());
    }

    let mut bundles_per_line: Vec<Vec<(usize, AnnotationBundle)>> = vec![Vec::new(); lines.len()];

    let mut keys: Vec<u64> = Vec::new();
    let mut key_line_idx: Vec<usize> = Vec::new();
    let mut key_alt_idx: Vec<usize> = Vec::new();

    keys.reserve(lines.len());
    key_line_idx.reserve(lines.len());
    key_alt_idx.reserve(lines.len());

    for (line_idx, line) in lines.iter().enumerate() {
        let Some(ref parsed) = line.parsed else {
            continue;
        };
        let Some(chr_id) = line.chr_id else {
            continue;
        };
        for (alt_idx, alt) in line.alt_alleles.iter().enumerate() {
            keys.push(make_key(chr_id, parsed.pos, &parsed.ref_allele, alt));
            key_line_idx.push(line_idx);
            key_alt_idx.push(alt_idx);
        }
    }

    let mut offset = 0usize;
    while offset < keys.len() {
        let end = (offset + gpu.batch_cap).min(keys.len());
        let idxs = gpu.run_batch(&keys[offset..end])?;
        for (i, idx) in idxs.iter().enumerate() {
            let global = offset + i;
            let line_idx = key_line_idx[global];
            let alt_idx = key_alt_idx[global];

            let Some(ref parsed) = lines[line_idx].parsed else {
                continue;
            };
            let Some(chr_id) = lines[line_idx].chr_id else {
                continue;
            };
            if *idx == u32::MAX {
                continue;
            }

            let bundle = ani.build_bundle_from_entry(&ani.entries[*idx as usize]);
            bundles_per_line[line_idx].push((alt_idx, bundle));
        }
        offset = end;
    }

    let outputs: Vec<String> = lines
        .par_iter()
        .enumerate()
        .map(|(i, line)| {
            if let Some(ref parsed) = line.parsed {
                annotate_record_with_bundles(
                    parsed,
                    &bundles_per_line[i],
                    field_meta,
                    column_modes,
                    sample_map,
                    info_overwrite_all,
                    format_overwrite_all,
                    false,
                )
            } else {
                line.raw.clone()
            }
        })
        .collect();

    for out in outputs {
        writeln!(writer, "{}", out)?;
    }

    Ok(())
}

fn make_key(chr_id: u8, pos: u32, ref_allele: &str, alt: &str) -> u64 {
    let mut h = (chr_id as u64) << 32 | pos as u64;
    h ^= fast_hash64(ref_allele.as_bytes());
    h ^= fast_hash64(alt.as_bytes());
    h
}

fn build_entry_keys(ani: &AniIndex) -> Vec<u64> {
    let mut keys = Vec::with_capacity(ani.entries.len());
    for entry in &ani.entries {
        let ref_str = ani.read_cstring(entry.ref_ofs as usize);
        let alt_str = ani.read_cstring(entry.alt_ofs as usize);
        let key = make_key(entry.chr_id, entry.pos, ref_str.as_ref(), alt_str.as_ref());
        keys.push(key);
    }
    keys
}
