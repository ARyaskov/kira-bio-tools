use std::fs::File;
use std::io::{BufWriter, Write as _};

use anyhow::{Context, Result};
use clap::Parser;
use noodles_sam::{self as sam, alignment::RecordBuf};

use crate::bam::BamReader;
use crate::bam::pileup::{LiveRead, build_live_from_cram};
use crate::cli::args::{MpileupArgs, SolidArgs};
use crate::cli::commands::mpileup::run_mpileup_from_bams;

use kira_ls_aligner::cli::commands::mem::{FusedAlignerParams, build_short_pe_aligner};
use kira_ls_aligner::index::{Index, IndexConfig};
use kira_ls_aligner::io::read_reference;

/// Fused align → sort/markdup → mpileup pipeline in a single process. Records
/// flow in memory between stages; only the output VCF touches disk.
pub fn cmd_solid(args: SolidArgs) -> Result<()> {
    eprintln!("[KIRA_SOLID] fused pipeline (in-memory): align -> sort/markdup -> mpileup");

    // ── 1. Build the aligner + resolve/load the index ───────────────────────
    // Accept literal `\t` in --aligner-rg (shells don't expand it).
    let read_group = args.aligner_rg.as_ref().map(|s| s.replace("\\t", "\t"));
    let paired = args.aligner_r2.is_some();
    let n_files = if paired { 2 } else { 1 };
    let params = FusedAlignerParams {
        reference: args.aligner_ref.clone(),
        index: args.aligner_index.clone(),
        threads: args.threads,
        num_p_threads: None,
        num_e_threads: None,
        batch_bases: args.aligner_batch,
        read_group: read_group.clone(),
        paired,
        interleaved: false,
        insert_size: args.aligner_insert_size.clone(),
        n_read_files: n_files,
    };
    let (aligner, idx_path) = build_short_pe_aligner(&params)?;
    let index = match idx_path {
        Some(p) => {
            eprintln!("[KIRA_SOLID] loading index {}", p.display());
            Index::load(&p).context("load index")?
        }
        None => {
            eprintln!("[KIRA_SOLID] building index in memory from {}", args.aligner_ref.display());
            let reference = read_reference(&args.aligner_ref).context("load reference")?;
            Index::build(reference, IndexConfig::default())
        }
    };
    let mut reads = vec![args.aligner_r1.clone()];
    if let Some(r2) = &args.aligner_r2 {
        reads.push(r2.clone());
    }

    let mp_args = build_mpileup_args(&args)?;
    let sample = read_group
        .as_deref()
        .and_then(|rg| rg.split('\t').find_map(|f| f.strip_prefix("SM:")))
        .map(|s| s.to_string())
        .unwrap_or_else(|| "sample".to_string());

    if args.window_mb > 0 {
        return run_windowed(&args, aligner, index, &reads, &mp_args, &sample);
    }

    // ── 2. Align → records, batch by batch, no SAM text ─────────────────────
    let t = std::time::Instant::now();
    let header = parse_header(&aligner.sam_header_bytes(&index).context("SAM header")?)?;
    let mut records: Vec<RecordBuf> = Vec::new();
    aligner
        .align_streaming(index, &reads, |batch| {
            super::solid_records::append_batch(batch, MAX_ALIGNMENTS, &mut records);
            Ok(())
        })
        .context("alignment")?;
    eprintln!(
        "[KIRA_SOLID] aligned {} records in {:.1}s (no SAM text)",
        records.len(),
        t.elapsed().as_secs_f64()
    );

    // ── 4. Coordinate sort (+markdup) in RAM ────────────────────────────────
    let t = std::time::Instant::now();
    let sorted = kira_bam::sort::sort_and_markdup_in_memory(records, args.bam_markdup)
        .context("sort/markdup")?;
    eprintln!(
        "[KIRA_SOLID] sorted{} {} records in {:.1}s",
        if args.bam_markdup { "+markdup" } else { "" },
        sorted.len(),
        t.elapsed().as_secs_f64()
    );

    // KIRA_SOLID_DUMP_BAM=path — dump the post-markdup BAM for offline caller iteration.
    if let Ok(dump) = std::env::var("KIRA_SOLID_DUMP_BAM") {
        use noodles_bam as bam;
        use noodles_sam::alignment::io::Write as _;
        let f = File::create(&dump).context("create dump BAM")?;
        let mut w = bam::io::Writer::new(f);
        w.write_header(&header).context("write dump BAM header")?;
        for rb in &sorted {
            w.write_alignment_record(&header, rb as &dyn sam::alignment::Record)
                .context("write dump BAM record")?;
        }
        w.try_finish().context("finish dump BAM")?;
        eprintln!("[KIRA_SOLID] dumped post-markdup BAM -> {}", dump);
    }

    // ── 5. RecordBuf → LiveRead (mpileup's record form) ─────────────────────
    // Consume `sorted` as we go: each `RecordBuf` is dropped right after its
    // `LiveRead` is built, so the two representations never both exist in full.
    // Iterating by reference instead keeps a complete second copy of every
    // alignment resident for the rest of the run.
    let live: Vec<LiveRead> = sorted
        .into_iter()
        .filter_map(|rb| build_live_from_cram(&rb, 0))
        .collect();

    // ── 6. In-memory BamReader ──────────────────────────────────────────────
    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|k| k.to_string())
        .collect();
    let bam = BamReader::from_parts(live, header, vec![sample], ref_names);

    // ── 7. mpileup → VCF (the only disk write) ──────────────────────────────
    let mut out = BufWriter::with_capacity(1 << 20, File::create(&args.output).context("create VCF")?);
    let t = std::time::Instant::now();
    run_mpileup_from_bams(vec![bam], &mp_args, &mut out).context("mpileup")?;
    eprintln!(
        "[KIRA_SOLID] mpileup -> {} in {:.1}s",
        args.output.display(),
        t.elapsed().as_secs_f64()
    );
    Ok(())
}

/// The fused aligner is configured for primary-only output.
const MAX_ALIGNMENTS: usize = 1;

/// Removes the window scratch directory when the windowed run ends, including on
/// the error paths.
struct ScratchDir(std::path::PathBuf);

impl Drop for ScratchDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn parse_header(bytes: &[u8]) -> Result<sam::Header> {
    sam::io::Reader::new(std::io::Cursor::new(bytes))
        .read_header()
        .context("parse SAM header")
}

/// Assemble the mpileup arguments shared by both execution paths.
fn build_mpileup_args(args: &SolidArgs) -> Result<MpileupArgs> {
    let mut argv: Vec<String> = vec![
        "solid-mpileup".to_string(),
        "--output".to_string(),
        args.output.display().to_string(),
        "--fasta-ref".to_string(),
        args.aligner_ref.display().to_string(),
        "--annotate".to_string(),
        args.mpileup_annotate.clone(),
        "--threads".to_string(),
        args.threads.to_string(),
        "--max-depth".to_string(),
        args.mpileup_max_depth.to_string(),
    ];
    if args.mpileup_variants_only {
        argv.push("--variants-only".to_string());
    }
    if args.mpileup_min_mq > 0 {
        argv.push("--min-MQ".to_string());
        argv.push(args.mpileup_min_mq.to_string());
    }
    argv.push("--min-BQ".to_string());
    argv.push(args.mpileup_min_bq.to_string());
    // run_mpileup_from_bams ignores `inputs`; satisfy the required positional.
    argv.push("<mem>".to_string());
    MpileupArgs::try_parse_from(&argv).context("assemble mpileup args")
}

/// Windowed execution: spill alignments to per-window temporary BAMs, then sort,
/// deduplicate and call one window at a time.
///
/// Peak memory becomes one window's depth rather than the whole run, at the cost
/// of writing the alignments to scratch once. See `solid_windows` for the
/// boundary handling — reads that straddle a window are written to both sides.
fn run_windowed(
    args: &SolidArgs,
    aligner: kira_ls_aligner::aligner_core::Aligner,
    index: Index,
    reads: &[std::path::PathBuf],
    mp_args: &MpileupArgs,
    sample: &str,
) -> Result<()> {
    use super::solid_windows::{WindowSpiller, resolve_tmpdir};

    let window_len: u32 = args
        .window_mb
        .saturating_mul(1_000_000)
        .max(1_000_000);
    let tmp_root = resolve_tmpdir(args.window_tmpdir.as_deref(), &args.output);
    let scratch = tmp_root.join(format!("kira-solid-win-{}", std::process::id()));
    std::fs::create_dir_all(&scratch)
        .with_context(|| format!("create window scratch {}", scratch.display()))?;
    // Removed on the way out; per-window files are unlinked as each is consumed.
    let _scratch_guard = ScratchDir(scratch.clone());
    eprintln!(
        "[KIRA_SOLID] windowed mode: {} Mb windows, scratch {}",
        window_len / 1_000_000,
        scratch.display()
    );

    // ── Phase 1: align and spill ────────────────────────────────────────────
    let t = std::time::Instant::now();
    let header = parse_header(&aligner.sam_header_bytes(&index).context("SAM header")?)?;
    let mut spiller = WindowSpiller::new(scratch.clone(), window_len);
    let mut batch_records: Vec<RecordBuf> = Vec::new();
    let mut total: u64 = 0;
    {
        let hdr = &header;
        let spiller_ref = &mut spiller;
        let batch_buf = &mut batch_records;
        let total_ref = &mut total;
        aligner
            .align_streaming(index, reads, move |batch| {
                // Spilled as produced; only one batch is ever resident.
                batch_buf.clear();
                super::solid_records::append_batch(batch, MAX_ALIGNMENTS, batch_buf);
                for rec in batch_buf.iter() {
                    spiller_ref.push(hdr, rec)?;
                }
                *total_ref += batch_buf.len() as u64;
                Ok(())
            })
            .context("alignment")?;
    }
    let unplaced = spiller.unplaced_records();
    let windows = spiller.finish().context("finish window spill")?;
    eprintln!(
        "[KIRA_SOLID] aligned {total} records into {} window(s) ({unplaced} unplaced) in {:.1}s",
        windows.len(),
        t.elapsed().as_secs_f64()
    );

    // ── Phase 2: one window at a time ───────────────────────────────────────
    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|k| k.to_string())
        .collect();
    let mut out =
        BufWriter::with_capacity(1 << 20, File::create(&args.output).context("create VCF")?);
    let mut header_written = false;
    let t = std::time::Instant::now();

    for w in &windows {
        let (_, records) = w.load().context("load window")?;
        let sorted = kira_bam::sort::sort_and_markdup_in_memory(records, args.bam_markdup)
            .context("sort/markdup")?;
        let live: Vec<LiveRead> = sorted
            .into_iter()
            .filter_map(|rb| build_live_from_cram(&rb, 0))
            .collect();
        let bam = BamReader::from_parts(
            live,
            header.clone(),
            vec![sample.to_string()],
            ref_names.clone(),
        );

        let mut vcf: Vec<u8> = Vec::new();
        run_mpileup_from_bams(vec![bam], mp_args, &mut vcf).context("mpileup")?;
        append_window_vcf(&mut out, &vcf, w, &header, &mut header_written)?;
        w.remove();
    }
    out.flush().context("flush VCF")?;
    eprintln!(
        "[KIRA_SOLID] called {} window(s) -> {} in {:.1}s",
        windows.len(),
        args.output.display(),
        t.elapsed().as_secs_f64()
    );
    Ok(())
}

/// Append one window's VCF to the output, writing the header once and keeping
/// only the calls that belong to this window.
///
/// A window's spill deliberately contains reads starting before it, so its pileup
/// also produces calls in the preceding window's territory. Those are dropped
/// here; the window that owns them emits them from its own complete pileup.
fn append_window_vcf<W: std::io::Write>(
    out: &mut W,
    vcf: &[u8],
    window: &super::solid_windows::SpilledWindow,
    header: &sam::Header,
    header_written: &mut bool,
) -> Result<()> {
    let (ref_id, widx) = window.id;
    let Some((name, seq)) = header.reference_sequences().get_index(ref_id) else {
        return Ok(());
    };
    let contig = String::from_utf8_lossy(name.as_ref()).into_owned();
    let len = usize::from(seq.length());
    let win_len = window.window_len();
    let lo = widx as usize * win_len + 1; // 1-based inclusive
    let hi = ((widx as usize + 1) * win_len).min(len); // 1-based inclusive

    for line in vcf.split(|&b| b == b'\n') {
        if line.is_empty() {
            continue;
        }
        if line[0] == b'#' {
            if !*header_written {
                out.write_all(line)?;
                out.write_all(b"\n")?;
            }
            continue;
        }
        let mut fields = line.splitn(3, |&b| b == b'\t');
        let (Some(chrom), Some(pos)) = (fields.next(), fields.next()) else {
            continue;
        };
        if chrom != contig.as_bytes() {
            continue;
        }
        let Ok(pos) = std::str::from_utf8(pos).unwrap_or("").parse::<usize>() else {
            continue;
        };
        if pos < lo || pos > hi {
            continue;
        }
        out.write_all(line)?;
        out.write_all(b"\n")?;
    }
    *header_written = true;
    Ok(())
}
