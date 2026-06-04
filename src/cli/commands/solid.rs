use std::fs::File;
use std::io::BufWriter;

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

    // ── 2. Align → SAM bytes in RAM ─────────────────────────────────────────
    let t = std::time::Instant::now();
    let sam_bytes = aligner.align_to_sam_bytes(index, &reads).context("alignment")?;
    eprintln!(
        "[KIRA_SOLID] aligned: {:.1} MB SAM in {:.1}s",
        sam_bytes.len() as f64 / 1e6,
        t.elapsed().as_secs_f64()
    );

    // ── 3. Parse SAM bytes → RecordBuf (noodles) ────────────────────────────
    let mut reader = sam::io::Reader::new(std::io::Cursor::new(&sam_bytes[..]));
    let header = reader.read_header().context("parse SAM header")?;
    let mut records: Vec<RecordBuf> = Vec::new();
    for rec in reader.record_bufs(&header) {
        records.push(rec.context("parse SAM record")?);
    }
    drop(sam_bytes);
    eprintln!("[KIRA_SOLID] parsed {} records", records.len());

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
    let live: Vec<LiveRead> = sorted
        .iter()
        .filter_map(|rb| build_live_from_cram(rb, 0))
        .collect();

    // ── 6. In-memory BamReader ──────────────────────────────────────────────
    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|k| k.to_string())
        .collect();
    let sample = read_group
        .as_deref()
        .and_then(|rg| rg.split('\t').find_map(|f| f.strip_prefix("SM:")))
        .map(|s| s.to_string())
        .unwrap_or_else(|| "sample".to_string());
    let bam = BamReader::from_parts(live, header, vec![sample], ref_names);

    // ── 7. mpileup → VCF (the only disk write) ──────────────────────────────
    let mut mp_argv: Vec<String> = vec![
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
        mp_argv.push("--variants-only".to_string());
    }
    if args.mpileup_min_mq > 0 {
        mp_argv.push("--min-MQ".to_string());
        mp_argv.push(args.mpileup_min_mq.to_string());
    }
    mp_argv.push("--min-BQ".to_string());
    mp_argv.push(args.mpileup_min_bq.to_string());
    // run_mpileup_from_bams ignores `inputs`; satisfy the required positional.
    mp_argv.push("<mem>".to_string());
    let mp_args = MpileupArgs::try_parse_from(&mp_argv).context("assemble mpileup args")?;

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
