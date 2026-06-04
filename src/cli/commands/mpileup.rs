use crate::bam::{
    AnnotateSpec, BamReader, ErrorModel, FlagFilters, LiveRead,
    PileupSite, PresetConfig, mpileup_engine_from_records, mpileup_engine_multi,
    parse_samples_filter,
};
use crate::call::GvcfBlocker;
use crate::cli::args::MpileupArgs;
use anyhow::{Context, Result, bail};
use fxhash::FxHashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

pub fn cmd_mpileup(args: MpileupArgs) -> Result<()> {
    if args.inputs.is_empty() { bail!("mpileup: need at least one BAM input"); }

    let preset = args.config.as_deref().map(PresetConfig::parse).transpose()?;
    let min_mq = preset.as_ref().and_then(|p| p.min_mq).unwrap_or(args.min_mq) as u8;
    let min_bq = preset.as_ref().and_then(|p| p.min_bq).unwrap_or(args.min_bq) as u8;
    let max_depth = preset.as_ref().and_then(|p| p.max_depth).unwrap_or(args.max_depth);
    let _indel_size = preset.as_ref().and_then(|p| p.indel_size).unwrap_or(args.indel_size);
    let baq_off = preset.as_ref().and_then(|p| p.no_baq).unwrap_or(args.no_baq);
    let skip_indels = args.skip_indels;
    let annotate = AnnotateSpec::parse(args.annotate.as_deref())?;
    let flag_filters = FlagFilters::from_full(args.ef, args.df, args.if_, args.nf,
        args.rf.as_deref(), args.ff.as_deref())?;
    let gap_frac = preset.as_ref().and_then(|p| p.gap_frac).unwrap_or(args.gap_frac);
    let min_ireads = args.min_ireads;
    let sample_filter = parse_samples_filter(args.samples.as_deref(), args.samples_file.as_deref())?;

    let fasta = args.fasta_ref.as_ref().map(load_fasta).transpose()?;
    let out_path = args.output.clone().unwrap_or_else(|| PathBuf::from("out.mpileup.vcf"));
    let mut out = BufWriter::with_capacity(1 << 20, File::create(&out_path).context("create output")?);

    let timing = std::env::var("KIRA_BT_TIMING").is_ok();
    let mut bams: Vec<BamReader> = Vec::with_capacity(args.inputs.len());
    // Single combined skip-list across all samples — pre-scan once per BAM, merged later.
    let mut combined_filter: Option<crate::bam::pos_filter::InterestingMap> = None;
    for p in &args.inputs {
        let t0 = std::time::Instant::now();
        let mut r = if let Some(reg) = &args.regions {
            BamReader::open_with_region(p, reg)?
        } else {
            BamReader::open(p)?
        };
        if timing { eprintln!("[KIRA_BT] BAM load: {:.1}s, {} records", t0.elapsed().as_secs_f64(), r.records_buf.len()); }
        if !flag_filters.exclude_flags == 0 || flag_filters.require_flags != 0 {
            r.records_buf.retain(|lr| flag_filters.passes(lr.flags));
        }
        let _ = max_depth;
        if !baq_off {
            if let Some(fa) = &fasta {
                let ref_names = r.ref_names.clone();
                let t_pre = std::time::Instant::now();
                let pre_arc = std::sync::Arc::new(std::mem::take(&mut r.records_buf));
                let pre = crate::bam::pos_filter::pre_scan(
                    std::slice::from_ref(&pre_arc),
                    fa,
                    &ref_names,
                    args.min_alt_reads.max(2),
                );
                r.records_buf = std::sync::Arc::try_unwrap(pre_arc)
                    .unwrap_or_else(|arc| (*arc).clone());
                if timing {
                    eprintln!(
                        "[KIRA_BT] pre-scan: {:.1}s, skip-list={}, baq-skip={}/{}",
                        t_pre.elapsed().as_secs_f64(),
                        pre.pos_filter.total(),
                        pre.skipped_baq,
                        pre.needs_baq.len(),
                    );
                }
                let t1 = std::time::Instant::now();
                crate::bam::reader::apply_hmm_baq_to_reads_masked(
                    &mut r.records_buf,
                    &ref_names,
                    fa,
                    Some(&pre.needs_baq),
                );
                if timing { eprintln!("[KIRA_BT] BAQ: {:.1}s", t1.elapsed().as_secs_f64()); }
                // First sample wins; subsequent samples would need merging — out of scope here.
                if combined_filter.is_none() {
                    combined_filter = Some(pre.pos_filter);
                }
            }
        }
        bams.push(r);
    }

    let mut samples: Vec<String> = Vec::with_capacity(bams.len());
    let mut sample_keep: Vec<bool> = Vec::with_capacity(bams.len());
    for (i, b) in bams.iter().enumerate() {
        let s = b.samples.first().cloned().unwrap_or_else(|| {
            args.inputs[i].file_stem().and_then(|s| s.to_str()).unwrap_or("sample").to_string()
        });
        let keep = sample_filter.as_ref().map_or(true, |list| list.iter().any(|n| n == &s));
        sample_keep.push(keep);
        samples.push(s);
    }
    let any_filter = sample_filter.is_some();

    let ref_names: Vec<String> = bams[0].ref_names.clone();

    write_vcf_header(&mut out, &annotate, &ref_names, &samples, &sample_keep, any_filter)?;

    let em = ErrorModel::new();
    let gvcf_opts: Option<Vec<u32>> = args.gvcf.as_ref().map(|s| {
        s.split(',').filter_map(|t| t.parse::<u32>().ok()).collect()
    });
    let mut gvcf_blocker: Option<GvcfBlocker> = gvcf_opts.map(GvcfBlocker::new);

    let n_chunks = args.threads.max(1);
    // Chunk-parallel walk: applies when threads > 1 AND no gvcf state to maintain across positions.
    let parallel_eligible = n_chunks > 1 && gvcf_blocker.is_none();

    let hp_indel = std::env::var("KIRA_HP_INDEL").map(|v| v != "0").unwrap_or(true);
    let ctx = EmitCtx {
        args: &args,
        em: &em,
        annotate: &annotate,
        fasta: fasta.as_ref(),
        ref_names: &ref_names,
        sample_keep: &sample_keep,
        any_filter,
        skip_indels,
        min_ireads,
        gap_frac,
        hp_indel,
    };

    let t_engine = std::time::Instant::now();
    if parallel_eligible {
        if timing { eprintln!("[KIRA_BT] engine: parallel, {} chunks", n_chunks); }
        run_parallel(bams, n_chunks, &ctx, min_mq, min_bq, skip_indels, combined_filter, &mut out)?;
    } else {
        if timing { eprintln!("[KIRA_BT] engine: sequential"); }
        run_sequential(bams, &ctx, min_mq, min_bq, skip_indels, &mut gvcf_blocker, &mut out)?;
        if let Some(mut g) = gvcf_blocker { g.flush(&mut out)?; }
    }
    if timing { eprintln!("[KIRA_BT] engine done in {:.1}s", t_engine.elapsed().as_secs_f64()); }

    out.flush()?;
    Ok(())
}

/// Run mpileup over already-built [`BamReader`]s (records in memory) and write
/// VCF to `out`. Same logic as [`cmd_mpileup`] minus the file-opening — for the
/// fused `solid` pipeline that hands sorted in-memory records straight here.
pub fn run_mpileup_from_bams(
    mut bams: Vec<BamReader>,
    args: &MpileupArgs,
    out: &mut BufWriter<File>,
) -> Result<()> {
    if bams.is_empty() {
        bail!("mpileup: no records");
    }
    let preset = args.config.as_deref().map(PresetConfig::parse).transpose()?;
    let min_mq = preset.as_ref().and_then(|p| p.min_mq).unwrap_or(args.min_mq) as u8;
    let min_bq = preset.as_ref().and_then(|p| p.min_bq).unwrap_or(args.min_bq) as u8;
    let max_depth = preset.as_ref().and_then(|p| p.max_depth).unwrap_or(args.max_depth);
    let baq_off = preset.as_ref().and_then(|p| p.no_baq).unwrap_or(args.no_baq);
    let skip_indels = args.skip_indels;
    let annotate = AnnotateSpec::parse(args.annotate.as_deref())?;
    let flag_filters = FlagFilters::from_full(
        args.ef, args.df, args.if_, args.nf, args.rf.as_deref(), args.ff.as_deref(),
    )?;
    let gap_frac = preset.as_ref().and_then(|p| p.gap_frac).unwrap_or(args.gap_frac);
    let min_ireads = args.min_ireads;
    let sample_filter = parse_samples_filter(args.samples.as_deref(), args.samples_file.as_deref())?;
    let fasta = args.fasta_ref.as_ref().map(load_fasta).transpose()?;
    let timing = std::env::var("KIRA_BT_TIMING").is_ok();

    let mut combined_filter: Option<crate::bam::pos_filter::InterestingMap> = None;
    for r in bams.iter_mut() {
        if !flag_filters.exclude_flags == 0 || flag_filters.require_flags != 0 {
            r.records_buf.retain(|lr| flag_filters.passes(lr.flags));
        }
        let _ = max_depth;
        if !baq_off {
            if let Some(fa) = &fasta {
                let ref_names = r.ref_names.clone();
                let t_pre = std::time::Instant::now();
                let pre_arc = std::sync::Arc::new(std::mem::take(&mut r.records_buf));
                let pre = crate::bam::pos_filter::pre_scan(
                    std::slice::from_ref(&pre_arc), fa, &ref_names, args.min_alt_reads.max(2),
                );
                r.records_buf = std::sync::Arc::try_unwrap(pre_arc).unwrap_or_else(|arc| (*arc).clone());
                if timing {
                    eprintln!("[KIRA_BT] pre-scan: {:.1}s, skip-list={}", t_pre.elapsed().as_secs_f64(), pre.pos_filter.total());
                }
                let t1 = std::time::Instant::now();
                crate::bam::reader::apply_hmm_baq_to_reads_masked(&mut r.records_buf, &ref_names, fa, Some(&pre.needs_baq));
                if timing { eprintln!("[KIRA_BT] BAQ: {:.1}s", t1.elapsed().as_secs_f64()); }
                if combined_filter.is_none() {
                    combined_filter = Some(pre.pos_filter);
                }
            }
        }
    }

    let mut samples: Vec<String> = Vec::with_capacity(bams.len());
    let mut sample_keep: Vec<bool> = Vec::with_capacity(bams.len());
    for (i, b) in bams.iter().enumerate() {
        let s = b.samples.first().cloned().unwrap_or_else(|| format!("sample{i}"));
        let keep = sample_filter.as_ref().map_or(true, |list| list.iter().any(|n| n == &s));
        sample_keep.push(keep);
        samples.push(s);
    }
    let any_filter = sample_filter.is_some();
    let ref_names: Vec<String> = bams[0].ref_names.clone();

    write_vcf_header(out, &annotate, &ref_names, &samples, &sample_keep, any_filter)?;

    let em = ErrorModel::new();
    let gvcf_opts: Option<Vec<u32>> = args
        .gvcf
        .as_ref()
        .map(|s| s.split(',').filter_map(|t| t.parse::<u32>().ok()).collect());
    let mut gvcf_blocker: Option<GvcfBlocker> = gvcf_opts.map(GvcfBlocker::new);

    let n_chunks = args.threads.max(1);
    let parallel_eligible = n_chunks > 1 && gvcf_blocker.is_none();

    let hp_indel = std::env::var("KIRA_HP_INDEL").map(|v| v != "0").unwrap_or(true);
    let ctx = EmitCtx {
        args,
        em: &em,
        annotate: &annotate,
        fasta: fasta.as_ref(),
        ref_names: &ref_names,
        sample_keep: &sample_keep,
        any_filter,
        skip_indels,
        min_ireads,
        gap_frac,
        hp_indel,
    };

    if parallel_eligible {
        run_parallel(bams, n_chunks, &ctx, min_mq, min_bq, skip_indels, combined_filter, out)?;
    } else {
        run_sequential(bams, &ctx, min_mq, min_bq, skip_indels, &mut gvcf_blocker, out)?;
        if let Some(mut g) = gvcf_blocker {
            g.flush(out)?;
        }
    }
    out.flush()?;
    Ok(())
}

fn write_vcf_header<W: Write>(
    out: &mut W,
    annotate: &AnnotateSpec,
    ref_names: &[String],
    samples: &[String],
    sample_keep: &[bool],
    any_filter: bool,
) -> Result<()> {
    writeln!(out, "##fileformat=VCFv4.2")?;
    writeln!(out, "##source=kira_bt_mpileup")?;
    writeln!(out, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">")?;
    if annotate.info_ad { writeln!(out, "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Total allelic depths\">")?; }
    if annotate.info_adf { writeln!(out, "##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Forward strand allelic depths\">")?; }
    if annotate.info_adr { writeln!(out, "##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Reverse strand allelic depths\">")?; }
    writeln!(out, "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"Mean MAPQ at site\">")?;
    writeln!(out, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alt allele counts\">")?;
    writeln!(out, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total alleles\">")?;
    writeln!(out, "##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Variant is an indel\">")?;
    writeln!(out, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")?;
    if annotate.fmt_dp { writeln!(out, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Sample depth\">")?; }
    if annotate.fmt_ad { writeln!(out, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">")?; }
    if annotate.fmt_adf { writeln!(out, "##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Forward-strand allelic depths\">")?; }
    if annotate.fmt_adr { writeln!(out, "##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Reverse-strand allelic depths\">")?; }
    if annotate.fmt_qs { writeln!(out, "##FORMAT=<ID=QS,Number=R,Type=Integer,Description=\"Quality sum per allele\">")?; }
    if annotate.fmt_sp { writeln!(out, "##FORMAT=<ID=SP,Number=1,Type=Integer,Description=\"Strand-bias P-value (Phred)\">")?; }
    if annotate.fmt_scr { writeln!(out, "##FORMAT=<ID=SCR,Number=1,Type=Integer,Description=\"Soft-clip read count\">")?; }
    if annotate.fmt_pl { writeln!(out, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled likelihoods\">")?; }
    if annotate.fmt_gq { writeln!(out, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">")?; }
    for rn in ref_names { writeln!(out, "##contig=<ID={}>", rn)?; }
    write!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    for (i, s) in samples.iter().enumerate() {
        if !any_filter || sample_keep[i] { write!(out, "\t{}", s)?; }
    }
    writeln!(out)?;
    Ok(())
}

struct EmitCtx<'a> {
    args: &'a MpileupArgs,
    em: &'a ErrorModel,
    annotate: &'a AnnotateSpec,
    fasta: Option<&'a Fasta>,
    ref_names: &'a [String],
    sample_keep: &'a [bool],
    any_filter: bool,
    skip_indels: bool,
    min_ireads: u32,
    gap_frac: f64,
    /// Homopolymer/STR-aware stringency for short indels (default on; KIRA_HP_INDEL=0 disables).
    hp_indel: bool,
}

fn run_sequential(
    mut bams: Vec<BamReader>,
    ctx: &EmitCtx,
    min_mq: u8,
    min_bq: u8,
    skip_indels: bool,
    gvcf_blocker: &mut Option<GvcfBlocker>,
    out: &mut BufWriter<File>,
) -> Result<()> {
    mpileup_engine_multi(&mut bams, min_mq, min_bq, skip_indels, &mut |site, overlapping| {
        let mut buf: Vec<u8> = Vec::with_capacity(256);
        if let Some(g) = gvcf_blocker.as_mut() {
            if !emit_site_or_gvcf(site, overlapping, ctx, &mut buf, Some(g), out) {
                return;
            }
        } else {
            emit_site(site, overlapping, ctx, &mut buf);
        }
        let _ = out.write_all(&buf);
    })?;
    Ok(())
}

fn run_parallel(
    mut bams: Vec<BamReader>,
    n_chunks: usize,
    ctx: &EmitCtx,
    min_mq: u8,
    min_bq: u8,
    skip_indels: bool,
    pos_filter: Option<crate::bam::pos_filter::InterestingMap>,
    out: &mut BufWriter<File>,
) -> Result<()> {
    use rayon::prelude::*;

    let n_samples = bams.len();
    let records_per_sample: Vec<std::sync::Arc<Vec<LiveRead>>> = bams
        .iter_mut()
        .enumerate()
        .map(|(i, b)| {
            let mut v = std::mem::take(&mut b.records_buf);
            for lr in v.iter_mut() { lr.sample_idx = i; }
            std::sync::Arc::new(v)
        })
        .collect();

    let chunks = compute_chunks(&records_per_sample, n_chunks);
    if chunks.is_empty() {
        return Ok(());
    }

    // pos_filter was built during pre-scan in cmd_mpileup — reuse it (no duplicate work).
    // Fallback: build now if it wasn't (e.g. BAQ-off mode).
    let pos_filter = pos_filter.or_else(|| {
        ctx.fasta.map(|fa| {
            let t = std::time::Instant::now();
            let m = crate::bam::pos_filter::build(&records_per_sample, fa, ctx.ref_names);
            if std::env::var("KIRA_BT_TIMING").is_ok() {
                eprintln!("[KIRA_BT] late skip-list: {} in {:.1}s", m.total(), t.elapsed().as_secs_f64());
            }
            m
        })
    });

    use indicatif::{ProgressBar, ProgressStyle};
    let pb = ProgressBar::new(chunks.len() as u64);
    pb.set_style(
        ProgressStyle::with_template("[ENGINE] {bar:40.magenta/blue} {pos}/{len} chunks ({per_sec}, ETA {eta})")
            .unwrap_or_else(|_| ProgressStyle::default_bar()),
    );
    pb.set_draw_target(indicatif::ProgressDrawTarget::stderr_with_hz(2));

    // BAM records are sorted by (ref_id, ref_start). For each chunk's [start, end) window we
    // need records whose ref_id matches AND that overlap the window. Using partition_point we
    // get a tight slice in O(log N) instead of scanning all 12M reads N times. The left
    // boundary backs up by MAX_READ_LEN so we don't miss reads that started just before the
    // chunk but extend into it.
    const MAX_READ_LEN_FLANK: u32 = 1024;
    let buffers: Vec<Vec<u8>> = chunks
        .par_iter()
        .map(|chunk| {
            let mut buf: Vec<u8> = Vec::with_capacity(1 << 20);
            let lo_pos = chunk.start.saturating_sub(MAX_READ_LEN_FLANK);
            let per_sample: Vec<Vec<LiveRead>> = records_per_sample
                .iter()
                .map(|all| {
                    let lo = all.partition_point(|lr| {
                        (lr.ref_id, lr.ref_start) < (chunk.ref_id, lo_pos)
                    });
                    let hi = all.partition_point(|lr| {
                        (lr.ref_id, lr.ref_start) < (chunk.ref_id, chunk.end)
                    });
                    all[lo..hi]
                        .iter()
                        .filter(|lr| lr.ref_id == chunk.ref_id && lr.ref_end_cached > chunk.start)
                        .cloned()
                        .collect()
                })
                .collect();
            // Engine skips ref-only positions via pos_filter, and we still bound emit to the chunk window.
            let _ = mpileup_engine_from_records(
                per_sample,
                min_mq,
                min_bq,
                skip_indels,
                pos_filter.as_ref(),
                &mut |site, overlapping| {
                    if site.ref_id != chunk.ref_id { return; }
                    if site.pos < chunk.start || site.pos >= chunk.end { return; }
                    emit_site(site, overlapping, ctx, &mut buf);
                },
            );
            let _ = n_samples;
            pb.inc(1);
            buf
        })
        .collect();
    pb.finish_with_message("ENGINE done");

    for buf in buffers {
        out.write_all(&buf)?;
    }
    Ok(())
}

#[derive(Clone, Copy, Debug)]
struct ChunkSpec {
    ref_id: usize,
    start: u32,
    end: u32,
}

fn compute_chunks(records_per_sample: &[std::sync::Arc<Vec<LiveRead>>], n_chunks: usize) -> Vec<ChunkSpec> {
    // Group positions per ref_id across all samples.
    let mut per_ref: FxHashMap<usize, (u32, u32)> = FxHashMap::default();
    for sample in records_per_sample {
        for lr in sample.iter() {
            let entry = per_ref.entry(lr.ref_id).or_insert((u32::MAX, 0));
            entry.0 = entry.0.min(lr.ref_start);
            entry.1 = entry.1.max(lr.ref_end());
        }
    }
    if per_ref.is_empty() {
        return Vec::new();
    }
    let mut ref_ids: Vec<usize> = per_ref.keys().copied().collect();
    ref_ids.sort_unstable();

    // Distribute chunk budget proportionally to each ref's position span.
    let total_span: u64 = per_ref.values().map(|(a, b)| (*b as u64).saturating_sub(*a as u64).max(1)).sum();
    let mut chunks: Vec<ChunkSpec> = Vec::new();
    for rid in ref_ids {
        let (lo, hi) = per_ref[&rid];
        let span = (hi as u64).saturating_sub(lo as u64).max(1);
        let share = ((span * n_chunks as u64) / total_span).max(1) as usize;
        let step = (span / share as u64).max(1);
        for i in 0..share {
            let s = lo as u64 + i as u64 * step;
            let e = if i + 1 == share { hi as u64 } else { s + step };
            chunks.push(ChunkSpec {
                ref_id: rid,
                start: s as u32,
                end: e as u32,
            });
        }
    }
    chunks
}

/// Per-site emit, writing the formatted VCF line into `out` (a Vec<u8> buffer
/// in the parallel runner; a BufWriter slice in the sequential one).
fn emit_site(
    site: &PileupSite,
    overlapping_reads: &[&LiveRead],
    ctx: &EmitCtx,
    out: &mut Vec<u8>,
) {
    let args = ctx.args;
    let annotate = ctx.annotate;
    let chr = ctx.ref_names.get(site.ref_id).map(|s| s.as_str()).unwrap_or(".");
    let pos1 = site.pos as u64 + 1;
    let ref_base = ctx.fasta.and_then(|fa| fa.base(chr, pos1 as u32)).unwrap_or(b'N');
    if ref_base == b'N' || ref_base == b'n' { return; }
    let ref_idx = match ref_base { b'A' | b'a' => 0, b'C' | b'c' => 1, b'G' | b'g' => 2, b'T' | b't' => 3, _ => 4 };
    let agg = site.aggregated();
    let total: u32 = agg.base_counts.iter().sum();
    if total == 0 { return; }

    let (alt_min_reads, alt_min_af) = if args.variants_only {
        (args.min_alt_reads, args.min_af)
    } else {
        (1u32, 0.0f64)
    };
    let total_f = total.max(1) as f64;
    let mut snv_alts: Vec<(u8, u32, u32)> = Vec::new();
    for i in 0..4 {
        if i == ref_idx { continue; }
        let c = agg.base_counts[i];
        let af = (c as f64) / total_f;
        if c >= alt_min_reads && af >= alt_min_af {
            snv_alts.push((b"ACGT"[i], c, agg.base_quals[i]));
        }
    }
    snv_alts.sort_by(|a, b| b.1.cmp(&a.1));

    let all_alts: Vec<String> = snv_alts.iter().map(|(b, _, _)| (*b as char).to_string()).collect();

    let alt_str = if all_alts.is_empty() { ".".into() } else { all_alts.join(",") };
    let ref_for_record = if ref_base == b'N' && !snv_alts.is_empty() { snv_alts[0].0 } else { ref_base };
    let n_alleles = all_alts.len() + 1;
    let mean_mq = if total > 0 { (agg.mq_sum as f64) / (total as f64) } else { 0.0 };

    let mut per_sample_cols: Vec<String> = Vec::new();
    let mut site_qual: u8 = 0;
    let mut any_variant = false;
    for (si, s) in site.per_sample.iter().enumerate() {
        if ctx.any_filter && !ctx.sample_keep.get(si).copied().unwrap_or(false) { continue; }
        if s.depth == 0 {
            let mut parts: Vec<String> = vec!["./.".into()];
            if annotate.fmt_dp { parts.push("0".into()); }
            if annotate.fmt_ad { parts.push(vec!["0"; n_alleles].join(",")); }
            if annotate.fmt_pl { parts.push("0".into()); }
            if annotate.fmt_gq { parts.push("0".into()); }
            per_sample_cols.push(parts.join(":"));
            continue;
        }
        let mut counts: Vec<u32> = vec![0; n_alleles];
        let mut quals: Vec<u32> = vec![0; n_alleles];
        if ref_idx < 4 { counts[0] = s.base_counts[ref_idx]; quals[0] = s.base_quals[ref_idx]; }
        for (i, (b, _, _)) in snv_alts.iter().enumerate() {
            let bi = match *b { b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3, _ => continue };
            counts[i + 1] = s.base_counts[bi];
            quals[i + 1] = s.base_quals[bi];
        }
        let gl_raw = ctx.em.likelihoods(n_alleles, &counts, &quals);
        let gl = gl_raw.with_prior(n_alleles, args.prior);
        let (gi, gj) = gl.most_likely_gt(n_alleles);
        if gi != 0 || gj != 0 { any_variant = true; }
        site_qual = site_qual.max(gl.qual());
        let gt = format!("{}/{}", gi, gj);
        let mut parts: Vec<String> = vec![gt];
        if annotate.fmt_dp { parts.push(s.depth.to_string()); }
        if annotate.fmt_ad { parts.push(counts.iter().map(u32::to_string).collect::<Vec<_>>().join(",")); }
        if annotate.fmt_qs { parts.push(quals.iter().map(u32::to_string).collect::<Vec<_>>().join(",")); }
        if annotate.fmt_pl { parts.push(gl.to_pl_string()); }
        if annotate.fmt_gq {
            let mut sorted = gl.pl.clone(); sorted.sort();
            let gq = sorted.get(1).copied().unwrap_or(0);
            parts.push(gq.to_string());
        }
        per_sample_cols.push(parts.join(":"));
    }
    let have_indel = !ctx.skip_indels && {
        let m = if args.variants_only { args.min_alt_reads.max(ctx.min_ireads) } else { ctx.min_ireads };
        agg.ins_alleles.iter().any(|(_, c)| *c >= m) || agg.del_alleles.iter().any(|(_, c)| *c >= m)
    };
    if args.variants_only && (!any_variant || site_qual < args.min_qual) && !have_indel {
        return;
    }

    let ac: u32 = snv_alts.iter().map(|(_, c, _)| *c).sum();
    let an = total + ac;
    let mut info = format!("DP={};MQ={:.1};AC={};AN={}", total, mean_mq, ac, an);
    if annotate.info_ad {
        let mut ad = vec![0u32; n_alleles];
        if ref_idx < 4 { ad[0] = agg.base_counts[ref_idx]; }
        for (i, (_, c, _)) in snv_alts.iter().enumerate() { ad[i + 1] = *c; }
        info.push_str(&format!(";AD={}", ad.iter().map(u32::to_string).collect::<Vec<_>>().join(",")));
    }

    let mut fmt_keys = vec!["GT"];
    if annotate.fmt_dp { fmt_keys.push("DP"); }
    if annotate.fmt_ad { fmt_keys.push("AD"); }
    if annotate.fmt_qs { fmt_keys.push("QS"); }
    if annotate.fmt_pl { fmt_keys.push("PL"); }
    if annotate.fmt_gq { fmt_keys.push("GQ"); }

    let snv_emit = !snv_alts.is_empty() || !args.variants_only;
    if snv_emit && !(args.variants_only && (!any_variant || site_qual < args.min_qual)) {
        let qual_str = if site_qual == 0 { ".".to_string() } else { site_qual.to_string() };
        let _ = write!(out, "{}\t{}\t.\t{}\t{}\t{}\t.\t{}\t{}",
            chr, pos1, ref_for_record as char, alt_str, qual_str, info, fmt_keys.join(":"));
        for col in &per_sample_cols { let _ = write!(out, "\t{}", col); }
        let _ = writeln!(out);
    }

    if !ctx.skip_indels {
        emit_indels(site, &agg, chr, pos1, ref_base, total, mean_mq, ctx, out, overlapping_reads);
    }

    // Local indel recovery: at non-reference sites, re-assemble the reads to find
    // indels the aligner hid as mismatches. Gated by KIRA_ASSEMBLE; sorts downstream.
    // Only where the pileup found NO indel of its own (`!have_indel`) — the
    // mismatch-modeled FN — so we add new recall, not duplicates of CIGAR indels.
    if assemble_mode() && !ctx.skip_indels && !have_indel {
        let cfg = asm_cfg();
        let nonref = total.saturating_sub(if ref_idx < 4 { agg.base_counts[ref_idx] } else { 0 });
        if nonref >= cfg.nonref {
            if let Some(fa) = ctx.fasta {
                let pos0 = site.pos as u32;
                let w_lo = pos0.saturating_sub(cfg.up);
                if let Some(refw) = fa.slice_bytes(chr, w_lo + 1, cfg.len) {
                    if let Some(call) =
                        crate::call::haplotype::assemble_indel(overlapping_reads, w_lo, refw, cfg.min_sup, cfg.max_mm)
                    {
                        let n_cols = site
                            .per_sample
                            .iter()
                            .enumerate()
                            .filter(|(si, _)| !ctx.any_filter || ctx.sample_keep.get(*si).copied().unwrap_or(false))
                            .count();
                        emit_assembled(call, n_cols, chr, ctx, out, annotate);
                    }
                }
            }
        }
    }
}

/// `KIRA_INDEL_REALIGN`: 0/off, 1/ins/on = insertions only, 2/all = insertions+deletions (default off).
fn realign_mode() -> u8 {
    static V: std::sync::OnceLock<u8> = std::sync::OnceLock::new();
    *V.get_or_init(|| match std::env::var("KIRA_INDEL_REALIGN") {
        Ok(v) if v == "2" || v.eq_ignore_ascii_case("all") => 2,
        Ok(v)
            if v == "1"
                || v.eq_ignore_ascii_case("ins")
                || v.eq_ignore_ascii_case("on")
                || v.eq_ignore_ascii_case("true") =>
        {
            1
        }
        _ => 0,
    })
}

/// Realignment flank window (bp each side). `KIRA_REALIGN_W`, default 20.
fn realign_window() -> u32 {
    static V: std::sync::OnceLock<u32> = std::sync::OnceLock::new();
    *V.get_or_init(|| {
        std::env::var("KIRA_REALIGN_W").ok().and_then(|s| s.parse().ok()).filter(|&w: &u32| w >= 4).unwrap_or(20)
    })
}

/// `KIRA_REALIGN_MARGIN` — a read supports the indel iff `d_alt + margin <= d_ref` (default 1).
fn realign_margin() -> i32 {
    static V: std::sync::OnceLock<i32> = std::sync::OnceLock::new();
    *V.get_or_init(|| {
        std::env::var("KIRA_REALIGN_MARGIN").ok().and_then(|s| s.parse().ok()).unwrap_or(1)
    })
}

/// Capped Levenshtein (edit) distance, case-insensitive. Returns `cap+1` once it provably exceeds `cap`.
fn edit_dist_capped(a: &[u8], b: &[u8], cap: u32) -> u32 {
    let (n, m) = (a.len(), b.len());
    if n == 0 {
        return (m as u32).min(cap + 1);
    }
    if m == 0 {
        return (n as u32).min(cap + 1);
    }
    if (n as i64 - m as i64).unsigned_abs() as u32 > cap {
        return cap + 1;
    }
    let mut prev: Vec<u32> = (0..=m as u32).collect();
    let mut cur = vec![0u32; m + 1];
    for i in 1..=n {
        cur[0] = i as u32;
        let mut row_min = cur[0];
        for j in 1..=m {
            let cost = if a[i - 1].eq_ignore_ascii_case(&b[j - 1]) { 0 } else { 1 };
            cur[j] = (prev[j] + 1).min(cur[j - 1] + 1).min(prev[j - 1] + cost);
            row_min = row_min.min(cur[j]);
        }
        if row_min > cap {
            return cap + 1;
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    prev[m]
}

/// Count reads whose window realigns better to ref+ins than ref; `u32::MAX` if realignment impossible.
fn realign_ins_support(reads: &[&LiveRead], chr: &str, pos1: u64, ins: &[u8], ctx: &EmitCtx) -> u32 {
    let w = realign_window();
    let margin = realign_margin();
    let Some(fa) = ctx.fasta else { return u32::MAX };
    let Some(anchor0) = (pos1 as u32).checked_sub(1) else { return u32::MAX };
    let lo = anchor0.saturating_sub(w);
    let hi = anchor0 + 1 + w; // [lo, hi) reference window, 0-based
    let wlen = (hi - lo) as usize;
    let Some(refbytes) = fa.slice_bytes(chr, lo + 1, wlen) else { return u32::MAX };
    if refbytes.len() < wlen {
        return u32::MAX;
    }
    let anchor_off = (anchor0 + 1 - lo) as usize; // index just past the anchor base
    if anchor_off > wlen {
        return u32::MAX;
    }
    let href: &[u8] = &refbytes[..wlen];
    let mut halt: Vec<u8> = Vec::with_capacity(wlen + ins.len());
    halt.extend_from_slice(&refbytes[..anchor_off]);
    halt.extend_from_slice(ins);
    halt.extend_from_slice(&refbytes[anchor_off..wlen]);
    let cap = 40u32;
    let mut genuine = 0u32;
    for lr in reads {
        let Some(sub) = lr.query_window(lo, hi) else { continue };
        if sub.is_empty() {
            continue;
        }
        let d_ref = edit_dist_capped(&sub, href, cap) as i32;
        let d_alt = edit_dist_capped(&sub, &halt, cap) as i32;
        if d_alt + margin <= d_ref {
            genuine += 1;
        }
    }
    genuine
}

/// Count reads whose window realigns better to ref-with-deletion than ref; `u32::MAX` if impossible.
fn realign_del_support(reads: &[&LiveRead], chr: &str, pos1: u64, del_len: u32, ctx: &EmitCtx) -> u32 {
    let w = realign_window();
    let margin = realign_margin();
    let Some(fa) = ctx.fasta else { return u32::MAX };
    let Some(anchor0) = (pos1 as u32).checked_sub(1) else { return u32::MAX };
    let lo = anchor0.saturating_sub(w);
    let ext_hi = anchor0 + 1 + w + del_len; // extend right by del_len for post-deletion context
    let wlen = (ext_hi - lo) as usize;
    let Some(refbytes) = fa.slice_bytes(chr, lo + 1, wlen) else { return u32::MAX };
    if refbytes.len() < wlen {
        return u32::MAX;
    }
    let anchor_off = (anchor0 + 1 - lo) as usize;
    let dl = del_len as usize;
    if anchor_off + dl > wlen {
        return u32::MAX;
    }
    let href: &[u8] = &refbytes[..wlen];
    // Halt = ref with `del_len` bases removed after the anchor (length wlen - dl).
    let mut halt: Vec<u8> = Vec::with_capacity(wlen - dl);
    halt.extend_from_slice(&refbytes[..anchor_off]);
    halt.extend_from_slice(&refbytes[anchor_off + dl..wlen]);
    let cap = 40u32;
    let mut genuine = 0u32;
    for lr in reads {
        let Some(sub) = lr.query_window(lo, ext_hi) else { continue };
        if sub.is_empty() {
            continue;
        }
        let d_ref = edit_dist_capped(&sub, href, cap) as i32;
        let d_alt = edit_dist_capped(&sub, &halt, cap) as i32;
        if d_alt + margin <= d_ref {
            genuine += 1;
        }
    }
    genuine
}

/// Emit one VCF record per indel allele (insertion or deletion). Each record
/// is fully normalized: REF includes the deleted bases (read from FASTA) for
/// deletions, REF is single-base for insertions.
fn emit_indels(
    site: &PileupSite,
    agg: &crate::bam::SampleSiteCounts,
    chr: &str,
    pos1: u64,
    ref_base: u8,
    total: u32,
    mean_mq: f64,
    ctx: &EmitCtx,
    out: &mut Vec<u8>,
    overlapping_reads: &[&LiveRead],
) {
    let args = ctx.args;
    let annotate = ctx.annotate;
    let (ind_min, ind_af) = if args.variants_only {
        (args.min_alt_reads.max(ctx.min_ireads), args.min_af.max(ctx.gap_frac))
    } else {
        (ctx.min_ireads, ctx.gap_frac)
    };
    let denom = total.max(1) as f64;

    // For each indel allele, emit a biallelic record (REF + this one ALT).
    for (seq, c) in &agg.ins_alleles {
        let (mn, af) = if ctx.hp_indel {
            indel_repeat_stringency(ctx.fasta, chr, pos1, seq.as_bytes(), ind_min, ind_af)
        } else { (ind_min, ind_af) };
        if *c < mn { continue; }
        let frac = (*c as f64) / denom;
        if frac < af { continue; }
        // Re-apply count/fraction thresholds to realignment-confirmed support.
        if realign_mode() >= 1 {
            let genuine = realign_ins_support(overlapping_reads, chr, pos1, seq.as_bytes(), ctx);
            if genuine != u32::MAX && (genuine < mn || (genuine as f64) / denom < af) {
                continue;
            }
        }
        let ref_str = (ref_base as char).to_string();
        let alt_str = format!("{}{}", ref_base as char, seq);
        emit_one_indel(site, agg, chr, pos1, &ref_str, &alt_str, *c, total, mean_mq, ctx, out, annotate);
    }
    for (l, c) in &agg.del_alleles {
        // REF = anchor + deleted bases from FASTA; ALT = anchor.
        let Some(fa) = ctx.fasta else { continue };
        let Some(deleted) = fa.slice_bytes(chr, pos1 as u32 + 1, *l as usize) else { continue };
        if deleted.iter().any(|&b| b == b'N' || b == b'n') { continue; }
        let (mn, af) = if ctx.hp_indel {
            indel_repeat_stringency(ctx.fasta, chr, pos1, deleted, ind_min, ind_af)
        } else { (ind_min, ind_af) };
        if *c < mn { continue; }
        let frac = (*c as f64) / denom;
        if frac < af { continue; }
        if realign_mode() >= 2 {
            let genuine = realign_del_support(overlapping_reads, chr, pos1, *l, ctx);
            if genuine != u32::MAX && (genuine < mn || (genuine as f64) / denom < af) {
                continue;
            }
        }
        let mut ref_str = String::with_capacity(1 + deleted.len());
        ref_str.push(ref_base as char);
        for &b in deleted { ref_str.push(b as char); }
        let alt_str = (ref_base as char).to_string();
        emit_one_indel(site, agg, chr, pos1, &ref_str, &alt_str, *c, total, mean_mq, ctx, out, annotate);
    }
}

fn emit_one_indel(
    site: &PileupSite,
    _agg: &crate::bam::SampleSiteCounts,
    chr: &str,
    pos1: u64,
    ref_str: &str,
    alt_str: &str,
    indel_count: u32,
    total: u32,
    mean_mq: f64,
    ctx: &EmitCtx,
    out: &mut Vec<u8>,
    annotate: &AnnotateSpec,
) {
    let args = ctx.args;
    // Aggregate per-sample for this specific indel allele. For each sample, build
    // a synthetic biallelic [ref, indel] counts/quals pair: ref = depth - indel,
    // indel = matched_count, qual proxy = 30 per supporting read.
    let n_alleles = 2usize;
    let mut per_sample_cols: Vec<String> = Vec::new();
    let mut site_qual: u8 = 0;
    let mut any_variant = false;

    // KIRA_INDEL_HPQUAL: off/0 = qual 30; "F,S" = floor F, slope S; 1/auto = 8,5 (default off).
    let indel_qual: u32 = {
        static CFG: std::sync::OnceLock<Option<(u32, u32)>> = std::sync::OnceLock::new();
        let cfg = *CFG.get_or_init(|| match std::env::var("KIRA_INDEL_HPQUAL") {
            Ok(v) if v == "0" || v.eq_ignore_ascii_case("off") => None,
            Ok(v) if v == "1" || v.eq_ignore_ascii_case("auto") => Some((8, 5)),
            Ok(v) => {
                let mut it = v.split(',').filter_map(|x| x.trim().parse::<u32>().ok());
                match (it.next(), it.next()) {
                    (Some(f), Some(s)) => Some((f, s)),
                    _ => None,
                }
            }
            Err(_) => None,
        });
        match cfg {
            None => 30,
            Some((floor, slope)) => {
                let unit: &[u8] = if alt_str.len() > ref_str.len() {
                    &alt_str.as_bytes()[1..]
                } else {
                    &ref_str.as_bytes()[1..]
                };
                let u = unit.len();
                let run_units = ctx
                    .fasta
                    .and_then(|fa| {
                        if u == 0 || u > 6 {
                            return None;
                        }
                        let c = fa.slice_bytes(chr, pos1 as u32 + 1, 64)?;
                        let mut units = 0u32;
                        let mut i = 0usize;
                        while i + u <= c.len() && c[i..i + u].eq_ignore_ascii_case(unit) {
                            units += 1;
                            i += u;
                        }
                        Some(units)
                    })
                    .unwrap_or(0);
                let penalty = slope.saturating_mul(run_units.saturating_sub(2));
                30u32.saturating_sub(penalty).max(floor)
            }
        }
    };

    for (si, s) in site.per_sample.iter().enumerate() {
        if ctx.any_filter && !ctx.sample_keep.get(si).copied().unwrap_or(false) { continue; }
        if s.depth == 0 {
            let mut parts: Vec<String> = vec!["./.".into()];
            if annotate.fmt_dp { parts.push("0".into()); }
            if annotate.fmt_ad { parts.push("0,0".into()); }
            if annotate.fmt_pl { parts.push("0".into()); }
            if annotate.fmt_gq { parts.push("0".into()); }
            per_sample_cols.push(parts.join(":"));
            continue;
        }
        // Indel count for THIS allele in THIS sample.
        let indel_c = if alt_str.len() > ref_str.len() {
            // insertion: seq = alt_str[1..]
            let seq = &alt_str[1..];
            s.ins_alleles.iter().find(|(k, _)| k == seq).map(|(_, c)| *c).unwrap_or(0)
        } else {
            // deletion: length = ref_str.len() - 1
            let dl = (ref_str.len() - 1) as u32;
            s.del_alleles.iter().find(|(l, _)| *l == dl).map(|(_, c)| *c).unwrap_or(0)
        };
        let ref_c = s.depth.saturating_sub(indel_c);
        let counts = vec![ref_c, indel_c];
        let quals = vec![ref_c * 30, indel_c * indel_qual];

        let gl_raw = ctx.em.likelihoods(n_alleles, &counts, &quals);
        let gl = gl_raw.with_prior(n_alleles, args.prior);
        let (gi, gj) = gl.most_likely_gt(n_alleles);
        if gi != 0 || gj != 0 { any_variant = true; }
        site_qual = site_qual.max(gl.qual());
        let gt = format!("{}/{}", gi, gj);
        let mut parts: Vec<String> = vec![gt];
        if annotate.fmt_dp { parts.push(s.depth.to_string()); }
        if annotate.fmt_ad { parts.push(format!("{},{}", ref_c, indel_c)); }
        if annotate.fmt_qs { parts.push(format!("{},{}", ref_c * 30, indel_c * 30)); }
        if annotate.fmt_pl { parts.push(gl.to_pl_string()); }
        if annotate.fmt_gq {
            let mut sorted = gl.pl.clone(); sorted.sort();
            let gq = sorted.get(1).copied().unwrap_or(0);
            parts.push(gq.to_string());
        }
        per_sample_cols.push(parts.join(":"));
    }
    if args.variants_only && (!any_variant || site_qual < args.min_qual) {
        return;
    }

    let an = total + indel_count;
    let mut info = format!("DP={};MQ={:.1};AC={};AN={};INDEL", total, mean_mq, indel_count, an);
    if annotate.info_ad {
        let ref_total = total.saturating_sub(indel_count);
        info.push_str(&format!(";AD={},{}", ref_total, indel_count));
    }

    let mut fmt_keys = vec!["GT"];
    if annotate.fmt_dp { fmt_keys.push("DP"); }
    if annotate.fmt_ad { fmt_keys.push("AD"); }
    if annotate.fmt_qs { fmt_keys.push("QS"); }
    if annotate.fmt_pl { fmt_keys.push("PL"); }
    if annotate.fmt_gq { fmt_keys.push("GQ"); }

    let qual_str = if site_qual == 0 { ".".to_string() } else { site_qual.to_string() };
    let _ = write!(out, "{}\t{}\t.\t{}\t{}\t{}\t.\t{}\t{}",
        chr, pos1, ref_str, alt_str, qual_str, info, fmt_keys.join(":"));
    for col in &per_sample_cols { let _ = write!(out, "\t{}", col); }
    let _ = writeln!(out);
}

/// `KIRA_ASSEMBLE=1` enables local indel recovery (mm-filter assembler) at
/// non-reference sites — recovers indels the aligner modelled as mismatches and
/// never placed in a CIGAR. Records may land slightly out of position order, so
/// the output must be sorted downstream.
fn assemble_mode() -> bool {
    static V: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *V.get_or_init(|| matches!(std::env::var("KIRA_ASSEMBLE"),
        Ok(v) if v == "1" || v.eq_ignore_ascii_case("on") || v.eq_ignore_ascii_case("yes")))
}

struct AsmCfg { up: u32, len: usize, nonref: u32, min_sup: u32, max_mm: u32 }
/// Assembler window/gate, env-tunable. The mismatch-modeled indel's anchor sits
/// UPSTREAM of the non-ref trigger, so `up` must reach back far enough to cover it.
fn asm_cfg() -> &'static AsmCfg {
    static C: std::sync::OnceLock<AsmCfg> = std::sync::OnceLock::new();
    C.get_or_init(|| {
        let g = |k: &str, d: u32| std::env::var(k).ok().and_then(|v| v.parse().ok()).unwrap_or(d);
        AsmCfg {
            up: g("KIRA_ASM_UP", 6),
            len: g("KIRA_ASM_LEN", 57) as usize,
            nonref: g("KIRA_ASM_NONREF", 3),
            min_sup: g("KIRA_ASM_SUP", 3),
            max_mm: g("KIRA_ASM_MM", 3),
        }
    })
}

thread_local! {
    static ASM_EMITTED: std::cell::RefCell<std::collections::HashSet<(u64, String, String)>> =
        std::cell::RefCell::new(std::collections::HashSet::new());
}

fn emit_assembled(
    call: crate::call::haplotype::AssembledCall,
    n_cols: usize,
    chr: &str,
    ctx: &EmitCtx,
    out: &mut Vec<u8>,
    annotate: &AnnotateSpec,
) {
    let args = ctx.args;
    let ref_c = call.total.saturating_sub(call.support);
    let alt_c = call.support;
    let gl = ctx.em.likelihoods(2, &[ref_c, alt_c], &[ref_c * 30, alt_c * 30]).with_prior(2, args.prior);
    let (gi, gj) = gl.most_likely_gt(2);
    let qual = gl.qual();
    if args.variants_only && ((gi == 0 && gj == 0) || qual < args.min_qual) { return; }
    let key = (call.pos1, call.ref_str.clone(), call.alt_str.clone());
    let dup = ASM_EMITTED.with(|e| {
        let mut s = e.borrow_mut();
        if s.contains(&key) { return true; }
        if s.len() > 200_000 { s.clear(); }
        s.insert(key);
        false
    });
    if dup { return; }
    let mut info = format!("DP={};AC={};AN={};INDEL;ASSEMBLED",
        call.total, if gi == 0 && gj == 0 { 0 } else { alt_c }, call.total);
    if annotate.info_ad { info.push_str(&format!(";AD={},{}", ref_c, alt_c)); }
    let mut fmt_keys = vec!["GT"];
    if annotate.fmt_dp { fmt_keys.push("DP"); }
    if annotate.fmt_ad { fmt_keys.push("AD"); }
    if annotate.fmt_qs { fmt_keys.push("QS"); }
    if annotate.fmt_pl { fmt_keys.push("PL"); }
    if annotate.fmt_gq { fmt_keys.push("GQ"); }
    let mut p0: Vec<String> = vec![format!("{}/{}", gi, gj)];
    if annotate.fmt_dp { p0.push(call.total.to_string()); }
    if annotate.fmt_ad { p0.push(format!("{},{}", ref_c, alt_c)); }
    if annotate.fmt_qs { p0.push(format!("{},{}", ref_c * 30, alt_c * 30)); }
    if annotate.fmt_pl { p0.push(gl.to_pl_string()); }
    if annotate.fmt_gq { let mut s = gl.pl.clone(); s.sort(); p0.push(s.get(1).copied().unwrap_or(0).to_string()); }
    let col0 = p0.join(":");
    let mut pe: Vec<String> = vec!["./.".into()];
    if annotate.fmt_dp { pe.push("0".into()); }
    if annotate.fmt_ad { pe.push("0,0".into()); }
    if annotate.fmt_qs { pe.push("0,0".into()); }
    if annotate.fmt_pl { pe.push("0".into()); }
    if annotate.fmt_gq { pe.push("0".into()); }
    let cole = pe.join(":");
    let qual_str = if qual == 0 { ".".to_string() } else { qual.to_string() };
    let _ = write!(out, "{}\t{}\t.\t{}\t{}\t{}\t.\t{}\t{}",
        chr, call.pos1, call.ref_str, call.alt_str, qual_str, info, fmt_keys.join(":"));
    for si in 0..n_cols.max(1) {
        let _ = write!(out, "\t{}", if si == 0 { &col0 } else { &cole });
    }
    let _ = writeln!(out);
}

/// Wrapper for the sequential path when gvcf blocking is active. Returns false
/// when the site was emitted as a gVCF ref-block (caller should not append the
/// SNV record).
fn emit_site_or_gvcf<W: Write>(
    site: &PileupSite,
    overlapping_reads: &[&LiveRead],
    ctx: &EmitCtx,
    buf: &mut Vec<u8>,
    gvcf: Option<&mut GvcfBlocker>,
    direct_out: &mut W,
) -> bool {
    let chr = ctx.ref_names.get(site.ref_id).map(|s| s.as_str()).unwrap_or(".");
    let pos1 = site.pos as u64 + 1;
    let ref_base = ctx.fasta.and_then(|fa| fa.base(chr, pos1 as u32)).unwrap_or(b'N');
    if ref_base == b'N' { return false; }
    let agg = site.aggregated();
    let total: u32 = agg.base_counts.iter().sum();
    if total == 0 { return false; }
    // For gVCF, branch on whether any alts were observed at all.
    let any_nonref = (0..4).filter(|&i| {
        let rb = match ref_base { b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3, _ => 4 };
        i != rb && agg.base_counts[i] > 0
    }).count() > 0;
    if let Some(g) = gvcf {
        if !any_nonref {
            let _ = g.add_ref_site(chr, pos1, &(ref_base as char).to_string(), total, 0.0, direct_out);
            return false;
        }
        let _ = g.flush(direct_out);
    }
    emit_site(site, overlapping_reads, ctx, buf);
    true
}

fn cap_depth(reads: &mut Vec<crate::bam::pileup::LiveRead>, max_depth: u32) {
    if max_depth == 0 || reads.len() <= max_depth as usize { return; }
    let stride = reads.len() / max_depth as usize;
    let kept: Vec<_> = reads.iter().step_by(stride.max(1)).cloned().collect();
    *reads = kept;
}

struct Fasta { seqs: FxHashMap<String, Vec<u8>> }
impl Fasta {
    fn base(&self, chr: &str, pos: u32) -> Option<u8> {
        self.seqs.get(chr).and_then(|s| s.get((pos as usize).saturating_sub(1)).copied())
    }
    fn slice_bytes(&self, chr: &str, pos: u32, len: usize) -> Option<&[u8]> {
        let s = self.seqs.get(chr)?;
        let start = (pos as usize).saturating_sub(1);
        let end = (start + len).min(s.len());
        if end <= start { return None; }
        Some(&s[start..end])
    }
}
impl crate::bam::reader::FastaLike for Fasta {
    fn slice(&self, chr: &str, pos: u32, len: usize) -> Option<&[u8]> { self.slice_bytes(chr, pos, len) }
}

/// Repeat-context stringency for short indels. Homopolymer/STR tracts are the dominant source
/// of false-positive indels (polymerase slippage), so for a short inserted/deleted `unit` that
/// tiles a long repeat run immediately 3' of the anchor `pos1`, require more supporting reads and
/// a higher alt fraction. Returns the (possibly tightened) (min_count, min_frac). No-op for
/// non-repeat or long indels. This is what a mature caller (bcftools) does via HMM realignment;
/// here we approximate it at emit time using the reference context.
fn indel_repeat_stringency(
    fa: Option<&Fasta>, chr: &str, pos1: u64, unit: &[u8], base_min: u32, base_af: f64,
) -> (u32, f64) {
    let u = unit.len();
    if u == 0 || u > 4 { return (base_min, base_af); }
    let fa = match fa { Some(f) => f, None => return (base_min, base_af) };
    let ctx = match fa.slice_bytes(chr, pos1 as u32 + 1, 64) { Some(c) => c, None => return (base_min, base_af) };
    let mut bases = 0u32;
    let mut i = 0usize;
    while i + u <= ctx.len() && ctx[i..i + u].eq_ignore_ascii_case(unit) { bases += u as u32; i += u; }
    let copies = bases / (u as u32);
    let strict = (u == 1 && bases >= 5) || (u >= 2 && copies >= 4);
    if strict { (base_min.max(4), base_af + 0.12) } else { (base_min, base_af) }
}

fn load_fasta(p: &PathBuf) -> Result<Fasta> {
    let mut seqs: FxHashMap<String, Vec<u8>> = FxHashMap::default();
    let data = std::fs::read(p).with_context(|| format!("open fasta {:?}", p))?;
    let mut name: Option<String> = None;
    let mut cur: Vec<u8> = Vec::new();
    for line in data.split(|&b| b == b'\n') {
        if line.is_empty() { continue; }
        if line[0] == b'>' {
            if let Some(n) = name.take() { seqs.insert(n, std::mem::take(&mut cur)); }
            let rest = &line[1..];
            let end = rest.iter().position(|&b| b == b' ' || b == b'\t' || b == b'\r').unwrap_or(rest.len());
            name = Some(std::str::from_utf8(&rest[..end])?.to_string());
        } else {
            for &b in line { if b != b'\r' { cur.push(b.to_ascii_uppercase()); } }
        }
    }
    if let Some(n) = name { seqs.insert(n, cur); }
    Ok(Fasta { seqs })
}
