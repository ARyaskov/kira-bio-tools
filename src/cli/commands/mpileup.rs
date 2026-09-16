use crate::bam::bcf_layout::{self, IndelSpec, LayoutCtx};
use crate::bam::pileup::{EngineOpts, IndelEvent, SiteHook, strand_bias_phred};
use noodles_sam::alignment::record::cigar::op::Kind;
use crate::bam::pos_filter::InterestingMap;
use crate::bam::samples::{SampleMap, SampleOpts};
use crate::bam::{
    AnnotateSpec, BamReader, ErrorModel, FlagFilters, LiveRead, PileupSite, PresetConfig,
    mpileup_engine_from_records_opts,
};
use crate::call::GvcfBlocker;
use crate::cli::args::MpileupArgs;
use crate::regions::RegionSet;
use crate::bam::errmod::{Drand48, pack_base};
use crate::call::haplotype::haplotype_pls;
use crate::call::mcall::{CallResult, CallSite, Caller, CallerOpts, PL_MISSING};
use anyhow::{Context, Result, bail};
use fxhash::FxHashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

/// Reference contigs are loaded on demand through the `.fai` index.
pub(crate) type Fasta = crate::fasta::IndexedFasta;

pub(crate) fn load_fasta(p: &PathBuf) -> Result<Fasta> {
    Fasta::open(p)
}

/// Format a site/genotype QUAL float the way bcftools prints it: `.` for 0,
/// an integer when whole, otherwise up to two decimals.
fn fmt_qual(q: f64) -> String {
    if q <= 0.0 {
        ".".to_string()
    } else if (q - q.round()).abs() < 1e-3 {
        format!("{:.0}", q)
    } else {
        format!("{:.2}", q)
    }
}

/// Per-BAM preparation shared by both entry points: flag filters, the
/// reference pre-scan (skip-list + BAQ mask), BAQ, and the BAQ-less capping
/// fallback. Returns the union skip-list across samples.
fn prepare_bams(
    bams: &mut [BamReader],
    args: &MpileupArgs,
    flag_filters: &FlagFilters,
    baq_off: bool,
    fasta: Option<&Fasta>,
    timing: bool,
) -> Option<InterestingMap> {
    let mut combined: Option<InterestingMap> = None;
    // Same threshold the emitter uses, so the sequential and parallel paths
    // visit exactly the sites that can produce output.
    let min_alt = if args.variants_only { args.min_alt_reads.max(1) } else { 1 };
    for r in bams.iter_mut() {
        if flag_filters.is_active() {
            r.records_buf.retain(|lr| flag_filters.passes(lr.flags));
        }
        // The bcftools layout realigns lazily per site (`mplp_realn`).
        if !baq_off && !args.bcftools_layout {
            match fasta {
                Some(fa) => {
                    let ref_names = r.ref_names.clone();
                    let t_pre = std::time::Instant::now();
                    let pre_arc = std::sync::Arc::new(std::mem::take(&mut r.records_buf));
                    let pre = crate::bam::pos_filter::pre_scan(std::slice::from_ref(&pre_arc), fa, &ref_names, min_alt);
                    r.records_buf = std::sync::Arc::try_unwrap(pre_arc).unwrap_or_else(|arc| (*arc).clone());
                    if timing {
                        eprintln!(
                            "[KIRA_BT] pre-scan: {:.1}s, skip-list={}, baq-skip={}/{}",
                            t_pre.elapsed().as_secs_f64(),
                            pre.pos_filter.total(),
                            pre.skipped_baq,
                            pre.needs_baq.len(),
                        );
                    }
                    if args.recal {
                        // Empirical qualities learned at non-candidate sites, before BAQ.
                        let t_rc = std::time::Instant::now();
                        let arc = std::sync::Arc::new(std::mem::take(&mut r.records_buf));
                        let table = crate::bam::pos_filter::RecalTable::build(std::slice::from_ref(&arc), fa, &ref_names, &pre.pos_filter);
                        r.records_buf = std::sync::Arc::try_unwrap(arc).unwrap_or_else(|a| (*a).clone());
                        table.apply(&mut r.records_buf, fa, &ref_names);
                        if timing {
                            eprintln!("[KIRA_BT] recal: {} cells calibrated in {:.1}s", table.n_calibrated(), t_rc.elapsed().as_secs_f64());
                        }
                    }
                    let t1 = std::time::Instant::now();
                    crate::bam::reader::apply_hmm_baq_to_reads_masked(&mut r.records_buf, &ref_names, fa, Some(&pre.needs_baq), &args.nm_weight, args.redo_baq);
                    if timing { eprintln!("[KIRA_BT] BAQ: {:.1}s", t1.elapsed().as_secs_f64()); }
                    combined = Some(match combined.take() {
                        None => pre.pos_filter,
                        Some(mut m) => {
                            m.merge(pre.pos_filter);
                            m
                        }
                    });
                }
                // No reference: no BAQ, as in bcftools.
                None => {}
            }
        }
    }
    combined
}

pub fn cmd_mpileup(args: MpileupArgs) -> Result<()> {
    if args.inputs.is_empty() { bail!("mpileup: need at least one BAM input"); }

    let preset = args.config.as_deref().map(PresetConfig::parse).transpose()?;
    let min_mq = preset.as_ref().and_then(|p| p.min_mq).unwrap_or(args.min_mq) as u8;
    let min_bq = resolve_min_bq(&args, preset.as_ref());
    let _max_depth = preset.as_ref().and_then(|p| p.max_depth).unwrap_or(args.max_depth);
    let _indel_size = preset.as_ref().and_then(|p| p.indel_size).unwrap_or(args.indel_size);
    let baq_off = preset.as_ref().and_then(|p| p.no_baq).unwrap_or(args.no_baq);
    let skip_indels = args.skip_indels;
    let annotate_spec = if args.annotate.is_empty() { None } else { Some(args.annotate.join(",")) };
    let annotate = parse_annotate(&args, annotate_spec.as_deref())?;
    if args.bcftools_layout {
        if args.variants_only {
            bail!("mpileup: --bcftools-layout writes genotype likelihoods only; -v/--variants-only needs the kira layout");
        }
        if args.gvcf.is_some() {
            bail!("mpileup: --gvcf is not available with --bcftools-layout");
        }
        if args.assemble {
            bail!("mpileup: --assemble is not available with --bcftools-layout");
        }
    }
    let flag_filters = flag_filters(&args)?;
    let gap_frac = preset.as_ref().and_then(|p| p.gap_frac).unwrap_or(args.gap_frac);
    let min_ireads = args.min_ireads;
    let sopts = SampleOpts::parse(args.samples.as_deref(), args.samples_file.as_deref(), args.read_groups.as_deref(), args.ignore_rg)?;
    let (smap, bam_ids) = build_sample_map(&args.inputs, args.fasta_ref.as_deref(), &sopts)?;
    if args.recal && args.stream {
        bail!("mpileup: --recal needs two passes over the reads and cannot be combined with --stream");
    }
    if args.recal && args.fasta_ref.is_none() {
        bail!("mpileup: --recal needs the reference (-f)");
    }
    if !matches!(args.indel_realign.as_str(), "off" | "ins" | "all") {
        bail!("mpileup: --indel-realign expects off, ins or all (got {:?})", args.indel_realign);
    }
    let regions = collect_regions(&args)?;

    let fasta = args.fasta_ref.as_ref().map(load_fasta).transpose()?;
    let out_path = args.output.clone().unwrap_or_else(|| PathBuf::from("out.mpileup.vcf"));
    let mut out = BufWriter::with_capacity(1 << 20, File::create(&out_path).context("create output")?);

    let timing = std::env::var("KIRA_BT_TIMING").is_ok();
    if args.stream {
        run_streaming(&args, &flag_filters, baq_off, fasta, &annotate, min_mq, min_bq, skip_indels, min_ireads, gap_frac, &smap, &bam_ids, regions.as_ref(), timing, &mut out)?;
        out.flush()?;
        return Ok(());
    }
    if args.bcftools_layout && (args.regions.is_some() || args.regions_file.is_some()) {
        // bcftools walks `-r` regions one at a time, each with its own
        // pileup, so overlapping regions repeat their records.
        let targets = collect_targets(&args)?;
        let mut header_written = false;
        for (c, b, e) in ordered_regions(&args)? {
            let mut bams: Vec<BamReader> = Vec::with_capacity(args.inputs.len());
            for (p, id) in args.inputs.iter().zip(&bam_ids) {
                let Some(id) = id else { continue };
                bams.push(open_input(p, Some(&[(c.clone(), b, e)]), args.fasta_ref.as_deref(), &smap, *id)?);
            }
            let combined_filter = prepare_bams(&mut bams, &args, &flag_filters, baq_off, fasta.as_ref(), timing);
            let mut one = RegionSet::default();
            one.add(&c, b, e);
            one.finalize();
            run_prepared(bams, &args, combined_filter, fasta.as_ref(), &annotate, min_mq, min_bq, skip_indels, min_ireads, gap_frac, &smap.names, Some(&one), targets.as_ref(), !header_written, timing, &mut out)?;
            header_written = true;
        }
        out.flush()?;
        return Ok(());
    }
    let region_list: Option<Vec<(String, u32, u32)>> =
        regions.as_ref().map(|rs| rs.iter().map(|(c, b, e)| (c.to_string(), b, e)).collect());
    let mut bams: Vec<BamReader> = Vec::with_capacity(args.inputs.len());
    for (p, id) in args.inputs.iter().zip(&bam_ids) {
        let Some(id) = id else { continue };
        let t0 = std::time::Instant::now();
        let r = open_input(p, region_list.as_deref(), args.fasta_ref.as_deref(), &smap, *id)?;
        if timing { eprintln!("[KIRA_BT] BAM load: {:.1}s, {} records", t0.elapsed().as_secs_f64(), r.records_buf.len()); }
        bams.push(r);
    }
    let combined_filter = prepare_bams(&mut bams, &args, &flag_filters, baq_off, fasta.as_ref(), timing);

    run_prepared(bams, &args, combined_filter, fasta.as_ref(), &annotate, min_mq, min_bq, skip_indels, min_ireads, gap_frac, &smap.names, regions.as_ref(), None, true, timing, &mut out)?;
    out.flush()?;
    Ok(())
}

/// `-Q`: the preset, the flag, then the layout's default (bcftools 1.19
/// skips only quality 0; kira keeps its historical 13).
fn resolve_min_bq(args: &MpileupArgs, preset: Option<&PresetConfig>) -> u8 {
    preset
        .and_then(|p| p.min_bq)
        .or(args.min_bq)
        .unwrap_or(if args.bcftools_layout { 1 } else { 13 })
        .min(255) as u8
}

/// Read-level flag filters. The bcftools layout uses bcftools' rules alone:
/// `--ff` is the only exclusion (its default drops unmapped, secondary,
/// QC-fail and duplicate reads, not supplementary ones) and `--rf` skips
/// reads with every listed bit unset.
fn flag_filters(args: &MpileupArgs) -> Result<FlagFilters> {
    if args.bcftools_layout {
        let mut f = FlagFilters::from_full(args.ef, args.df, args.if_, args.nf, None, args.ff.as_deref())?
            .with_skip_all(args.skip_all_set.as_deref(), args.skip_all_unset.as_deref())?;
        if let Some(rf) = &args.rf {
            f.skip_all_unset |= crate::bam::mpileup_opts::parse_sam_flags(rf)?;
        }
        return Ok(f);
    }
    FlagFilters::from_full(args.ef, args.df, args.if_, args.nf, args.rf.as_deref(), args.ff.as_deref())?
        .with_skip_all(args.skip_all_set.as_deref(), args.skip_all_unset.as_deref())
}

/// Samples from the input headers under `-s/-S/-G` (bcftools' `bam_smpl`);
/// per input, its id in the map or `None` when bcftools would drop the file.
fn build_sample_map(inputs: &[PathBuf], reference: Option<&Path>, sopts: &SampleOpts) -> Result<(SampleMap, Vec<Option<usize>>)> {
    let mut smap = SampleMap::default();
    let mut ids = Vec::with_capacity(inputs.len());
    for p in inputs {
        let h = BamReader::read_header_only(p, reference).with_context(|| format!("open {}", p.display()))?;
        let id = smap.add_bam(&h, &p.to_string_lossy(), sopts);
        if id.is_none() {
            eprintln!("[mpileup] warning: no selected read group in {}, skipping the file", p.display());
        }
        ids.push(id);
    }
    if smap.names.is_empty() {
        bail!("mpileup: failed to find a file header with usable read groups");
    }
    Ok((smap, ids))
}

/// Load one input with every read tagged by its sample.
fn open_input(p: &Path, regions: Option<&[(String, u32, u32)]>, reference: Option<&Path>, smap: &SampleMap, bam_id: usize) -> Result<BamReader> {
    let resolve = |rg: Option<&[u8]>| smap.sample_id(bam_id, rg);
    match regions {
        Some(r) => BamReader::open_with_regions_resolved(p, r, reference, Some(&resolve)),
        None => BamReader::open_with_reference_resolved(p, reference, Some(&resolve)),
    }
}

/// One `BamReader` per sample, reads merged across inputs in position order.
fn regroup_by_sample(bams: Vec<BamReader>, n: usize) -> Vec<BamReader> {
    let header = bams[0].header.clone();
    let ref_names = bams[0].ref_names.clone();
    let mut per: Vec<Vec<LiveRead>> = (0..n).map(|_| Vec::new()).collect();
    let single_source = bams.len() == 1;
    for b in bams {
        for lr in b.records_buf {
            if lr.sample_idx < n {
                per[lr.sample_idx].push(lr);
            }
        }
    }
    per.into_iter()
        .map(|mut v| {
            if !single_source {
                v.sort_by_key(|lr| (lr.ref_id, lr.ref_start));
            }
            BamReader::from_parts(v, header.clone(), Vec::new(), ref_names.clone())
        })
        .collect()
}

/// `-a` under the requested layout: kira's tag set, or bcftools' defaults
/// and tag names with `--bcftools-layout`.
fn parse_annotate(args: &MpileupArgs, spec: Option<&str>) -> Result<AnnotateSpec> {
    if args.bcftools_layout { AnnotateSpec::parse_bcftools(spec) } else { AnnotateSpec::parse(spec) }
}

/// Read and base filters of the pileup pass. The bcftools layout applies
/// bcftools' defaults (`--delta-BQ 30`, `--max-BQ 60`, orphan pairs dropped
/// unless `-A`); the kira layout applies the two quality knobs only when
/// given explicitly.
fn engine_opts(args: &MpileupArgs, min_mq: u8, min_bq: u8, skip_indels: bool) -> EngineOpts {
    let mut o = EngineOpts::new(min_mq, min_bq, skip_indels, !args.ignore_overlaps);
    let clamp = |v: u32| v.min(255) as u8;
    if args.bcftools_layout {
        o.max_bq = Some(clamp(args.max_bq.unwrap_or(60)));
        o.delta_bq = Some(clamp(args.delta_bq.unwrap_or(30)));
        o.skip_orphans = !args.count_orphans;
        o.keep_n = true;
        o.emit_covered = true;
        o.max_depth = Some(args.max_depth);
        o.builtin_flag_skip = false;
    } else {
        o.max_bq = args.max_bq.map(clamp);
        o.delta_bq = args.delta_bq.map(clamp);
    }
    o
}

/// `-r/-R` regions in bcftools' iteration order: contigs by first
/// appearance, regions by start within a contig, overlaps kept.
fn ordered_regions(args: &MpileupArgs) -> Result<Vec<(String, u32, u32)>> {
    let mut regs: Vec<(String, u32, u32)> = Vec::new();
    if let Some(s) = &args.regions {
        for tok in s.split(',') {
            let t = tok.trim();
            if t.is_empty() {
                continue;
            }
            regs.push(crate::regions::parse_region_spec(t)?);
        }
    }
    if let Some(f) = &args.regions_file {
        regs.extend(crate::cli::commands::tabix::read_region_file(f)?);
    }
    let mut order: Vec<String> = Vec::new();
    for (c, _, _) in &regs {
        if !order.contains(c) {
            order.push(c.clone());
        }
    }
    regs.sort_by_key(|(c, _, _)| order.iter().position(|x| x == c).unwrap_or(usize::MAX));
    Ok(regs)
}

/// bcftools' `mplp_realn` (1.19, `MPLP_REALN_PARTIAL`): at a site where some
/// pileup read carries an indel, BAQ is applied once to the reads that have
/// not been realigned yet, skipping reads whose indel sits far from both
/// ends when the pileup is deep enough. Positions before the first
/// qualifying site keep the reads' original qualities.
fn bcftools_realn_hook<'a>(fasta: &'a Fasta, ref_names: &'a [String], kept: Vec<bool>, max_read_len: usize, redo: bool) -> Box<SiteHook<'a>> {
    Box::new(move |live: &mut [LiveRead], _ref_id: usize, pos: u32| {
        let (mut nt, mut has_indel, mut has_clip) = (0i64, 0i64, 0i64);
        let (mut min_indel, mut max_indel) = (i64::MAX, i64::MIN);
        for lr in live.iter_mut() {
            if !kept.get(lr.sample_idx).copied().unwrap_or(false) {
                continue;
            }
            nt += 1;
            lr.seek(pos);
            let indel = match lr.indel_after(pos) {
                Some(IndelEvent::Ins(s)) => s.len() as i64,
                Some(IndelEvent::Del(l)) => -(l as i64),
                None => 0,
            };
            let cigar_indel = lr.cigar_pairs.iter().any(|(k, _)| matches!(k, Kind::Insertion | Kind::Deletion | Kind::Skip));
            if cigar_indel || indel != 0 {
                has_indel += 1;
            }
            if lr.has_softclip {
                has_clip += 1;
            }
            max_indel = max_indel.max(indel);
            min_indel = min_indel.min(indel);
        }
        if has_indel == 0
            || ((has_clip as f64) < 0.2 * nt as f64 && max_indel == min_indel && ((has_indel as f64) < 0.1 * nt as f64 || has_indel == 1))
        {
            return;
        }
        let realn_dist = 40 + 10 * i64::from(nt < 40) + 10 * i64::from(nt < 20);
        for lr in live.iter_mut() {
            if !kept.get(lr.sample_idx).copied().unwrap_or(false) || lr.realigned {
                continue;
            }
            lr.realigned = true;
            let l_qseq = lr.seq().len();
            if l_qseq > max_read_len {
                continue;
            }
            let cig = &lr.cigar_pairs;
            if nt > 15 && cig.len() > 1 {
                let long = l_qseq > 500;
                let is_match = |k: Kind| matches!(k, Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch);
                let (mut lm, mut rm, mut nm) = (0i64, 0i64, 0usize);
                for &(k, l) in cig.iter() {
                    if long && matches!(k, Kind::HardClip | Kind::SoftClip) {
                        continue;
                    }
                    if is_match(k) {
                        lm += l as i64;
                        nm += 1;
                    } else {
                        break;
                    }
                }
                if nm != cig.len() {
                    for &(k, l) in cig.iter().rev() {
                        if long && matches!(k, Kind::HardClip | Kind::SoftClip) {
                            continue;
                        }
                        if is_match(k) {
                            rm += l as i64;
                        } else {
                            break;
                        }
                    }
                    if lm >= realn_dist * 4 && rm >= realn_dist * 4 {
                        continue;
                    }
                    if lm >= realn_dist && rm >= realn_dist && (has_clip as f64) < (0.15 + 0.05 * f64::from(nt > 20)) * nt as f64 {
                        continue;
                    }
                }
            }
            if l_qseq > 500 {
                let rl: i64 = cig.iter().filter(|(k, _)| matches!(k, Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion | Kind::Skip)).map(|(_, l)| *l as i64).sum();
                if (rl - l_qseq as i64).abs() * l_qseq as i64 >= 500_000 {
                    continue;
                }
            }
            crate::bam::reader::baq_read(lr, ref_names, fasta, redo);
        }
    })
}

/// `-t/-T` targets, applied as a position filter (never as an iterator).
fn collect_targets(args: &MpileupArgs) -> Result<Option<RegionSet>> {
    let mut set = RegionSet::default();
    let mut any = false;
    if let Some(s) = &args.targets {
        for (c, b, e) in RegionSet::from_cli(s)?.iter() {
            set.add(c, b, e);
        }
        any = true;
    }
    if let Some(f) = &args.targets_file {
        for (c, b, e) in RegionSet::from_file(f)?.iter() {
            set.add(c, b, e);
        }
        any = true;
    }
    if !any {
        return Ok(None);
    }
    set.finalize();
    Ok(Some(set))
}

/// `-r/-R/-t/-T`: comma lists and region files (BED or `chr:beg-end`),
/// merged per contig so overlapping regions fetch every read once.
fn collect_regions(args: &MpileupArgs) -> Result<Option<RegionSet>> {
    let mut set = RegionSet::default();
    let mut any = false;
    for s in [&args.regions, &args.targets].into_iter().flatten() {
        for (c, b, e) in RegionSet::from_cli(s)?.iter() {
            set.add(c, b, e);
        }
        any = true;
    }
    for f in [&args.regions_file, &args.targets_file].into_iter().flatten() {
        for (c, b, e) in RegionSet::from_file(f)?.iter() {
            set.add(c, b, e);
        }
        any = true;
    }
    if !any {
        return Ok(None);
    }
    set.finalize();
    Ok(Some(set))
}

/// `--stream`: decode each BAM on its own thread and pile up straight from
/// the channels, so no BAM is held in memory. Flag filters and BAQ run per
/// read on the decoder threads; without the reference pre-scan every
/// covered site is visited.
#[allow(clippy::too_many_arguments)]
fn run_streaming<W: Write>(
    args: &MpileupArgs,
    flag_filters: &FlagFilters,
    baq_off: bool,
    fasta: Option<Fasta>,
    annotate: &AnnotateSpec,
    min_mq: u8,
    min_bq: u8,
    skip_indels: bool,
    min_ireads: u32,
    gap_frac: f64,
    smap: &SampleMap,
    bam_ids: &[Option<usize>],
    regions: Option<&RegionSet>,
    timing: bool,
    out: &mut W,
) -> Result<()> {
    use crate::bam::{ReadHook, StreamingBam};
    let fasta = fasta.map(std::sync::Arc::new);
    let workers = (crate::threads::bam_workers() / args.inputs.len().max(1)).max(1);
    let bcftools_layout = args.bcftools_layout;
    let redo_baq = args.redo_baq;
    // One decoder per input feeds one engine source, so each input must map
    // to exactly its own sample, in input order.
    let smap_arc = std::sync::Arc::new(smap.clone());
    let mut stream_inputs: Vec<(&PathBuf, usize)> = Vec::new();
    for (p, id) in args.inputs.iter().zip(bam_ids) {
        let Some(id) = id else { continue };
        let ss = smap.file_samples(*id);
        if ss.len() != 1 || ss[0] != stream_inputs.len() {
            bail!("mpileup --stream: each input must carry exactly one sample of its own; {} does not (use the default mode)", p.display());
        }
        stream_inputs.push((p, *id));
    }
    let mut streams = Vec::with_capacity(stream_inputs.len());
    for (i, (p, bam_id)) in stream_inputs.iter().enumerate() {
        let bam_id = *bam_id;
        let flags = flag_filters.clone();
        let fa = fasta.clone();
        let nm_weight = args.nm_weight.clone();
        let mut nm_ready = false;
        let hook: ReadHook = Box::new(move |lr: &mut LiveRead, ref_names: &[String]| -> bool {
            if flags.is_active() && !flags.passes(lr.flags) {
                return false;
            }
            if baq_off || bcftools_layout {
                return true;
            }
            // Without a reference there is no BAQ, as in bcftools.
            if let Some(fa) = fa.as_deref() {
                if !nm_ready {
                    crate::bam::baq::init_nm_profile(&nm_weight, lr.seq().len());
                    nm_ready = true;
                }
                crate::bam::reader::baq_read(lr, ref_names, fa, redo_baq);
            }
            true
        });
        let sm = std::sync::Arc::clone(&smap_arc);
        let resolver: std::sync::Arc<dyn Fn(Option<&[u8]>) -> Option<usize> + Send + Sync> =
            std::sync::Arc::new(move |rg| sm.sample_id(bam_id, rg));
        let _ = i;
        let s = StreamingBam::open_with_resolver(p, resolver, workers, Some(hook)).with_context(|| format!("open {}", p.display()))?;
        streams.push(s);
    }

    let samples: Vec<String> = smap.names.clone();
    let sample_keep: Vec<bool> = vec![true; samples.len()];
    let any_filter = false;
    let ref_names = streams[0].ref_names.clone();
    let ref_lengths = streams[0].ref_lengths.clone();
    if args.bcftools_layout {
        write_bcftools_header(out, annotate, &ref_names, &ref_lengths, &samples, &sample_keep, any_filter, args.fasta_ref.as_deref())?;
    } else {
        write_vcf_header(out, annotate, &ref_names, &ref_lengths, &samples, &sample_keep, any_filter)?;
    }

    let em = ErrorModel::new();
    let mut gvcf_blocker: Option<GvcfBlocker> = args
        .gvcf
        .as_ref()
        .map(|s| GvcfBlocker::new(s.split(',').filter_map(|t| t.parse::<u32>().ok()).collect()));
    let caller = site_caller(args, &sample_keep, any_filter);
    let eo = engine_opts(args, min_mq, min_bq, skip_indels);
    let ctx = EmitCtx {
        args,
        em: &em,
        seed: args.seed,
        rng: std::sync::Mutex::new(Drand48::new(args.seed as u32)),
        caller: &caller,
        annotate,
        fasta: fasta.as_deref(),
        ref_names: &ref_names,
        sample_keep: &sample_keep,
        any_filter,
        skip_indels,
        min_ireads,
        gap_frac,
        hp_indel: args.hp_indel,
        realign: realign_level(args),
        assemble: args.assemble,
        regions,
        targets: None,
        engine: eo,
    };

    if timing { eprintln!("[KIRA_BT] engine: streaming, {} decoder thread(s)", streams.len()); }
    let mut err: Option<anyhow::Error> = None;
    let mut buf: Vec<u8> = Vec::with_capacity(256);
    let mut hook: Option<Box<SiteHook<'_>>> = match (args.bcftools_layout && !baq_off, fasta.as_deref()) {
        (true, Some(fa)) => Some(bcftools_realn_hook(fa, &ref_names, sample_keep.clone(), args.max_read_len, args.redo_baq)),
        _ => None,
    };
    crate::bam::mpileup_engine_streaming_opts(streams, &eo, hook.as_deref_mut(), &mut |site, reads| {
        if err.is_none() {
            if let Err(e) = write_site(site, reads, &ctx, &mut gvcf_blocker, &mut buf, out) {
                err = Some(e);
            }
        }
    })?;
    if let Some(e) = err {
        return Err(e);
    }
    if let Some(mut g) = gvcf_blocker {
        g.flush(out)?;
    }
    Ok(())
}

/// Render one site (or fold it into the open gVCF block) and write it.
fn write_site<W: Write>(
    site: &PileupSite,
    reads: &[LiveRead],
    ctx: &EmitCtx,
    gvcf: &mut Option<GvcfBlocker>,
    buf: &mut Vec<u8>,
    out: &mut W,
) -> Result<()> {
    buf.clear();
    if let Some(g) = gvcf.as_mut() {
        if !emit_site_or_gvcf(site, reads, ctx, buf, Some(g), out)? {
            return Ok(());
        }
    } else {
        emit_site(site, reads, ctx, buf);
    }
    out.write_all(buf).context("write mpileup output")
}

/// Run mpileup over already-built [`BamReader`]s (records in memory) and write
/// VCF to `out`. Same logic as [`cmd_mpileup`] minus the file-opening — for the
/// fused `solid` pipeline that hands sorted in-memory records straight here.
pub fn run_mpileup_from_bams<W: Write>(
    mut bams: Vec<BamReader>,
    args: &MpileupArgs,
    out: &mut W,
) -> Result<()> {
    if bams.is_empty() {
        bail!("mpileup: no records");
    }
    let preset = args.config.as_deref().map(PresetConfig::parse).transpose()?;
    let min_mq = preset.as_ref().and_then(|p| p.min_mq).unwrap_or(args.min_mq) as u8;
    let min_bq = resolve_min_bq(args, preset.as_ref());
    let baq_off = preset.as_ref().and_then(|p| p.no_baq).unwrap_or(args.no_baq);
    let skip_indels = args.skip_indels;
    let annotate_spec = if args.annotate.is_empty() { None } else { Some(args.annotate.join(",")) };
    let annotate = parse_annotate(args, annotate_spec.as_deref())?;
    let flag_filters = flag_filters(args)?;
    let gap_frac = preset.as_ref().and_then(|p| p.gap_frac).unwrap_or(args.gap_frac);
    let min_ireads = args.min_ireads;
    let fasta = args.fasta_ref.as_ref().map(load_fasta).transpose()?;
    let timing = std::env::var("KIRA_BT_TIMING").is_ok();

    // In-memory inputs: one sample per reader, named from its read groups.
    let mut samples: Vec<String> = Vec::with_capacity(bams.len());
    for (i, b) in bams.iter_mut().enumerate() {
        for lr in b.records_buf.iter_mut() {
            lr.sample_idx = i;
        }
        samples.push(b.samples.first().cloned().unwrap_or_else(|| format!("sample{i}")));
    }
    let combined_filter = prepare_bams(&mut bams, args, &flag_filters, baq_off, fasta.as_ref(), timing);
    run_prepared(bams, args, combined_filter, fasta.as_ref(), &annotate, min_mq, min_bq, skip_indels, min_ireads, gap_frac, &samples, None, None, true, timing, out)?;
    out.flush()?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_prepared<W: Write>(
    bams: Vec<BamReader>,
    args: &MpileupArgs,
    combined_filter: Option<InterestingMap>,
    fasta: Option<&Fasta>,
    annotate: &AnnotateSpec,
    min_mq: u8,
    min_bq: u8,
    skip_indels: bool,
    min_ireads: u32,
    gap_frac: f64,
    sample_names: &[String],
    regions: Option<&RegionSet>,
    targets: Option<&RegionSet>,
    write_header: bool,
    timing: bool,
    out: &mut W,
) -> Result<()> {
    if bams.is_empty() {
        bail!("mpileup: no input left after sample selection");
    }
    let samples: Vec<String> = sample_names.to_vec();
    let sample_keep: Vec<bool> = vec![true; samples.len()];
    let any_filter = false;
    let bams = regroup_by_sample(bams, samples.len());

    let ref_names: Vec<String> = bams[0].ref_names.clone();
    let ref_lengths: Vec<u64> = bams[0].ref_lengths.clone();

    if write_header {
        if args.bcftools_layout {
            write_bcftools_header(out, annotate, &ref_names, &ref_lengths, &samples, &sample_keep, any_filter, args.fasta_ref.as_deref())?;
        } else {
            write_vcf_header(out, annotate, &ref_names, &ref_lengths, &samples, &sample_keep, any_filter)?;
        }
    }

    let em = ErrorModel::new();
    let gvcf_opts: Option<Vec<u32>> = args.gvcf.as_ref().map(|s| {
        s.split(',').filter_map(|t| t.parse::<u32>().ok()).collect()
    });
    let mut gvcf_blocker: Option<GvcfBlocker> = gvcf_opts.map(GvcfBlocker::new);

    let n_chunks = args.threads.max(1);
    // Lazy per-site BAQ makes the bcftools layout order-dependent, so it
    // always runs sequentially.
    let parallel_eligible = n_chunks > 1 && gvcf_blocker.is_none() && !args.bcftools_layout;
    // The skip-list drops reference-only sites; that is only valid when such
    // sites are never emitted. It is then applied on every path so the output
    // does not depend on the thread count.
    let pos_filter = if args.variants_only && gvcf_blocker.is_none() { combined_filter } else { None };

    let caller = site_caller(args, &sample_keep, any_filter);
    let ctx = EmitCtx {
        args,
        em: &em,
        seed: args.seed,
        rng: std::sync::Mutex::new(Drand48::new(args.seed as u32)),
        caller: &caller,
        annotate,
        fasta,
        ref_names: &ref_names,
        sample_keep: &sample_keep,
        any_filter,
        skip_indels,
        min_ireads,
        gap_frac,
        hp_indel: args.hp_indel,
        realign: realign_level(args),
        assemble: args.assemble,
        regions,
        targets,
        engine: engine_opts(args, min_mq, min_bq, skip_indels),
    };

    let t_engine = std::time::Instant::now();
    if parallel_eligible {
        if timing { eprintln!("[KIRA_BT] engine: parallel, {} chunks", n_chunks); }
        run_parallel(bams, n_chunks, &ctx, min_mq, min_bq, skip_indels, pos_filter, out)?;
    } else {
        if timing { eprintln!("[KIRA_BT] engine: sequential"); }
        let baq_off = args.config.as_deref().map(PresetConfig::parse).transpose()?.and_then(|p| p.no_baq).unwrap_or(args.no_baq);
        let mut hook: Option<Box<SiteHook<'_>>> = match (args.bcftools_layout && !baq_off, fasta) {
            (true, Some(fa)) => Some(bcftools_realn_hook(fa, &ref_names, sample_keep.clone(), args.max_read_len, args.redo_baq)),
            _ => None,
        };
        run_sequential(bams, &ctx, min_mq, min_bq, skip_indels, pos_filter.as_ref(), hook.as_deref_mut(), &mut gvcf_blocker, out)?;
        if let Some(mut g) = gvcf_blocker { g.flush(out)?; }
    }
    if timing { eprintln!("[KIRA_BT] engine done in {:.1}s", t_engine.elapsed().as_secs_f64()); }
    Ok(())
}

fn write_vcf_header<W: Write>(
    out: &mut W,
    annotate: &AnnotateSpec,
    ref_names: &[String],
    ref_lengths: &[u64],
    samples: &[String],
    sample_keep: &[bool],
    any_filter: bool,
) -> Result<()> {
    writeln!(out, "##fileformat=VCFv4.2")?;
    writeln!(out, "##FILTER=<ID=PASS,Description=\"All filters passed\">")?;
    writeln!(out, "##source=kira_bt_mpileup")?;
    writeln!(out, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">")?;
    writeln!(out, "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"Mean MAPQ at site\">")?;
    writeln!(out, "##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">")?;
    if annotate.info_ad { writeln!(out, "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Total allelic depths (high-quality bases)\">")?; }
    if annotate.info_adf { writeln!(out, "##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Total allelic depths on the forward strand (high-quality bases)\">")?; }
    if annotate.info_adr { writeln!(out, "##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Total allelic depths on the reverse strand (high-quality bases)\">")?; }
    if annotate.info_scr { writeln!(out, "##INFO=<ID=SCR,Number=1,Type=Integer,Description=\"Number of soft-clipped reads (at high-quality bases)\">")?; }
    if annotate.info_dp4 { writeln!(out, "##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Number of high-quality ref-forward, ref-reverse, alt-forward and alt-reverse bases\">")?; }
    writeln!(out, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")?;
    if annotate.fmt_pl { writeln!(out, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">")?; }
    if annotate.fmt_dp { writeln!(out, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of high-quality bases\">")?; }
    if annotate.fmt_ad { writeln!(out, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths (high-quality bases)\">")?; }
    if annotate.fmt_adf { writeln!(out, "##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths on the forward strand (high-quality bases)\">")?; }
    if annotate.fmt_adr { writeln!(out, "##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths on the reverse strand (high-quality bases)\">")?; }
    if annotate.fmt_dv { writeln!(out, "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Number of high-quality non-reference bases\">")?; }
    if annotate.fmt_dp4 { writeln!(out, "##FORMAT=<ID=DP4,Number=4,Type=Integer,Description=\"Number of high-quality ref-forward, ref-reverse, alt-forward and alt-reverse bases\">")?; }
    if annotate.fmt_qs { writeln!(out, "##FORMAT=<ID=QS,Number=R,Type=Integer,Description=\"Quality sum per allele\">")?; }
    if annotate.fmt_sp { writeln!(out, "##FORMAT=<ID=SP,Number=1,Type=Integer,Description=\"Phred-scaled strand bias P-value\">")?; }
    if annotate.fmt_scr { writeln!(out, "##FORMAT=<ID=SCR,Number=1,Type=Integer,Description=\"Number of soft-clipped reads (at high-quality bases)\">")?; }
    if annotate.fmt_gq { writeln!(out, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">")?; }
    for (i, rn) in ref_names.iter().enumerate() {
        match ref_lengths.get(i) {
            Some(len) if *len > 0 => writeln!(out, "##contig=<ID={},length={}>", rn, len)?,
            _ => writeln!(out, "##contig=<ID={}>", rn)?,
        }
    }
    write!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    for (i, s) in samples.iter().enumerate() {
        if !any_filter || sample_keep[i] { write!(out, "\t{}", s)?; }
    }
    writeln!(out)?;
    Ok(())
}

/// bcftools `mpileup` header: reference, contigs, then the tag lines in
/// bcftools' order.
#[allow(clippy::too_many_arguments)]
fn write_bcftools_header<W: Write>(
    out: &mut W,
    annotate: &AnnotateSpec,
    ref_names: &[String],
    ref_lengths: &[u64],
    samples: &[String],
    sample_keep: &[bool],
    any_filter: bool,
    fasta: Option<&Path>,
) -> Result<()> {
    writeln!(out, "##fileformat=VCFv4.2")?;
    writeln!(out, "##FILTER=<ID=PASS,Description=\"All filters passed\">")?;
    if let Some(f) = fasta {
        writeln!(out, "##reference=file://{}", f.display())?;
    }
    for (i, rn) in ref_names.iter().enumerate() {
        writeln!(out, "##contig=<ID={},length={}>", rn, ref_lengths.get(i).copied().unwrap_or(0))?;
    }
    for l in bcf_layout::header_lines(annotate) {
        writeln!(out, "{l}")?;
    }
    write!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    for (i, s) in samples.iter().enumerate() {
        if !any_filter || sample_keep[i] { write!(out, "\t{}", s)?; }
    }
    writeln!(out)?;
    Ok(())
}

struct EmitCtx<'a> {
    args: &'a MpileupArgs,
    /// Dependent-error genotype likelihoods (htslib errmod).
    em: &'a ErrorModel,
    /// `--seed`: the bcftools layout draws from one `hts_drand48` stream,
    /// the kira layout seeds one per site (see [`Drand48::for_site`]).
    seed: u64,
    rng: std::sync::Mutex<Drand48>,
    /// Multiallelic caller (bcftools `call -m`) applied to every site.
    caller: &'a Caller,
    annotate: &'a AnnotateSpec,
    fasta: Option<&'a Fasta>,
    ref_names: &'a [String],
    sample_keep: &'a [bool],
    any_filter: bool,
    skip_indels: bool,
    min_ireads: u32,
    gap_frac: f64,
    /// `--hp-indel`: homopolymer/STR-aware stringency and quality for short indels.
    hp_indel: bool,
    /// `--indel-realign`: 0 off, 1 insertions, 2 insertions and deletions.
    realign: u8,
    /// `--assemble`: local indel recovery in active regions.
    assemble: bool,
    /// Requested regions: only sites inside them are emitted.
    regions: Option<&'a RegionSet>,
    /// `-t/-T` targets (bcftools layout), a second position filter.
    targets: Option<&'a RegionSet>,
    /// The pileup pass' filters, for the bcftools layout's quality rules.
    engine: EngineOpts,
}

fn realign_level(args: &MpileupArgs) -> u8 {
    match args.indel_realign.as_str() {
        "ins" => 1,
        "all" => 2,
        _ => 0,
    }
}

/// The caller shared by every site: θ from `--prior`, Watterson-scaled over
/// the samples in the output.
fn site_caller(args: &MpileupArgs, sample_keep: &[bool], any_filter: bool) -> Caller {
    let n = if any_filter { sample_keep.iter().filter(|k| **k).count() } else { sample_keep.len() };
    Caller::new(CallerOpts { theta: args.prior, ploidy: 2, ..CallerOpts::default() }, n.max(1))
}

impl EmitCtx<'_> {
    fn in_regions(&self, chr: &str, pos1: u64) -> bool {
        let p = pos1.min(u32::MAX as u64) as u32;
        self.regions.is_none_or(|rs| rs.contains(chr, p)) && self.targets.is_none_or(|ts| ts.contains(chr, p))
    }

    /// Samples in output order (`-s/-S` applied).
    fn kept_samples(&self, n: usize) -> Vec<usize> {
        (0..n).filter(|&si| !self.any_filter || self.sample_keep.get(si).copied().unwrap_or(false)).collect()
    }

    fn layout(&self) -> LayoutCtx<'_> {
        LayoutCtx { tags: self.annotate, em: self.em, rng: &self.rng, delta_bq: self.engine.delta_bq, max_bq: self.engine.max_bq }
    }
}

#[allow(clippy::too_many_arguments)]
fn run_sequential<W: Write>(
    mut bams: Vec<BamReader>,
    ctx: &EmitCtx,
    min_mq: u8,
    min_bq: u8,
    skip_indels: bool,
    pos_filter: Option<&InterestingMap>,
    hook: Option<&mut SiteHook<'_>>,
    gvcf_blocker: &mut Option<GvcfBlocker>,
    out: &mut W,
) -> Result<()> {
    let records_per_sample: Vec<Vec<LiveRead>> = bams
        .iter_mut()
        .enumerate()
        .map(|(i, b)| {
            let mut v = std::mem::take(&mut b.records_buf);
            for lr in v.iter_mut() { lr.sample_idx = i; }
            v
        })
        .collect();
    let mut err: Option<anyhow::Error> = None;
    let mut buf: Vec<u8> = Vec::with_capacity(256);
    let eo = engine_opts(ctx.args, min_mq, min_bq, skip_indels);
    mpileup_engine_from_records_opts(records_per_sample, &eo, pos_filter, hook, &mut |site, overlapping| {
        if err.is_none() {
            if let Err(e) = write_site(site, overlapping, ctx, gvcf_blocker, &mut buf, out) {
                err = Some(e);
            }
        }
    })?;
    match err {
        Some(e) => Err(e),
        None => Ok(()),
    }
}

#[allow(clippy::too_many_arguments)]
fn run_parallel<W: Write>(
    mut bams: Vec<BamReader>,
    n_chunks: usize,
    ctx: &EmitCtx,
    min_mq: u8,
    min_bq: u8,
    skip_indels: bool,
    pos_filter: Option<InterestingMap>,
    out: &mut W,
) -> Result<()> {
    use rayon::prelude::*;

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

    // The skip-list normally comes from the pre-scan; build it here only when
    // BAQ was off and reference-only sites are dropped anyway.
    let pos_filter = pos_filter.or_else(|| {
        if !ctx.args.variants_only { return None; }
        ctx.fasta.map(|fa| {
            let t = std::time::Instant::now();
            let m = crate::bam::pos_filter::pre_scan(&records_per_sample, fa, ctx.ref_names, ctx.args.min_alt_reads.max(1)).pos_filter;
            if std::env::var("KIRA_BT_TIMING").is_ok() {
                eprintln!("[KIRA_BT] late skip-list: {} in {:.1}s", m.total(), t.elapsed().as_secs_f64());
            }
            m
        })
    });

    // Reads that start before a chunk but reach into it must be included, so
    // the left boundary backs up by the longest reference span seen in the data.
    let max_span: u32 = records_per_sample
        .iter()
        .flat_map(|v| v.iter())
        .map(|lr| lr.ref_end_cached.saturating_sub(lr.ref_start))
        .max()
        .unwrap_or(0)
        .max(1);

    use indicatif::{ProgressBar, ProgressStyle};
    let pb = ProgressBar::new(chunks.len() as u64);
    pb.set_style(
        ProgressStyle::with_template("[ENGINE] {bar:40.magenta/blue} {pos}/{len} chunks ({per_sec}, ETA {eta})")
            .unwrap_or_else(|_| ProgressStyle::default_bar()),
    );
    let no_pb = std::env::var("KIRA_BT_NO_PROGRESS").is_ok();
    pb.set_draw_target(if no_pb { indicatif::ProgressDrawTarget::hidden() } else { indicatif::ProgressDrawTarget::stderr_with_hz(2) });

    let buffers: Vec<Vec<u8>> = chunks
        .par_iter()
        .map(|chunk| {
            let mut buf: Vec<u8> = Vec::with_capacity(1 << 20);
            let lo_pos = chunk.start.saturating_sub(max_span);
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
            let _ = mpileup_engine_from_records_opts(
                per_sample,
                &engine_opts(ctx.args, min_mq, min_bq, skip_indels),
                pos_filter.as_ref(),
                None,
                &mut |site, overlapping| {
                    if site.ref_id != chunk.ref_id { return; }
                    if site.pos < chunk.start || site.pos >= chunk.end { return; }
                    emit_site(site, overlapping, ctx, &mut buf);
                },
            );
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

/// FORMAT keys in output order for the requested annotations.
fn format_keys(annotate: &AnnotateSpec) -> Vec<&'static str> {
    let mut keys = vec!["GT"];
    if annotate.fmt_dp { keys.push("DP"); }
    if annotate.fmt_ad { keys.push("AD"); }
    if annotate.fmt_adf { keys.push("ADF"); }
    if annotate.fmt_adr { keys.push("ADR"); }
    if annotate.fmt_dv { keys.push("DV"); }
    if annotate.fmt_dp4 { keys.push("DP4"); }
    if annotate.fmt_qs { keys.push("QS"); }
    if annotate.fmt_sp { keys.push("SP"); }
    if annotate.fmt_scr { keys.push("SCR"); }
    if annotate.fmt_pl { keys.push("PL"); }
    if annotate.fmt_gq { keys.push("GQ"); }
    keys
}

/// Sample column for a sample with no reads at the site.
fn empty_sample_col(annotate: &AnnotateSpec, n_alleles: usize) -> String {
    let zeros = vec!["0"; n_alleles].join(",");
    let mut parts: Vec<String> = vec!["./.".into()];
    if annotate.fmt_dp { parts.push("0".into()); }
    if annotate.fmt_ad { parts.push(zeros.clone()); }
    if annotate.fmt_adf { parts.push(zeros.clone()); }
    if annotate.fmt_adr { parts.push(zeros.clone()); }
    if annotate.fmt_dv { parts.push("0".into()); }
    if annotate.fmt_dp4 { parts.push("0,0,0,0".into()); }
    if annotate.fmt_qs { parts.push(zeros); }
    if annotate.fmt_sp { parts.push("0".into()); }
    if annotate.fmt_scr { parts.push("0".into()); }
    if annotate.fmt_pl { parts.push("0".into()); }
    if annotate.fmt_gq { parts.push("0".into()); }
    parts.join(":")
}

fn join_u32(v: &[u32]) -> String {
    v.iter().map(u32::to_string).collect::<Vec<_>>().join(",")
}

/// ref-fwd, ref-rev, alt-fwd, alt-rev from per-allele counts and forward counts.
fn dp4_of(counts: &[u32], fwd: &[u32]) -> String {
    let rf = fwd.first().copied().unwrap_or(0);
    let rr = counts.first().copied().unwrap_or(0).saturating_sub(rf);
    let af: u32 = fwd.iter().skip(1).sum();
    let ar: u32 = counts.iter().skip(1).zip(fwd.iter().skip(1)).map(|(c, f)| c.saturating_sub(*f)).sum();
    format!("{rf},{rr},{af},{ar}")
}

/// Per-sample FORMAT values shared by SNV and indel records.
#[allow(clippy::too_many_arguments)]
fn sample_parts(annotate: &AnnotateSpec, gt: String, depth: u32, counts: &[u32], fwd: &[u32], quals: &[u32], n_softclip: u32, pl: &str, gq: u32) -> Vec<String> {
    let rev: Vec<u32> = counts.iter().zip(fwd.iter()).map(|(c, f)| c.saturating_sub(*f)).collect();
    let alt_fwd: u32 = fwd.iter().skip(1).sum();
    let alt_rev: u32 = rev.iter().skip(1).sum();
    let mut parts: Vec<String> = vec![gt];
    if annotate.fmt_dp { parts.push(depth.to_string()); }
    if annotate.fmt_ad { parts.push(join_u32(counts)); }
    if annotate.fmt_adf { parts.push(join_u32(fwd)); }
    if annotate.fmt_adr { parts.push(join_u32(&rev)); }
    if annotate.fmt_dv { parts.push(counts.iter().skip(1).sum::<u32>().to_string()); }
    if annotate.fmt_dp4 { parts.push(dp4_of(counts, fwd)); }
    if annotate.fmt_qs { parts.push(join_u32(quals)); }
    if annotate.fmt_sp { parts.push(strand_bias_phred(fwd.first().copied().unwrap_or(0), rev.first().copied().unwrap_or(0), alt_fwd, alt_rev).to_string()); }
    if annotate.fmt_scr { parts.push(n_softclip.to_string()); }
    if annotate.fmt_pl { parts.push(pl.to_string()); }
    if annotate.fmt_gq { parts.push(gq.to_string()); }
    parts
}

/// Site-level INFO: depth, mean MAPQ, optional allelic depths by strand and soft-clip count.
fn site_info(total: u32, mean_mq: f64, ad: &[u32], adf: &[u32], scr: u32, annotate: &AnnotateSpec, indel: bool) -> String {
    let mut info = format!("DP={};MQ={:.1}", total, mean_mq);
    if indel { info.push_str(";INDEL"); }
    if annotate.info_ad { info.push_str(&format!(";AD={}", join_u32(ad))); }
    if annotate.info_adf { info.push_str(&format!(";ADF={}", join_u32(adf))); }
    if annotate.info_adr {
        let adr: Vec<u32> = ad.iter().zip(adf.iter()).map(|(a, f)| a.saturating_sub(*f)).collect();
        info.push_str(&format!(";ADR={}", join_u32(&adr)));
    }
    if annotate.info_dp4 { info.push_str(&format!(";DP4={}", dp4_of(ad, adf))); }
    if annotate.info_scr { info.push_str(&format!(";SCR={}", scr)); }
    info
}

/// Per-site emit, writing the formatted VCF line into `out` (a Vec<u8> buffer
/// in the parallel runner; a BufWriter slice in the sequential one).
fn emit_site(
    site: &PileupSite,
    overlapping_reads: &[LiveRead],
    ctx: &EmitCtx,
    out: &mut Vec<u8>,
) {
    let args = ctx.args;
    let annotate = ctx.annotate;
    let chr = ctx.ref_names.get(site.ref_id).map(|s| s.as_str()).unwrap_or(".");
    let pos1 = site.pos as u64 + 1;
    if !ctx.in_regions(chr, pos1) { return; }
    if args.bcftools_layout {
        emit_site_bcftools(site, overlapping_reads, ctx, chr, pos1, out);
        return;
    }
    let ref_base = ctx.fasta.and_then(|fa| fa.base(chr, pos1 as u32)).unwrap_or(b'N');
    if ref_base == b'N' || ref_base == b'n' { return; }
    let ref_idx = match ref_base { b'A' | b'a' => 0, b'C' | b'c' => 1, b'G' | b'g' => 2, b'T' | b't' => 3, _ => return };
    let agg = site.aggregated();
    let total: u32 = agg.base_counts.iter().sum();
    if total == 0 { return; }

    let (alt_min_reads, alt_min_af) = if args.variants_only {
        (args.min_alt_reads, args.min_af)
    } else {
        (1u32, 0.0f64)
    };
    let total_f = total.max(1) as f64;
    // Candidate ALT bases as (ACGT index, count), most frequent first.
    let mut snv_alts: Vec<(usize, u32)> = Vec::new();
    for i in 0..4 {
        if i == ref_idx { continue; }
        let c = agg.base_counts[i];
        if c >= alt_min_reads && (c as f64) / total_f >= alt_min_af {
            snv_alts.push((i, c));
        }
    }
    snv_alts.sort_by(|a, b| b.1.cmp(&a.1));
    // Candidate alleles in REF-first order, as ACGT indices, plus the least
    // supported other base standing in for bcftools' unseen `<*>` allele:
    // it is never emitted, but weighs the reference hypothesis against
    // "anything else" so reference sites get a real QUAL.
    let mut cand: Vec<usize> = std::iter::once(ref_idx).chain(snv_alts.iter().map(|a| a.0)).collect();
    let unseen = (0..4usize).filter(|i| !cand.contains(i)).min_by_key(|&i| agg.base_counts[i]).map(|i| {
        cand.push(i);
        cand.len() - 1
    });
    let n_cand = cand.len();
    let n_gt = n_cand * (n_cand + 1) / 2;
    let codes: Vec<u8> = cand.iter().map(|&c| c as u8).collect();
    let mean_mq = (agg.mq_sum as f64) / (total as f64);

    // Per-sample PLs from the dependent-error model, then the multiallelic caller.
    let mut kept_samples: Vec<usize> = Vec::new();
    let mut pls_all: Vec<i32> = Vec::new();
    let mut bases: Vec<u16> = Vec::new();
    // INFO/QS as `bcf_call_combine` builds it: per-sample quality sums per
    // allele, normalised within the sample, summed over samples.
    let mut qs = vec![0.0f64; n_cand];
    for (si, s) in site.per_sample.iter().enumerate() {
        if ctx.any_filter && !ctx.sample_keep.get(si).copied().unwrap_or(false) { continue; }
        kept_samples.push(si);
        if s.depth == 0 {
            pls_all.extend(std::iter::repeat_n(PL_MISSING, n_gt));
            continue;
        }
        bases.clear();
        let mut qsum = [0u64; 16];
        for &(code, q) in &s.obs {
            qsum[(code & 0xf) as usize] += q as u64;
            bases.push(pack_base(q, code & 0x10 != 0, code & 0xf));
        }
        let qtotal: u64 = qsum.iter().sum();
        if qtotal > 0 {
            for (k, &c) in codes.iter().enumerate() {
                qs[k] += qsum[c as usize] as f64 / qtotal as f64;
            }
        }
        let pl = ctx.em.pls(&mut bases, &codes, &mut Drand48::for_site(ctx.seed, site.ref_id, site.pos, si));
        pls_all.extend(pl.iter().map(|&v| v as i32));
    }
    if let Some(u) = unseen {
        qs[u] = 0.0;
    }
    let mut call_site = CallSite { n_samples: kept_samples.len(), n_alleles: n_cand, pls: pls_all, is_indel: false, depths: None, qs: Some(qs), unseen, sample_af: None, prior_an_ac: None };
    let CallResult::Called { alleles_kept, qual: site_qual, gts, gqs, pls: pls_sub, .. } = ctx.caller.call_site(&mut call_site) else { return };
    let site_qual = site_qual.unwrap_or(0.0);
    let any_variant = alleles_kept.len() > 1;

    let have_indel = !ctx.skip_indels && {
        let m = if args.variants_only { args.min_alt_reads.max(ctx.min_ireads) } else { ctx.min_ireads };
        agg.ins_alleles.iter().any(|(_, c, _)| *c >= m) || agg.del_alleles.iter().any(|(_, c, _)| *c >= m)
    };
    if args.variants_only && (!any_variant || site_qual < args.min_qual as f64) && !have_indel {
        return;
    }

    let n_alleles = alleles_kept.len();
    let alt_str = if n_alleles == 1 {
        ".".to_string()
    } else {
        alleles_kept[1..].iter().map(|&a| (b"ACGT"[cand[a as usize]] as char).to_string()).collect::<Vec<_>>().join(",")
    };
    let mut per_sample_cols: Vec<String> = Vec::with_capacity(kept_samples.len());
    for (k, &si) in kept_samples.iter().enumerate() {
        let s = &site.per_sample[si];
        if s.depth == 0 {
            per_sample_cols.push(empty_sample_col(annotate, n_alleles));
            continue;
        }
        let counts: Vec<u32> = alleles_kept.iter().map(|&a| s.base_counts[cand[a as usize]]).collect();
        let fwd: Vec<u32> = alleles_kept.iter().map(|&a| s.fwd_counts[cand[a as usize]]).collect();
        let quals: Vec<u32> = alleles_kept.iter().map(|&a| s.base_quals[cand[a as usize]]).collect();
        let gt = gt_string(&alleles_kept, gts[k]);
        let pl = pl_string_or_ref(&pls_sub[k], &call_site.pls[k * n_gt..(k + 1) * n_gt]);
        let parts = sample_parts(annotate, gt, s.depth, &counts, &fwd, &quals, s.n_softclip, &pl, gqs[k]);
        per_sample_cols.push(parts.join(":"));
    }

    let ad: Vec<u32> = alleles_kept.iter().map(|&a| agg.base_counts[cand[a as usize]]).collect();
    let adf: Vec<u32> = alleles_kept.iter().map(|&a| agg.fwd_counts[cand[a as usize]]).collect();
    let info = site_info(total, mean_mq, &ad, &adf, agg.n_softclip, annotate, false);
    let fmt_keys = format_keys(annotate);

    let snv_emit = any_variant || !args.variants_only;
    if snv_emit && !(args.variants_only && site_qual < args.min_qual as f64) {
        let qual_str = fmt_qual(site_qual);
        let _ = write!(out, "{}\t{}\t.\t{}\t{}\t{}\t.\t{}\t{}",
            chr, pos1, ref_base as char, alt_str, qual_str, info, fmt_keys.join(":"));
        for col in &per_sample_cols { let _ = write!(out, "\t{}", col); }
        let _ = writeln!(out);
    }

    if !ctx.skip_indels {
        emit_indels(site, &agg, chr, pos1, ref_base, total, mean_mq, ctx, out, overlapping_reads);
    }

    // `--assemble`: at non-reference sites without a CIGAR indel, realign the
    // reads to find indels the aligner modelled as mismatches. Records may
    // land out of position order; sort downstream.
    if ctx.assemble && !ctx.skip_indels && !have_indel {
        let cfg = asm_cfg();
        let nonref = total.saturating_sub(agg.base_counts[ref_idx]);
        if nonref >= cfg.nonref {
            if let Some(fa) = ctx.fasta {
                let w_lo = site.pos.saturating_sub(cfg.up);
                if let Some(refw) = fa.slice_bytes(chr, w_lo + 1, cfg.len) {
                    if let Some(call) = crate::call::haplotype::assemble_indel(overlapping_reads, w_lo, refw, cfg.min_sup, cfg.max_mm) {
                        emit_assembled(call, site, overlapping_reads, chr, ctx, out, annotate);
                    }
                }
            }
        }
    }
}

/// bcftools layout: one SNP record per covered site (reference `N` when the
/// base is unknown), then kira's indel candidates as biallelic INDEL records.
fn emit_site_bcftools(site: &PileupSite, reads: &[LiveRead], ctx: &EmitCtx, chr: &str, pos1: u64, out: &mut Vec<u8>) {
    let kept = ctx.kept_samples(site.per_sample.len());
    if !reads.iter().any(|r| kept.contains(&r.sample_idx)) {
        return;
    }
    let ref_base = ctx.fasta.and_then(|fa| fa.base(chr, pos1 as u32)).unwrap_or(b'N');
    let lctx = ctx.layout();
    bcf_layout::render_snp(site, reads, &kept, ref_base, chr, pos1, &lctx, out);
    if ctx.skip_indels || matches!(ref_base, b'N' | b'n') {
        return;
    }
    let agg = site.aggregated();
    let total: u32 = agg.base_counts.iter().sum();
    if total == 0 {
        return;
    }
    let mean_mq = (agg.mq_sum as f64) / (total as f64);
    emit_indels(site, &agg, chr, pos1, ref_base, total, mean_mq, ctx, out, reads);
}

/// `a/b` over the kept-allele indices; `./.` for a missing genotype.
fn gt_string(alleles_kept: &[u32], gt: Option<(u32, u32)>) -> String {
    let Some(gt) = gt else { return "./.".to_string() };
    let gi = alleles_kept.iter().position(|&x| x == gt.0).unwrap_or(0);
    let gj = alleles_kept.iter().position(|&x| x == gt.1).unwrap_or(0);
    format!("{}/{}", gi.min(gj), gi.max(gj))
}

fn pl_string(pl: &[i32]) -> String {
    pl.iter().map(|&v| if v == PL_MISSING { ".".to_string() } else { v.to_string() }).collect::<Vec<_>>().join(",")
}

/// The caller drops PL at reference-only sites; mpileup output keeps the
/// homozygous-reference value so FORMAT/PL is always present.
fn pl_string_or_ref(pl: &[i32], raw: &[i32]) -> String {
    if pl.is_empty() { pl_string(&raw[..1.min(raw.len())]) } else { pl_string(pl) }
}

/// Like [`pl_string_or_ref`], but a record that prints its ALT anyway (an
/// indel candidate called 0/0) keeps the full likelihood vector.
fn pl_string_or_raw(pl: &[i32], raw: &[i32]) -> String {
    if pl.is_empty() { pl_string(raw) } else { pl_string(pl) }
}

/// Realignment flank (bp each side of the anchor).
const REALIGN_WINDOW: u32 = 20;
/// A read supports the indel when the alt haplotype is this many nats more
/// likely than the reference window (one decade).
const REALIGN_MARGIN: f64 = 2.3;

/// Reads whose window fits `halt` at least [`REALIGN_MARGIN`] better than `href`.
fn glocal_support(reads: &[LiveRead], lo: u32, hi: u32, href: &[u8], halt: &[u8]) -> u32 {
    let mut genuine = 0u32;
    for lr in reads {
        let Some((bases, quals)) = lr.query_window_qual(lo, hi) else { continue };
        if bases.is_empty() {
            continue;
        }
        let ll_ref = crate::call::pairhmm::read_vs_hap_loglik(&bases, &quals, href);
        let ll_alt = crate::call::pairhmm::read_vs_hap_loglik(&bases, &quals, halt);
        if ll_alt - ll_ref >= REALIGN_MARGIN {
            genuine += 1;
        }
    }
    genuine
}

/// Count reads whose window realigns better to ref+ins than ref; `u32::MAX` if realignment impossible.
fn realign_ins_support(reads: &[LiveRead], chr: &str, pos1: u64, ins: &[u8], ctx: &EmitCtx) -> u32 {
    let w = REALIGN_WINDOW;
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
    glocal_support(reads, lo, hi, href, &halt)
}

/// Count reads whose window realigns better to ref-with-deletion than ref; `u32::MAX` if impossible.
fn realign_del_support(reads: &[LiveRead], chr: &str, pos1: u64, del_len: u32, ctx: &EmitCtx) -> u32 {
    let w = REALIGN_WINDOW;
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
    glocal_support(reads, lo, ext_hi, href, &halt)
}

/// Emit one VCF record per indel allele (insertion or deletion). Each record
/// is fully normalized: REF includes the deleted bases (read from FASTA) for
/// deletions, REF is single-base for insertions.
#[allow(clippy::too_many_arguments)]
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
    overlapping_reads: &[LiveRead],
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
    for (seq, c, f) in &agg.ins_alleles {
        let (mn, af) = if ctx.hp_indel {
            indel_repeat_stringency(ctx.fasta, chr, pos1, seq.as_bytes(), ind_min, ind_af)
        } else { (ind_min, ind_af) };
        if *c < mn { continue; }
        let frac = (*c as f64) / denom;
        if frac < af { continue; }
        // Re-apply count/fraction thresholds to realignment-confirmed support.
        if ctx.realign >= 1 {
            let genuine = realign_ins_support(overlapping_reads, chr, pos1, seq.as_bytes(), ctx);
            if genuine != u32::MAX && (genuine < mn || (genuine as f64) / denom < af) {
                continue;
            }
        }
        let ref_str = (ref_base as char).to_string();
        let alt_str = format!("{}{}", ref_base as char, seq);
        if args.bcftools_layout {
            let iq = if ctx.hp_indel { hp_indel_qual(ctx, chr, pos1, &ref_str, &alt_str) } else { 30 };
            let spec = IndelSpec { ref_str: &ref_str, alt_str: &alt_str, is_ins: true, ins_seq: seq.as_bytes(), del_len: 0, indel_qual: iq };
            bcf_layout::render_indel(site, overlapping_reads, &ctx.kept_samples(site.per_sample.len()), chr, pos1, &spec, &ctx.layout(), out);
            continue;
        }
        emit_one_indel(site, agg, chr, pos1, &ref_str, &alt_str, *c, *f, total, mean_mq, ctx, out, annotate);
    }
    for (l, c, f) in &agg.del_alleles {
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
        if ctx.realign >= 2 {
            let genuine = realign_del_support(overlapping_reads, chr, pos1, *l, ctx);
            if genuine != u32::MAX && (genuine < mn || (genuine as f64) / denom < af) {
                continue;
            }
        }
        let mut ref_str = String::with_capacity(1 + deleted.len());
        ref_str.push(ref_base as char);
        for &b in deleted { ref_str.push(b as char); }
        let alt_str = (ref_base as char).to_string();
        if args.bcftools_layout {
            let iq = if ctx.hp_indel { hp_indel_qual(ctx, chr, pos1, &ref_str, &alt_str) } else { 30 };
            let spec = IndelSpec { ref_str: &ref_str, alt_str: &alt_str, is_ins: false, ins_seq: &[], del_len: *l, indel_qual: iq };
            bcf_layout::render_indel(site, overlapping_reads, &ctx.kept_samples(site.per_sample.len()), chr, pos1, &spec, &ctx.layout(), out);
            continue;
        }
        emit_one_indel(site, agg, chr, pos1, &ref_str, &alt_str, *c, *f, total, mean_mq, ctx, out, annotate);
    }
}

#[allow(clippy::too_many_arguments)]
fn emit_one_indel(
    site: &PileupSite,
    agg: &crate::bam::SampleSiteCounts,
    chr: &str,
    pos1: u64,
    ref_str: &str,
    alt_str: &str,
    indel_count: u32,
    indel_fwd: u32,
    total: u32,
    mean_mq: f64,
    ctx: &EmitCtx,
    out: &mut Vec<u8>,
    annotate: &AnnotateSpec,
) {
    let args = ctx.args;
    // Per sample, a synthetic biallelic pileup: reference reads at quality 30,
    // indel reads at `indel_qual`, strands from the forward counts; the
    // dependent-error model then gives PLs and the caller the site.
    let n_alleles = 2usize;
    let indel_qual: u8 = if ctx.hp_indel { hp_indel_qual(ctx, chr, pos1, ref_str, alt_str) } else { 30 };
    let is_ins = alt_str.len() > ref_str.len();
    let mut kept_samples: Vec<usize> = Vec::new();
    let mut per_counts: Vec<(u32, u32, u32, u32)> = Vec::new();
    let mut pls_all: Vec<i32> = Vec::new();
    let mut bases: Vec<u16> = Vec::new();
    let mut qs = [0.0f64; 2];
    for (si, s) in site.per_sample.iter().enumerate() {
        if ctx.any_filter && !ctx.sample_keep.get(si).copied().unwrap_or(false) { continue; }
        kept_samples.push(si);
        if s.depth == 0 {
            per_counts.push((0, 0, 0, 0));
            pls_all.extend([PL_MISSING; 3]);
            continue;
        }
        // Indel count (and forward-strand count) for THIS allele in THIS sample.
        let (indel_c, indel_f) = if is_ins {
            let seq = &alt_str[1..];
            s.ins_alleles.iter().find(|(k, _, _)| k == seq).map(|(_, c, f)| (*c, *f)).unwrap_or((0, 0))
        } else {
            let dl = (ref_str.len() - 1) as u32;
            s.del_alleles.iter().find(|(l, _, _)| *l == dl).map(|(_, c, f)| (*c, *f)).unwrap_or((0, 0))
        };
        let ref_c = s.depth.saturating_sub(indel_c);
        let ref_f_total: u32 = s.fwd_counts.iter().sum();
        let ref_f = ref_f_total.saturating_sub(indel_f).min(ref_c);
        bases.clear();
        bases.extend((0..ref_c).map(|i| pack_base(30, i >= ref_f, 0)));
        bases.extend((0..indel_c).map(|i| pack_base(indel_qual, i >= indel_f, 1)));
        let pl = ctx.em.pls(&mut bases, &[0, 1], &mut Drand48::for_site(ctx.seed, site.ref_id, site.pos, si));
        pls_all.extend(pl.iter().map(|&v| v as i32));
        per_counts.push((ref_c, indel_c, ref_f, indel_f));
        let (qr, qa) = (ref_c as f64 * 30.0, indel_c as f64 * indel_qual as f64);
        if qr + qa > 0.0 {
            qs[0] += qr / (qr + qa);
            qs[1] += qa / (qr + qa);
        }
    }
    let mut call_site = CallSite { n_samples: kept_samples.len(), n_alleles, pls: pls_all, is_indel: true, depths: None, qs: Some(qs.to_vec()), unseen: None, sample_af: None, prior_an_ac: None };
    let CallResult::Called { alleles_kept, qual: site_qual, gts, gqs, pls: pls_sub, .. } = ctx.caller.call_site(&mut call_site) else { return };
    let site_qual = site_qual.unwrap_or(0.0);
    let any_variant = alleles_kept.len() > 1;
    if args.variants_only && (!any_variant || site_qual < args.min_qual as f64) {
        return;
    }
    let mut per_sample_cols: Vec<String> = Vec::with_capacity(kept_samples.len());
    for (k, &si) in kept_samples.iter().enumerate() {
        let s = &site.per_sample[si];
        if s.depth == 0 {
            per_sample_cols.push(empty_sample_col(annotate, n_alleles));
            continue;
        }
        let (ref_c, indel_c, ref_f, indel_f) = per_counts[k];
        let counts = vec![ref_c, indel_c];
        let fwd = vec![ref_f, indel_f];
        let quals = vec![ref_c * 30, indel_c * indel_qual as u32];
        let gt = gt_string(&alleles_kept, gts[k]);
        let pl = pl_string_or_raw(&pls_sub[k], &call_site.pls[k * 3..(k + 1) * 3]);
        let parts = sample_parts(annotate, gt, s.depth, &counts, &fwd, &quals, s.n_softclip, &pl, gqs[k]);
        per_sample_cols.push(parts.join(":"));
    }

    let ref_total = total.saturating_sub(indel_count);
    let ref_fwd_total: u32 = agg.fwd_counts.iter().sum::<u32>().saturating_sub(indel_fwd).min(ref_total);
    let ad = vec![ref_total, indel_count];
    let adf = vec![ref_fwd_total, indel_fwd];
    let info = site_info(total, mean_mq, &ad, &adf, agg.n_softclip, annotate, true);
    let fmt_keys = format_keys(annotate);

    let qual_str = fmt_qual(site_qual);
    let _ = write!(out, "{}\t{}\t.\t{}\t{}\t{}\t.\t{}\t{}",
        chr, pos1, ref_str, alt_str, qual_str, info, fmt_keys.join(":"));
    for col in &per_sample_cols { let _ = write!(out, "\t{}", col); }
    let _ = writeln!(out);
}

/// `--hp-indel` quality for an indel allele: 30, minus 5 per repeat unit of
/// the indel sequence beyond two in the reference downstream of the anchor,
/// floored at 8 (polymerase slippage makes long tracts unreliable).
fn hp_indel_qual(ctx: &EmitCtx, chr: &str, pos1: u64, ref_str: &str, alt_str: &str) -> u8 {
    let unit: &[u8] = if alt_str.len() > ref_str.len() { &alt_str.as_bytes()[1..] } else { &ref_str.as_bytes()[1..] };
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
    let penalty = 5u32.saturating_mul(run_units.saturating_sub(2));
    30u32.saturating_sub(penalty).max(8) as u8
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

/// Emit an assembled indel: per-sample PLs come from read-vs-haplotype
/// likelihoods on the pair-HMM, the site from the multiallelic caller.
fn emit_assembled(
    call: crate::call::haplotype::AssembledCall,
    site: &PileupSite,
    reads: &[LiveRead],
    chr: &str,
    ctx: &EmitCtx,
    out: &mut Vec<u8>,
    annotate: &AnnotateSpec,
) {
    let args = ctx.args;
    let hs = haplotype_pls(reads, site.per_sample.len(), &call);
    let mut kept_samples: Vec<usize> = Vec::new();
    let mut pls_all: Vec<i32> = Vec::new();
    let mut qs = [0.0f64; 2];
    for si in 0..site.per_sample.len() {
        if ctx.any_filter && !ctx.sample_keep.get(si).copied().unwrap_or(false) { continue; }
        kept_samples.push(si);
        let h = hs[si];
        let depth = h.n_ref + h.n_alt;
        if depth == 0 {
            pls_all.extend([PL_MISSING; 3]);
        } else {
            pls_all.extend(h.pl.iter().map(|&v| v as i32));
            qs[0] += h.n_ref as f64 / depth as f64;
            qs[1] += h.n_alt as f64 / depth as f64;
        }
    }
    let mut call_site = CallSite { n_samples: kept_samples.len(), n_alleles: 2, pls: pls_all, is_indel: true, depths: None, qs: Some(qs.to_vec()), unseen: None, sample_af: None, prior_an_ac: None };
    let CallResult::Called { alleles_kept, qual, gts, gqs, pls: pls_sub, .. } = ctx.caller.call_site(&mut call_site) else { return };
    let qual = qual.unwrap_or(0.0);
    // Assembly only reports variant calls.
    if alleles_kept.len() < 2 || (args.variants_only && qual < args.min_qual as f64) { return; }
    let key = (call.pos1, call.ref_str.clone(), call.alt_str.clone());
    let dup = ASM_EMITTED.with(|e| {
        let mut s = e.borrow_mut();
        if s.contains(&key) { return true; }
        if s.len() > 200_000 { s.clear(); }
        s.insert(key);
        false
    });
    if dup { return; }
    let ref_c = call.total.saturating_sub(call.support);
    let alt_c = call.support;
    let mut info = format!("DP={};INDEL;ASSEMBLED", call.total);
    if annotate.info_ad { info.push_str(&format!(";AD={},{}", ref_c, alt_c)); }
    let fmt_keys = format_keys(annotate);
    let qual_str = fmt_qual(qual);
    let _ = write!(out, "{}\t{}\t.\t{}\t{}\t{}\t.\t{}\t{}",
        chr, call.pos1, call.ref_str, call.alt_str, qual_str, info, fmt_keys.join(":"));
    for (k, &si) in kept_samples.iter().enumerate() {
        let h = hs[si];
        let depth = h.n_ref + h.n_alt;
        if depth == 0 {
            let _ = write!(out, "\t{}", empty_sample_col(annotate, 2));
            continue;
        }
        // Strand counts are not tracked per haplotype; ADF/ADR/DP4/SP read as 0.
        let counts = vec![h.n_ref, h.n_alt];
        let fwd = vec![0u32, 0u32];
        let quals = vec![h.n_ref * 30, h.n_alt * 30];
        let gt = gt_string(&alleles_kept, gts[k]);
        let parts = sample_parts(annotate, gt, depth, &counts, &fwd, &quals, 0, &pl_string(&pls_sub[k]), gqs[k]);
        let _ = write!(out, "\t{}", parts.join(":"));
    }
    let _ = writeln!(out);
}

/// Wrapper for the sequential path when gvcf blocking is active. Returns false
/// when the site was emitted as a gVCF ref-block (caller should not append the
/// SNV record).
fn emit_site_or_gvcf<W: Write>(
    site: &PileupSite,
    overlapping_reads: &[LiveRead],
    ctx: &EmitCtx,
    buf: &mut Vec<u8>,
    gvcf: Option<&mut GvcfBlocker>,
    direct_out: &mut W,
) -> Result<bool> {
    let chr = ctx.ref_names.get(site.ref_id).map(|s| s.as_str()).unwrap_or(".");
    let pos1 = site.pos as u64 + 1;
    if !ctx.in_regions(chr, pos1) { return Ok(false); }
    let ref_base = ctx.fasta.and_then(|fa| fa.base(chr, pos1 as u32)).unwrap_or(b'N');
    if ref_base == b'N' { return Ok(false); }
    let agg = site.aggregated();
    let total: u32 = agg.base_counts.iter().sum();
    if total == 0 { return Ok(false); }
    // For gVCF, branch on whether any alts were observed at all.
    let any_nonref = (0..4).filter(|&i| {
        let rb = match ref_base { b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3, _ => 4 };
        i != rb && agg.base_counts[i] > 0
    }).count() > 0;
    if let Some(g) = gvcf {
        if !any_nonref {
            g.add_ref_site(chr, pos1, &(ref_base as char).to_string(), total, 0.0, direct_out)?;
            return Ok(false);
        }
        g.flush(direct_out)?;
    }
    emit_site(site, overlapping_reads, ctx, buf);
    Ok(true)
}

/// Repeat-context stringency for short indels. Homopolymer/STR tracts are the dominant source
/// of false-positive indels (polymerase slippage), so for a short inserted/deleted `unit` that
/// tiles a long repeat run immediately 3' of the anchor `pos1`, require more supporting reads and
/// a higher alt fraction. Returns the (possibly tightened) (min_count, min_frac). No-op for
/// non-repeat or long indels.
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
