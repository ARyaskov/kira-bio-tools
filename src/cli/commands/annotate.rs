use crate::annotate;
#[cfg(feature = "gpu")]
use crate::annotate::cpu_v2::{
    ColumnSpec, build_sample_map, extract_samples_from_headers, iter_ani_header_lines,
    merge_annotation_headers,
};
use crate::cli::args::{AnnotateArgs, AnnotateIndexArgs, AnnotateServeArgs, DbBuildArgs};
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
#[cfg(feature = "gpu")]
use std::time::Instant;

pub fn cmd_annotate(args: AnnotateArgs) -> Result<()> {
    let mut input_path = resolve_input_path(&args.input, args.cache_plain)?;
    let mut extra_header_lines = read_header_lines(args.header_lines.as_deref())?;
    for h in &args.header_line {
        let trimmed = h.trim();
        if trimmed.is_empty() { continue; }
        if !trimmed.starts_with("##") {
            anyhow::bail!("--header-line must start with ##: {}", trimmed);
        }
        extra_header_lines.push(trimmed.to_string());
    }

    let region_filter = build_region_filter(args.regions.as_deref(), args.regions_file.as_deref())?;
    if let Some(rf) = &region_filter {
        let filtered = pre_filter_regions(&input_path, rf)
            .context("pre-filtering by --regions")?;
        input_path = filtered;
    }

    let mut postproc = build_postproc(&args, &extra_header_lines)?;

    let out = args.output.clone().unwrap_or_else(|| {
        let mut p = args.input.clone();
        p.set_extension("annot.vcf");
        p
    });

    let output_kind = match args.output_type.as_deref() {
        Some(s) => Some(annotate::postproc::parse_output_type(s)?),
        None => None,
    };
    if matches!(output_kind, Some(annotate::postproc::OutputKind::Bcf(_))) {
        anyhow::bail!("-O u|b (BCF) not yet supported — use v or z");
    }
    let needs_bgzf = matches!(output_kind, Some(annotate::postproc::OutputKind::VcfGz(_))) || args.bgzf_after;
    let effective_bgzf_level = match output_kind {
        Some(annotate::postproc::OutputKind::VcfGz(lvl)) => Some(lvl),
        _ => args.bgzf_level,
    };
    let needs_postproc = !postproc_is_noop(&postproc)
        || args.include.is_some()
        || args.exclude.is_some()
        || args.samples.is_some()
        || args.samples_file.is_some();
    let intermediate_path: PathBuf = if needs_postproc || needs_bgzf {
        let mut tmp = out.clone();
        tmp.set_extension("kira-bt-tmp.vcf");
        tmp
    } else {
        out.clone()
    };
    let out_for_annotate = intermediate_path.clone();
    let bgzf_after: Option<PathBuf> = if needs_bgzf { Some(intermediate_path.clone()) } else { None };

    // ---- B2 phase 2 auto-detect / auto-build of .ktile sidecar ----
    //
    // Default behaviour: look for `<input>.ktile` next to the input VCF.
    // If fresh → use it (mmap, no BGZF decode). If missing or stale →
    // smart auto-build with a size guard.
    //
    // **Size guard**: the ktile reader uses mmap, which only beats BGZF
    // decode when the decompressed file fits in the OS page cache. For
    // multi-sample 1000 G-style inputs (1.3 GB compressed → 62 GB ktile
    // on chr1) the sidecar exceeds RAM and degrades to page-cache
    // thrashing — the build is pure overhead. The auto-skip heuristic
    // turns off auto-build when the source is larger than the threshold
    // (compressed > 500 MB by default, override with
    // `KIRA_BT_KTILE_MAX_INPUT_MB`).
    //
    // Opt-outs / overrides:
    //   --no-ktile           — disable ktile entirely (always BGZF+parse)
    //   --no-build-ktile     — use sidecar if present, don't build
    //   --force-build-ktile  — override the size guard, always build
    //
    // `UnifiedVcfReader::open` dispatches on the `.ktile` extension to
    // the mmap-based `KtileSourceReader`, so CPU / GPU / OpenCL backends
    // all benefit transparently — no per-backend wiring needed.
    let input_path = resolve_ktile_or_fallback(
        &input_path,
        args.no_ktile,
        !args.no_build_ktile,
        args.force_build_ktile,
    );

    let annotations = args.annotations.as_ref();
    let ani_path = if let Some(ani) = &args.ani {
        ani.clone()
    } else if let Some(ann) = annotations {
        if ann.extension().unwrap_or_default() == "ani" {
            ann.clone()
        } else {
            let mut p = ann.clone();
            p.set_extension("ani");
            p
        }
    } else {
        anyhow::bail!("Either --annotations or --ani must be provided");
    };

    if args.ani.is_none() {
        let Some(ann) = annotations else {
            anyhow::bail!("--annotations is required to build ANI index");
        };
        let ext = ann.extension().and_then(|e| e.to_str());
        if ext != Some("ani") {
            if should_rebuild_ani(ann, &ani_path, args.columns.as_deref())? {
                eprintln!("[annotate] Building ANI index from source...");
                if ext == Some("tab") {
                    annotate::build_ani_index_from_tab(ann, &ani_path, args.columns.as_deref())?;
                } else {
                    annotate::build_ani_index_auto_v2(ann, &ani_path)?;
                }
                write_ani_meta(ann, &ani_path, args.columns.as_deref())?;
            } else {
                eprintln!("[annotate] Using existing ANI index.");
            }
        } else if !ani_path.exists() {
            anyhow::bail!("ANI file not found: {:?}", ani_path);
        }
    } else if !ani_path.exists() {
        anyhow::bail!("ANI file not found: {:?}", ani_path);
    }

    eprintln!("[annotate] ANI = {:?}", ani_path);
    eprintln!("[annotate] Input = {:?}", input_path);
    eprintln!("[annotate] Output = {:?}", out);

    let mut columns: Vec<String> = if let Some(cols_str) = &args.columns {
        cols_str.split(',').map(|s| s.to_string()).collect()
    } else {
        Vec::new()
    };
    if let Some(cf) = &args.columns_file {
        let (file_cols, _types) = annotate::postproc::read_columns_file(cf)?;
        columns.extend(file_cols);
    }
    postproc.match_tags = columns.iter()
        .filter_map(|c| {
            let s = c.trim().trim_start_matches(|ch| ch == '+' || ch == '-' || ch == '=');
            let key = s.strip_prefix("INFO/").unwrap_or(s);
            if key.is_empty() || key == "CHROM" || key == "POS" || key == "ID" || key == "REF" || key == "ALT" || key == "QUAL" || key == "FILTER" || key == "~ID" { None }
            else { Some(key.to_string()) }
        })
        .collect();
    let bgzf_level = effective_bgzf_level;
    let mmap_output = args.mmap_output;
    let mmap_no_flush = args.mmap_no_flush;
    let ram_output = args.ram_output;
    let ram_max_mb = args.ram_max_mb;

    let mut ran = false;

    #[cfg(feature = "gpu")]
    if args.gpu && extra_header_lines.is_empty() {
        eprintln!("[annotate] Using CUDA GPU backend...");
        let ani = annotate::AniIndex::open(&ani_path)?;
        annotate::cuda::annotate_vcf_ani_gpu(
            &ani,
            &input_path,
            &out_for_annotate,
            &columns,
            bgzf_level,
            mmap_output,
            mmap_no_flush,
            ram_output,
            ram_max_mb,
        )?;
        ran = true;
    }

    if !ran {
        annotate::cpu_v2::annotate_vcf_ani_v2_with_extra_headers(
            &ani_path,
            &input_path,
            &out_for_annotate,
            &columns,
            bgzf_level,
            mmap_output,
            mmap_no_flush,
            ram_output,
            ram_max_mb,
            &extra_header_lines,
        )?;
        ran = true;
    }

    if ran {
        if needs_postproc {
            let final_or_bgzf_tmp: PathBuf = if needs_bgzf {
                let mut p = out.clone();
                p.set_extension("kira-bt-pp.vcf");
                p
            } else { out.clone() };
            apply_postproc_stream_with_args(&out_for_annotate, &final_or_bgzf_tmp, &args, postproc)?;
            let _ = std::fs::remove_file(&out_for_annotate);

            if needs_bgzf {
                compress_plain_to_bgzf(&final_or_bgzf_tmp, &out, effective_bgzf_level)?;
                let _ = std::fs::remove_file(&final_or_bgzf_tmp);
            }
        } else if let Some(tmp) = bgzf_after {
            compress_plain_to_bgzf(&tmp, &out, effective_bgzf_level)?;
            let _ = std::fs::remove_file(tmp);
        }

        if let Some(idx_fmt) = &args.write_index {
            write_output_index(&out, idx_fmt)?;
        }
    }
    Ok(())
}

fn postproc_is_noop(pp: &annotate::postproc::PostProcessor) -> bool {
    pp.remove.is_none() && pp.set_id.is_none() && pp.mark_sites.is_none()
        && pp.include.is_none() && pp.exclude.is_none() && !pp.keep_sites
        && pp.rename_chrs.is_none() && pp.rename_annots.is_none()
        && pp.samples_keep.is_none() && pp.extra_header_lines.is_empty()
        && !pp.no_version
}

fn build_postproc(
    args: &AnnotateArgs,
    extra_header_lines: &[String],
) -> Result<annotate::postproc::PostProcessor> {
    use annotate::postproc::PostProcessor;
    let mut pp = PostProcessor::default();
    pp.no_version = args.no_version;
    pp.extra_header_lines = extra_header_lines.to_vec();
    pp.keep_sites = args.keep_sites;

    if let Some(spec) = &args.remove {
        pp.remove = Some(PostProcessor::parse_remove(spec)?);
    }
    if let Some(spec) = &args.set_id {
        pp.set_id = Some(PostProcessor::parse_set_id(spec)?);
    }
    if let Some(spec) = &args.mark_sites {
        pp.mark_sites = Some(PostProcessor::parse_mark_sites(spec)?);
    }
    if let Some(p) = &args.rename_chrs {
        pp.rename_chrs = Some(PostProcessor::read_rename_chrs(p)?);
    }
    if let Some(p) = &args.rename_annots {
        pp.rename_annots = Some(PostProcessor::read_rename_annots(p)?);
    }
    if args.include.is_some() && args.exclude.is_some() {
        anyhow::bail!("-i/--include and -e/--exclude are mutually exclusive");
    }
    Ok(pp)
}

fn build_predicates_from_headers(
    args: &AnnotateArgs,
    headers: &[String],
    pp: &mut annotate::postproc::PostProcessor,
) -> Result<()> {
    use crate::filter::FilterEngine;
    if let Some(expr) = &args.include {
        let engine = FilterEngine::new(headers, Some(expr.as_str()), false)
            .context("-i/--include expression")?;
        pp.include = Some(annotate::postproc::Predicate { engine });
    }
    if let Some(expr) = &args.exclude {
        let engine = FilterEngine::new(headers, Some(expr.as_str()), false)
            .context("-e/--exclude expression")?;
        pp.exclude = Some(annotate::postproc::Predicate { engine });
    }
    Ok(())
}

fn build_region_filter(
    regions: Option<&str>,
    regions_file: Option<&Path>,
) -> Result<Option<annotate::postproc::RegionFilter>> {
    use annotate::postproc::RegionFilter;
    if let Some(s) = regions {
        return Ok(Some(RegionFilter::from_cli(s)?));
    }
    if let Some(p) = regions_file {
        return Ok(Some(RegionFilter::from_file(p)?));
    }
    Ok(None)
}

fn pre_filter_regions(input: &Path, rf: &annotate::postproc::RegionFilter) -> Result<PathBuf> {
    use crate::vcf::UnifiedVcfReader;
    let mut tmp = std::env::temp_dir();
    tmp.push(format!(
        "kira-bt-regions-{}-{}.vcf",
        std::process::id(),
        input.file_name().and_then(|n| n.to_str()).unwrap_or("input")
    ));
    let mut reader = UnifiedVcfReader::open(input).context("open input for region pre-filter")?;
    let headers = reader.header()?;
    let mut out = std::io::BufWriter::with_capacity(8 * 1024 * 1024, File::create(&tmp)?);
    for h in &headers {
        out.write_all(h.as_bytes())?;
        out.write_all(b"\n")?;
    }
    while let Some(line) = reader.read_line()? {
        if rf.line_passes(&line) {
            out.write_all(line.as_bytes())?;
            out.write_all(b"\n")?;
        }
    }
    out.flush()?;
    Ok(tmp)
}

fn apply_postproc_stream_with_args(
    input: &Path,
    output: &Path,
    args: &AnnotateArgs,
    mut pp: annotate::postproc::PostProcessor,
) -> Result<()> {
    use annotate::postproc::{HeaderOptions, apply_to_header, version_header_line};

    let reader = BufReader::with_capacity(8 * 1024 * 1024, File::open(input)?);
    let mut writer = std::io::BufWriter::with_capacity(8 * 1024 * 1024, File::create(output)?);

    let mut headers: Vec<String> = Vec::new();
    let mut iter = reader.lines();
    let mut first_record: Option<String> = None;
    for line in iter.by_ref() {
        let line = line?;
        if line.starts_with('#') {
            headers.push(line);
        } else {
            first_record = Some(line);
            break;
        }
    }

    build_predicates_from_headers(args, &headers, &mut pp)?;
    let input_samples = extract_samples_from_chrom_line(&headers);
    if !input_samples.is_empty() {
        if let Some(spec) = &args.samples {
            let (names, inverse) = annotate::postproc::parse_samples_cli(spec);
            pp.samples_keep = Some(annotate::postproc::resolve_samples_keep(&input_samples, &names, inverse));
        } else if let Some(path) = &args.samples_file {
            let (names, inverse) = annotate::postproc::read_samples_file(path)?;
            pp.samples_keep = Some(annotate::postproc::resolve_samples_keep(&input_samples, &names, inverse));
        }
    }

    let version = version_header_line();
    let opts = HeaderOptions {
        no_version: pp.no_version,
        extra_header_lines: &pp.extra_header_lines,
        remove: pp.remove.as_ref(),
        rename_chrs: pp.rename_chrs.as_ref(),
        rename_annots: pp.rename_annots.as_ref(),
        mark_sites: pp.mark_sites.as_ref(),
        set_id: pp.set_id.is_some(),
        samples_keep: pp.samples_keep.as_deref(),
        version_line: Some(&version),
    };
    let out_headers = apply_to_header(headers, &opts);
    for h in &out_headers {
        writer.write_all(h.as_bytes())?;
        writer.write_all(b"\n")?;
    }

    if let Some(line) = first_record {
        write_processed_line(&mut writer, &line, &pp)?;
    }
    for line in iter {
        let line = line?;
        write_processed_line(&mut writer, &line, &pp)?;
    }
    writer.flush()?;
    Ok(())
}

fn write_processed_line<W: Write>(
    writer: &mut W,
    line: &str,
    pp: &annotate::postproc::PostProcessor,
) -> Result<()> {
    use annotate::postproc::{LineAction, process_record_line};
    match process_record_line(line, pp, false) {
        LineAction::Keep => { writer.write_all(line.as_bytes())?; writer.write_all(b"\n")?; }
        LineAction::Replace(s) => { writer.write_all(s.as_bytes())?; writer.write_all(b"\n")?; }
        LineAction::Drop => {}
    }
    Ok(())
}

fn extract_samples_from_chrom_line(headers: &[String]) -> Vec<String> {
    for h in headers {
        if h.starts_with("#CHROM") {
            let cols: Vec<&str> = h.split('\t').collect();
            if cols.len() > 9 {
                return cols[9..].iter().map(|s| s.to_string()).collect();
            }
        }
    }
    Vec::new()
}

fn write_output_index(out: &Path, fmt: &str) -> Result<()> {
    let ext = if fmt.eq_ignore_ascii_case("tbi") { "tbi" } else { "csi" };
    let idx = std::path::PathBuf::from(format!("{}.{}", out.display(), ext));
    match crate::csi::build_csi_index(out, &idx) {
        Ok(_) => Ok(()),
        Err(e) => {
            eprintln!("[annotate] -W {}: {}", fmt, e);
            Ok(())
        }
    }
}

/// Default ceiling (compressed input MB) for the .ktile auto-build size
/// guard. With the new sliding-window pool reader, the previous
/// "RAM thrashing" reason for the guard is gone — the guard is now
/// purely a *disk-space* protection: a 60 GB sidecar for a single
/// 1.3 GB compressed input is wasteful unless you re-annotate often.
/// Default 2 GB compressed (≈ 60-80 GB decompressed sidecar). Override
/// with `KIRA_BT_KTILE_MAX_INPUT_MB`.
const KTILE_AUTOBUILD_MAX_INPUT_MB_DEFAULT: u64 = 2 * 1024;

/// Returns the auto-build size threshold in bytes. Reads
/// `KIRA_BT_KTILE_MAX_INPUT_MB` env var (parsed as MB) or falls back to
/// the default 500 MB. Used by [`resolve_ktile_or_fallback`] to decide
/// whether a compressed input is small enough that an auto-built
/// `.ktile` sidecar would actually fit in OS page cache.
fn ktile_autobuild_threshold_bytes() -> u64 {
    let mb = std::env::var("KIRA_BT_KTILE_MAX_INPUT_MB")
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(KTILE_AUTOBUILD_MAX_INPUT_MB_DEFAULT);
    mb.saturating_mul(1024 * 1024)
}

/// B2 phase 2 helper: decide whether to feed the annotate pipeline from
/// the original VCF or from a `<input>.ktile` sidecar.
///
/// Returns the path that `UnifiedVcfReader::open` should be called with
/// (which dispatches on the `.ktile` extension internally). Behaviour:
///   - `no_ktile = true`: always returns original `input_path`.
///   - sidecar absent + `build_if_missing` + small input: builds + uses it.
///   - sidecar absent + `build_if_missing` + huge input (above size
///     guard) + !force_build: skip auto-build with informative message.
///   - sidecar absent + `force_build`: always build regardless of size.
///   - sidecar absent + `!build_if_missing`: original `input_path`.
///   - sidecar present + fresh: ktile path.
///   - sidecar present + stale + `build_if_missing` + small: rebuild.
///   - sidecar present + stale + `build_if_missing` + huge + !force: VCF.
///   - sidecar present + stale + `!build_if_missing`: original `input_path`.
///   - sidecar present + freshness-check Err: original `input_path`.
fn resolve_ktile_or_fallback(
    input_path: &std::path::Path,
    no_ktile: bool,
    build_if_missing: bool,
    force_build: bool,
) -> std::path::PathBuf {
    use crate::annotate::ktile::{
        KtileFreshness, check_ktile_freshness, ktile_path_for, write_ktile_from_vcf,
    };

    if no_ktile {
        return input_path.to_path_buf();
    }

    let ktile_path = ktile_path_for(input_path);

    // Pre-build size guard (relaxed). The ktile reader now adapts —
    // small ktiles get whole-file mmap (fastest), big ktiles use a
    // sliding-window pool reader that bounds RAM at the chunk size
    // (default 2 GB). So even multi-GB ktiles process cleanly without
    // page-cache thrashing.
    //
    // The size guard remains as a *disk-space* guard: if the user's
    // input is, say, 1.3 GB compressed → 60 GB ktile, building it
    // eats 60 GB of disk for a sidecar that may save 30 % per run.
    // Default ceiling 2 GB compressed (≈ 80 GB decompressed); above
    // this the auto-build is skipped to protect disk space. Force
    // with `--force-build-ktile` or raise `KIRA_BT_KTILE_MAX_INPUT_MB`.
    let should_autobuild = |reason: &str| -> bool {
        if force_build {
            eprintln!(
                "[annotate] --force-build-ktile: building .ktile despite size guard ({reason})"
            );
            return true;
        }
        let threshold = ktile_autobuild_threshold_bytes();
        let input_size = std::fs::metadata(input_path).map(|m| m.len()).unwrap_or(0);
        if input_size > threshold {
            let mb = input_size / (1024 * 1024);
            let limit_mb = threshold / (1024 * 1024);
            eprintln!(
                "[annotate] Skipping .ktile auto-build: input is {} MB compressed (> {} MB guard); \
                 sidecar would consume substantial disk space (~30-50× compressed size). \
                 Reader-side memory is fine — the sliding-window pool keeps RAM bounded. \
                 To build anyway: --force-build-ktile, or raise KIRA_BT_KTILE_MAX_INPUT_MB to {}. \
                 Using VCF path for this run.",
                mb,
                limit_mb,
                mb + 1
            );
            return false;
        }
        true
    };

    // Helper: actually build the sidecar, log progress + outcome, and
    // return the ktile_path on success or the VCF path on failure.
    let do_build = |reason: &str| -> std::path::PathBuf {
        eprintln!(
            "[annotate] Building .ktile sidecar at {:?} ({}) ...",
            ktile_path, reason
        );
        match write_ktile_from_vcf(input_path, &ktile_path) {
            Ok(stats) => {
                let mb = stats.bytes_written as f64 / (1024.0 * 1024.0);
                eprintln!(
                    "[annotate] Built .ktile with {} records ({:.1} MB) — \
                     subsequent annotate runs against this input will skip BGZF decode.",
                    stats.n_records, mb
                );
                ktile_path.clone()
            }
            Err(e) => {
                eprintln!(
                    "[annotate] Warning: .ktile build failed ({}), falling back to VCF",
                    e
                );
                input_path.to_path_buf()
            }
        }
    };

    if !ktile_path.exists() {
        if build_if_missing && should_autobuild("missing") {
            return do_build("missing, auto-building");
        }
        if !build_if_missing {
            eprintln!(
                "[annotate] No .ktile sidecar (auto-build disabled); using VCF."
            );
        }
        return input_path.to_path_buf();
    }

    match check_ktile_freshness(&ktile_path, input_path) {
        Ok(KtileFreshness::Fresh) => {
            eprintln!(
                "[annotate] Using .ktile sidecar (fresh): {:?}",
                ktile_path
            );
            ktile_path
        }
        Ok(KtileFreshness::StaleSize { ktile, source }) => {
            if build_if_missing && should_autobuild("stale: size mismatch") {
                do_build(&format!(
                    "stale: size changed {} → {}, auto-rebuilding",
                    ktile, source
                ))
            } else {
                eprintln!(
                    "[annotate] .ktile is stale (size {} → {}); using VCF. \
                     Drop --no-build-ktile to auto-rebuild.",
                    ktile, source
                );
                input_path.to_path_buf()
            }
        }
        Ok(KtileFreshness::StaleMtime { ktile, source }) => {
            if build_if_missing && should_autobuild("stale: mtime mismatch") {
                do_build(&format!(
                    "stale: mtime changed {} → {}, auto-rebuilding",
                    ktile, source
                ))
            } else {
                eprintln!(
                    "[annotate] .ktile is stale (mtime {} → {}); using VCF. \
                     Drop --no-build-ktile to auto-rebuild.",
                    ktile, source
                );
                input_path.to_path_buf()
            }
        }
        Ok(KtileFreshness::Unknown) => {
            if build_if_missing && should_autobuild("source metadata missing") {
                do_build("source metadata missing, auto-rebuilding")
            } else {
                eprintln!(
                    "[annotate] .ktile lacks source metadata; using VCF (drop --no-build-ktile to rebuild)."
                );
                input_path.to_path_buf()
            }
        }
        Err(e) => {
            eprintln!(
                "[annotate] .ktile freshness check failed ({}); using VCF.",
                e
            );
            input_path.to_path_buf()
        }
    }
}

pub fn cmd_annotate_serve(args: AnnotateServeArgs) -> Result<()> {
    let use_gpu = args.gpu;
    let use_cpu = !use_gpu;

    #[cfg(not(feature = "gpu"))]
    if use_gpu {
        anyhow::bail!("GPU feature not enabled");
    }

    let annotations = args.annotations.as_ref();
    let ani_path = if let Some(ani) = &args.ani {
        ani.clone()
    } else if let Some(ann) = annotations {
        if ann.extension().unwrap_or_default() == "ani" {
            ann.clone()
        } else {
            let mut p = ann.clone();
            p.set_extension("ani");
            p
        }
    } else {
        anyhow::bail!("Either --annotations or --ani must be provided");
    };

    if args.ani.is_none() {
        let Some(ann) = annotations else {
            anyhow::bail!("--annotations is required to build ANI index");
        };
        let ext = ann.extension().and_then(|e| e.to_str());
        if ext != Some("ani") {
            if should_rebuild_ani(ann, &ani_path, args.columns.as_deref())? {
                eprintln!("[annotate] Building ANI index from source...");
                if ext == Some("tab") {
                    annotate::build_ani_index_from_tab(ann, &ani_path, args.columns.as_deref())?;
                } else {
                    annotate::build_ani_index_auto_v2(ann, &ani_path)?;
                }
                write_ani_meta(ann, &ani_path, args.columns.as_deref())?;
            } else {
                eprintln!("[annotate] Using existing ANI index.");
            }
        } else if !ani_path.exists() {
            anyhow::bail!("ANI file not found: {:?}", ani_path);
        }
    } else if !ani_path.exists() {
        anyhow::bail!("ANI file not found: {:?}", ani_path);
    }

    #[cfg(feature = "gpu")]
    let ani = annotate::AniIndex::open(&ani_path)?;

    #[cfg(feature = "gpu")]
    let mut cuda_state: Option<annotate::cuda::GpuAnnotator> = if use_gpu {
        let start = Instant::now();
        let state = annotate::cuda::GpuAnnotator::new(&ani)?;
        eprintln!("[gpu] warmup: {:.3}s", start.elapsed().as_secs_f64());
        Some(state)
    } else {
        None
    };

    let default_columns = parse_columns(args.columns.as_deref());
    let mut stdout = std::io::stdout();
    let stdin = std::io::stdin();

    for line in stdin.lock().lines() {
        let line = line?;
        let line = line.trim_end();
        if line.is_empty() {
            continue;
        }
        let cmd = line.trim();
        if cmd.eq_ignore_ascii_case("quit") || cmd.eq_ignore_ascii_case("exit") {
            break;
        }

        let tokens = if line.contains('\t') {
            line.split('\t').collect::<Vec<_>>()
        } else {
            line.split_whitespace().collect::<Vec<_>>()
        };
        if tokens.len() < 2 {
            writeln!(stdout, "ERR\tmissing input/output")?;
            stdout.flush()?;
            continue;
        }

        let mut kv_start = tokens.len();
        for (i, t) in tokens.iter().enumerate() {
            if t.contains('=') {
                kv_start = i;
                break;
            }
        }
        let (path_tokens, kv_tokens) = tokens.split_at(kv_start);
        if path_tokens.len() < 2 || path_tokens.len() % 2 != 0 {
            writeln!(stdout, "ERR\tneed input/output pairs")?;
            stdout.flush()?;
            continue;
        }

        let mut columns = default_columns.clone();
        let mut bgzf_level = args.bgzf_level;
        let mut cache_plain = args.cache_plain;
        let mut bgzf_after = args.bgzf_after;
        let mut mmap_output = args.mmap_output;
        let mut mmap_no_flush = args.mmap_no_flush;
        let mut ram_output = args.ram_output;
        let mut ram_max_mb = args.ram_max_mb;

        for kv in kv_tokens {
            if let Some((k, v)) = kv.split_once('=') {
                let key = k.trim();
                let val = v.trim();
                match key {
                    "c" => {
                        columns = parse_columns(Some(val));
                    }
                    "columns" => {
                        columns = parse_columns(Some(val));
                    }
                    "bgzf_level" => {
                        if let Ok(n) = val.parse::<u32>() {
                            bgzf_level = Some(n);
                        }
                    }
                    "cache_plain" => {
                        if let Some(b) = parse_bool(val) {
                            cache_plain = b;
                        }
                    }
                    "bgzf_after" => {
                        if let Some(b) = parse_bool(val) {
                            bgzf_after = b;
                        }
                    }
                    "mmap_output" => {
                        if let Some(b) = parse_bool(val) {
                            mmap_output = b;
                        }
                    }
                    "mmap_no_flush" => {
                        if let Some(b) = parse_bool(val) {
                            mmap_no_flush = b;
                        }
                    }
                    "ram_output" => {
                        if let Some(b) = parse_bool(val) {
                            ram_output = b;
                        }
                    }
                    "ram_max_mb" => {
                        if let Ok(n) = val.parse::<u32>() {
                            ram_max_mb = n;
                        }
                    }
                    _ => {}
                }
            }
        }

        let mut input_outputs = Vec::with_capacity(path_tokens.len() / 2);
        let mut it = path_tokens.iter();
        while let (Some(input), Some(output)) = (it.next(), it.next()) {
            input_outputs.push((PathBuf::from(*input), PathBuf::from(*output)));
        }

        if use_gpu && input_outputs.len() > 1 {
            #[cfg(feature = "gpu")]
            {
                let mut jobs = Vec::with_capacity(input_outputs.len());
                let mut bgzf_after_jobs = Vec::with_capacity(input_outputs.len());

                for (input, output) in &input_outputs {
                    let input_path = match resolve_input_path(input, cache_plain) {
                        Ok(v) => v,
                        Err(e) => {
                            writeln!(stdout, "ERR\t{}", e)?;
                            stdout.flush()?;
                            jobs.clear();
                            break;
                        }
                    };
                    let (out_for_annotate, bgzf_after_tmp) = match resolve_output_path(output, true)
                    {
                        Ok(v) => v,
                        Err(e) => {
                            writeln!(stdout, "ERR\t{}", e)?;
                            stdout.flush()?;
                            jobs.clear();
                            break;
                        }
                    };

                    let use_bgzf = false;

                    let input_reader = annotate::VcfAnnotationReader::open(&input_path)?;
                    let streaming_reader = annotate::StreamingVcfReader::new(input_reader);
                    let (headers, _reader) = streaming_reader.into_headers_and_self()?;
                    let ani_headers = iter_ani_header_lines(&ani);
                    let merged_headers = merge_annotation_headers(
                        &headers,
                        &ani_headers,
                        &ColumnSpec::parse_all(&columns),
                    )?;
                    let input_samples = extract_samples_from_headers(&headers);
                    let db_samples = extract_samples_from_headers(&ani_headers);
                    let sample_map = build_sample_map(&input_samples, &db_samples);

                    jobs.push(annotate::cuda::GpuJob {
                        input: input_path,
                        output: out_for_annotate,
                        use_bgzf,
                        headers: merged_headers,
                        sample_map: std::sync::Arc::new(sample_map),
                    });
                    bgzf_after_jobs.push((bgzf_after_tmp, output.clone()));
                }

                if !jobs.is_empty() {
                    let result = annotate::cuda::annotate_vcf_ani_gpu_multi_with_state(
                        cuda_state.as_mut().unwrap(),
                        &ani,
                        jobs,
                        &columns,
                        bgzf_level,
                        mmap_output,
                        mmap_no_flush,
                        ram_output,
                        ram_max_mb,
                        2,
                    );

                    let result = match result {
                        Ok(()) => {
                            for (tmp, output) in bgzf_after_jobs {
                                if let Some(tmp) = tmp {
                                    compress_plain_to_bgzf(&tmp, &output, bgzf_level)?;
                                    let _ = std::fs::remove_file(tmp);
                                }
                            }
                            Ok(())
                        }
                        Err(e) => Err(e),
                    };

                    match result {
                        Ok(()) => {
                            writeln!(stdout, "OK\tmulti")?;
                        }
                        Err(e) => {
                            writeln!(stdout, "ERR\t{}", e)?;
                        }
                    }
                    stdout.flush()?;
                }
            }
            #[cfg(not(feature = "gpu"))]
            {
                writeln!(stdout, "ERR\tGPU feature not enabled")?;
                stdout.flush()?;
            }
            continue;
        }

        let (input, output) = (&input_outputs[0].0, &input_outputs[0].1);
        let input_path = match resolve_input_path(input, cache_plain) {
            Ok(v) => v,
            Err(e) => {
                writeln!(stdout, "ERR\t{}", e)?;
                stdout.flush()?;
                continue;
            }
        };
        let (out_for_annotate, bgzf_after_tmp) = match resolve_output_path(output, bgzf_after) {
            Ok(v) => v,
            Err(e) => {
                writeln!(stdout, "ERR\t{}", e)?;
                stdout.flush()?;
                continue;
            }
        };

        let result = if use_gpu {
            #[cfg(feature = "gpu")]
            {
                annotate::cuda::annotate_vcf_ani_gpu_with_state(
                    cuda_state.as_mut().unwrap(),
                    &ani,
                    &input_path,
                    &out_for_annotate,
                    &columns,
                    bgzf_level,
                    mmap_output,
                    mmap_no_flush,
                    ram_output,
                    ram_max_mb,
                )
            }
            #[cfg(not(feature = "gpu"))]
            {
                anyhow::bail!("GPU feature not enabled")
            }
        } else if use_cpu {
            annotate::cpu_v2::annotate_vcf_ani_v2(
                &ani_path,
                &input_path,
                &out_for_annotate,
                &columns,
                bgzf_level,
                mmap_output,
                mmap_no_flush,
                ram_output,
                ram_max_mb,
            )
        } else {
            anyhow::bail!("No backend selected")
        };

        let result = match result {
            Ok(()) => {
                if let Some(tmp) = bgzf_after_tmp {
                    if let Err(e) = compress_plain_to_bgzf(&tmp, output, bgzf_level) {
                        Err(e)
                    } else {
                        let _ = std::fs::remove_file(tmp);
                        Ok(())
                    }
                } else {
                    Ok(())
                }
            }
            Err(e) => Err(e),
        };

        match result {
            Ok(()) => {
                writeln!(stdout, "OK\t{}", output.display())?;
            }
            Err(e) => {
                writeln!(stdout, "ERR\t{}", e)?;
            }
        }
        stdout.flush()?;
    }
    Ok(())
}

fn resolve_input_path(input: &Path, cache_plain: bool) -> Result<PathBuf> {
    if !cache_plain {
        return Ok(input.to_path_buf());
    }

    let ext = input.extension().and_then(|e| e.to_str()).unwrap_or("");
    if !matches!(ext, "gz" | "bgz" | "bgzf") {
        return Ok(input.to_path_buf());
    }

    let stem = input
        .file_stem()
        .and_then(|s| s.to_str())
        .context("Input file has no stem")?;
    let mut cached = std::env::current_dir().context("Failed to get CWD")?;
    cached.push(stem);

    if cached.exists() {
        return Ok(cached);
    }

    let input_file = File::open(input).with_context(|| format!("Failed to open {:?}", input))?;
    let mut decoder = flate2::read::MultiGzDecoder::new(input_file);
    let mut out =
        File::create(&cached).with_context(|| format!("Failed to create {:?}", cached))?;
    std::io::copy(&mut decoder, &mut out).context("Failed to decompress input")?;
    out.flush().ok();

    Ok(cached)
}

fn parse_columns(s: Option<&str>) -> Vec<String> {
    if let Some(cols_str) = s {
        if cols_str.trim().is_empty() {
            Vec::new()
        } else {
            cols_str.split(',').map(|v| v.to_string()).collect()
        }
    } else {
        Vec::new()
    }
}

fn read_header_lines(path: Option<&Path>) -> Result<Vec<String>> {
    let Some(path) = path else {
        return Ok(Vec::new());
    };
    let file = File::open(path).with_context(|| format!("Failed to open {:?}", path))?;
    let reader = BufReader::new(file);
    let mut lines = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if !trimmed.starts_with("##") {
            anyhow::bail!("Header line must start with ##: {}", trimmed);
        }
        lines.push(trimmed.to_string());
    }
    Ok(lines)
}

fn should_rebuild_ani(source: &Path, ani: &Path, columns: Option<&str>) -> Result<bool> {
    if !ani.exists() {
        return Ok(true);
    }
    let source_meta = std::fs::metadata(source)
        .with_context(|| format!("Failed to stat annotation source {:?}", source))?;
    let ani_meta =
        std::fs::metadata(ani).with_context(|| format!("Failed to stat ANI index {:?}", ani))?;
    if let (Ok(src_time), Ok(ani_time)) = (source_meta.modified(), ani_meta.modified()) {
        if src_time > ani_time {
            return Ok(true);
        }
    }
    if source.extension().and_then(|e| e.to_str()) == Some("tab") {
        return Ok(!ani_meta_matches(source, ani, columns, &source_meta)?);
    }
    Ok(false)
}

fn ani_meta_matches(
    source: &Path,
    ani: &Path,
    columns: Option<&str>,
    source_meta: &std::fs::Metadata,
) -> Result<bool> {
    let meta_path = ani.with_extension("ani.meta");
    if !meta_path.exists() {
        return Ok(false);
    }
    let text = std::fs::read_to_string(&meta_path)
        .with_context(|| format!("Failed to read {:?}", meta_path))?;
    let expected_path = source.canonicalize().unwrap_or_else(|_| source.to_path_buf());
    let expected_len = source_meta.len().to_string();
    let expected_columns = columns.unwrap_or("");
    let mut path_ok = false;
    let mut len_ok = false;
    let mut columns_ok = false;
    for line in text.lines() {
        if let Some(v) = line.strip_prefix("source=") {
            path_ok = v == expected_path.to_string_lossy();
        } else if let Some(v) = line.strip_prefix("source_len=") {
            len_ok = v == expected_len;
        } else if let Some(v) = line.strip_prefix("columns=") {
            columns_ok = v == expected_columns;
        }
    }
    Ok(path_ok && len_ok && columns_ok)
}

fn write_ani_meta(source: &Path, ani: &Path, columns: Option<&str>) -> Result<()> {
    let source_meta = std::fs::metadata(source)
        .with_context(|| format!("Failed to stat annotation source {:?}", source))?;
    let source_path = source.canonicalize().unwrap_or_else(|_| source.to_path_buf());
    let meta_path = ani.with_extension("ani.meta");
    let text = format!(
        "source={}\nsource_len={}\ncolumns={}\n",
        source_path.to_string_lossy(),
        source_meta.len(),
        columns.unwrap_or("")
    );
    std::fs::write(&meta_path, text).with_context(|| format!("Failed to write {:?}", meta_path))
}

fn parse_bool(s: &str) -> Option<bool> {
    match s.trim().to_ascii_lowercase().as_str() {
        "1" | "true" | "yes" | "y" => Some(true),
        "0" | "false" | "no" | "n" => Some(false),
        _ => None,
    }
}

fn resolve_output_path(output: &Path, bgzf_after: bool) -> Result<(PathBuf, Option<PathBuf>)> {
    if !bgzf_after {
        return Ok((output.to_path_buf(), None));
    }
    let ext = output.extension().and_then(|e| e.to_str()).unwrap_or("");
    if !matches!(ext, "gz" | "bgz" | "bgzf") {
        return Ok((output.to_path_buf(), None));
    }
    let stem = output
        .file_stem()
        .and_then(|s| s.to_str())
        .context("Output file has no stem")?;
    let mut tmp = output
        .parent()
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| PathBuf::from("."));
    tmp.push(format!("{}.tmp.vcf", stem));
    let tmp_clone = tmp.clone();
    Ok((tmp, Some(tmp_clone)))
}

fn compress_plain_to_bgzf(input: &Path, output: &Path, level: Option<u32>) -> Result<()> {
    let timing = std::env::var("KIRA_BT_TIMING").is_ok();
    let start = std::time::Instant::now();
    if timing {
        eprintln!("[bgzf-after] start: input={:?}, output={:?}", input, output);
    }
    let file = File::open(input).with_context(|| format!("Failed to open {:?}", input))?;
    let mut reader = BufReader::with_capacity(8 * 1024 * 1024, file);
    let lvl = level.unwrap_or(1).min(9);
    if timing {
        eprintln!("[bgzf-after] level: {}", lvl);
    }
    let mut writer =
        crate::bgzf::BgzfWriter::with_compression(output, flate2::Compression::new(lvl))?;
    std::io::copy(&mut reader, &mut writer).context("Failed to write BGZF")?;
    writer.finish().context("Failed to finalize BGZF")?;
    if timing {
        let input_len = std::fs::metadata(input).map(|m| m.len()).unwrap_or(0);
        let output_len = std::fs::metadata(output).map(|m| m.len()).unwrap_or(0);
        eprintln!(
            "[bgzf-after] done: {:.3}s, input={} bytes, output={} bytes",
            start.elapsed().as_secs_f64(),
            input_len,
            output_len
        );
    }
    Ok(())
}

pub fn cmd_annotate_index(args: AnnotateIndexArgs) -> Result<()> {
    let out = args.output.clone().unwrap_or_else(|| {
        let mut p = args.input.clone();
        p.set_extension("ani");
        p
    });

    eprintln!("[annotate-index] Input  = {:?}", args.input);
    eprintln!("[annotate-index] Output = {:?}", out);

    let ext = args.input.extension().and_then(|e| e.to_str());

    if ext == Some("tab") {
        annotate::build_ani_index_from_tab(&args.input, &out, None)?;
    } else {
        annotate::build_ani_index_auto_v2(&args.input, &out)?;
    }
    write_ani_meta(&args.input, &out, None)?;

    Ok(())
}

pub fn cmd_db_build(args: DbBuildArgs) -> Result<()> {
    let out = args.output.clone().unwrap_or_else(|| {
        let mut p = args.input.clone();
        p.set_extension("ani");
        p
    });

    eprintln!("[db-build] Input: {:?}", args.input);
    eprintln!("[db-build] Output: {:?}", out);

    let ext = args.input.extension().and_then(|e| e.to_str());

    if ext == Some("tab") {
        annotate::build_ani_index_from_tab(&args.input, &out, None)?;
    } else {
        annotate::build_ani_index_auto_v2(&args.input, &out)?;
    }
    write_ani_meta(&args.input, &out, None)?;

    eprintln!("[db-build] Done");
    Ok(())
}
