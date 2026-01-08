use crate::annotate;
use crate::annotate::cpu_v2::{
    build_sample_map, extract_samples_from_headers, iter_ani_header_lines,
    merge_annotation_headers, ColumnSpec,
};
use crate::cli::args::{AnnotateArgs, AnnotateIndexArgs, AnnotateServeArgs, DbBuildArgs};
use crate::detect_format;
use crate::VcfFormat;
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

pub fn cmd_annotate(args: AnnotateArgs) -> Result<()> {
    let input_path = resolve_input_path(&args.input, args.cache_plain)?;
    let out = args.output.clone().unwrap_or_else(|| {
        let mut p = args.input.clone();
        p.set_extension("annot.vcf");
        p
    });
    let (out_for_annotate, bgzf_after) = resolve_output_path(&out, args.bgzf_after)?;

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

    if args.ani.is_none() && !ani_path.exists() {
        eprintln!("[annotate] ANI index not found, building from source...");

        let Some(ann) = annotations else {
            anyhow::bail!("--annotations is required to build ANI index");
        };
        let ext = ann.extension().and_then(|e| e.to_str());

        if ext == Some("tab") {
            annotate::build_ani_index_from_tab(ann, &ani_path, args.columns.as_deref())?;
        } else {
            annotate::build_ani_index_auto_v2(ann, &ani_path)?;
        }
    } else if args.ani.is_some() && !ani_path.exists() {
        anyhow::bail!("ANI file not found: {:?}", ani_path);
    }

    eprintln!("[annotate] ANI = {:?}", ani_path);
    eprintln!("[annotate] Input = {:?}", input_path);
    eprintln!("[annotate] Output = {:?}", out);

    let columns: Vec<String> = if let Some(cols_str) = &args.columns {
        cols_str.split(',').map(|s| s.to_string()).collect()
    } else {
        Vec::new()
    };
    let bgzf_level = args.bgzf_level;
    let mmap_output = args.mmap_output;
    let mmap_no_flush = args.mmap_no_flush;
    let ram_output = args.ram_output;
    let ram_max_mb = args.ram_max_mb;

    let mut ran = false;

    #[cfg(feature = "gpu")]
    if args.gpu {
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

    #[cfg(feature = "opencl")]
    if !ran && args.opencl {
        eprintln!("[annotate] Using OpenCL backend...");
        let ani = annotate::AniIndex::open(&ani_path)?;
        annotate::opencl::annotate_vcf_opencl_v2(
            &ani,
            &input_path,
            &out_for_annotate,
            &columns,
            bgzf_level,
            mmap_output,
            mmap_no_flush,
            ram_output,
            ram_max_mb,
            200_000,
        )?;
        ran = true;
    }

    if !ran {
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
        )?;
        ran = true;
    }

    if ran {
        if let Some(tmp) = bgzf_after {
            compress_plain_to_bgzf(&tmp, &out, bgzf_level)?;
            let _ = std::fs::remove_file(tmp);
        }
    }
    Ok(())
}

pub fn cmd_annotate_serve(args: AnnotateServeArgs) -> Result<()> {
    let use_gpu = args.gpu;
    let use_opencl = args.opencl || !args.gpu;

    if use_gpu && use_opencl && args.gpu && args.opencl {
        anyhow::bail!("Select only one backend: --gpu or --opencl");
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

    if args.ani.is_none() && !ani_path.exists() {
        eprintln!("[annotate] ANI index not found, building from source...");
        let Some(ann) = annotations else {
            anyhow::bail!("--annotations is required to build ANI index");
        };
        let ext = ann.extension().and_then(|e| e.to_str());
        if ext == Some("tab") {
            annotate::build_ani_index_from_tab(ann, &ani_path, args.columns.as_deref())?;
        } else {
            annotate::build_ani_index_auto_v2(ann, &ani_path)?;
        }
    } else if args.ani.is_some() && !ani_path.exists() {
        anyhow::bail!("ANI file not found: {:?}", ani_path);
    }

    let ani = annotate::AniIndex::open(&ani_path)?;

    let mut opencl_gpu = {
        #[cfg(feature = "opencl")]
        {
            if use_opencl {
                Some(annotate::opencl::OpenCLv2::new(&ani, 200_000)?)
            } else {
                None
            }
        }
        #[cfg(not(feature = "opencl"))]
        {
            if use_opencl {
                anyhow::bail!("OpenCL feature not enabled");
            }
            None
        }
    };

    let mut cuda_state = {
        #[cfg(feature = "gpu")]
        {
            if use_gpu {
                let start = Instant::now();
                let state = annotate::cuda::GpuAnnotator::new(&ani)?;
                eprintln!("[gpu] warmup: {:.3}s", start.elapsed().as_secs_f64());
                Some(state)
            } else {
                None
            }
        }
        #[cfg(not(feature = "gpu"))]
        {
            if use_gpu {
                anyhow::bail!("GPU feature not enabled");
            }
            None
        }
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

                    let _input_format = detect_format(&input_path)?;
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
                        sample_map: Arc::new(sample_map),
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
        } else {
            #[cfg(feature = "opencl")]
            {
                annotate::opencl::annotate_vcf_opencl_v2_with_gpu(
                    opencl_gpu.as_mut().unwrap(),
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
            #[cfg(not(feature = "opencl"))]
            {
                anyhow::bail!("OpenCL feature not enabled")
            }
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

    eprintln!("[db-build] Done");
    Ok(())
}
