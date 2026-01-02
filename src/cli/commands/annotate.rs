use crate::annotate;
use crate::cli::args::{AnnotateArgs, AnnotateIndexArgs, DbBuildArgs};
use anyhow::Result;

pub fn cmd_annotate(args: AnnotateArgs) -> Result<()> {
    let out = args.output.clone().unwrap_or_else(|| {
        let mut p = args.input.clone();
        p.set_extension("annot.vcf");
        p
    });

    let ani_path = if args.annotations.extension().unwrap_or_default() == "ani" {
        args.annotations.clone()
    } else {
        let mut p = args.annotations.clone();
        p.set_extension("ani");
        p
    };

    if !ani_path.exists() {
        eprintln!("[annotate] ANI index not found, building from source...");

        let ext = args.annotations.extension().and_then(|e| e.to_str());

        if ext == Some("tab") {
            annotate::build_ani_index_from_tab(
                &args.annotations,
                &ani_path,
                args.columns.as_deref(),
            )?;
        } else {
            annotate::build_ani_index_auto_v2(&args.annotations, &ani_path)?;
        }
    }

    eprintln!("[annotate] ANI = {:?}", ani_path);
    eprintln!("[annotate] Input = {:?}", args.input);
    eprintln!("[annotate] Output = {:?}", out);

    #[cfg(feature = "gpu")]
    if args.gpu {
        eprintln!("[annotate] Using CUDA GPU backend…");
        let ani = annotate::AniIndex::open(&ani_path)?;
        let gpu = annotate::cuda::GpuAni::load(&ani)?;
        annotate::cuda::annotate_vcf_ani_gpu(&gpu, &ani, &args.input, &out)?;
        return Ok(());
    }

    let columns: Vec<String> = if let Some(cols_str) = &args.columns {
        cols_str.split(',').map(|s| s.to_string()).collect()
    } else {
        Vec::new()
    };

    #[cfg(feature = "opencl")]
    if args.opencl {
        eprintln!("[annotate] Using OpenCL backend…");
        let ani = annotate::AniIndex::open(&ani_path)?;
        let gpu = annotate::opencl::OpenCLv2::new(&ani, 200_000)?;
        annotate::opencl::annotate_vcf_opencl_v2(&gpu, &ani, &args.input, &out, &columns)?;
        return Ok(());
    }

    annotate::cpu_v2::annotate_vcf_ani_v2(&ani_path, &args.input, &out, &columns)?;
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
