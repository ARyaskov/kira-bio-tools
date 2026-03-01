use anyhow::Result;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{Cursor, Read};

use crate::cli::args::HeadArgs;

pub fn cmd_head(args: HeadArgs) -> Result<()> {
    let cfg = parse_args(&args.bcftools_args)?;
    let content = read_input(&cfg)?;
    let lines: Vec<&str> = content.lines().collect();

    let mut header = Vec::new();
    let mut records = Vec::new();
    let mut chrom_line = None;

    for line in lines {
        if line.starts_with('#') {
            header.push(line);
            if line.starts_with("#CHROM") {
                chrom_line = Some(line);
            }
        } else {
            records.push(line);
        }
    }

    let mut out = Vec::new();
    let h_count = cfg.h.unwrap_or(usize::MAX);
    if cfg.h.is_some() {
        for line in header.iter().take(h_count) {
            out.push(*line);
        }
    } else if cfg.s.is_some() {
    } else {
        for line in &header {
            out.push(*line);
        }
    }

    let mut record_count = cfg.n.unwrap_or(0);
    if let Some(s) = cfg.s {
        record_count = s;
        if let Some(chrom) = chrom_line {
            if !out.iter().any(|l| *l == chrom) {
                out.push(chrom);
            }
        }
    }

    if record_count > 0 {
        for line in records.iter().take(record_count) {
            out.push(*line);
        }
    }

    if !out.is_empty() {
        println!("{}", out.join("\n"));
    }

    Ok(())
}

#[derive(Default)]
struct HeadConfig {
    input: Option<String>,
    h: Option<usize>,
    n: Option<usize>,
    s: Option<usize>,
}

fn parse_args(args: &[String]) -> Result<HeadConfig> {
    let mut cfg = HeadConfig::default();
    let mut i = 0usize;
    while i < args.len() {
        let a = &args[i];
        if a == "-h" {
            i += 1;
            let v = args
                .get(i)
                .ok_or_else(|| anyhow::anyhow!("missing value for -h"))?;
            cfg.h = Some(v.parse::<usize>()?);
        } else if let Some(v) = a.strip_prefix("-h") {
            if !v.is_empty() {
                cfg.h = Some(v.parse::<usize>()?);
            }
        } else if a == "-n" {
            i += 1;
            let v = args
                .get(i)
                .ok_or_else(|| anyhow::anyhow!("missing value for -n"))?;
            cfg.n = Some(v.parse::<usize>()?);
        } else if let Some(v) = a.strip_prefix("-n") {
            if !v.is_empty() {
                cfg.n = Some(v.parse::<usize>()?);
            }
        } else if a == "-s" {
            i += 1;
            let v = args
                .get(i)
                .ok_or_else(|| anyhow::anyhow!("missing value for -s"))?;
            cfg.s = Some(v.parse::<usize>()?);
        } else if let Some(v) = a.strip_prefix("-s") {
            if !v.is_empty() {
                cfg.s = Some(v.parse::<usize>()?);
            }
        } else if !a.starts_with('-') {
            cfg.input = Some(a.clone());
        }
        i += 1;
    }
    Ok(cfg)
}

fn read_input(cfg: &HeadConfig) -> Result<String> {
    match cfg.input.as_deref() {
        Some(path) if path != "-" => {
            let mut bytes = Vec::new();
            File::open(path)?.read_to_end(&mut bytes)?;
            decode_bytes(bytes, path)
        }
        _ => {
            let mut bytes = Vec::new();
            std::io::stdin().read_to_end(&mut bytes)?;
            decode_bytes(bytes, "")
        }
    }
}

fn decode_bytes(bytes: Vec<u8>, path_hint: &str) -> Result<String> {
    let is_gzip_magic = bytes.len() >= 2 && bytes[0] == 0x1f && bytes[1] == 0x8b;
    let is_gzip_ext = matches!(
        std::path::Path::new(path_hint)
            .extension()
            .and_then(|x| x.to_str()),
        Some("gz" | "bgz" | "bgzf")
    );

    if is_gzip_magic || is_gzip_ext {
        let mut decoded = String::new();
        MultiGzDecoder::new(Cursor::new(bytes)).read_to_string(&mut decoded)?;
        return Ok(decoded);
    }

    String::from_utf8(bytes)
        .map_err(|_| anyhow::anyhow!("non-UTF8 input is not supported by native head yet"))
}
