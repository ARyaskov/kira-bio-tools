use crate::annotate::postproc::{OutputKind, parse_output_type, version_header_line};
use crate::cli::args::ConcatArgs;
use crate::vcf::UnifiedVcfReader;
use anyhow::{Context, Result, bail};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

pub fn cmd_concat(args: ConcatArgs) -> Result<()> {
    let mut inputs: Vec<PathBuf> = args.inputs.clone();
    if let Some(fl) = &args.file_list {
        for line in BufReader::new(File::open(fl)?).lines() {
            let l = line?;
            let t = l.trim();
            if t.is_empty() || t.starts_with('#') { continue; }
            inputs.push(PathBuf::from(t));
        }
    }
    if inputs.is_empty() { bail!("concat: no inputs"); }

    let kind = args.output_type.as_deref().map(parse_output_type).transpose()?.unwrap_or(OutputKind::Vcf);
    let mut sink = open_sink(args.output.as_deref(), kind)?;

    if args.naive || args.naive_force {
        return concat_naive(&inputs, &mut sink);
    }
    if args.ligate || args.ligate_force || args.ligate_warn {
        return concat_ligate(&inputs, args, &mut sink);
    }
    concat_streaming(&inputs, args, &mut sink)
}

fn concat_ligate(inputs: &[PathBuf], args: ConcatArgs, sink: &mut Box<dyn Write>) -> Result<()> {
    let mut first = UnifiedVcfReader::open(&inputs[0]).context("open first input")?;
    let headers = first.header()?;
    let version = version_header_line();
    let mut wrote_chrom = false;
    for h in &headers {
        if h.starts_with("#CHROM") {
            if !args.no_version {
                sink.write_all(version.as_bytes())?; sink.write_all(b"\n")?;
            }
            wrote_chrom = true;
        }
        sink.write_all(h.as_bytes())?; sink.write_all(b"\n")?;
    }
    if !wrote_chrom && !args.no_version {
        sink.write_all(version.as_bytes())?; sink.write_all(b"\n")?;
    }

    let mut last_lines: std::collections::HashMap<String, (String, i32)> = std::collections::HashMap::new();
    let mut phase_swap: std::collections::HashMap<String, bool> = std::collections::HashMap::new();
    let mut samples: Vec<String> = Vec::new();
    for h in &headers {
        if h.starts_with("#CHROM") {
            let cols: Vec<&str> = h.split('\t').collect();
            if cols.len() > 9 { samples = cols[9..].iter().map(|s| s.to_string()).collect(); }
            break;
        }
    }

    let min_pq = args.min_pq as i32;
    let n_inputs = inputs.len();
    let readers: Vec<Box<dyn Iterator<Item = Result<String>>>> = Vec::new();
    let _ = readers;

    // Iterate file-by-file maintaining boundary phase state
    let mut process = |r: &mut UnifiedVcfReader, is_first: bool, sink: &mut Box<dyn Write>| -> Result<()> {
        let mut local_swap: std::collections::HashMap<String, bool> = phase_swap.clone();
        while let Some(line) = r.read_line()? {
            if line.is_empty() || line.as_bytes()[0] == b'#' { continue; }
            let key = dedup_key(&line);
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 8 { sink.write_all(line.as_bytes())?; sink.write_all(b"\n")?; continue; }
            let chrom_pos = format!("{}:{}", cols[0], cols[1]);
            let _ = chrom_pos;
            if let Some((_prev, _prev_pq)) = last_lines.get(&key) {
                if !is_first {
                    if let Some((new_line, swap_map)) = ligate_phase(&line, &samples, min_pq) {
                        for (s, sw) in swap_map.iter() {
                            *local_swap.entry(s.clone()).or_insert(false) ^= *sw;
                        }
                        let _ = new_line;
                    }
                    continue;
                }
            }
            let final_line = apply_phase_swap(&line, &samples, &local_swap);
            sink.write_all(final_line.as_bytes())?; sink.write_all(b"\n")?;
            last_lines.insert(key, (line, 0));
            if last_lines.len() > 4096 {
                let drop_keys: Vec<String> = last_lines.keys().take(2048).cloned().collect();
                for k in drop_keys { last_lines.remove(&k); }
            }
        }
        phase_swap = local_swap;
        Ok(())
    };
    process(&mut first, true, sink)?;
    for (i, p) in inputs[1..].iter().enumerate() {
        let mut r = UnifiedVcfReader::open(p).with_context(|| format!("open {:?}", p))?;
        let _ = r.header()?;
        process(&mut r, i + 1 == 0, sink)?;
    }
    let _ = n_inputs;
    sink.flush()?;
    Ok(())
}

fn ligate_phase(line: &str, _samples: &[String], min_pq: i32) -> Option<(String, std::collections::HashMap<String, bool>)> {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 10 { return None; }
    let fmt: Vec<&str> = cols[8].split(':').collect();
    let gt_idx = fmt.iter().position(|k| *k == "GT")?;
    let pq_idx = fmt.iter().position(|k| *k == "PQ");
    let mut swap_map: std::collections::HashMap<String, bool> = std::collections::HashMap::new();
    for (i, samp) in cols[9..].iter().enumerate() {
        let parts: Vec<&str> = samp.split(':').collect();
        if let (Some(gt), Some(_)) = (parts.get(gt_idx), Some(())) {
            if let Some(pq_i) = pq_idx {
                let pq: i32 = parts.get(pq_i).and_then(|s| s.parse().ok()).unwrap_or(0);
                if pq < min_pq { continue; }
            }
            let parts2: Vec<&str> = gt.split('|').collect();
            if parts2.len() == 2 && parts2[0] != parts2[1] {
                swap_map.insert(format!("sample_{}", i), parts2[0] > parts2[1]);
            }
        }
    }
    Some((line.to_string(), swap_map))
}

fn apply_phase_swap(line: &str, _samples: &[String], swap_map: &std::collections::HashMap<String, bool>) -> String {
    if swap_map.is_empty() { return line.to_string(); }
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 10 { return line.to_string(); }
    let fmt: Vec<&str> = cols[8].split(':').collect();
    let Some(gt_idx) = fmt.iter().position(|k| *k == "GT") else { return line.to_string(); };
    let mut new_cols: Vec<String> = cols[..9].iter().map(|s| s.to_string()).collect();
    for (i, samp) in cols[9..].iter().enumerate() {
        let key = format!("sample_{}", i);
        let swap = swap_map.get(&key).copied().unwrap_or(false);
        if !swap { new_cols.push(samp.to_string()); continue; }
        let parts: Vec<&str> = samp.split(':').collect();
        let mut nps: Vec<String> = parts.iter().map(|s| s.to_string()).collect();
        if let Some(gt) = parts.get(gt_idx) {
            if gt.contains('|') {
                let p: Vec<&str> = gt.split('|').collect();
                if p.len() == 2 { nps[gt_idx] = format!("{}|{}", p[1], p[0]); }
            }
        }
        new_cols.push(nps.join(":"));
    }
    new_cols.join("\t")
}

fn concat_streaming(inputs: &[PathBuf], args: ConcatArgs, sink: &mut Box<dyn Write>) -> Result<()> {
    let mut first = UnifiedVcfReader::open(&inputs[0]).context("open first input")?;
    let headers = first.header()?;

    let out_headers = if args.drop_genotypes { drop_format_headers(headers) } else { headers };
    let version = version_header_line();
    let mut wrote_chrom = false;
    for h in &out_headers {
        if h.starts_with("#CHROM") {
            if !args.no_version {
                sink.write_all(version.as_bytes())?; sink.write_all(b"\n")?;
            }
            wrote_chrom = true;
        }
        let line = if args.drop_genotypes && h.starts_with("#CHROM") { trim_chrom_samples(h) } else { h.clone() };
        sink.write_all(line.as_bytes())?; sink.write_all(b"\n")?;
    }
    if !wrote_chrom && !args.no_version {
        sink.write_all(version.as_bytes())?; sink.write_all(b"\n")?;
    }

    let mut seen: HashSet<String> = HashSet::new();
    write_records(&mut first, &args, sink, &mut seen)?;
    for p in &inputs[1..] {
        let mut r = UnifiedVcfReader::open(p).with_context(|| format!("open {:?}", p))?;
        let _ = r.header()?;
        write_records(&mut r, &args, sink, &mut seen)?;
    }
    sink.flush()?;
    Ok(())
}

fn write_records(
    r: &mut UnifiedVcfReader,
    args: &ConcatArgs,
    sink: &mut Box<dyn Write>,
    seen: &mut HashSet<String>,
) -> Result<()> {
    while let Some(line) = r.read_line()? {
        if line.is_empty() || line.as_bytes()[0] == b'#' { continue; }
        if args.remove_duplicates || args.rm_dups.is_some() {
            let key = dedup_key(&line);
            if !seen.insert(key) { continue; }
        }
        let out_line = if args.drop_genotypes { drop_genotypes_line(&line) } else { line };
        sink.write_all(out_line.as_bytes())?; sink.write_all(b"\n")?;
    }
    Ok(())
}

fn concat_naive(inputs: &[PathBuf], sink: &mut Box<dyn Write>) -> Result<()> {
    for (i, p) in inputs.iter().enumerate() {
        let mut r = UnifiedVcfReader::open(p)?;
        let h = r.header()?;
        if i == 0 {
            for line in &h { sink.write_all(line.as_bytes())?; sink.write_all(b"\n")?; }
        }
        while let Some(line) = r.read_line()? {
            if line.starts_with('#') { continue; }
            sink.write_all(line.as_bytes())?; sink.write_all(b"\n")?;
        }
    }
    sink.flush()?;
    Ok(())
}

fn open_sink(path: Option<&Path>, kind: OutputKind) -> Result<Box<dyn Write>> {
    match (path, kind) {
        (None, OutputKind::Vcf) => Ok(Box::new(BufWriter::with_capacity(1 << 20, std::io::stdout()))),
        (None, OutputKind::VcfGz(_)) => bail!("-O z requires -o FILE"),
        (Some(p), OutputKind::Vcf) => Ok(Box::new(BufWriter::with_capacity(1 << 20, File::create(p)?))),
        (Some(p), OutputKind::VcfGz(lvl)) => {
            let w = crate::bgzf::BgzfWriter::with_compression(p, flate2::Compression::new(lvl))?;
            Ok(Box::new(w))
        }
        (_, OutputKind::Bcf(_)) => bail!("-O u|b (BCF) not yet supported in concat"),
    }
}

fn dedup_key(line: &str) -> String {
    let cols: Vec<&str> = line.splitn(6, '\t').collect();
    if cols.len() < 5 { return line.to_string(); }
    format!("{}\t{}\t{}\t{}", cols[0], cols[1], cols[3], cols[4])
}

fn drop_format_headers(h: Vec<String>) -> Vec<String> {
    h.into_iter().filter(|l| !l.starts_with("##FORMAT=")).collect()
}

fn trim_chrom_samples(h: &str) -> String {
    let cols: Vec<&str> = h.split('\t').collect();
    cols.iter().take(8).copied().collect::<Vec<_>>().join("\t")
}

fn drop_genotypes_line(line: &str) -> String {
    let cols: Vec<&str> = line.split('\t').collect();
    cols.iter().take(8).copied().collect::<Vec<_>>().join("\t")
}

#[allow(dead_code)]
fn naive_bgzf_concat(inputs: &[PathBuf], output: &Path) -> Result<()> {
    let mut out = BufWriter::new(File::create(output)?);
    for p in inputs {
        let mut f = File::open(p)?;
        let mut buf = vec![0u8; 1 << 20];
        loop {
            let n = f.read(&mut buf)?;
            if n == 0 { break; }
            out.write_all(&buf[..n])?;
        }
    }
    out.flush()?;
    Ok(())
}
