use anyhow::{Context, Result, bail};
use std::collections::{BinaryHeap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::cli::args::SortArgs;
use crate::csi::{IndexKind, build_index};
use crate::vcf::sink::{OutputKind, parse_output_type};
use crate::vcf::{UnifiedVcfReader, VcfSink};

pub fn cmd_sort(args: SortArgs) -> Result<()> {
    let mut cfg = SortCfg {
        output: args.output.clone(),
        output_type: args.output_type.clone(),
        max_mem: args.max_mem.clone(),
        temp_dir: args.temp_dir.clone(),
        write_index: args.write_index.clone(),
    };
    apply_passthrough(&mut cfg, &args.passthrough)?;

    let kind = cfg.output_type.as_deref().map(parse_output_type).transpose()?.unwrap_or(OutputKind::Vcf);
    let max_mem_bytes = parse_max_mem(&cfg.max_mem)?;
    let temp_dir = cfg.temp_dir.clone().unwrap_or_else(std::env::temp_dir);

    let mut reader = UnifiedVcfReader::open(&args.input).context("open input")?;
    let header = reader.header()?;
    let header_with_pass = ensure_pass_filter(header);

    // Contig ranks: header order first, then contigs in order of appearance.
    let mut contig_order: HashMap<String, usize> = build_contig_order(&header_with_pass);

    let mut sink = VcfSink::open(cfg.output.as_deref(), kind, &header_with_pass)?;
    sink.write_header(&header_with_pass)?;

    let chunk_paths = chunk_and_sort(&mut reader, &temp_dir, max_mem_bytes, &mut contig_order)?;

    if chunk_paths.len() <= 1 {
        if let Some(p) = chunk_paths.first() {
            let reader = BufReader::new(File::open(p)?);
            for line in reader.lines() {
                sink.write_line(&line?)?;
            }
        }
    } else {
        kway_merge_sink(&chunk_paths, &mut sink, &contig_order)?;
    }

    for p in &chunk_paths {
        let _ = std::fs::remove_file(p);
    }

    sink.finish()?;

    if let (Some(kind_s), Some(out)) = (cfg.write_index.as_deref(), cfg.output.as_deref()) {
        if matches!(kind, OutputKind::VcfGz(_) | OutputKind::Bcf(_)) && out != Path::new("-") {
            let (ik, ext) = if kind_s == "tbi" { (IndexKind::Tbi, "tbi") } else { (IndexKind::Csi, "csi") };
            let idx_path = PathBuf::from(format!("{}.{}", out.display(), ext));
            build_index(out, &idx_path, ik, None).with_context(|| format!("-W: write {}", idx_path.display()))?;
        } else {
            eprintln!("[sort] -W: index requires BGZF/BCF output to a file; skipping");
        }
    }
    Ok(())
}

struct SortCfg {
    output: Option<PathBuf>,
    output_type: Option<String>,
    max_mem: String,
    temp_dir: Option<PathBuf>,
    write_index: Option<String>,
}

fn apply_passthrough(cfg: &mut SortCfg, args: &[String]) -> Result<()> {
    let mut i = 0usize;
    while i < args.len() {
        let a = args[i].as_str();
        let mut take = |cfg_field: &mut dyn FnMut(&str) -> Result<()>| -> Result<()> {
            i += 1;
            let v = args.get(i).ok_or_else(|| anyhow::anyhow!("missing value for {a}"))?;
            cfg_field(v)
        };
        match a {
            "-o" | "--output" => take(&mut |v| { cfg.output = Some(PathBuf::from(v)); Ok(()) })?,
            "-O" | "--output-type" => take(&mut |v| { cfg.output_type = Some(v.to_string()); Ok(()) })?,
            "-m" | "--max-mem" => take(&mut |v| { cfg.max_mem = v.to_string(); Ok(()) })?,
            "-T" | "--temp-dir" => take(&mut |v| { cfg.temp_dir = Some(PathBuf::from(v)); Ok(()) })?,
            "-W" | "--write-index" => cfg.write_index = Some("csi".into()),
            "--no-version" | "--verbosity" | "-v" => {}
            other => {
                if let Some(v) = other.strip_prefix("-m") {
                    cfg.max_mem = v.to_string();
                } else if let Some(v) = other.strip_prefix("-O") {
                    cfg.output_type = Some(v.to_string());
                } else if let Some(v) = other.strip_prefix("-W=") {
                    cfg.write_index = Some(v.to_string());
                } else {
                    bail!("sort: unknown option {other:?}");
                }
            }
        }
        i += 1;
    }
    Ok(())
}

fn kway_merge_sink(paths: &[PathBuf], sink: &mut VcfSink, contig_order: &HashMap<String, usize>) -> Result<()> {
    let mut readers: Vec<std::io::Lines<BufReader<File>>> = paths
        .iter()
        .map(|p| Ok::<_, std::io::Error>(BufReader::new(File::open(p)?).lines()))
        .collect::<Result<Vec<_>, _>>()?;
    let mut heap: BinaryHeap<HeapEntry> = BinaryHeap::new();
    for (i, r) in readers.iter_mut().enumerate() {
        if let Some(Ok(line)) = r.next() {
            heap.push(HeapEntry { key: sort_key(&line, 0, contig_order), line, src: i });
        }
    }
    while let Some(e) = heap.pop() {
        sink.write_line(&e.line)?;
        if let Some(Ok(line)) = readers[e.src].next() {
            heap.push(HeapEntry { key: sort_key(&line, 0, contig_order), line, src: e.src });
        }
    }
    Ok(())
}

fn parse_max_mem(s: &str) -> Result<usize> {
    let s = s.trim();
    let (num, unit) = s.split_at(s.find(|c: char| !c.is_ascii_digit() && c != '.').unwrap_or(s.len()));
    let n: f64 = num.parse().with_context(|| format!("parse --max-mem {s:?}"))?;
    let mult: f64 = match unit.to_ascii_uppercase().as_str() {
        "" | "B" => 1.0,
        "K" | "KB" => 1024.0,
        "M" | "MB" => 1024.0 * 1024.0,
        "G" | "GB" => 1024.0 * 1024.0 * 1024.0,
        u => bail!("unknown --max-mem unit {u:?}"),
    };
    Ok((n * mult) as usize)
}

fn build_contig_order(header: &[String]) -> HashMap<String, usize> {
    let mut m = HashMap::new();
    for h in header {
        if let Some((id, _)) = crate::vcf::header::parse_contig_line(h) {
            let n = m.len();
            m.entry(id).or_insert(n);
        }
    }
    m
}

fn ensure_pass_filter(mut header: Vec<String>) -> Vec<String> {
    if header.iter().any(|h| h.starts_with("##FILTER=<ID=PASS,")) {
        return header;
    }
    let insert_at = header.iter().position(|h| h.starts_with("##reference=")).unwrap_or(1).min(header.len());
    header.insert(insert_at, "##FILTER=<ID=PASS,Description=\"All filters passed\">".to_string());
    header
}

/// Sort key: contig rank (header order, then first appearance), position,
/// then REF/ALT and input order for stability.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
struct SortKey {
    rank: usize,
    pos: u32,
    refa: String,
    alt: String,
    idx: usize,
}

fn sort_key(line: &str, idx: usize, contig_order: &HashMap<String, usize>) -> SortKey {
    let mut parts = line.split('\t');
    let chrom = parts.next().unwrap_or("");
    let pos = parts.next().and_then(|x| x.parse::<u32>().ok()).unwrap_or(u32::MAX);
    let _id = parts.next();
    let refa = parts.next().unwrap_or("").to_string();
    let alt = parts.next().unwrap_or("").to_string();
    let rank = contig_order.get(chrom).copied().unwrap_or(usize::MAX);
    SortKey { rank, pos, refa, alt, idx }
}

fn chunk_and_sort(
    reader: &mut UnifiedVcfReader,
    temp_dir: &Path,
    max_bytes: usize,
    contig_order: &mut HashMap<String, usize>,
) -> Result<Vec<PathBuf>> {
    let mut paths: Vec<PathBuf> = Vec::new();
    let mut buf: Vec<(SortKey, String)> = Vec::new();
    let mut buf_bytes: usize = 0;
    let mut idx = 0usize;
    while let Some(line) = reader.read_line()? {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let chrom = line.split('\t').next().unwrap_or("");
        if !contig_order.contains_key(chrom) {
            let n = contig_order.len();
            contig_order.insert(chrom.to_string(), n);
        }
        buf_bytes += line.len() + 1;
        buf.push((sort_key(&line, idx, contig_order), line));
        idx += 1;
        if buf_bytes >= max_bytes {
            flush_chunk(&mut buf, temp_dir, paths.len(), &mut paths)?;
            buf.clear();
            buf_bytes = 0;
        }
    }
    if !buf.is_empty() {
        flush_chunk(&mut buf, temp_dir, paths.len(), &mut paths)?;
    }
    Ok(paths)
}

fn flush_chunk(buf: &mut [(SortKey, String)], temp_dir: &Path, idx: usize, paths: &mut Vec<PathBuf>) -> Result<()> {
    buf.sort_by(|a, b| a.0.cmp(&b.0));
    let p = temp_dir.join(format!("kira-sort-{}-{}.vcf", std::process::id(), idx));
    let mut w = BufWriter::with_capacity(1 << 20, File::create(&p).with_context(|| format!("create temp {}", p.display()))?);
    for (_, l) in buf.iter() {
        writeln!(w, "{}", l)?;
    }
    w.flush()?;
    paths.push(p);
    Ok(())
}

struct HeapEntry {
    key: SortKey,
    line: String,
    src: usize,
}
impl Eq for HeapEntry {}
impl PartialEq for HeapEntry {
    fn eq(&self, o: &Self) -> bool {
        self.cmp(o) == std::cmp::Ordering::Equal
    }
}
impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for HeapEntry {
    // Reversed so the BinaryHeap pops the smallest key first.
    fn cmp(&self, o: &Self) -> std::cmp::Ordering {
        o.key.cmp(&self.key).then_with(|| o.src.cmp(&self.src))
    }
}

#[cfg(test)]
#[path = "../../../tests/unit/cli_commands_sort.rs"]
mod tests;
