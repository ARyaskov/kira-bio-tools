use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use std::collections::{BinaryHeap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use crate::cli::args::SortArgs;
use crate::annotate::postproc::{OutputKind, parse_output_type};

pub fn cmd_sort(args: SortArgs) -> Result<()> {
    let out_path = args.output.clone().unwrap_or_else(|| PathBuf::from("out.sort.vcf"));
    let kind = args.output_type.as_deref().map(parse_output_type).transpose()?.unwrap_or(OutputKind::Vcf);
    let max_mem_bytes = parse_max_mem(&args.max_mem)?;
    let temp_dir = args.temp_dir.clone().unwrap_or_else(std::env::temp_dir);

    let header = read_header(&args.input)?;
    let contig_order = build_contig_order(&header);
    let header_with_pass = ensure_pass_filter(header);

    let mut sink = SortSink::open(&out_path, kind, &header_with_pass)?;
    sink.write_headers(&header_with_pass)?;

    let (n_chunks, chunk_paths) = chunk_and_sort(&args.input, &temp_dir, max_mem_bytes, &contig_order)?;

    if n_chunks == 1 {
        let reader = BufReader::new(File::open(&chunk_paths[0])?);
        for line in reader.lines() {
            let line = line?;
            sink.write_line(&line)?;
        }
    } else {
        kway_merge_sink(&chunk_paths, &mut sink, &contig_order)?;
    }

    for p in &chunk_paths { let _ = std::fs::remove_file(p); }

    sink.finish()?;

    if let Some(kind) = &args.write_index {
        if matches!(parse_output_type(args.output_type.as_deref().unwrap_or("v"))?, OutputKind::VcfGz(_) | OutputKind::Bcf(_)) {
            let idx_path = std::path::PathBuf::from(format!("{}.{}", out_path.display(), if kind == "tbi" { "tbi" } else { "csi" }));
            eprintln!("[sort] -W {}: writing index to {:?}", kind, idx_path);
            let _ = crate::csi::build_csi_index(&out_path, &idx_path);
        } else {
            eprintln!("[sort] -W: index requires BGZF output; skipping");
        }
    }
    Ok(())
}

struct SortSink {
    inner: Box<dyn Write>,
    bcf: Option<crate::bcf::BcfWriter>,
}

impl SortSink {
    fn open(p: &Path, kind: OutputKind, headers: &[String]) -> Result<Self> {
        match kind {
            OutputKind::Vcf => Ok(Self { inner: Box::new(BufWriter::with_capacity(1 << 20, File::create(p)?)), bcf: None }),
            OutputKind::VcfGz(lvl) => Ok(Self { inner: Box::new(crate::bgzf::BgzfWriter::with_compression(p, flate2::Compression::new(lvl))?), bcf: None }),
            OutputKind::Bcf(lvl) => {
                let compressed = lvl > 0;
                let w = crate::bcf::BcfWriter::create(p, compressed, lvl, headers)?;
                Ok(Self { inner: Box::new(std::io::sink()), bcf: Some(w) })
            }
        }
    }
    fn write_headers(&mut self, headers: &[String]) -> Result<()> {
        if self.bcf.is_some() { return Ok(()); }
        for h in headers {
            self.inner.write_all(h.as_bytes())?;
            self.inner.write_all(b"\n")?;
        }
        Ok(())
    }
    fn write_line(&mut self, line: &str) -> Result<()> {
        if line.is_empty() { return Ok(()); }
        if let Some(bcf) = self.bcf.as_mut() {
            if !line.starts_with('#') { bcf.write_vcf_line(line)?; }
            return Ok(());
        }
        self.inner.write_all(line.as_bytes())?;
        self.inner.write_all(b"\n")?;
        Ok(())
    }
    fn finish(mut self) -> Result<()> {
        if let Some(bcf) = self.bcf.take() { bcf.finish()?; return Ok(()); }
        self.inner.flush()?;
        Ok(())
    }
}

fn kway_merge_sink(paths: &[PathBuf], sink: &mut SortSink, contig_order: &HashMap<String, usize>) -> Result<()> {
    let mut readers: Vec<std::io::Lines<BufReader<File>>> = paths.iter()
        .map(|p| Ok::<_, std::io::Error>(BufReader::new(File::open(p)?).lines())).collect::<Result<Vec<_>, _>>()?;
    let mut heap: BinaryHeap<HeapEntry> = BinaryHeap::new();
    for (i, r) in readers.iter_mut().enumerate() {
        if let Some(Ok(line)) = r.next() {
            let _ = contig_order;
            heap.push(HeapEntry { key: sort_key(&line, 0), line, src: i });
        }
    }
    while let Some(e) = heap.pop() {
        sink.write_line(&e.line)?;
        if let Some(Ok(line)) = readers[e.src].next() {
            heap.push(HeapEntry { key: sort_key(&line, 0), line, src: e.src });
        }
    }
    Ok(())
}

fn parse_max_mem(s: &str) -> Result<usize> {
    let s = s.trim();
    let (num, unit) = s.split_at(s.find(|c: char| !c.is_ascii_digit() && c != '.').unwrap_or(s.len()));
    let n: f64 = num.parse().with_context(|| format!("parse --max-mem {s:?}"))?;
    let mult: f64 = match unit.to_ascii_uppercase().as_str() {
        "" | "B" => 1.0, "K" | "KB" => 1024.0, "M" | "MB" => 1024.0 * 1024.0, "G" | "GB" => 1024.0 * 1024.0 * 1024.0,
        u => bail!("unknown --max-mem unit {u:?}"),
    };
    Ok((n * mult) as usize)
}

fn read_header(input: &Path) -> Result<Vec<String>> {
    let reader: Box<dyn Read> = if matches!(input.extension().and_then(|x| x.to_str()), Some("gz" | "bgz" | "bgzf")) {
        Box::new(MultiGzDecoder::new(File::open(input)?))
    } else { Box::new(File::open(input)?) };
    let mut header = Vec::new();
    for line in BufReader::new(reader).lines() {
        let line = line?;
        if line.starts_with('#') { header.push(line); } else { break; }
    }
    Ok(header)
}

fn build_contig_order(header: &[String]) -> HashMap<String, usize> {
    let mut m = HashMap::new();
    for h in header {
        if let Some(id) = parse_contig_id(h) {
            let n = m.len();
            m.entry(id).or_insert(n);
        }
    }
    m
}

fn ensure_pass_filter(mut header: Vec<String>) -> Vec<String> {
    if header.iter().any(|h| h.starts_with("##FILTER=<ID=PASS,")) { return header; }
    let insert_at = header.iter().position(|h| h.starts_with("##reference=")).unwrap_or(1);
    header.insert(insert_at, "##FILTER=<ID=PASS,Description=\"All filters passed\">".to_string());
    header
}

fn parse_contig_id(line: &str) -> Option<String> {
    if !line.starts_with("##contig=<") { return None; }
    let body = line.trim_start_matches("##contig=<").trim_end_matches('>');
    for kv in body.split(',') {
        if let Some(id) = kv.strip_prefix("ID=") { return Some(id.to_string()); }
    }
    None
}

#[derive(Clone)]
struct SortKey { chrom: String, pos: u32, refa: String, alt: String, idx: usize }

fn sort_key(line: &str, idx: usize) -> SortKey {
    let mut parts = line.split('\t');
    let chrom = parts.next().unwrap_or("").to_string();
    let pos = parts.next().and_then(|x| x.parse::<u32>().ok()).unwrap_or(u32::MAX);
    let _id = parts.next();
    let refa = parts.next().unwrap_or("").to_string();
    let alt = parts.next().unwrap_or("").to_string();
    SortKey { chrom, pos, refa, alt, idx }
}

fn cmp_keys(a: &SortKey, b: &SortKey, contig_order: &HashMap<String, usize>) -> std::cmp::Ordering {
    let ao = contig_order.get(&a.chrom).copied().unwrap_or(usize::MAX);
    let bo = contig_order.get(&b.chrom).copied().unwrap_or(usize::MAX);
    ao.cmp(&bo).then_with(|| a.chrom.cmp(&b.chrom)).then_with(|| a.pos.cmp(&b.pos))
        .then_with(|| a.refa.cmp(&b.refa)).then_with(|| a.alt.cmp(&b.alt))
        .then_with(|| a.idx.cmp(&b.idx))
}

fn chunk_and_sort(input: &Path, temp_dir: &Path, max_bytes: usize, contig_order: &HashMap<String, usize>) -> Result<(usize, Vec<PathBuf>)> {
    let reader: Box<dyn Read> = if matches!(input.extension().and_then(|x| x.to_str()), Some("gz" | "bgz" | "bgzf")) {
        Box::new(MultiGzDecoder::new(File::open(input)?))
    } else { Box::new(File::open(input)?) };
    let mut paths: Vec<PathBuf> = Vec::new();
    let mut buf: Vec<(SortKey, String)> = Vec::new();
    let mut buf_bytes: usize = 0;
    let mut idx = 0usize;
    for line in BufReader::new(reader).lines() {
        let line = line?;
        if line.starts_with('#') { continue; }
        buf_bytes += line.len() + 1;
        buf.push((sort_key(&line, idx), line));
        idx += 1;
        if buf_bytes >= max_bytes {
            flush_chunk(&mut buf, temp_dir, paths.len(), contig_order, &mut paths)?;
            buf.clear();
            buf_bytes = 0;
        }
    }
    if !buf.is_empty() { flush_chunk(&mut buf, temp_dir, paths.len(), contig_order, &mut paths)?; }
    let n = paths.len();
    Ok((n.max(1), if paths.is_empty() {
        let p = temp_dir.join(format!("kira-sort-empty-{}.vcf", std::process::id()));
        File::create(&p)?;
        vec![p]
    } else { paths }))
}

fn flush_chunk(buf: &mut Vec<(SortKey, String)>, temp_dir: &Path, idx: usize, contig_order: &HashMap<String, usize>, paths: &mut Vec<PathBuf>) -> Result<()> {
    buf.sort_by(|a, b| cmp_keys(&a.0, &b.0, contig_order));
    let p = temp_dir.join(format!("kira-sort-{}-{}.vcf", std::process::id(), idx));
    let mut w = BufWriter::with_capacity(1 << 20, File::create(&p)?);
    for (_, l) in buf.iter() { writeln!(w, "{}", l)?; }
    w.flush()?;
    paths.push(p);
    Ok(())
}

struct HeapEntry { key: SortKey, line: String, src: usize }
impl Eq for HeapEntry {}
impl PartialEq for HeapEntry { fn eq(&self, o: &Self) -> bool { self.cmp(o) == std::cmp::Ordering::Equal } }
impl PartialOrd for HeapEntry { fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> { Some(self.cmp(o)) } }
impl Ord for HeapEntry {
    fn cmp(&self, o: &Self) -> std::cmp::Ordering {
        o.key.chrom.cmp(&self.key.chrom).then_with(|| o.key.pos.cmp(&self.key.pos))
            .then_with(|| o.key.refa.cmp(&self.key.refa)).then_with(|| o.key.alt.cmp(&self.key.alt))
            .then_with(|| o.src.cmp(&self.src))
    }
}

