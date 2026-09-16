//! `sort`: external merge sort by (contig rank, POS), ties broken like
//! bcftools' `cmp_bcf_pos` (alleles compared case-insensitively, one by
//! one), then input order. Records are kept as text lines or, for a BCF
//! input, as their encoded blocks, so nothing is decoded or re-encoded
//! unless the output type differs; chunks that do not fit in `-m` go to
//! BGZF-compressed temporary files and come back through a k-way merge.

use anyhow::{Context, Result, bail};
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap};
use std::fs::File;
use std::io::{BufRead, Read};
use std::path::{Path, PathBuf};

use crate::bcf::record::{decode_record_to_vcf, record_meta, shared_alleles};
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

    let kind = cfg
        .output_type
        .as_deref()
        .map(parse_output_type)
        .transpose()?
        .unwrap_or(OutputKind::Vcf);
    let max_mem_bytes = parse_max_mem(&cfg.max_mem)?;
    let temp_dir = cfg.temp_dir.clone().unwrap_or_else(std::env::temp_dir);

    let mut reader = UnifiedVcfReader::open(&args.input).context("open input")?;
    let header = reader.header()?;
    let header_with_pass = ensure_pass_filter(header);

    // Contig ranks: header order first, then contigs in order of appearance
    // (a BCF's rid already is the header order).
    let mut contig_order: HashMap<String, usize> = build_contig_order(&header_with_pass);

    let mut sink = VcfSink::open(cfg.output.as_deref(), kind, &header_with_pass)?;
    sink.write_header(&header_with_pass)?;

    let bcf_in = reader.is_bcf();
    let mut sorter = Sorter::new(bcf_in, temp_dir, max_mem_bytes);
    if bcf_in {
        while let Some((shared, indiv)) = reader.read_raw_bcf()? {
            let Some(m) = record_meta(&shared) else {
                bail!("corrupt BCF record")
            };
            let key = pack_key(m.rid.max(0) as u64, m.pos.max(0) as u32 + 1);
            sorter.push(key, Data::bcf(shared, indiv))?;
        }
    } else {
        while let Some(line) = reader.read_line()? {
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let chrom = line.split('\t').next().unwrap_or("");
            if !contig_order.contains_key(chrom) {
                let n = contig_order.len();
                contig_order.insert(chrom.to_string(), n);
            }
            let key = sort_key(&line, &contig_order);
            sorter.push(key, Data::Text(line))?;
        }
    }

    // Raw BCF records go straight to a BCF sink; anything else is decoded.
    let dict = reader.bcf_dict().cloned();
    let raw_out = bcf_in && sink.is_bcf();
    let mut emit = |data: &Data| -> Result<()> {
        match data {
            Data::Text(line) => sink.write_line(line),
            Data::Bcf(framed) => {
                if raw_out {
                    let (shared, indiv) = split_framed(framed)?;
                    sink.write_raw_bcf(shared, indiv)
                } else {
                    let dict = dict.as_ref().context("BCF dictionary")?;
                    match decode_record_to_vcf(&mut std::io::Cursor::new(framed.as_slice()), dict)?
                    {
                        Some(line) => sink.write_line(&line),
                        None => Ok(()),
                    }
                }
            }
        }
    };
    let key_of = |d: &Data| -> u64 {
        match d {
            Data::Text(l) => sort_key(l, &contig_order),
            Data::Bcf(_) => bcf_key(d),
        }
    };
    sorter.finish(&key_of, &mut emit)?;
    sink.finish()?;

    if let (Some(kind_s), Some(out)) = (cfg.write_index.as_deref(), cfg.output.as_deref()) {
        if matches!(kind, OutputKind::VcfGz(_) | OutputKind::Bcf(_)) && out != Path::new("-") {
            let (ik, ext) = if kind_s == "tbi" {
                (IndexKind::Tbi, "tbi")
            } else {
                (IndexKind::Csi, "csi")
            };
            let idx_path = PathBuf::from(format!("{}.{}", out.display(), ext));
            build_index(out, &idx_path, ik, None)
                .with_context(|| format!("-W: write {}", idx_path.display()))?;
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
            let v = args
                .get(i)
                .ok_or_else(|| anyhow::anyhow!("missing value for {a}"))?;
            cfg_field(v)
        };
        match a {
            "-o" | "--output" => take(&mut |v| {
                cfg.output = Some(PathBuf::from(v));
                Ok(())
            })?,
            "-O" | "--output-type" => take(&mut |v| {
                cfg.output_type = Some(v.to_string());
                Ok(())
            })?,
            "-m" | "--max-mem" => take(&mut |v| {
                cfg.max_mem = v.to_string();
                Ok(())
            })?,
            "-T" | "--temp-dir" => take(&mut |v| {
                cfg.temp_dir = Some(PathBuf::from(v));
                Ok(())
            })?,
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

fn parse_max_mem(s: &str) -> Result<usize> {
    let s = s.trim();
    let (num, unit) = s.split_at(
        s.find(|c: char| !c.is_ascii_digit() && c != '.')
            .unwrap_or(s.len()),
    );
    let n: f64 = num
        .parse()
        .with_context(|| format!("parse --max-mem {s:?}"))?;
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
    let insert_at = header
        .iter()
        .position(|h| h.starts_with("##reference="))
        .unwrap_or(1)
        .min(header.len());
    header.insert(
        insert_at,
        "##FILTER=<ID=PASS,Description=\"All filters passed\">".to_string(),
    );
    header
}

/// `rank << 32 | pos`: one integer compare orders records by contig and
/// position; unknown contigs sort last.
fn pack_key(rank: u64, pos: u32) -> u64 {
    (rank.min(u32::MAX as u64) << 32) | pos as u64
}

/// The packed key of a text record.
fn sort_key(line: &str, contig_order: &HashMap<String, usize>) -> u64 {
    let mut parts = line.split('\t');
    let chrom = parts.next().unwrap_or("");
    let pos = parts
        .next()
        .and_then(|x| x.parse::<u32>().ok())
        .unwrap_or(u32::MAX);
    let rank = contig_order
        .get(chrom)
        .map(|&r| r as u64)
        .unwrap_or(u32::MAX as u64);
    pack_key(rank, pos)
}

/// One record: its text line, or an encoded BCF record framed as
/// `l_shared, l_indiv, shared, indiv` (the file layout).
enum Data {
    Text(String),
    Bcf(Vec<u8>),
}

impl Data {
    fn bcf(shared: Vec<u8>, indiv: Vec<u8>) -> Self {
        let mut framed = Vec::with_capacity(8 + shared.len() + indiv.len());
        framed.extend_from_slice(&(shared.len() as u32).to_le_bytes());
        framed.extend_from_slice(&(indiv.len() as u32).to_le_bytes());
        framed.extend_from_slice(&shared);
        framed.extend_from_slice(&indiv);
        Self::Bcf(framed)
    }

    fn bytes(&self) -> usize {
        match self {
            Self::Text(l) => l.len() + 1,
            Self::Bcf(b) => b.len(),
        }
    }

    /// bcftools' `cmp_bcf_pos` tie-break: alleles compared case-insensitively
    /// one by one, a shorter allele list first.
    fn cmp_alleles(&self, other: &Self) -> Ordering {
        match (self, other) {
            (Self::Text(a), Self::Text(b)) => cmp_allele_lists(text_alleles(a), text_alleles(b)),
            (Self::Bcf(a), Self::Bcf(b)) => {
                let (sa, sb) = (split_framed(a).map(|x| x.0), split_framed(b).map(|x| x.0));
                match (sa, sb) {
                    (Ok(sa), Ok(sb)) => cmp_allele_lists(
                        shared_alleles(sa).unwrap_or_default().into_iter(),
                        shared_alleles(sb).unwrap_or_default().into_iter(),
                    ),
                    _ => Ordering::Equal,
                }
            }
            _ => Ordering::Equal,
        }
    }
}

fn text_alleles(line: &str) -> impl Iterator<Item = &[u8]> {
    let mut cols = line.split('\t');
    let refa = cols.nth(3).unwrap_or("");
    let alt = cols.next().unwrap_or("");
    std::iter::once(refa.as_bytes()).chain(alt.split(',').map(str::as_bytes))
}

fn cmp_allele_lists<'a, A, B>(mut a: A, mut b: B) -> Ordering
where
    A: Iterator<Item = &'a [u8]>,
    B: Iterator<Item = &'a [u8]>,
{
    loop {
        match (a.next(), b.next()) {
            (None, None) => return Ordering::Equal,
            (None, Some(_)) => return Ordering::Less,
            (Some(_), None) => return Ordering::Greater,
            (Some(x), Some(y)) => {
                let o = x
                    .iter()
                    .map(u8::to_ascii_lowercase)
                    .cmp(y.iter().map(u8::to_ascii_lowercase));
                if o != Ordering::Equal {
                    return o;
                }
            }
        }
    }
}

fn split_framed(framed: &[u8]) -> Result<(&[u8], &[u8])> {
    if framed.len() < 8 {
        bail!("corrupt buffered BCF record");
    }
    let ls = u32::from_le_bytes(framed[0..4].try_into().unwrap()) as usize;
    let li = u32::from_le_bytes(framed[4..8].try_into().unwrap()) as usize;
    if framed.len() != 8 + ls + li {
        bail!("corrupt buffered BCF record");
    }
    Ok((&framed[8..8 + ls], &framed[8 + ls..]))
}

struct Item {
    key: u64,
    idx: u32,
    data: Data,
}

fn cmp_items(a: &Item, b: &Item) -> Ordering {
    a.key
        .cmp(&b.key)
        .then_with(|| a.data.cmp_alleles(&b.data))
        .then(a.idx.cmp(&b.idx))
}

/// Buffers records up to `max_bytes`, spilling sorted chunks to BGZF
/// temporary files, and replays everything in order.
struct Sorter {
    bcf: bool,
    temp_dir: PathBuf,
    max_bytes: usize,
    buf: Vec<Item>,
    buf_bytes: usize,
    idx: u32,
    chunks: Vec<PathBuf>,
}

impl Sorter {
    fn new(bcf: bool, temp_dir: PathBuf, max_bytes: usize) -> Self {
        Self {
            bcf,
            temp_dir,
            max_bytes,
            buf: Vec::new(),
            buf_bytes: 0,
            idx: 0,
            chunks: Vec::new(),
        }
    }

    fn push(&mut self, key: u64, data: Data) -> Result<()> {
        self.buf_bytes += data.bytes() + std::mem::size_of::<Item>();
        self.buf.push(Item {
            key,
            idx: self.idx,
            data,
        });
        self.idx = self.idx.wrapping_add(1);
        if self.buf_bytes >= self.max_bytes {
            self.spill()?;
        }
        Ok(())
    }

    fn spill(&mut self) -> Result<()> {
        if self.buf.is_empty() {
            return Ok(());
        }
        self.buf.sort_by(cmp_items);
        let ext = if self.bcf { "bcf" } else { "vcf" };
        let p = self.temp_dir.join(format!(
            "kira-sort-{}-{}.{ext}.gz",
            std::process::id(),
            self.chunks.len()
        ));
        let mut w = crate::bgzf::BgzfWriter::with_compression(&p, flate2::Compression::fast())
            .with_context(|| format!("create temp {}", p.display()))?;
        for item in &self.buf {
            match &item.data {
                Data::Text(l) => {
                    w.write_all(l.as_bytes())?;
                    w.write_all(b"\n")?;
                }
                Data::Bcf(b) => w.write_all(b)?,
            }
        }
        w.finish()
            .with_context(|| format!("write temp {}", p.display()))?;
        self.chunks.push(p);
        self.buf.clear();
        self.buf_bytes = 0;
        Ok(())
    }

    /// Emit every record in order; the in-memory buffer joins the k-way merge
    /// as one more source, so a single chunk never touches the disk twice.
    /// `key_of` recomputes the position key of a record read back from a
    /// chunk.
    fn finish(
        mut self,
        key_of: &dyn Fn(&Data) -> u64,
        emit: &mut dyn FnMut(&Data) -> Result<()>,
    ) -> Result<()> {
        self.buf.sort_by(cmp_items);
        if self.chunks.is_empty() {
            for item in &self.buf {
                emit(&item.data)?;
            }
            return Ok(());
        }
        let mut sources: Vec<Source> = Vec::with_capacity(self.chunks.len() + 1);
        for p in &self.chunks {
            let f = File::open(p).with_context(|| format!("open temp {}", p.display()))?;
            sources.push(Source::Chunk(noodles_bgzf::io::Reader::new(f), self.bcf));
        }
        let mem: Vec<Item> = std::mem::take(&mut self.buf);
        sources.push(Source::Memory(mem.into_iter()));
        let mut heap: BinaryHeap<HeapEntry> = BinaryHeap::new();
        for (i, s) in sources.iter_mut().enumerate() {
            if let Some(item) = s.next_item(key_of)? {
                heap.push(HeapEntry { item, src: i });
            }
        }
        while let Some(e) = heap.pop() {
            emit(&e.item.data)?;
            if let Some(item) = sources[e.src].next_item(key_of)? {
                heap.push(HeapEntry { item, src: e.src });
            }
        }
        for p in &self.chunks {
            let _ = std::fs::remove_file(p);
        }
        Ok(())
    }
}

enum Source {
    Chunk(noodles_bgzf::io::Reader<File>, bool),
    Memory(std::vec::IntoIter<Item>),
}

impl Source {
    /// Chunks are sorted, so records come back in key order; the key is
    /// re-derived for the heap and `idx` restarts per chunk.
    fn next_item(&mut self, key_of: &dyn Fn(&Data) -> u64) -> Result<Option<Item>> {
        let data = match self {
            Source::Memory(it) => return Ok(it.next()),
            Source::Chunk(r, false) => {
                let mut line = String::new();
                if r.read_line(&mut line)? == 0 {
                    return Ok(None);
                }
                Data::Text(line.trim_end_matches(['\r', '\n']).to_string())
            }
            Source::Chunk(r, true) => {
                let mut hdr = [0u8; 8];
                match r.read_exact(&mut hdr) {
                    Ok(()) => {}
                    Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => return Ok(None),
                    Err(e) => return Err(e.into()),
                }
                let ls = u32::from_le_bytes(hdr[0..4].try_into().unwrap()) as usize;
                let li = u32::from_le_bytes(hdr[4..8].try_into().unwrap()) as usize;
                let mut framed = Vec::with_capacity(8 + ls + li);
                framed.extend_from_slice(&hdr);
                framed.resize(8 + ls + li, 0);
                r.read_exact(&mut framed[8..])?;
                Data::Bcf(framed)
            }
        };
        Ok(Some(Item {
            key: key_of(&data),
            idx: 0,
            data,
        }))
    }
}

/// The packed key of an encoded record, from its `shared` block.
fn bcf_key(data: &Data) -> u64 {
    match data {
        Data::Bcf(framed) => split_framed(framed)
            .ok()
            .and_then(|(s, _)| record_meta(s))
            .map(|m| pack_key(m.rid.max(0) as u64, m.pos.max(0) as u32 + 1))
            .unwrap_or(u64::MAX),
        Data::Text(_) => u64::MAX,
    }
}

/// `src` keeps the merge stable between sources.
struct HeapEntry {
    item: Item,
    src: usize,
}

impl HeapEntry {
    fn ordering(&self, o: &Self) -> Ordering {
        cmp_items(&self.item, &o.item).then(self.src.cmp(&o.src))
    }
}

impl Eq for HeapEntry {}
impl PartialEq for HeapEntry {
    fn eq(&self, o: &Self) -> bool {
        self.ordering(o) == Ordering::Equal
    }
}
impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for HeapEntry {
    // Reversed so the BinaryHeap pops the smallest key first.
    fn cmp(&self, o: &Self) -> Ordering {
        o.ordering(self)
    }
}

#[cfg(test)]
#[path = "../../../tests/unit/cli_commands_sort.rs"]
mod tests;
