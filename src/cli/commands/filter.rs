use anyhow::Result;
use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, mpsc};
use std::time::{Duration, Instant};

use crate::bgzf::BgzfWriter;
use crate::cli::args::FilterArgs;
use crate::filter::FilterEngine;
use crate::filter_arch;
use crate::kbi::KbiIndex;
use crate::util::{Region, chr_name_to_id};
use crate::vcf::VcfParser;
use crate::vcf::{VcfReader, VcfRecord};

pub fn cmd_filter(args: &FilterArgs) -> Result<()> {
    let total_start = Instant::now();
    eprintln!(
        "filter start backend={} gpu_requested={} opencl_requested={} gpu_active={} opencl_active={}",
        filter_arch::arch_name(),
        args.gpu,
        args.opencl,
        false,
        false
    );
    let (expr, exclude) = resolve_expr(args);

    let mut reader = VcfReader::open(&args.input)?;
    let mut headers = reader.header()?;
    ensure_filter_header(&mut headers, "PASS", "All filters passed");

    if (args.mask.is_some() || args.mask_file.is_some()) && args.soft_filter.is_none() {
        anyhow::bail!("--mask and --mask-file require --soft-filter");
    }

    let engine = FilterEngine::new(&headers, expr.as_deref(), exclude)?;
    let filter_name = prepare_filter_header(&mut headers, args, expr.as_deref())?;
    let snp_gap = parse_snp_gap(args.snp_gap.as_deref())?;
    if let Some(ref cfg) = snp_gap {
        let desc = format!("SNP within {} bp of {}", cfg.gap, cfg.type_str);
        ensure_filter_header(&mut headers, "SnpGap", &desc);
    }
    if let Some(gap) = args.indel_gap {
        let desc = format!("Indel within {} bp of an indel", gap);
        ensure_filter_header(&mut headers, "IndelGap", &desc);
    }

    let mask = load_mask(args)?;
    let mode = parse_mode(args.mode.as_deref());
    let set_gts = parse_set_gts(args.set_gts.as_deref())?;

    let header_start = Instant::now();
    let mut writer = open_writer(args.output.as_ref(), args.output_type.as_deref())?;
    for h in &headers {
        writer.write_all(h.as_bytes())?;
        writer.write_all(b"\n")?;
    }
    let header_time = header_start.elapsed();

    let use_parallel =
        args.threads.unwrap_or(0) > 1 && snp_gap.is_none() && args.indel_gap.is_none();
    if use_parallel {
        let stats = run_parallel_filter(
            &mut reader,
            &mut writer,
            Arc::new(engine),
            mask.map(Arc::new),
            mode,
            filter_name.clone(),
            set_gts,
            args.threads.unwrap_or(0),
        )?;
        writer.finish()?;
        let total_time = total_start.elapsed();
        eprintln!(
            "filter done records={} dropped={} header_ms={} eval_ms={} mask_ms={} mode_ms={} setgt_ms={} gap_ms=0 flush_ms=0 write_ms={} total_ms={}",
            stats.records(),
            stats.dropped(),
            duration_ms(header_time),
            duration_ms(Duration::from_nanos(stats.eval_ns())),
            duration_ms(Duration::from_nanos(stats.mask_ns())),
            duration_ms(Duration::from_nanos(stats.mode_ns())),
            duration_ms(Duration::from_nanos(stats.setgt_ns())),
            duration_ms(Duration::from_nanos(stats.write_ns())),
            duration_ms(total_time)
        );
        return Ok(());
    }

    let mut stats = FilterStats::default();
    let mut buffer: VecDeque<BufferedRecord> = VecDeque::new();
    while let Some(rec) = reader.next_record()? {
        stats.records += 1;
        let eval_start = Instant::now();
        let res = engine.eval(&rec)?;
        stats.eval_time += eval_start.elapsed();
        let mut pass = res.pass_site;
        if let Some(ref mask_cfg) = mask {
            let mask_start = Instant::now();
            if !mask_pass(mask_cfg, &rec) {
                pass = false;
            }
            stats.mask_time += mask_start.elapsed();
        }

        if !pass && filter_name.is_none() && set_gts.is_none() {
            stats.filtered += 1;
            continue;
        }

        let mode_start = Instant::now();
        let mut out_rec = rec;
        match apply_filter_mode(&out_rec.filter, pass, filter_name.as_deref(), &mode) {
            Some(f) => out_rec.filter = f,
            None => continue,
        }
        stats.filter_mode_time += mode_start.elapsed();

        if let Some(mode) = set_gts {
            let set_start = Instant::now();
            set_genotypes(&mut out_rec, pass, res.pass_samples.as_ref(), mode);
            stats.set_gts_time += set_start.elapsed();
        }

        if snp_gap.is_some() || args.indel_gap.is_some() {
            let gap_start = Instant::now();
            let var_type = variant_type_flags(&out_rec);
            let ref_len = out_rec.ref_allele.len() as u32;
            if let Some(last) = buffer.back() {
                if last.rec.chrom != out_rec.chrom {
                    let flush_start = Instant::now();
                    flush_buffer(&mut buffer, &mut writer, filter_name.is_some())?;
                    stats.flush_time += flush_start.elapsed();
                }
            }
            buffer.push_back(BufferedRecord {
                rec: out_rec,
                var_type,
                ref_len,
                snp_gap_set: false,
                indel_gap_set: false,
                indel_gap_filtered: false,
            });
            buffered_filters(
                &mut buffer,
                snp_gap.as_ref(),
                args.indel_gap,
                filter_name.is_some(),
                false,
                &mut writer,
            )?;
            stats.gap_time += gap_start.elapsed();
        } else {
            let write_start = Instant::now();
            write_record(&mut writer, &out_rec)?;
            stats.write_time += write_start.elapsed();
        }
    }

    if snp_gap.is_some() || args.indel_gap.is_some() {
        let gap_start = Instant::now();
        buffered_filters(
            &mut buffer,
            snp_gap.as_ref(),
            args.indel_gap,
            filter_name.is_some(),
            true,
            &mut writer,
        )?;
        stats.gap_time += gap_start.elapsed();
        let flush_start = Instant::now();
        flush_buffer(&mut buffer, &mut writer, filter_name.is_some())?;
        stats.flush_time += flush_start.elapsed();
    }

    writer.finish()?;
    let total_time = total_start.elapsed();
    eprintln!(
        "filter done records={} dropped={} header_ms={} eval_ms={} mask_ms={} mode_ms={} setgt_ms={} gap_ms={} flush_ms={} write_ms={} total_ms={}",
        stats.records,
        stats.filtered,
        duration_ms(header_time),
        duration_ms(stats.eval_time),
        duration_ms(stats.mask_time),
        duration_ms(stats.filter_mode_time),
        duration_ms(stats.set_gts_time),
        duration_ms(stats.gap_time),
        duration_ms(stats.flush_time),
        duration_ms(stats.write_time),
        duration_ms(total_time)
    );
    Ok(())
}

fn resolve_expr(args: &FilterArgs) -> (Option<String>, bool) {
    if let Some(expr) = &args.include {
        return (Some(expr.clone()), false);
    }
    if let Some(expr) = &args.exclude {
        return (Some(expr.clone()), true);
    }
    if let Some(expr) = &args.expr {
        return (Some(expr.clone()), false);
    }
    (None, false)
}

#[derive(Clone, Copy)]
struct FilterMode {
    add: bool,
    reset: bool,
}

fn parse_mode(mode: Option<&str>) -> FilterMode {
    let mut m = FilterMode {
        add: false,
        reset: false,
    };
    if let Some(s) = mode {
        if s.contains('+') {
            m.add = true;
        }
        if s.contains('x') {
            m.reset = true;
        }
    }
    m
}

fn apply_filter_mode(
    current: &str,
    pass: bool,
    soft_filter: Option<&str>,
    mode: &FilterMode,
) -> Option<String> {
    let no_filter = current == "." || current == "PASS" || current.is_empty();
    if pass {
        if mode.reset || no_filter {
            return Some("PASS".to_string());
        }
        return Some(current.to_string());
    }
    let Some(name) = soft_filter else {
        return Some(current.to_string());
    };
    if mode.add {
        if no_filter {
            return Some(name.to_string());
        }
        let mut parts: Vec<&str> = current.split(';').filter(|s| !s.is_empty()).collect();
        if !parts.iter().any(|s| *s == name) {
            parts.push(name);
        }
        return Some(parts.join(";"));
    }
    Some(name.to_string())
}

fn prepare_filter_header(
    headers: &mut Vec<String>,
    args: &FilterArgs,
    expr: Option<&str>,
) -> Result<Option<String>> {
    let Some(name) = args.soft_filter.clone() else {
        return Ok(None);
    };

    let mut filter_name = name;
    if filter_name == "+" {
        filter_name = unique_filter_name(headers);
    }

    let desc = if let Some(expr) = expr {
        format!(
            "Set if {}true: {}",
            if args.exclude.is_some() { "" } else { "not " },
            expr
        )
    } else if args.mask.is_some() || args.mask_file.is_some() {
        "Record masked by region".to_string()
    } else {
        "Set if filter expression is true".to_string()
    };

    ensure_filter_header(headers, &filter_name, &desc);
    Ok(Some(filter_name))
}

fn unique_filter_name(headers: &[String]) -> String {
    let mut used = HashSet::new();
    for h in headers {
        if let Some(id) = extract_filter_id(h) {
            used.insert(id);
        }
    }
    let mut i = 1;
    loop {
        let name = format!("Filter{}", i);
        if !used.contains(&name) {
            return name;
        }
        i += 1;
    }
}

fn extract_filter_id(line: &str) -> Option<String> {
    if !line.starts_with("##FILTER=<") {
        return None;
    }
    let id_pos = line.find("ID=")? + 3;
    let rest = &line[id_pos..];
    let end = rest.find(',').or_else(|| rest.find('>'))?;
    Some(rest[..end].to_string())
}

fn ensure_filter_header(headers: &mut Vec<String>, id: &str, desc: &str) {
    for h in headers.iter() {
        if let Some(existing) = extract_filter_id(h) {
            if existing == id {
                return;
            }
        }
    }

    let line = format!(
        "##FILTER=<ID={},Description=\"{}\">",
        id,
        desc.replace('"', "\\\"")
    );
    let mut inserted = false;
    for i in 0..headers.len() {
        if headers[i].starts_with("#CHROM") {
            headers.insert(i, line.clone());
            inserted = true;
            break;
        }
    }
    if !inserted {
        headers.push(line);
    }
}

enum OutputWriter {
    Stdout(BufWriter<io::Stdout>),
    File(BufWriter<File>),
    Bgzf(BgzfWriter),
}

impl OutputWriter {
    fn write_all(&mut self, data: &[u8]) -> io::Result<()> {
        match self {
            OutputWriter::Stdout(w) => w.write_all(data),
            OutputWriter::File(w) => w.write_all(data),
            OutputWriter::Bgzf(w) => w.write_all(data),
        }
    }

    fn finish(self) -> io::Result<()> {
        match self {
            OutputWriter::Stdout(mut w) => w.flush(),
            OutputWriter::File(mut w) => w.flush(),
            OutputWriter::Bgzf(w) => w.finish(),
        }
    }
}

fn open_writer(path: Option<&PathBuf>, output_type: Option<&str>) -> Result<OutputWriter> {
    let out_type = output_type.unwrap_or("v");
    let use_bgzf = out_type.starts_with('z');
    if out_type.starts_with('b') || out_type.starts_with('u') {
        anyhow::bail!("BCF output not supported");
    }
    if use_bgzf {
        let Some(path) = path else {
            anyhow::bail!("BGZF output requires -o");
        };
        let writer = BgzfWriter::create(path)?;
        return Ok(OutputWriter::Bgzf(writer));
    }
    match path {
        Some(p) => Ok(OutputWriter::File(BufWriter::new(File::create(p)?))),
        None => Ok(OutputWriter::Stdout(BufWriter::new(io::stdout()))),
    }
}

fn write_record(writer: &mut OutputWriter, rec: &VcfRecord) -> Result<()> {
    let line = record_line(rec);
    writer.write_all(line.as_bytes())?;
    Ok(())
}

fn record_line(rec: &VcfRecord) -> String {
    let mut line = String::new();
    line.push_str(&rec.chrom);
    line.push('\t');
    line.push_str(&rec.pos.to_string());
    line.push('\t');
    line.push_str(&rec.id);
    line.push('\t');
    line.push_str(&rec.ref_allele);
    line.push('\t');
    line.push_str(&rec.alt);
    line.push('\t');
    line.push_str(&rec.qual);
    line.push('\t');
    line.push_str(&rec.filter);
    line.push('\t');
    line.push_str(&rec.info);
    if let Some(fmt) = &rec.format {
        line.push('\t');
        line.push_str(fmt);
        for sample in &rec.samples {
            line.push('\t');
            line.push_str(sample);
        }
    }
    line.push('\n');
    line
}

struct ParallelStats {
    records: AtomicU64,
    dropped: AtomicU64,
    eval_ns: AtomicU64,
    mask_ns: AtomicU64,
    mode_ns: AtomicU64,
    setgt_ns: AtomicU64,
    write_ns: AtomicU64,
}

impl ParallelStats {
    fn new() -> Self {
        Self {
            records: AtomicU64::new(0),
            dropped: AtomicU64::new(0),
            eval_ns: AtomicU64::new(0),
            mask_ns: AtomicU64::new(0),
            mode_ns: AtomicU64::new(0),
            setgt_ns: AtomicU64::new(0),
            write_ns: AtomicU64::new(0),
        }
    }

    fn records(&self) -> u64 {
        self.records.load(Ordering::Relaxed)
    }

    fn dropped(&self) -> u64 {
        self.dropped.load(Ordering::Relaxed)
    }

    fn eval_ns(&self) -> u64 {
        self.eval_ns.load(Ordering::Relaxed)
    }

    fn mask_ns(&self) -> u64 {
        self.mask_ns.load(Ordering::Relaxed)
    }

    fn mode_ns(&self) -> u64 {
        self.mode_ns.load(Ordering::Relaxed)
    }

    fn setgt_ns(&self) -> u64 {
        self.setgt_ns.load(Ordering::Relaxed)
    }

    fn write_ns(&self) -> u64 {
        self.write_ns.load(Ordering::Relaxed)
    }
}

fn run_parallel_filter(
    reader: &mut VcfReader,
    writer: &mut OutputWriter,
    engine: Arc<FilterEngine>,
    mask: Option<Arc<MaskSource>>,
    mode: FilterMode,
    filter_name: Option<String>,
    set_gts: Option<SetGtsMode>,
    threads: usize,
) -> Result<Arc<ParallelStats>> {
    let (res_tx, res_rx) = mpsc::channel::<Result<(u64, Option<String>)>>();
    let stats = Arc::new(ParallelStats::new());

    let mut handles = Vec::new();
    let mut worker_txs = Vec::with_capacity(threads);
    for _ in 0..threads {
        let (tx, rx) = mpsc::channel::<(u64, VcfRecord)>();
        worker_txs.push(tx);
        let engine = Arc::clone(&engine);
        let res_tx = res_tx.clone();
        let stats = Arc::clone(&stats);
        let mask = mask.clone();
        let filter_name = filter_name.clone();
        let mode = mode;
        let set_gts = set_gts;
        let handle = std::thread::spawn(move || {
            loop {
                let (idx, rec) = match rx.recv() {
                    Ok(v) => v,
                    Err(_) => break,
                };
                stats.records.fetch_add(1, Ordering::Relaxed);
                let eval_start = Instant::now();
                let res = match engine.eval(&rec) {
                    Ok(v) => v,
                    Err(e) => {
                        let _ = res_tx.send(Err(e));
                        continue;
                    }
                };
                stats
                    .eval_ns
                    .fetch_add(eval_start.elapsed().as_nanos() as u64, Ordering::Relaxed);
                let mut pass = res.pass_site;
                if let Some(ref mask_cfg) = mask {
                    let mask_start = Instant::now();
                    if !mask_pass(mask_cfg, &rec) {
                        pass = false;
                    }
                    stats
                        .mask_ns
                        .fetch_add(mask_start.elapsed().as_nanos() as u64, Ordering::Relaxed);
                }
                if !pass && filter_name.is_none() && set_gts.is_none() {
                    stats.dropped.fetch_add(1, Ordering::Relaxed);
                    let _ = res_tx.send(Ok((idx, None)));
                    continue;
                }
                let mode_start = Instant::now();
                let mut out_rec = rec;
                match apply_filter_mode(&out_rec.filter, pass, filter_name.as_deref(), &mode) {
                    Some(f) => out_rec.filter = f,
                    None => {
                        stats.dropped.fetch_add(1, Ordering::Relaxed);
                        let _ = res_tx.send(Ok((idx, None)));
                        continue;
                    }
                }
                stats
                    .mode_ns
                    .fetch_add(mode_start.elapsed().as_nanos() as u64, Ordering::Relaxed);
                if let Some(mode) = set_gts {
                    let set_start = Instant::now();
                    set_genotypes(&mut out_rec, pass, res.pass_samples.as_ref(), mode);
                    stats
                        .setgt_ns
                        .fetch_add(set_start.elapsed().as_nanos() as u64, Ordering::Relaxed);
                }
                let line = record_line(&out_rec);
                let _ = res_tx.send(Ok((idx, Some(line))));
            }
        });
        handles.push(handle);
    }
    drop(res_tx);

    let mut sent = 0u64;
    while let Some(rec) = reader.next_record()? {
        let idx = (sent as usize) % threads;
        worker_txs[idx].send((sent, rec))?;
        sent += 1;
    }
    drop(worker_txs);

    let mut next = 0u64;
    let mut pending: HashMap<u64, Option<String>> = HashMap::new();
    let mut received = 0u64;
    let mut write_ns = 0u64;
    while received < sent {
        let msg = match res_rx.recv() {
            Ok(v) => v,
            Err(_) => break,
        };
        let (idx, line) = msg?;
        received += 1;
        pending.insert(idx, line);
        while let Some(line) = pending.remove(&next) {
            if let Some(s) = line {
                let write_start = Instant::now();
                writer.write_all(s.as_bytes())?;
                write_ns += write_start.elapsed().as_nanos() as u64;
            }
            next += 1;
        }
    }

    for h in handles {
        let _ = h.join();
    }
    stats.write_ns.fetch_add(write_ns, Ordering::Relaxed);
    Ok(stats)
}

#[derive(Default)]
struct FilterStats {
    records: u64,
    filtered: u64,
    eval_time: Duration,
    mask_time: Duration,
    filter_mode_time: Duration,
    set_gts_time: Duration,
    gap_time: Duration,
    flush_time: Duration,
    write_time: Duration,
}

fn duration_ms(d: Duration) -> u64 {
    (d.as_secs_f64() * 1000.0) as u64
}

#[derive(Clone, Copy)]
enum SetGtsMode {
    Missing,
    Ref,
}

fn parse_set_gts(val: Option<&str>) -> Result<Option<SetGtsMode>> {
    let Some(v) = val else {
        return Ok(None);
    };
    match v {
        "." => Ok(Some(SetGtsMode::Missing)),
        "0" => Ok(Some(SetGtsMode::Ref)),
        _ => anyhow::bail!("The argument to -S not recognised: {}", v),
    }
}

struct MaskConfig {
    regions: HashMap<String, Vec<(u32, u32)>>,
    overlap: u8,
    negate: bool,
}

struct KbiMask {
    index: KbiIndex,
    path: PathBuf,
    overlap: u8,
    negate: bool,
    has_data: bool,
}

enum MaskSource {
    Regions(MaskConfig),
    Kbi(KbiMask),
}

fn load_mask(args: &FilterArgs) -> Result<Option<MaskSource>> {
    let mut negate = false;
    let mut regions = HashMap::new();
    if let Some(mask) = &args.mask {
        let mut raw = mask.as_str();
        if let Some(rest) = raw.strip_prefix('^') {
            negate = true;
            raw = rest;
        }
        for part in raw.split(',') {
            if let Some((chr, start, end)) = parse_region_str(part) {
                regions
                    .entry(chr)
                    .or_insert_with(Vec::new)
                    .push((start, end));
            }
        }
    }
    if let Some(path) = &args.mask_file {
        let mut path_str = path.to_string_lossy().to_string();
        if let Some(rest) = path_str.strip_prefix('^') {
            negate = true;
            path_str = rest.to_string();
        }
        let mask_path = PathBuf::from(&path_str);
        if let Some(kbi_mask) = try_load_kbi_mask(&mask_path, args.mask_overlap, negate)? {
            return Ok(Some(MaskSource::Kbi(kbi_mask)));
        }
        let file = File::open(&path_str)?;
        let reader = io::BufReader::new(file);
        for line in reader.lines() {
            let line = line?;
            if let Some((chr, start, end)) = parse_mask_line(&line) {
                regions
                    .entry(chr)
                    .or_insert_with(Vec::new)
                    .push((start, end));
            }
        }
    }
    if regions.is_empty() {
        return Ok(None);
    }
    let overlap = args.mask_overlap.unwrap_or(1);
    Ok(Some(MaskSource::Regions(MaskConfig {
        regions,
        overlap,
        negate,
    })))
}

fn parse_region_str(s: &str) -> Option<(String, u32, u32)> {
    let trimmed = s.trim();
    if trimmed.is_empty() || trimmed.starts_with('#') {
        return None;
    }
    if let Some(region) = Region::parse(trimmed) {
        let start = region.start.map(|v| v.saturating_sub(1)).unwrap_or(0);
        let end = region.end.map(|v| v.saturating_sub(1)).unwrap_or(u32::MAX);
        return Some((region.chr, start, end));
    }
    None
}

fn parse_mask_line(line: &str) -> Option<(String, u32, u32)> {
    let tab_count = line.as_bytes().iter().filter(|b| **b == b'\t').count();
    if tab_count >= 7 {
        let mut parser = VcfParser::new(line);
        let fields = parser.parse_standard_fields()?;
        let pos = fields.pos.parse::<u32>().ok()?;
        let pos0 = pos.saturating_sub(1);
        let end = pos0.saturating_add(fields.ref_allele.len().saturating_sub(1) as u32);
        return Some((fields.chrom.to_string(), pos0, end));
    }
    parse_region_str(line)
}

fn mask_pass(mask: &MaskSource, rec: &VcfRecord) -> bool {
    match mask {
        MaskSource::Regions(cfg) => mask_pass_regions(cfg, rec),
        MaskSource::Kbi(cfg) => mask_pass_kbi(cfg, rec),
    }
}

fn mask_pass_regions(mask: &MaskConfig, rec: &VcfRecord) -> bool {
    let (beg, end) = mask_bounds(rec, mask.overlap);
    let mut hit = false;
    if let Some(list) = mask.regions.get(&rec.chrom) {
        for (s, e) in list {
            if *s <= end && *e >= beg {
                hit = true;
                break;
            }
        }
    }
    if mask.negate { hit } else { !hit }
}

fn mask_pass_kbi(mask: &KbiMask, rec: &VcfRecord) -> bool {
    let (beg, end) = mask_bounds(rec, mask.overlap);
    let Some(chr_id) = chr_name_to_id(&rec.chrom) else {
        return true;
    };
    let mut hit = false;
    if !mask.has_data {
        hit = mask.index.has_range(chr_id, beg, end);
    } else {
        let hits = mask.index.range(chr_id, beg, end);
        for (_pos, offset) in hits {
            if let Ok(line) = crate::fetch_line(&mask.path, offset) {
                if let Some((mchr, mbeg, mend)) = parse_mask_vcf_bounds(&line, mask.overlap) {
                    if mchr == chr_id && mbeg <= end && mend >= beg {
                        hit = true;
                        break;
                    }
                }
            }
        }
    }
    if mask.negate { hit } else { !hit }
}

fn variant_bounds(rec: &VcfRecord, pos0: u32) -> (u32, u32) {
    let ref_bytes = rec.ref_allele.as_bytes();
    let mut off = ref_bytes.len();
    for alt in rec.alt.split(',') {
        let alt_bytes = alt.as_bytes();
        let mut j = 0usize;
        while j < ref_bytes.len() && j < alt_bytes.len() && ref_bytes[j] == alt_bytes[j] {
            j += 1;
        }
        if j < off {
            off = j;
        }
    }
    let beg = pos0.saturating_add(off as u32);
    let end = pos0.saturating_add(rec.ref_allele.len().saturating_sub(1) as u32);
    (beg, end)
}

fn mask_bounds(rec: &VcfRecord, overlap: u8) -> (u32, u32) {
    let pos0 = rec.pos.saturating_sub(1);
    match overlap {
        0 => (pos0, pos0),
        2 => variant_bounds(rec, pos0),
        _ => {
            let end = pos0.saturating_add(rec.ref_allele.len().saturating_sub(1) as u32);
            (pos0, end)
        }
    }
}

fn parse_mask_vcf_bounds(line: &str, overlap: u8) -> Option<(u8, u32, u32)> {
    let mut parser = VcfParser::new(line);
    let fields = parser.parse_standard_fields()?;
    let chr_id = chr_name_to_id(fields.chrom)?;
    let pos = fields.pos.parse::<u32>().ok()?;
    let pos0 = pos.saturating_sub(1);
    let (beg, end) = match overlap {
        0 => (pos0, pos0),
        2 => variant_bounds_fields(fields.ref_allele, fields.alt, pos0),
        _ => {
            let end = pos0.saturating_add(fields.ref_allele.len().saturating_sub(1) as u32);
            (pos0, end)
        }
    };
    Some((chr_id, beg, end))
}

fn variant_bounds_fields(ref_allele: &str, alt: &str, pos0: u32) -> (u32, u32) {
    let ref_bytes = ref_allele.as_bytes();
    let mut off = ref_bytes.len();
    for alt_allele in alt.split(',') {
        let alt_bytes = alt_allele.as_bytes();
        let mut j = 0usize;
        while j < ref_bytes.len() && j < alt_bytes.len() && ref_bytes[j] == alt_bytes[j] {
            j += 1;
        }
        if j < off {
            off = j;
        }
    }
    let beg = pos0.saturating_add(off as u32);
    let end = pos0.saturating_add(ref_bytes.len().saturating_sub(1) as u32);
    (beg, end)
}

fn try_load_kbi_mask(
    mask_path: &PathBuf,
    overlap: Option<u8>,
    negate: bool,
) -> Result<Option<KbiMask>> {
    if mask_path.extension().map(|e| e == "kbi").unwrap_or(false) {
        let index = KbiIndex::load(mask_path)?;
        let (data_path, has_data) = resolve_kbi_data_path(mask_path);
        return Ok(Some(KbiMask {
            index,
            path: data_path,
            overlap: overlap.unwrap_or(1),
            negate,
            has_data,
        }));
    }
    let kbi_path = mask_path.with_extension("kbi");
    if kbi_path.exists() {
        let index = KbiIndex::load(&kbi_path)?;
        return Ok(Some(KbiMask {
            index,
            path: mask_path.clone(),
            overlap: overlap.unwrap_or(1),
            negate,
            has_data: true,
        }));
    }
    Ok(None)
}

fn resolve_kbi_data_path(kbi_path: &PathBuf) -> (PathBuf, bool) {
    let base = kbi_path.with_extension("");
    if base.exists() {
        return (base, true);
    }
    if base.extension().map(|e| e == "vcf").unwrap_or(false) {
        let gz = base.with_extension("vcf.gz");
        if gz.exists() {
            return (gz, true);
        }
    }
    (kbi_path.clone(), false)
}

const VCF_SNP: u32 = 1;
const VCF_MNP: u32 = 2;
const VCF_INDEL: u32 = 4;
const VCF_OTHER: u32 = 8;
const VCF_BND: u32 = 16;
const VCF_OVERLAP: u32 = 32;

struct SnpGapConfig {
    gap: u32,
    type_mask: u32,
    type_str: String,
}

fn parse_snp_gap(val: Option<&str>) -> Result<Option<SnpGapConfig>> {
    let Some(raw) = val else {
        return Ok(None);
    };
    let mut parts = raw.splitn(2, ':');
    let gap_str = parts.next().unwrap_or("");
    let gap: u32 = gap_str.parse()?;
    let mut type_mask = 0u32;
    let type_str = if let Some(rest) = parts.next() {
        let mut names = Vec::new();
        for part in rest.split(',') {
            let p = part.trim().to_ascii_lowercase();
            if p.is_empty() {
                continue;
            }
            match p.as_str() {
                "indel" => type_mask |= VCF_INDEL,
                "mnp" => type_mask |= VCF_MNP,
                "bnd" => type_mask |= VCF_BND,
                "other" => type_mask |= VCF_OTHER,
                "overlap" => type_mask |= VCF_OVERLAP,
                _ => anyhow::bail!("Could not parse \"{}\" in \"--SnpGap {}\"", part, raw),
            }
            names.push(p);
        }
        if names.is_empty() {
            anyhow::bail!("Could not parse argument: --SnpGap {}", raw);
        }
        names.join(",")
    } else {
        type_mask = VCF_INDEL;
        "indel".to_string()
    };
    Ok(Some(SnpGapConfig {
        gap,
        type_mask,
        type_str,
    }))
}

struct BufferedRecord {
    rec: VcfRecord,
    var_type: u32,
    ref_len: u32,
    snp_gap_set: bool,
    indel_gap_set: bool,
    indel_gap_filtered: bool,
}

fn buffered_filters(
    buffer: &mut VecDeque<BufferedRecord>,
    snp_gap: Option<&SnpGapConfig>,
    indel_gap: Option<u32>,
    soft_filter: bool,
    finalize: bool,
    writer: &mut OutputWriter,
) -> Result<()> {
    if buffer.is_empty() {
        return Ok(());
    }
    let mut k_flush = 1usize;
    if let Some(gap) = indel_gap {
        k_flush = 0;
        let mut last_to: Option<u32> = None;
        let mut idx = 0usize;
        for rec in buffer.iter_mut() {
            if let Some(last) = last_to {
                if last < rec.rec.pos.saturating_sub(1) {
                    break;
                }
            }
            k_flush += 1;
            if rec.var_type & VCF_INDEL == 0 {
                idx += 1;
                continue;
            }
            rec.indel_gap_set = true;
            let to = rec
                .rec
                .pos
                .saturating_sub(1)
                .saturating_add(gap)
                .saturating_add(rec.ref_len.saturating_sub(1));
            last_to = Some(to);
            idx += 1;
        }
        if idx == buffer.len() && last_to.is_some() {
            k_flush = 0;
        }
        if k_flush > 0 {
            let mut max_qual = f64::NEG_INFINITY;
            let mut max_ac = i64::MIN;
            let mut best_idx: Option<usize> = None;
            for (i, rec) in buffer.iter().take(k_flush).enumerate() {
                if !rec.indel_gap_set {
                    continue;
                }
                let qual = rec.rec.qual.parse::<f64>().unwrap_or(f64::NEG_INFINITY);
                if qual > 0.0 {
                    if qual > max_qual {
                        max_qual = qual;
                        best_idx = Some(i);
                    }
                }
                if max_qual == f64::NEG_INFINITY {
                    let ac = first_alt_ac(&rec.rec).unwrap_or(i64::MIN);
                    if ac > max_ac {
                        max_ac = ac;
                        best_idx = Some(i);
                    }
                }
                if best_idx.is_none() {
                    best_idx = Some(i);
                }
            }
            for (i, rec) in buffer.iter_mut().take(k_flush).enumerate() {
                if !rec.indel_gap_set {
                    continue;
                }
                if Some(i) != best_idx {
                    rec.indel_gap_filtered = true;
                    rec.rec.filter = add_filter_tag(&rec.rec.filter, "IndelGap");
                }
            }
        }
    }
    if finalize {
        flush_n_records(buffer, buffer.len(), soft_filter, writer)?;
        return Ok(());
    }
    if let Some(cfg) = snp_gap {
        let mut j_flush = 0usize;
        let last_from = buffer.back().unwrap().rec.pos.saturating_sub(1);
        let var_type = buffer.back().unwrap().var_type;
        for rec in buffer.iter_mut() {
            let rec_to = rec
                .rec
                .pos
                .saturating_sub(1)
                .saturating_add(rec.ref_len.saturating_sub(1));
            if rec_to.saturating_add(cfg.gap) < last_from {
                j_flush += 1;
            } else if (var_type & cfg.type_mask) != 0
                && (rec.var_type & VCF_SNP) != 0
                && !rec.snp_gap_set
            {
                rec.snp_gap_set = true;
                rec.rec.filter = add_filter_tag(&rec.rec.filter, "SnpGap");
            } else if (var_type & VCF_SNP) != 0 && (rec.var_type & cfg.type_mask) != 0 {
                if let Some(last) = buffer.back_mut() {
                    last.snp_gap_set = true;
                    last.rec.filter = add_filter_tag(&last.rec.filter, "SnpGap");
                }
                break;
            }
        }
        let flush_n = if k_flush < j_flush { k_flush } else { j_flush };
        flush_n_records(buffer, flush_n, soft_filter, writer)?;
        return Ok(());
    }
    flush_n_records(buffer, k_flush, soft_filter, writer)?;
    Ok(())
}

fn flush_buffer(
    buffer: &mut VecDeque<BufferedRecord>,
    writer: &mut OutputWriter,
    soft_filter: bool,
) -> Result<()> {
    let n = buffer.len();
    flush_n_records(buffer, n, soft_filter, writer)
}

fn flush_n_records(
    buffer: &mut VecDeque<BufferedRecord>,
    n: usize,
    soft_filter: bool,
    writer: &mut OutputWriter,
) -> Result<()> {
    for _ in 0..n {
        if let Some(rec) = buffer.pop_front() {
            if !soft_filter && (rec.snp_gap_set || rec.indel_gap_filtered) {
                continue;
            }
            write_record(writer, &rec.rec)?;
        }
    }
    Ok(())
}

fn add_filter_tag(current: &str, tag: &str) -> String {
    if current.is_empty() || current == "." || current == "PASS" {
        return tag.to_string();
    }
    let mut parts: Vec<&str> = current.split(';').filter(|s| !s.is_empty()).collect();
    if !parts.iter().any(|s| *s == tag) {
        parts.push(tag);
    }
    parts.join(";")
}

fn variant_type_flags(rec: &VcfRecord) -> u32 {
    let ref_bytes = rec.ref_allele.as_bytes();
    let mut mask = 0u32;
    for alt in rec.alt.split(',') {
        let t = variant_type_for_alt(ref_bytes, alt.as_bytes());
        mask |= t;
    }
    mask
}

fn variant_type_for_alt(ref_bytes: &[u8], alt_bytes: &[u8]) -> u32 {
    if alt_bytes == b"*" {
        return VCF_OVERLAP;
    }
    if ref_bytes.len() == 1 && alt_bytes.len() == 1 {
        let a = alt_bytes[0];
        let r = ref_bytes[0];
        if a == b'.' || eq_icase(a, r) {
            return 0;
        }
        if a == b'X' || a == b'x' {
            return 0;
        }
        return VCF_SNP;
    }
    if alt_bytes.first() == Some(&b'<') {
        if alt_bytes.len() >= 3 {
            if (alt_bytes[1] == b'X' || alt_bytes[1] == b'x') && alt_bytes[2] == b'>' {
                return 0;
            }
            if alt_bytes[1] == b'*' && alt_bytes[2] == b'>' {
                return 0;
            }
            if alt_bytes.len() >= 9 && alt_bytes[1..].starts_with(b"NON_REF>") {
                return 0;
            }
        }
        return VCF_OTHER;
    }
    if alt_bytes[0] == b']' || alt_bytes[0] == b'[' {
        return VCF_BND;
    }

    let mut r_i = 0usize;
    let mut a_i = 0usize;
    while r_i < ref_bytes.len() && a_i < alt_bytes.len() && eq_icase(ref_bytes[r_i], alt_bytes[a_i])
    {
        r_i += 1;
        a_i += 1;
    }
    if a_i < alt_bytes.len() && r_i == ref_bytes.len() {
        if alt_bytes[alt_bytes.len() - 1] == b']' || alt_bytes[alt_bytes.len() - 1] == b'[' {
            return VCF_BND;
        }
        return VCF_INDEL;
    } else if r_i < ref_bytes.len() && a_i == alt_bytes.len() {
        return VCF_INDEL;
    } else if r_i == ref_bytes.len() && a_i == alt_bytes.len() {
        return 0;
    }

    let mut re = ref_bytes.len() - 1;
    let mut ae = alt_bytes.len() - 1;
    if alt_bytes[ae] == b']' || alt_bytes[ae] == b'[' {
        return VCF_BND;
    }
    while re > r_i && ae > a_i && eq_icase(ref_bytes[re], alt_bytes[ae]) {
        re -= 1;
        ae -= 1;
    }
    if ae == a_i {
        if re == r_i {
            return VCF_SNP;
        }
        if eq_icase(ref_bytes[re], alt_bytes[ae]) {
            return VCF_INDEL;
        }
        return VCF_OTHER;
    }
    if re == r_i {
        if eq_icase(ref_bytes[re], alt_bytes[ae]) {
            return VCF_INDEL;
        }
        return VCF_OTHER;
    }
    if (re as isize - r_i as isize) == (ae as isize - a_i as isize) {
        VCF_MNP
    } else {
        VCF_OTHER
    }
}

fn eq_icase(a: u8, b: u8) -> bool {
    a.to_ascii_uppercase() == b.to_ascii_uppercase()
}

fn set_genotypes(
    rec: &mut VcfRecord,
    pass: bool,
    pass_samples: Option<&Vec<bool>>,
    mode: SetGtsMode,
) {
    let Some(fmt) = &rec.format else {
        return;
    };
    let gt_index = fmt.split(':').position(|f| f == "GT");
    let Some(gt_index) = gt_index else {
        return;
    };
    let nsamples = rec.samples.len();
    let mut targets = vec![false; nsamples];
    if let Some(pass_vec) = pass_samples {
        let mut any_fail = false;
        for i in 0..nsamples.min(pass_vec.len()) {
            if !pass_vec[i] {
                targets[i] = true;
                any_fail = true;
            }
        }
        if !any_fail {
            return;
        }
    } else if pass {
        return;
    } else {
        for i in 0..nsamples {
            targets[i] = true;
        }
    }
    for i in 0..nsamples {
        if !targets[i] {
            continue;
        }
        let mut parts: Vec<String> = rec.samples[i].split(':').map(|s| s.to_string()).collect();
        if gt_index >= parts.len() {
            continue;
        }
        let gt = parts[gt_index].as_str();
        let sep = if gt.contains('|') { '|' } else { '/' };
        let ploidy = gt_ploidy(gt);
        let new_gt = match mode {
            SetGtsMode::Missing => build_gt(ploidy, sep, "."),
            SetGtsMode::Ref => build_gt(ploidy, sep, "0"),
        };
        parts[gt_index] = new_gt;
        rec.samples[i] = parts.join(":");
    }
    update_ac_an(rec);
}

fn gt_ploidy(gt: &str) -> usize {
    if gt.is_empty() {
        return 0;
    }
    if gt == "." {
        return 1;
    }
    gt.split(|c| c == '/' || c == '|').count()
}

fn build_gt(ploidy: usize, sep: char, allele: &str) -> String {
    if ploidy <= 1 {
        return allele.to_string();
    }
    let mut out = String::new();
    for i in 0..ploidy {
        if i > 0 {
            out.push(sep);
        }
        out.push_str(allele);
    }
    out
}

fn update_ac_an(rec: &mut VcfRecord) {
    let mut fields = parse_info(rec.info.as_str());
    let mut has_ac = false;
    let mut has_an = false;
    for (k, _) in &fields {
        if k == "AC" {
            has_ac = true;
        } else if k == "AN" {
            has_an = true;
        }
    }
    if !has_ac && !has_an {
        return;
    }
    let alt_count = alt_count(rec);
    let (an, ac) = compute_ac_an_from_gts(rec, alt_count);
    for (k, v) in fields.iter_mut() {
        if k == "AN" {
            *v = Some(an.to_string());
        } else if k == "AC" {
            if ac.is_empty() {
                *v = Some("0".to_string());
            } else {
                let vals: Vec<String> = ac.iter().map(|v| v.to_string()).collect();
                *v = Some(vals.join(","));
            }
        }
    }
    rec.info = rebuild_info(fields);
}

fn parse_info(info: &str) -> Vec<(String, Option<String>)> {
    if info.is_empty() || info == "." {
        return Vec::new();
    }
    let mut out = Vec::new();
    for part in info.split(';') {
        if part.is_empty() {
            continue;
        }
        if let Some((k, v)) = part.split_once('=') {
            out.push((k.to_string(), Some(v.to_string())));
        } else {
            out.push((part.to_string(), None));
        }
    }
    out
}

fn rebuild_info(fields: Vec<(String, Option<String>)>) -> String {
    if fields.is_empty() {
        return ".".to_string();
    }
    let mut parts = Vec::with_capacity(fields.len());
    for (k, v) in fields {
        if let Some(val) = v {
            parts.push(format!("{}={}", k, val));
        } else {
            parts.push(k);
        }
    }
    parts.join(";")
}

fn compute_ac_an_from_gts(rec: &VcfRecord, alt_count: usize) -> (i64, Vec<i64>) {
    let Some(fmt) = &rec.format else {
        return (0, vec![0; alt_count]);
    };
    let gt_index = match fmt.split(':').position(|f| f == "GT") {
        Some(i) => i,
        None => return (0, vec![0; alt_count]),
    };
    let mut an = 0i64;
    let mut ac = vec![0i64; alt_count];
    for sample in &rec.samples {
        let parts: Vec<&str> = sample.split(':').collect();
        if gt_index >= parts.len() {
            continue;
        }
        let gt = parts[gt_index];
        for allele in gt.split(|c| c == '/' || c == '|') {
            if allele == "." || allele.is_empty() {
                continue;
            }
            if let Ok(idx) = allele.parse::<i64>() {
                if idx >= 0 {
                    an += 1;
                }
                if idx > 0 {
                    let a = (idx - 1) as usize;
                    if a < ac.len() {
                        ac[a] += 1;
                    }
                }
            }
        }
    }
    (an, ac)
}

fn first_alt_ac(rec: &VcfRecord) -> Option<i64> {
    let alt_count = alt_count(rec);
    if alt_count == 0 {
        return None;
    }
    let (_, ac) = compute_ac_an_from_gts(rec, alt_count);
    ac.first().copied()
}

fn alt_count(rec: &VcfRecord) -> usize {
    if rec.alt.is_empty() || rec.alt == "." {
        0
    } else {
        rec.alt.split(',').count()
    }
}
