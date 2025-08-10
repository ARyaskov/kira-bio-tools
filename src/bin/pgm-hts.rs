use std::{
    cmp,
    fs::File,
    io::{self, Read},
    path::{Path, PathBuf},
    time::{Duration, Instant},
};

use clap::{ArgAction, Parser, Subcommand};

#[cfg(feature = "mmap")]
use memmap2::{Mmap, MmapOptions};

use memchr::{memchr, memchr_iter};
use rayon::prelude::*;

use pgm_index::persistence as persist;

// ---------- mmap helpers ----------
#[cfg(feature = "mmap")]
fn mmap_read_file(path: &Path) -> io::Result<Mmap> {
    let f = File::open(path)?;
    let mmap = unsafe { MmapOptions::new().map(&f)? };
    Ok(mmap)
}
#[cfg(not(feature = "mmap"))]
fn mmap_read_file(path: &Path) -> io::Result<Vec<u8>> {
    let mut f = File::open(path)?;
    let mut buf = Vec::new();
    f.read_to_end(&mut buf)?;
    Ok(buf)
}

// ---------- small utils ----------
#[inline]
fn parse_u64_fast(bytes: &[u8]) -> u64 {
    let mut v = 0u64;
    for &c in bytes {
        if c < b'0' || c > b'9' {
            break;
        }
        v = v * 10 + (c - b'0') as u64;
    }
    v
}

#[inline]
fn find_tab(bytes: &[u8]) -> Option<usize> {
    memchr(b'\t', bytes)
}

#[inline]
fn parse_line_chrom_pos_inline(line: &[u8]) -> Option<(usize, u64)> {
    if line.is_empty() || line[0] == b'#' {
        return None;
    }
    let t1 = find_tab(line)?;
    let t2 = find_tab(&line[t1 + 1..])?;
    let pos_s = &line[t1 + 1..t1 + 1 + t2];
    let pos = parse_u64_fast(pos_s);
    if pos == 0 {
        return None;
    }
    Some((t1, pos))
}

#[derive(Clone, Copy)]
struct LineRec {
    chrom_start: usize,
    chrom_len: usize,
    pos: u64,
    offset: u64,
}

fn parse_vcf_collect_sorted(
    bytes: &[u8],
) -> (
    Vec<u64>,
    Vec<u64>,
    Vec<persist::ChrIndexEntry>,
    Vec<u8>,
    u32,
    Vec<persist::ChrIdMapEntry>,
) {
    let n_threads = rayon::current_num_threads().max(1);
    let target_chunk = 48usize * 1024 * 1024;
    let mut ranges = Vec::new();
    let mut start = 0usize;

    while start < bytes.len() {
        let mut end = start.saturating_add(target_chunk);
        if end >= bytes.len() {
            end = bytes.len();
        } else {
            if let Some(off) = memchr(b'\n', &bytes[end..]) {
                end += off + 1;
            } else {
                end = bytes.len();
            }
        }
        ranges.push((start, end));
        start = end;
    }

    while ranges.len() < n_threads && ranges.iter().any(|&(s, e)| e - s > target_chunk) {
        if let Some((idx, &(s, e))) = ranges.iter().enumerate().max_by_key(|(_, r)| r.1 - r.0) {
            let mid_guess = s + (e - s) / 2;
            if let Some(off) = memchr(b'\n', &bytes[mid_guess..e]) {
                let mid = mid_guess + off + 1;
                let right = (mid, e);
                ranges[idx].1 = mid;
                ranges.push(right);
            } else {
                break;
            }
        } else {
            break;
        }
    }

    #[derive(Default)]
    struct ChunkOut {
        recs: Vec<LineRec>,
        names: Vec<Vec<u8>>,
    }

    let chunk_out: Vec<ChunkOut> = ranges
        .into_par_iter()
        .map(|(s, e)| {
            let mut out = ChunkOut::default();
            let buf = &bytes[s..e];

            let mut line_start = 0usize;
            use std::collections::HashSet;
            let mut local_set: HashSet<Vec<u8>> = HashSet::with_capacity(64);

            for nl_rel in memchr_iter(b'\n', buf) {
                let line_end = nl_rel;
                let line = &buf[line_start..line_end];
                if let Some((chrom_len, pos)) = parse_line_chrom_pos_inline(line) {
                    let chrom = &line[..chrom_len];
                    if local_set.insert(chrom.to_vec()) {
                        out.names.push(chrom.to_vec());
                    }
                    out.recs.push(LineRec {
                        chrom_start: s + line_start,
                        chrom_len,
                        pos,
                        offset: (s + line_start) as u64,
                    });
                }
                line_start = nl_rel + 1;
            }

            if line_start < buf.len() {
                let line = &buf[line_start..];
                if let Some((chrom_len, pos)) = parse_line_chrom_pos_inline(line) {
                    let chrom = &line[..chrom_len];
                    if local_set.insert(chrom.to_vec()) {
                        out.names.push(chrom.to_vec());
                    }
                    out.recs.push(LineRec {
                        chrom_start: s + line_start,
                        chrom_len,
                        pos,
                        offset: (s + line_start) as u64,
                    });
                }
            }
            out
        })
        .collect();

    let mut all_names: Vec<Vec<u8>> = Vec::new();
    for c in &chunk_out {
        all_names.extend_from_slice(&c.names);
    }
    all_names.sort_unstable();
    all_names.dedup();

    let mut chr_index: Vec<persist::ChrIndexEntry> = Vec::with_capacity(all_names.len());
    let mut names_blob: Vec<u8> = Vec::new();
    use std::collections::HashMap;
    let mut name2id: HashMap<Vec<u8>, u32> = HashMap::with_capacity(all_names.len());
    let mut idmap: Vec<persist::ChrIdMapEntry> = Vec::with_capacity(all_names.len() + 1);
    idmap.push(persist::ChrIdMapEntry {
        name_off: 0,
        name_len: 0,
        _pad: 0,
    });
    for (i, name) in all_names.into_iter().enumerate() {
        let id = (i as u32) + 1;
        let off = names_blob.len() as u64;
        let len = name.len() as u32;
        names_blob.extend_from_slice(&name);
        chr_index.push(persist::ChrIndexEntry {
            name_off: off,
            name_len: len,
            id,
        });
        name2id.insert(name, id);
        idmap.push(persist::ChrIdMapEntry {
            name_off: off,
            name_len: len,
            _pad: 0,
        });
    }
    let max_id = (idmap.len() as u32).saturating_sub(1);

    let mut pairs: Vec<(u64, u64)> =
        Vec::with_capacity(chunk_out.iter().map(|c| c.recs.len()).sum());
    for c in chunk_out {
        for r in c.recs {
            let name = &bytes[r.chrom_start..r.chrom_start + r.chrom_len];
            let id = *name2id.get(name).unwrap() as u64;
            let key = (id << 32) | r.pos;
            pairs.push((key, r.offset));
        }
    }

    let t_sort = Instant::now();
    let was_sorted = is_already_sorted(&pairs);
    if was_sorted {
        eprintln!("[sort] skipped (already sorted)");
    } else {
        eprintln!("[sort] radix-sort (u64 LSD, 8 passes)...");
        radix_sort_pairs_u64(&mut pairs);
    }
    let sort_ms = t_sort.elapsed();
    eprintln!("[sort] time={}", ms(sort_ms));

    let n = pairs.len();
    let mut keys = Vec::with_capacity(n);
    let mut offs = Vec::with_capacity(n);
    for (k, o) in pairs {
        keys.push(k);
        offs.push(o);
    }

    (keys, offs, chr_index, names_blob, max_id, idmap)
}

#[inline]
fn is_already_sorted(pairs: &[(u64, u64)]) -> bool {
    if pairs.len() < 2 {
        return true;
    }
    let n = pairs.len();
    let step = (n / 1024).max(1);
    let mut last = pairs[0].0;
    let mut i = 1usize;
    while i < n {
        let k = pairs[i].0;
        if k < last {
            return false;
        }
        last = k;
        i += step;
    }
    let start = n.saturating_sub(4096);
    let mut last = pairs[start].0;
    for j in start + 1..n {
        let k = pairs[j].0;
        if k < last {
            return false;
        }
        last = k;
    }
    true
}

fn radix_sort_pairs_u64(pairs: &mut [(u64, u64)]) {
    let n = pairs.len();
    if n <= 1 {
        return;
    }

    let mut tmp: Vec<(u64, u64)> = vec![(0, 0); n];
    let mut read_from_pairs = true; // pass 0: read pairs, write tmp

    for shift in (0..8).map(|b| b * 8) {
        let mut count = [0usize; 256];
        if read_from_pairs {
            for &(k, _) in pairs.iter() {
                count[((k >> shift) & 0xFF) as usize] += 1;
            }
        } else {
            for &(k, _) in tmp.iter() {
                count[((k >> shift) & 0xFF) as usize] += 1;
            }
        }

        let mut pos = [0usize; 256];
        let mut sum = 0usize;
        for i in 0..256 {
            pos[i] = sum;
            sum += count[i];
        }

        if read_from_pairs {
            for &(k, v) in pairs.iter() {
                let b = ((k >> shift) & 0xFF) as usize;
                let p = pos[b];
                tmp[p] = (k, v);
                pos[b] = p + 1;
            }
        } else {
            for &(k, v) in tmp.iter() {
                let b = ((k >> shift) & 0xFF) as usize;
                let p = pos[b];
                pairs[p] = (k, v);
                pos[b] = p + 1;
            }
        }

        read_from_pairs = !read_from_pairs;
    }
}

// ---------- PGM segments (f32 slope/intercept в persist::Segment) ----------
fn build_segments_from_keys(keys: &[u64], epsilon: usize) -> (Vec<persist::Segment>, Vec<u64>) {
    let n = keys.len();
    if n == 0 {
        return (Vec::new(), Vec::new());
    }
    let mut segs: Vec<persist::Segment> = Vec::new();
    let mut start = 0usize;

    while start < n {
        let end = (start + epsilon.saturating_mul(2)).min(n) - 1;
        let x0 = keys[start] as f64;
        let y0 = start as f64;
        let x1 = keys[end] as f64;
        let y1 = end as f64;

        let slope_f64 = if (x1 - x0).abs() > f64::EPSILON {
            (y1 - y0) / (x1 - x0)
        } else {
            0.0
        };
        let intercept_f64 = y0 - slope_f64 * x0;

        let seg = persist::Segment {
            key_lo: keys[start],
            key_hi: keys[end],
            slope: slope_f64 as f32,
            intercept: intercept_f64 as f32,
            base_rank: start as u64,
        };
        segs.push(seg);

        start = end + 1;
    }
    (segs, Vec::new())
}

// ---------- Varint (LEB128u) encode/decode ----------
#[inline]
fn put_varint_u64(mut v: u64, out: &mut Vec<u8>) {
    while v >= 0x80 {
        out.push(((v as u8) & 0x7F) | 0x80);
        v >>= 7;
    }
    out.push(v as u8);
}
#[inline]
fn get_varint_u64(mut p: usize, buf: &[u8]) -> (u64, usize) {
    let mut shift = 0;
    let mut x: u64 = 0;
    loop {
        let b = buf[p];
        p += 1;
        x |= ((b & 0x7F) as u64) << shift;
        if (b & 0x80) == 0 {
            break;
        }
        shift += 7;
    }
    (x, p)
}

fn build_varint_offsets(offsets: &[u64]) -> (Vec<u8>, Vec<persist::OffsCheckpoint>) {
    let n = offsets.len();
    if n == 0 {
        return (Vec::new(), Vec::new());
    }

    let mut blob = Vec::with_capacity(n * 3);
    let mut ckpts: Vec<persist::OffsCheckpoint> = Vec::new();

    ckpts.push(persist::OffsCheckpoint {
        index: 0,
        abs_offset: offsets[0],
        blob_pos: 0,
    });

    let mut prev = offsets[0];
    for i in 1..n {
        if i % persist::CKPT_STRIDE == 0 {
            ckpts.push(persist::OffsCheckpoint {
                index: i as u64,
                abs_offset: prev,
                blob_pos: blob.len() as u64,
            });
        }
        let delta = offsets[i] - prev;
        put_varint_u64(delta, &mut blob);
        prev = offsets[i];
    }

    (blob, ckpts)
}

fn decode_offsets_window(idx: &persist::PgmiIndex, lo: usize, hi: usize) -> Vec<u64> {
    if lo >= hi {
        return Vec::new();
    }
    let ckpts = idx.ckpts();
    let blob = idx.offsets_comp();
    if ckpts.is_empty() {
        return Vec::new();
    }

    let ci = match ckpts.binary_search_by(|c| (c.index as usize).cmp(&lo)) {
        Ok(i) => i,
        Err(0) => 0,
        Err(i) => i - 1,
    };
    let mut i = ckpts[ci].index as usize;
    let mut pos = ckpts[ci].blob_pos as usize;
    let mut off = ckpts[ci].abs_offset;

    while i < lo {
        let (d, p2) = get_varint_u64(pos, blob);
        pos = p2;
        off = off.wrapping_add(d);
        i += 1;
    }

    let need = hi - lo;
    let mut out = Vec::with_capacity(need);
    out.push(off);
    for _ in 1..need {
        if pos >= blob.len() {
            break;
        }
        let (d, p2) = get_varint_u64(pos, blob);
        pos = p2;
        off = off.wrapping_add(d);
        out.push(off);
    }
    out
}

#[derive(Clone)]
struct SplitMix64 {
    state: u64,
}
impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }
    #[inline]
    fn next_u64(&mut self) -> u64 {
        let mut z = {
            self.state = self.state.wrapping_add(0x9E3779B97F4A7C15);
            self.state
        };
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    #[inline]
    fn next_f64(&mut self) -> f64 {
        let u = self.next_u64() >> 11; // 53 bits
        (u as f64) * (1.0 / ((1u64 << 53) as f64))
    }
    #[inline]
    fn range_u64(&mut self, lo: u64, hi_exclusive: u64) -> u64 {
        if hi_exclusive <= lo {
            return lo;
        }
        lo + (self.next_u64() % (hi_exclusive - lo))
    }
    #[inline]
    fn range_usize(&mut self, hi_exclusive: usize) -> usize {
        if hi_exclusive == 0 {
            return 0;
        }
        (self.next_u64() as usize) % hi_exclusive
    }
}

// ---------- CLI ----------
#[derive(Parser, Debug)]
#[command(name = "pgm-hts")]
#[command(about = "PGM learned index tools (VCF/PGMI) with mmap I/O + chrmap + idmap + compressed offsets + batch bench", long_about = None)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand, Debug)]
enum Cmd {
    IndexVcf {
        #[arg(value_name = "FILE.VCF")]
        file: PathBuf,
        #[arg(long, default_value_t = 64)]
        epsilon: usize,
        #[arg(long, value_name = "OUT.PGMI")]
        out: Option<PathBuf>,
        #[arg(long = "mmap-save", action = ArgAction::SetTrue)]
        mmap_save: bool,
    },

    StatPgmi {
        #[arg(value_name = "FILE.PGMI")]
        file: PathBuf,
        #[arg(long, action = ArgAction::SetTrue)]
        chrmap: bool,
    },

    QueryPgmi {
        #[arg(value_name = "FILE.PGMI")]
        file: PathBuf,

        #[arg(long, value_name = "U64")]
        key: Option<u64>,

        #[arg(long, value_name = "CHR")]
        chr: Option<String>,
        #[arg(long, value_name = "START")]
        start: Option<u64>,
        #[arg(long, value_name = "END")]
        end: Option<u64>,

        #[arg(long, value_name = "FILE.VCF")]
        vcf: Option<PathBuf>,

        #[arg(long, default_value_t = 50)]
        max_print: usize,
    },

    BenchBatch {
        #[arg(value_name = "FILE.PGMI")]
        pgmi: PathBuf,

        #[arg(long, default_value = "random")]
        mode: String,

        #[arg(long, default_value_t = 100000)]
        n: usize,

        #[arg(long, default_value_t = 500)]
        width: u64,

        #[arg(long, default_value_t = 50)]
        clusters: usize,

        #[arg(long, action = ArgAction::SetTrue)]
        with_fetch: bool,

        #[arg(long, value_name = "FILE.VCF")]
        vcf: Option<PathBuf>,

        #[arg(long, default_value_t = 0)]
        threads: usize,
    },
}

fn main() -> Result<(), String> {
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::IndexVcf {
            file,
            epsilon,
            out,
            mmap_save,
        } => cmd_index_vcf(&file, epsilon, out.as_deref(), mmap_save).map_err(errs),
        Cmd::StatPgmi { file, chrmap } => cmd_stat_pgmi(&file, chrmap).map_err(errs),
        Cmd::QueryPgmi {
            file,
            key,
            chr,
            start,
            end,
            vcf,
            max_print,
        } => cmd_query_pgmi(&file, key, chr, start, end, vcf.as_deref(), max_print).map_err(errs),
        Cmd::BenchBatch {
            pgmi,
            mode,
            n,
            width,
            clusters,
            with_fetch,
            vcf,
            threads,
        } => cmd_bench_batch(
            &pgmi,
            &mode,
            n,
            width,
            clusters,
            with_fetch,
            vcf.as_deref(),
            threads,
        )
        .map_err(errs),
    }
}

fn errs(e: impl std::fmt::Debug) -> String {
    format!("Error: {:?}", e)
}

// ---------- commands ----------

fn cmd_index_vcf(
    file: &Path,
    epsilon: usize,
    out: Option<&Path>,
    mmap_save: bool,
) -> io::Result<()> {
    let t_all = Instant::now();

    let t0 = Instant::now();
    let bytes = mmap_read_file(file)?;
    #[cfg(feature = "mmap")]
    eprintln!("[vcf] mmap-read {}", file.display());
    let read_ms = t0.elapsed();

    let t1 = Instant::now();
    let (keys, offsets, chr_index, chr_names, max_id, idmap) = {
        #[cfg(feature = "mmap")]
        {
            parse_vcf_collect_sorted(&bytes[..])
        }
        #[cfg(not(feature = "mmap"))]
        {
            parse_vcf_collect_sorted(&bytes)
        }
    };
    let parse_ms = t1.elapsed();
    eprintln!(
        "[vcf] parsed {} records (chr={}, max_id={})",
        keys.len(),
        chr_index.len(),
        max_id
    );

    let t2 = Instant::now();
    let (segments, anchors) = build_segments_from_keys(&keys, epsilon);
    let build_ms = t2.elapsed();
    eprintln!("[pgm] segments: {}", segments.len());

    let t3 = Instant::now();
    let (offsets_comp, ckpts) = build_varint_offsets(&offsets);
    let comp_ms = t3.elapsed();
    let comp_ratio = if offsets.is_empty() {
        1.0
    } else {
        (offsets_comp.len() as f64) / ((offsets.len() as f64) * 8.0)
    };
    eprintln!(
        "[comp] offsets: {} → {} bytes (ratio {:.3}) ckpts={}",
        offsets.len() * 8,
        offsets_comp.len(),
        comp_ratio,
        ckpts.len()
    );

    let out_path = out.map(|p| p.to_path_buf()).unwrap_or_else(|| {
        let mut p = PathBuf::from(file);
        p.set_extension("vcf.pgmi");
        p
    });

    let t4 = Instant::now();
    #[cfg(feature = "mmap")]
    if mmap_save {
        persist::save_pgmi_mmap(
            &out_path,
            epsilon as u32,
            &segments,
            &anchors,
            &chr_index,
            &chr_names,
            &idmap,
            &offsets_comp,
            &ckpts,
        )?;
        eprintln!("[pgmi] saved via MmapMut in {:?}", t4.elapsed());
    } else {
        persist::save_pgmi_v24(
            &out_path,
            epsilon as u32,
            &segments,
            &anchors,
            &chr_index,
            &chr_names,
            &idmap,
            &offsets_comp,
            &ckpts,
        )?;
        eprintln!("[pgmi] saved (flat) in {:?}", t4.elapsed());
    }
    #[cfg(not(feature = "mmap"))]
    {
        persist::save_pgmi_v24(
            &out_path,
            epsilon as u32,
            &segments,
            &anchors,
            &chr_index,
            &chr_names,
            &idmap,
            &offsets_comp,
            &ckpts,
        )?;
        eprintln!("[pgmi] saved (flat) in {:?}", t4.elapsed());
    }

    eprintln!(
        "[time] read={} parse={} build={} compress={} save={} total={}",
        ms(read_ms),
        ms(parse_ms),
        ms(build_ms),
        ms(comp_ms),
        ms(t4.elapsed()),
        ms(t_all.elapsed())
    );
    Ok(())
}

fn cmd_stat_pgmi(path: &Path, dump_chrmap: bool) -> io::Result<()> {
    let t0 = Instant::now();
    let idx = persist::load_pgmi(path)?;
    let load_ms = t0.elapsed();

    eprintln!(
        "[pgmi] loaded via {} in {}",
        if cfg!(feature = "mmap") {
            "mmap"
        } else {
            "owned"
        },
        ms(load_ms)
    );
    eprintln!("epsilon       = {}", idx.epsilon());
    eprintln!("segments      = {}", idx.segments().len());
    eprintln!("anchors       = {}", idx.anchors().len());
    eprintln!("chr entries   = {}", idx.chr_index().len());
    eprintln!("idmap entries = {}", idx.idmap().len());
    eprintln!("ckpt stride   = {}", idx.header().ckpt_stride);
    eprintln!("ckpts         = {}", idx.ckpts().len());
    eprintln!("offs blob     = {} bytes", idx.offsets_comp().len());

    if dump_chrmap {
        let names = idx.chr_names();
        for (i, e) in idx.idmap().iter().enumerate().skip(1) {
            if e.name_len == 0 {
                continue;
            }
            let name = &names[e.name_off as usize..(e.name_off as usize + e.name_len as usize)];
            eprintln!("  id={:<3} name={}", i, String::from_utf8_lossy(name));
        }
    }
    Ok(())
}

fn chr_id_by_name(idx: &persist::PgmiIndex, name: &str) -> Option<u32> {
    let blob = idx.chr_names();
    for e in idx.chr_index() {
        let s = &blob[e.name_off as usize..(e.name_off as usize + e.name_len as usize)];
        if s == name.as_bytes() {
            return Some(e.id);
        }
    }
    None
}

fn cmd_query_pgmi(
    pgmi_path: &Path,
    key: Option<u64>,
    chr: Option<String>,
    start: Option<u64>,
    end: Option<u64>,
    vcf_path_opt: Option<&Path>,
    max_print: usize,
) -> io::Result<()> {
    let t_all = Instant::now();
    let t0 = Instant::now();
    let idx = persist::load_pgmi(pgmi_path)?;
    let load_ms = t0.elapsed();

    let eps = idx.epsilon() as usize;
    let segs = idx.segments();
    if segs.is_empty() {
        eprintln!("Empty index");
        return Ok(());
    }

    let t1 = Instant::now();
    let (k1, k2, chr_name, s_pos, e_pos) = if let Some(k) = key {
        (k, k, None, 0, 0)
    } else {
        let (chr, s, e) = match (chr, start, end) {
            (Some(c), Some(s), Some(e)) if s <= e => (c, s, e),
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "specify either --key <U64> OR --chr <NAME> --start <POS> --end <POS>",
                ))
            }
        };
        let id = chr_id_by_name(&idx, &chr).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("Chromosome not found: {}", chr),
            )
        })?;
        let k1 = ((id as u64) << 32) | s;
        let k2 = ((id as u64) << 32) | e;
        (k1, k2, Some(chr), s, e)
    };
    let chr_name_str: Option<&str> = chr_name.as_deref();
    let translate_ms = t1.elapsed();

    // predict+refine windows
    let t2 = Instant::now();
    let (l1, r1) = pgmf_predict_range(segs, eps, k1);
    let (l2, r2) = pgmf_predict_range(segs, eps, k2);
    let lo = cmp::min(l1, l2);
    let hi = cmp::max(r1, r2);
    let refine_ms = t2.elapsed();

    eprintln!(
        "[query] {} → keys [{:#x}, {:#x}]",
        if let Some(cn) = chr_name_str {
            format!("{}:{}-{}", cn, s_pos, e_pos)
        } else {
            format!("key={:#x}", k1)
        },
        k1,
        k2
    );
    eprintln!("[pgmf] union refine window = [{}, {})", lo, hi);

    let t3 = Instant::now();
    let mut printed = 0usize;
    if let Some(q_chr) = chr_name_str {
        let vcf_path = match vcf_path_opt {
            Some(p) => p.to_path_buf(),
            None => infer_vcf_from_pgmi(pgmi_path),
        };
        let vcf_mm = mmap_read_file(&vcf_path)?;
        #[cfg(feature = "mmap")]
        let vbytes: &[u8] = &vcf_mm[..];
        #[cfg(not(feature = "mmap"))]
        let vbytes: &[u8] = &vcf_mm;

        let offs = decode_offsets_window(&idx, lo, hi);
        eprintln!(
            "[fetch] source: {} (decode {} offsets)",
            vcf_path.display(),
            offs.len()
        );

        let (q_s, q_e) = (s_pos, e_pos);
        for off in offs {
            let o = off as usize;
            if o >= vbytes.len() {
                continue;
            }
            // [o .. next '\n']
            let mut j = o;
            if let Some(nl) = memchr(b'\n', &vbytes[o..]) {
                j = o + nl;
            } else {
                while j < vbytes.len() && vbytes[j] != b'\n' {
                    j += 1;
                }
            }
            let line = &vbytes[o..j];

            if let Some(t1) = memchr(b'\t', line) {
                if let Some(t2_rel) = memchr(b'\t', &line[t1 + 1..]) {
                    let chrom = &line[..t1];
                    let pos_s = &line[t1 + 1..t1 + 1 + t2_rel];
                    let pos = parse_u64_fast(pos_s);
                    if chrom == q_chr.as_bytes() && pos >= q_s && pos <= q_e {
                        println!("{}", String::from_utf8_lossy(line));
                        printed += 1;
                        if printed >= max_print {
                            eprintln!("[fetch] max_print reached ({} lines)", max_print);
                            break;
                        }
                    }
                }
            }
        }
    }
    let fetch_ms = t3.elapsed();

    eprintln!(
        "[time] load={} translate={} refine={} fetch={} total={}",
        ms(load_ms),
        ms(translate_ms),
        ms(refine_ms),
        ms(fetch_ms),
        ms(t_all.elapsed())
    );
    if printed == 0 && chr_name_str.is_some() {
        eprintln!("[fetch] no records printed (maybe region empty?)");
    }
    Ok(())
}

// ---------- bench-batch ----------

fn cmd_bench_batch(
    pgmi_path: &Path,
    mode: &str,
    n: usize,
    width: u64,
    clusters: usize,
    with_fetch: bool,
    vcf_path_opt: Option<&Path>,
    threads: usize,
) -> io::Result<()> {
    if threads > 0 {
        std::env::set_var("RAYON_NUM_THREADS", threads.to_string());
    }

    let t_load = Instant::now();
    let idx = persist::load_pgmi(pgmi_path)?;
    let load_ms = t_load.elapsed();

    eprintln!("[bench] loaded index in {}", ms(load_ms));
    let segs = idx.segments();
    if segs.is_empty() {
        return Ok(());
    }
    let eps = idx.epsilon() as usize;

    let mut max_pos_by_chr: Vec<u64> = vec![0; idx.idmap().len().max(1)];
    for s in segs {
        let id = (s.key_hi >> 32) as usize;
        let pos = (s.key_hi & 0xFFFF_FFFF) as u64;
        if id < max_pos_by_chr.len() {
            if pos > max_pos_by_chr[id] {
                max_pos_by_chr[id] = pos;
            }
        }
    }
    let max_id: u32 = (max_pos_by_chr.len() as u32).saturating_sub(1);

    #[cfg(feature = "mmap")]
    let (vbytes_opt, vcf_path_print): (Option<Vec<u8>>, String) = if with_fetch {
        let vcf_path = vcf_path_opt
            .map(|p| p.to_path_buf())
            .unwrap_or_else(|| infer_vcf_from_pgmi(pgmi_path));
        let mm = mmap_read_file(&vcf_path)?;
        (Some(mm[..].to_vec()), vcf_path.display().to_string())
    } else {
        (None, String::new())
    };
    #[cfg(not(feature = "mmap"))]
    let (vbytes_opt, vcf_path_print): (Option<Vec<u8>>, String) = if with_fetch {
        let vcf_path = vcf_path_opt
            .map(|p| p.to_path_buf())
            .unwrap_or_else(|| infer_vcf_from_pgmi(pgmi_path));
        let mm = mmap_read_file(&vcf_path)?;
        (Some(mm), vcf_path.display().to_string())
    } else {
        (None, String::new())
    };

    if with_fetch {
        eprintln!("[bench] with-fetch from {}", vcf_path_print);
    }

    let t_gen = Instant::now();
    let queries = gen_queries(mode, n, width, &segs, max_id, &max_pos_by_chr);
    let gen_ms = t_gen.elapsed();
    eprintln!(
        "[gen] mode={} n={} width={} gen_time={}",
        mode,
        n,
        width,
        ms(gen_ms)
    );

    let mut t_predict = Duration::ZERO;
    let mut t_decode = Duration::ZERO;
    let mut t_fetch = Duration::ZERO;

    let mut decoded_per_q: Vec<usize> = Vec::with_capacity(n);
    let mut refine_span_per_q: Vec<usize> = Vec::with_capacity(n);

    let t_all = Instant::now();
    for &(k1, k2) in &queries {
        let t1 = Instant::now();
        let (l1, r1) = pgmf_predict_range(segs, eps, k1);
        let (l2, r2) = pgmf_predict_range(segs, eps, k2);
        let lo = cmp::min(l1, l2);
        let hi = cmp::max(r1, r2);
        t_predict += t1.elapsed();

        let t2 = Instant::now();
        let offs = decode_offsets_window(&idx, lo, hi);
        t_decode += t2.elapsed();
        decoded_per_q.push(offs.len());
        refine_span_per_q.push(hi.saturating_sub(lo));

        if with_fetch {
            let t3 = Instant::now();
            let vbytes = vbytes_opt.as_ref().unwrap();

            let id = (k1 >> 32) as u32;
            let s_pos = (k1 & 0xFFFF_FFFF) as u64;
            let e_pos = (k2 & 0xFFFF_FFFF) as u64;
            let chr_blob = idx.chr_names();

            let chr_name = {
                let e = &idx.idmap()[id as usize];
                &chr_blob[e.name_off as usize..(e.name_off as usize + e.name_len as usize)]
            };
            for off in offs {
                let o = off as usize;
                if o >= vbytes.len() {
                    continue;
                }
                let mut j = o;
                if let Some(nl) = memchr(b'\n', &vbytes[o..]) {
                    j = o + nl;
                } else {
                    while j < vbytes.len() && vbytes[j] != b'\n' {
                        j += 1;
                    }
                }
                let line = &vbytes[o..j];

                if let Some(t1) = memchr(b'\t', line) {
                    if let Some(t2_rel) = memchr(b'\t', &line[t1 + 1..]) {
                        let chrom = &line[..t1];
                        let pos_s = &line[t1 + 1..t1 + 1 + t2_rel];
                        let pos = parse_u64_fast(pos_s);
                        if chrom == chr_name && pos >= s_pos && pos <= e_pos {
                            // no-op
                        }
                    }
                }
            }
            t_fetch += t3.elapsed();
        }
    }
    let total_ms = t_all.elapsed();

    let p50_span = percentile_usize(&mut refine_span_per_q.clone(), 50);
    let p95_span = percentile_usize(&mut refine_span_per_q, 95);
    let p50_dec = percentile_usize(&mut decoded_per_q.clone(), 50);
    let p95_dec = percentile_usize(&mut decoded_per_q, 95);

    let qps = (queries.len() as f64) / (total_ms.as_secs_f64());
    eprintln!(
        "[bench] predict+refine={} decode={} fetch={} total={}",
        ms(t_predict),
        ms(t_decode),
        ms(t_fetch),
        ms(total_ms)
    );
    eprintln!(
        "[bench] refine window size: p50={} p95={}",
        p50_span, p95_span
    );
    eprintln!(
        "[bench] decoded offsets per query: p50={} p95={}",
        p50_dec, p95_dec
    );
    eprintln!("[bench] throughput: {:.0} queries/sec", qps);
    Ok(())
}

fn percentile_usize(v: &mut [usize], p: u32) -> usize {
    if v.is_empty() {
        return 0;
    }
    let idx = ((v.len() as u128) * (p as u128) / 100u128).min((v.len() - 1) as u128) as usize;
    v.select_nth_unstable(idx).1.clone()
}

fn gen_queries(
    mode: &str,
    n: usize,
    width: u64,
    segs: &[persist::Segment],
    max_id: u32,
    max_pos_by_chr: &[u64],
) -> Vec<(u64, u64)> {
    match mode {
        "browser" => gen_browser(n, width, max_id, max_pos_by_chr),
        "clustered" => gen_clustered(n, width, segs, max_id),
        _ => gen_random(n, width, segs),
    }
}

fn gen_random(n: usize, width: u64, segs: &[persist::Segment]) -> Vec<(u64, u64)> {
    let mut rng = SplitMix64::new(0x1234_5678_9abc_def0);
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let si = rng.range_usize(segs.len());
        let s = &segs[si];
        let id = (s.key_lo >> 32) as u64;
        let lo = s.key_lo & 0xFFFF_FFFF;
        let hi = s.key_hi & 0xFFFF_FFFF;
        let pos = if hi > lo {
            rng.range_u64(lo, hi + 1)
        } else {
            lo
        };
        let k1 = (id << 32) | pos;
        let k2 = (id << 32) | pos.saturating_add(width);
        out.push((k1, k2));
    }
    out
}

fn gen_browser(n: usize, width: u64, max_id: u32, max_pos_by_chr: &[u64]) -> Vec<(u64, u64)> {
    let mut out = Vec::with_capacity(n);
    if max_id == 0 {
        return out;
    }

    let per_chr = (n as u64 + max_id as u64) / (max_id as u64);
    for id in 1..=max_id {
        if out.len() >= n {
            break;
        }
        let max_pos = max_pos_by_chr.get(id as usize).copied().unwrap_or(0);
        if max_pos == 0 {
            continue;
        }

        let mut pos: u64 = 1;
        let mut cnt = 0usize;
        while pos < max_pos && cnt < per_chr as usize && out.len() < n {
            let k1 = ((id as u64) << 32) | pos;
            let k2 = ((id as u64) << 32) | pos.saturating_add(width);
            out.push((k1, k2));
            pos = pos.saturating_add(width);
            cnt += 1;
        }
    }
    out
}

fn gen_clustered(n: usize, width: u64, segs: &[persist::Segment], max_id: u32) -> Vec<(u64, u64)> {
    let mut rng = SplitMix64::new(0x9e37_79b9_7f4a_7c15);
    let k = (max_id.max(1) as usize).min(64).max(1);
    let k = k.min(n.max(1));

    let mut centers: Vec<(u32, u64)> = Vec::with_capacity(k);
    for _ in 0..k {
        let si = rng.range_usize(segs.len());
        let s = &segs[si];
        let id = (s.key_lo >> 32) as u32;
        let lo = (s.key_lo & 0xFFFF_FFFF) as u64;
        let hi = (s.key_hi & 0xFFFF_FFFF) as u64;
        let pos = if hi > lo {
            rng.range_u64(lo, hi + 1)
        } else {
            lo
        };
        centers.push((id, pos));
    }

    let lambda = 1.0 / (width as f64 * 4.0);
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let (id, cpos) = centers[rng.range_usize(centers.len())];
        let u = rng.next_f64().max(1e-12);
        let offs = (-u.ln() / lambda) as u64;
        let dir = (rng.next_u64() & 1) == 0;
        let pos = if dir {
            cpos.saturating_add(offs)
        } else {
            cpos.saturating_sub(offs)
        };
        let k1 = ((id as u64) << 32) | pos;
        let k2 = ((id as u64) << 32) | pos.saturating_add(width);
        out.push((k1, k2));
    }
    out
}

fn pgmf_predict_range(segs: &[persist::Segment], eps: usize, key: u64) -> (usize, usize) {
    // locate segment
    let si = match segs.binary_search_by(|s| {
        if key < s.key_lo {
            std::cmp::Ordering::Greater
        } else if key > s.key_hi {
            std::cmp::Ordering::Less
        } else {
            std::cmp::Ordering::Equal
        }
    }) {
        Ok(i) => i,
        Err(_) => {
            let mut lo = 0usize;
            let mut hi = segs.len() - 1;
            while lo < hi {
                let mid = (lo + hi) >> 1;
                if key <= segs[mid].key_hi {
                    hi = mid;
                } else {
                    lo = mid + 1;
                }
            }
            lo
        }
    };
    let s = &segs[si];
    let kf = key as f32;
    let mut pred = (s.slope.mul_add(kf, s.intercept).floor() as i64); // f32 FMA
    if pred < 0 {
        pred = 0;
    }
    let base = s.base_rank as i64;
    let lo = (pred - eps as i64).max(base) as usize;
    let hi = (pred + eps as i64 + 1) as usize;

    //     eprintln!(
    //         "[pgmf] seg={} key={:#x} in [{}, {}], base_rank={}, slope={:.6e}, intercept={:.6e} → pred≈{} refine=[{}, {})",
    //         si, key, s.key_lo, s.key_hi, s.base_rank, s.slope, s.intercept, pred, lo, hi
    //     );
    (lo, hi)
}

fn infer_vcf_from_pgmi(pgmi: &Path) -> PathBuf {
    let s = pgmi.to_string_lossy();
    if s.ends_with(".vcf.pgmi") {
        let mut t = s.to_string();
        t.truncate(t.len() - ".pgmi".len());
        PathBuf::from(t)
    } else if s.ends_with(".pgmi") {
        let mut t = s.to_string();
        t.truncate(t.len() - ".pgmi".len());
        PathBuf::from(t)
    } else {
        pgmi.with_extension("")
    }
}

fn ms(d: std::time::Duration) -> String {
    let us = d.as_micros();
    if us < 1_000 {
        return format!("{}µs", us);
    }
    let ms = us as f64 / 1000.0;
    if ms < 1000.0 {
        return format!("{:.2}ms", ms);
    }
    format!("{:.3}s", ms / 1000.0)
}
