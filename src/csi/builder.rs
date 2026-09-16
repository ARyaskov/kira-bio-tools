//! Index construction (htslib `hts_idx_push` / `hts_idx_finish` semantics) and
//! the VCF/BCF drivers that feed it.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, Read};
use std::path::Path;

use anyhow::{Context, Result, bail};

use super::binning::{BinEntry, BinIndex, IndexKind, RefIndex, RefMeta, TabixHeader};
use super::utils::{bin_bot, bin_first, bin_parent, depth_for, max_pos, n_bins, reg2bin, tabix_depth_for};
use crate::bcf::{BCF_MAGIC, BcfReader};
use crate::bcf::record::{record_end0, record_meta};
use crate::bgzf::is_bgzf;
use crate::vcf::header::HeaderInfo;

/// Bins whose chunks span less than this many compressed bytes are folded
/// into their parent (`HTS_MIN_MARKER_DIST`).
const MIN_MARKER_DIST: u64 = 0x10000;

struct RefBuild {
    bins: BTreeMap<u32, Vec<(u64, u64)>>,
    lidx: Vec<u64>,
    meta: RefMeta,
    has_records: bool,
}

impl Default for RefBuild {
    fn default() -> Self {
        Self {
            bins: BTreeMap::new(),
            lidx: Vec::new(),
            meta: RefMeta { beg_off: u64::MAX, end_off: 0, n_mapped: 0, n_unmapped: 0 },
            has_records: false,
        }
    }
}

struct Pending {
    tid: usize,
    bin: u32,
    start: u64,
    end: u64,
}

pub struct IndexBuilder {
    kind: IndexKind,
    min_shift: u8,
    depth: u8,
    refs: Vec<RefBuild>,
    pending: Option<Pending>,
    last_tid: Option<usize>,
    last_beg: u64,
}

impl IndexBuilder {
    pub fn new(kind: IndexKind, min_shift: u8, depth: u8) -> Self {
        Self { kind, min_shift, depth, refs: Vec::new(), pending: None, last_tid: None, last_beg: 0 }
    }

    /// Register a record on contig `tid` covering `[beg0, end0)` stored at
    /// virtual positions `[start, end)`. Records must arrive in file order.
    pub fn add(&mut self, tid: usize, beg0: u64, mut end0: u64, start: u64, end: u64) -> Result<()> {
        if end0 <= beg0 {
            end0 = beg0 + 1;
        }
        if end0 > max_pos(self.min_shift, self.depth) {
            bail!(
                "position {} exceeds the index range ({} levels of 2^{} bp)",
                end0,
                self.depth,
                self.min_shift
            );
        }
        if let Some(lt) = self.last_tid {
            if lt == tid && beg0 < self.last_beg {
                bail!("input is not sorted by position (contig #{tid}: {} after {})", beg0 + 1, self.last_beg + 1);
            }
            if lt > tid {
                bail!("input is not sorted: contig #{tid} appears after contig #{lt}");
            }
        }
        self.last_tid = Some(tid);
        self.last_beg = beg0;

        if self.refs.len() <= tid {
            self.refs.resize_with(tid + 1, RefBuild::default);
        }
        let bin = reg2bin(beg0, end0, self.min_shift, self.depth);
        let r = &mut self.refs[tid];
        let bot = (beg0 >> self.min_shift) as usize;
        let top = ((end0 - 1) >> self.min_shift) as usize;
        if r.lidx.len() <= top {
            r.lidx.resize(top + 1, u64::MAX);
        }
        for i in bot..=top {
            if r.lidx[i] == u64::MAX || start < r.lidx[i] {
                r.lidx[i] = start;
            }
        }
        r.has_records = true;
        r.meta.n_mapped += 1;
        r.meta.beg_off = r.meta.beg_off.min(start);
        r.meta.end_off = r.meta.end_off.max(end);

        match &mut self.pending {
            Some(p) if p.tid == tid && p.bin == bin => {
                p.end = end;
            }
            _ => {
                self.flush_pending();
                self.pending = Some(Pending { tid, bin, start, end });
            }
        }
        Ok(())
    }

    fn flush_pending(&mut self) {
        let Some(p) = self.pending.take() else { return };
        let list = self.refs[p.tid].bins.entry(p.bin).or_default();
        if let Some(last) = list.last_mut()
            && last.1 >> 16 == p.start >> 16 {
                last.1 = p.end;
                return;
            }
        list.push((p.start, p.end));
    }

    /// Finish the index. `n_refs` pads the reference list (contigs without
    /// records still get an empty entry, as htslib does when they are named).
    pub fn finish(mut self, header: Option<TabixHeader>, n_refs: usize) -> BinIndex {
        self.flush_pending();
        if self.refs.len() < n_refs {
            self.refs.resize_with(n_refs, RefBuild::default);
        }
        let depth = self.depth;
        let kind = self.kind;
        let n_regular = n_bins(depth);
        let mut refs = Vec::with_capacity(self.refs.len());
        for mut r in self.refs {
            // update_loff
            let offset0 = if r.has_records { r.meta.beg_off } else { 0 };
            let mut l = 0usize;
            while l < r.lidx.len() && r.lidx[l] == u64::MAX {
                r.lidx[l] = offset0;
                l += 1;
            }
            if l == r.lidx.len() {
                r.lidx.clear();
            }
            for i in 1..r.lidx.len() {
                if r.lidx[i] == u64::MAX {
                    r.lidx[i] = r.lidx[i - 1];
                }
            }

            // compress_binning: fold small bins into parents, deepest level first.
            for level in (1..=depth).rev() {
                let lo = bin_first(level);
                let hi = bin_first(level + 1);
                let ids: Vec<u32> = r.bins.range(lo..hi).map(|(&id, _)| id).collect();
                for id in ids {
                    let list = r.bins.get_mut(&id).unwrap();
                    if list.len() > 1 {
                        list.sort_unstable();
                    }
                    let first = list.first().map(|c| c.0 >> 16).unwrap_or(0);
                    let last = list.last().map(|c| c.1 >> 16).unwrap_or(0);
                    // htslib folds a small bin only into a parent that
                    // already has records of its own.
                    let parent = bin_parent(id);
                    if last.saturating_sub(first) < MIN_MARKER_DIST && r.bins.contains_key(&parent) {
                        let moved = r.bins.remove(&id).unwrap();
                        r.bins.get_mut(&parent).unwrap().extend(moved);
                    }
                }
            }
            if let Some(list) = r.bins.get_mut(&0) {
                list.sort_unstable();
            }
            // merge chunks that touch the same BGZF block
            let mut bins = BTreeMap::new();
            for (id, mut list) in r.bins {
                if id >= n_regular {
                    continue;
                }
                list.sort_unstable();
                let mut merged: Vec<(u64, u64)> = Vec::with_capacity(list.len());
                for (s, e) in list {
                    if let Some(m) = merged.last_mut()
                        && m.1 >> 16 >= s >> 16 {
                            if e > m.1 {
                                m.1 = e;
                            }
                            continue;
                        }
                    merged.push((s, e));
                }
                let bot = bin_bot(id, depth) as usize;
                let loffset = if bot < r.lidx.len() { r.lidx[bot] } else { 0 };
                bins.insert(id, BinEntry { loffset, chunks: merged });
            }
            let meta = if r.has_records { Some(r.meta) } else { None };
            let linear = if kind == IndexKind::Tbi { r.lidx } else { Vec::new() };
            refs.push(RefIndex { bins, linear, meta });
        }
        BinIndex {
            kind,
            min_shift: self.min_shift,
            depth,
            header,
            refs,
            n_no_coor: Some(0),
        }
    }
}

/// `(beg0, end0)` of a VCF text record for indexing: `[POS-1, POS-1+len(REF))`,
/// or `[POS-1, INFO/END)` when END is present (`tbx_parse1`, VCF preset).
pub fn vcf_line_interval(line: &str) -> Option<(&str, u64, u64)> {
    let mut it = line.split('\t');
    let chrom = it.next()?;
    let pos: u64 = it.next()?.parse().ok()?;
    let _id = it.next()?;
    let ref_len = it.next()?.len() as u64;
    let _alt = it.next();
    let _qual = it.next();
    let _filter = it.next();
    let beg0 = pos.saturating_sub(1);
    let mut end0 = beg0 + ref_len.max(1);
    if let Some(info) = it.next() {
        for kv in info.split(';') {
            if let Some(v) = kv.strip_prefix("END=") {
                if let Ok(e) = v.parse::<u64>() {
                    end0 = e;
                }
                break;
            }
        }
    }
    Some((chrom, beg0, end0))
}

pub const TBX_GENERIC: i32 = 0;
pub const TBX_SAM: i32 = 1;
pub const TBX_VCF: i32 = 2;
pub const TBX_UCSC: i32 = 0x10000;

/// htslib `tbx_conf_t`: which columns carry the sequence name and the
/// coordinates (1-based column numbers), the comment character and the
/// number of leading lines to skip.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TabixConf {
    pub preset: i32,
    pub sc: i32,
    pub bc: i32,
    pub ec: i32,
    pub meta: u8,
    pub skip: i32,
}

impl TabixConf {
    pub const GFF: Self = Self { preset: TBX_GENERIC, sc: 1, bc: 4, ec: 5, meta: b'#', skip: 0 };
    pub const BED: Self = Self { preset: TBX_UCSC, sc: 1, bc: 2, ec: 3, meta: b'#', skip: 0 };
    pub const PSLTBL: Self = Self { preset: TBX_UCSC, sc: 15, bc: 17, ec: 18, meta: b'#', skip: 0 };
    pub const SAM: Self = Self { preset: TBX_SAM, sc: 3, bc: 4, ec: 0, meta: b'@', skip: 0 };
    pub const VCF: Self = Self { preset: TBX_VCF, sc: 1, bc: 2, ec: 0, meta: b'#', skip: 0 };

    /// `-p gff|bed|sam|vcf` (`psltbl` too).
    pub fn preset(name: &str) -> Option<Self> {
        match name.to_ascii_lowercase().as_str() {
            "gff" => Some(Self::GFF),
            "bed" => Some(Self::BED),
            "sam" => Some(Self::SAM),
            "vcf" => Some(Self::VCF),
            "psltbl" => Some(Self::PSLTBL),
            _ => None,
        }
    }

    /// Guess the preset from the file name, `.gz`/`.bgz` stripped.
    pub fn from_extension(path: &Path) -> Option<Self> {
        let name = path.file_name()?.to_str()?.to_ascii_lowercase();
        let base = name.strip_suffix(".gz").or_else(|| name.strip_suffix(".bgz")).unwrap_or(&name);
        let ext = base.rsplit('.').next()?;
        match ext {
            "vcf" => Some(Self::VCF),
            "bed" => Some(Self::BED),
            "gff" | "gff3" | "gtf" => Some(Self::GFF),
            "sam" => Some(Self::SAM),
            "psl" | "psltbl" => Some(Self::PSLTBL),
            _ => None,
        }
    }

    pub fn from_header(h: &TabixHeader) -> Self {
        Self { preset: h.format, sc: h.col_seq, bc: h.col_beg, ec: h.col_end, meta: h.meta as u8, skip: h.skip }
    }

    pub fn to_header(self, names: Vec<String>) -> TabixHeader {
        TabixHeader { format: self.preset, col_seq: self.sc, col_beg: self.bc, col_end: self.ec, meta: self.meta as i32, skip: self.skip, names }
    }

    /// The VCF preset on the standard columns, where [`vcf_line_interval`]
    /// applies.
    pub fn is_vcf(&self) -> bool {
        self.preset & 0xffff == TBX_VCF && self.sc == 1 && self.bc == 2
    }
}

/// `strtoll(s, .., 0)` without the hex/octal cases: leading blanks, a sign
/// and decimal digits.
fn leading_i64(s: &str) -> Option<i64> {
    let t = s.trim_start();
    let (neg, digits) = match t.strip_prefix('-') {
        Some(r) => (true, r),
        None => (false, t.strip_prefix('+').unwrap_or(t)),
    };
    let n = digits.bytes().take_while(u8::is_ascii_digit).count();
    if n == 0 {
        return None;
    }
    let v: i64 = digits[..n].parse().ok()?;
    Some(if neg { -v } else { v })
}

/// `tbx_parse1`: `(name, beg0, end0)` of a text line under `conf`.
pub fn tabix_line_interval<'a>(conf: &TabixConf, line: &'a str) -> Option<(&'a str, u64, u64)> {
    let mut name: Option<&str> = None;
    let (mut beg, mut end): (i64, i64) = (-1, -1);
    for (i, field) in line.split('\t').enumerate() {
        let id = i as i32 + 1;
        if id == conf.sc {
            name = Some(field);
        } else if id == conf.bc {
            beg = leading_i64(field)?;
            if conf.bc <= conf.ec {
                end = beg;
            }
            if conf.preset & TBX_UCSC == 0 {
                beg -= 1;
            } else if conf.bc <= conf.ec {
                end += 1;
            }
            if beg < 0 {
                eprintln!("Warning: Coordinate <= 0 detected. Did you forget to use the -0 option?");
                beg = 0;
            }
            if end < 1 {
                end = 1;
            }
        } else {
            match conf.preset & 0xffff {
                TBX_GENERIC => {
                    if id == conf.ec {
                        end = leading_i64(field)?;
                    }
                }
                TBX_SAM => {
                    if id == 6 {
                        // Reference span of the CIGAR: M, D and N consume reference.
                        let mut l: i64 = 0;
                        let mut num: i64 = 0;
                        for c in field.bytes() {
                            if c.is_ascii_digit() {
                                num = num * 10 + (c - b'0') as i64;
                            } else {
                                if matches!(c.to_ascii_uppercase(), b'M' | b'D' | b'N') {
                                    l += num;
                                }
                                num = 0;
                            }
                        }
                        if l == 0 {
                            l = 1;
                        }
                        end = beg + l;
                    }
                }
                TBX_VCF => {
                    if id == 4 {
                        if !field.is_empty() {
                            end = beg + field.len() as i64;
                        }
                    } else if id == 8 {
                        let v = field.strip_prefix("END=").or_else(|| field.find(";END=").map(|k| &field[k + 5..]));
                        if let Some(v) = v.filter(|v| !v.starts_with('.'))
                            && let Some(e) = leading_i64(v)
                                && e > beg {
                                    end = e;
                                }
                    }
                }
                _ => {}
            }
        }
    }
    let name = name?;
    if beg < 0 || end < 0 {
        return None;
    }
    Some((name, beg as u64, end as u64))
}

/// `adjust_n_lvls`: enough levels for `max_len` bases at `min_shift`.
fn adjust_n_lvls(min_shift: u8, mut n_lvls: u8, max_len: i64) -> u8 {
    let mut s: i64 = 1i64 << (min_shift as i64 + n_lvls as i64 * 3);
    while max_len > s {
        n_lvls += 1;
        s <<= 3;
    }
    n_lvls
}

/// The longest sequence a header line declares (`##contig` length, `@SQ LN`).
fn header_ref_len(conf: &TabixConf, line: &str) -> Option<i64> {
    match conf.preset & 0xffff {
        TBX_VCF => {
            let body = line.strip_prefix("##contig=<")?;
            let v = body.split(',').find_map(|kv| kv.strip_prefix("length="))?;
            leading_i64(v.trim_end_matches('>'))
        }
        TBX_SAM => {
            if !line.starts_with("@SQ") {
                return None;
            }
            line.split('\t').find_map(|f| f.strip_prefix("LN:")).and_then(leading_i64)
        }
        _ => None,
    }
}

/// Index a BGZF text file under a tabix column layout (`tbx_index`): the
/// first `conf.skip` lines and lines starting with `conf.meta` are header,
/// every other line contributes `[beg0, end0)` on its sequence. A CSI gets
/// as many levels as the longest declared sequence needs (100 Gbp when the
/// header declares none), as htslib does.
pub fn build_tabix_index(input: &Path, conf: &TabixConf, kind: IndexKind, min_shift: Option<u8>) -> Result<BinIndex> {
    if !is_bgzf(input).with_context(|| format!("open {}", input.display()))? {
        bail!("{}: not BGZF-compressed; only bgzip-compressed files can be indexed", input.display());
    }
    let (min_shift, mut depth) = match kind {
        IndexKind::Tbi => (14u8, 5u8),
        IndexKind::Csi => {
            let ms = min_shift.unwrap_or(14);
            (ms, tabix_depth_for(ms))
        }
    };
    let f = File::open(input).with_context(|| format!("open {}", input.display()))?;
    let mut reader = noodles_bgzf::io::Reader::new(f);
    let mut b: Option<IndexBuilder> = None;
    let mut names: Vec<String> = Vec::new();
    let mut name_ids: crate::fx::FxHashMap<String, usize> = crate::fx::FxHashMap::default();
    let mut line = String::new();
    let mut lineno = 0i64;
    let mut max_ref_len = 0i64;
    loop {
        let start = u64::from(reader.virtual_position());
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            break;
        }
        lineno += 1;
        let end = u64::from(reader.virtual_position());
        let is_meta = line.as_bytes().first() == Some(&conf.meta);
        if is_meta && kind == IndexKind::Csi
            && let Some(l) = header_ref_len(conf, line.trim_end_matches(['\r', '\n'])) {
                max_ref_len = max_ref_len.max(l);
            }
        if lineno <= conf.skip as i64 || is_meta || line.trim_end().is_empty() {
            continue;
        }
        let b = b.get_or_insert_with(|| {
            if kind == IndexKind::Csi {
                if max_ref_len == 0 {
                    max_ref_len = 100 * 1024 * 1024 * 1024;
                }
                depth = adjust_n_lvls(min_shift, depth, max_ref_len);
            }
            IndexBuilder::new(kind, min_shift, depth)
        });
        let l = line.trim_end_matches(['\r', '\n']);
        let interval = if conf.is_vcf() { vcf_line_interval(l) } else { tabix_line_interval(conf, l) };
        let Some((chrom, beg0, end0)) = interval else {
            // htslib skips lines it cannot parse.
            eprintln!("Warning: failed to parse line {lineno} of {}, was the wrong -p [type] used?\n{}", input.display(), l);
            continue;
        };
        let tid = match name_ids.get(chrom) {
            Some(&t) => t,
            None => {
                let t = names.len();
                names.push(chrom.to_string());
                name_ids.insert(chrom.to_string(), t);
                t
            }
        };
        b.add(tid, beg0, end0, start, end)?;
    }
    let b = b.unwrap_or_else(|| IndexBuilder::new(kind, min_shift, depth));
    let n = names.len();
    Ok(b.finish(Some(conf.to_header(names)), n))
}

fn is_bcf_bgzf(path: &Path) -> Result<bool> {
    let f = File::open(path)?;
    let mut r = noodles_bgzf::io::Reader::new(f);
    let mut magic = [0u8; 5];
    let n = r.read(&mut magic)?;
    Ok(n == 5 && magic == BCF_MAGIC)
}

/// Build a CSI or TBI index for a BGZF VCF or BCF, in memory.
pub fn build_index_in_memory(input: &Path, kind: IndexKind, min_shift: Option<u8>) -> Result<BinIndex> {
    if !is_bgzf(input).with_context(|| format!("open {}", input.display()))? {
        bail!(
            "{}: not BGZF-compressed; only bgzip-compressed VCF and BCF can be indexed",
            input.display()
        );
    }
    if is_bcf_bgzf(input)? {
        if kind == IndexKind::Tbi {
            bail!("TBI indices are only for VCF; use CSI for BCF");
        }
        return build_bcf_index(input, min_shift.unwrap_or(14));
    }
    build_vcf_index(input, kind, min_shift)
}

fn build_bcf_index(input: &Path, min_shift: u8) -> Result<BinIndex> {
    let mut reader = BcfReader::open(input)?;
    let info = HeaderInfo::parse(&reader.header_lines);
    let max_rid = reader.dict.contig_idx.values().copied().max().map(|m| m as usize + 1).unwrap_or(0);
    let max_len = info.contigs.iter().filter_map(|(_, _, l)| l).max();
    let depth = depth_for(min_shift, max_len.unwrap_or((1u64 << 31) - 1));
    let mut b = IndexBuilder::new(IndexKind::Csi, min_shift, depth);
    loop {
        let start = reader.virtual_position().unwrap_or(0);
        let Some((shared, _indiv)) = reader.read_record_raw()? else { break };
        let end = reader.virtual_position().unwrap_or(start);
        let Some(meta) = record_meta(&shared) else { bail!("corrupt BCF record") };
        if meta.rid < 0 {
            continue;
        }
        b.add(meta.rid as usize, meta.pos.max(0) as u64, record_end0(&meta), start, end)?;
    }
    Ok(b.finish(None, max_rid))
}

fn build_vcf_index(input: &Path, kind: IndexKind, min_shift: Option<u8>) -> Result<BinIndex> {
    build_tabix_index(input, &TabixConf::VCF, kind, min_shift)
}

/// Build and write an index next to `input` (or at `output`).
pub fn build_index(input: &Path, output: &Path, kind: IndexKind, min_shift: Option<u8>) -> Result<BinIndex> {
    let idx = build_index_in_memory(input, kind, min_shift)?;
    idx.save(output).with_context(|| format!("write {}", output.display()))?;
    Ok(idx)
}

/// Compatibility entry point: CSI with default parameters.
pub fn build_csi_index<P: AsRef<Path>>(vcf_path: P, output_path: P) -> Result<()> {
    build_index(vcf_path.as_ref(), output_path.as_ref(), IndexKind::Csi, None).map(|_| ())
}

/// `<path>.csi` or `<path>.tbi`, whichever exists (CSI preferred, like htslib).
pub fn find_index_for(path: &Path) -> Option<std::path::PathBuf> {
    let s = path.as_os_str().to_os_string();
    for ext in [".csi", ".tbi"] {
        let mut c = s.clone();
        c.push(ext);
        let p = std::path::PathBuf::from(c);
        if p.exists() {
            return Some(p);
        }
    }
    None
}

#[cfg(test)]
#[path = "../../tests/unit/csi_builder.rs"]
mod tests;
