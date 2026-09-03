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
        if let Some(last) = list.last_mut() {
            if last.1 >> 16 == p.start >> 16 {
                last.1 = p.end;
                return;
            }
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
                    if last.saturating_sub(first) < MIN_MARKER_DIST {
                        let moved = r.bins.remove(&id).unwrap();
                        r.bins.entry(bin_parent(id)).or_default().extend(moved);
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
                    if let Some(m) = merged.last_mut() {
                        if m.1 >> 16 >= s >> 16 {
                            if e > m.1 {
                                m.1 = e;
                            }
                            continue;
                        }
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
    let (min_shift, depth) = match kind {
        IndexKind::Tbi => (14u8, 5u8),
        IndexKind::Csi => {
            let ms = min_shift.unwrap_or(14);
            (ms, tabix_depth_for(ms))
        }
    };
    let f = File::open(input).with_context(|| format!("open {}", input.display()))?;
    let mut reader = noodles_bgzf::io::Reader::new(f);
    let mut b = IndexBuilder::new(kind, min_shift, depth);
    let mut names: Vec<String> = Vec::new();
    let mut name_ids: fxhash::FxHashMap<String, usize> = fxhash::FxHashMap::default();
    let mut line = String::new();
    loop {
        let start = u64::from(reader.virtual_position());
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            break;
        }
        let end = u64::from(reader.virtual_position());
        if line.starts_with('#') || line.trim_end().is_empty() {
            continue;
        }
        let Some((chrom, beg0, end0)) = vcf_line_interval(line.trim_end_matches(['\r', '\n'])) else {
            bail!("malformed VCF line: {}", line.trim_end());
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
    let n = names.len();
    Ok(b.finish(Some(TabixHeader::vcf(names)), n))
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
