//! Binning index model with htslib-compatible CSI and TBI codecs.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, Read, Write};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use super::utils::{bin_first, bin_parent, max_pos, metadata_bin, reg2bins};

pub const CSI_MAGIC: [u8; 4] = *b"CSI\x01";
pub const TBI_MAGIC: [u8; 4] = *b"TBI\x01";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IndexKind {
    Csi,
    Tbi,
}

/// Tabix header (`tbx_conf_t` plus names): the TBI header, or the CSI aux block.
#[derive(Clone, Debug)]
pub struct TabixHeader {
    pub format: i32,
    pub col_seq: i32,
    pub col_beg: i32,
    pub col_end: i32,
    pub meta: i32,
    pub skip: i32,
    pub names: Vec<String>,
}

impl TabixHeader {
    /// `tbx_conf_vcf`.
    pub fn vcf(names: Vec<String>) -> Self {
        Self { format: 2, col_seq: 1, col_beg: 2, col_end: 0, meta: b'#' as i32, skip: 0, names }
    }

    fn write_to(&self, out: &mut Vec<u8>) {
        for v in [self.format, self.col_seq, self.col_beg, self.col_end, self.meta, self.skip] {
            out.extend_from_slice(&v.to_le_bytes());
        }
        let l_nm: usize = self.names.iter().map(|n| n.len() + 1).sum();
        out.extend_from_slice(&(l_nm as i32).to_le_bytes());
        for n in &self.names {
            out.extend_from_slice(n.as_bytes());
            out.push(0);
        }
    }

    fn parse(c: &mut Cursor<'_>) -> io::Result<Self> {
        let format = c.i32()?;
        let col_seq = c.i32()?;
        let col_beg = c.i32()?;
        let col_end = c.i32()?;
        let meta = c.i32()?;
        let skip = c.i32()?;
        let l_nm = c.i32()?;
        if l_nm < 0 {
            return Err(bad("negative name block length"));
        }
        let block = c.bytes(l_nm as usize)?;
        let names = block
            .split(|&b| b == 0)
            .filter(|s| !s.is_empty())
            .map(|s| String::from_utf8_lossy(s).into_owned())
            .collect();
        Ok(Self { format, col_seq, col_beg, col_end, meta, skip, names })
    }
}

#[derive(Clone, Debug, Default)]
pub struct BinEntry {
    /// Virtual offset of the first record overlapping the bin's interval (CSI).
    pub loffset: u64,
    pub chunks: Vec<(u64, u64)>,
}

#[derive(Clone, Copy, Debug)]
pub struct RefMeta {
    pub beg_off: u64,
    pub end_off: u64,
    pub n_mapped: u64,
    pub n_unmapped: u64,
}

#[derive(Clone, Debug, Default)]
pub struct RefIndex {
    pub bins: BTreeMap<u32, BinEntry>,
    /// 16 kb linear index (TBI only).
    pub linear: Vec<u64>,
    pub meta: Option<RefMeta>,
}

#[derive(Clone, Debug)]
pub struct BinIndex {
    pub kind: IndexKind,
    pub min_shift: u8,
    pub depth: u8,
    pub header: Option<TabixHeader>,
    pub refs: Vec<RefIndex>,
    pub n_no_coor: Option<u64>,
}

fn bad(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.to_string())
}

struct Cursor<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn bytes(&mut self, n: usize) -> io::Result<&'a [u8]> {
        if self.pos + n > self.data.len() {
            return Err(bad("truncated index"));
        }
        let s = &self.data[self.pos..self.pos + n];
        self.pos += n;
        Ok(s)
    }
    fn i32(&mut self) -> io::Result<i32> {
        let b = self.bytes(4)?;
        Ok(i32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }
    fn u32(&mut self) -> io::Result<u32> {
        let b = self.bytes(4)?;
        Ok(u32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }
    fn u64(&mut self) -> io::Result<u64> {
        let b = self.bytes(8)?;
        Ok(u64::from_le_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]]))
    }
    fn remaining(&self) -> usize {
        self.data.len() - self.pos
    }
}

fn count(n: i32, what: &str) -> io::Result<usize> {
    if n < 0 {
        return Err(bad(&format!("negative {what}")));
    }
    Ok(n as usize)
}

impl BinIndex {
    pub fn load<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let mut f = File::open(path.as_ref())?;
        let mut raw = Vec::new();
        let mut head = [0u8; 2];
        let n = f.read(&mut head)?;
        f = File::open(path.as_ref())?;
        if n == 2 && head == [0x1f, 0x8b] {
            MultiGzDecoder::new(f).read_to_end(&mut raw)?;
        } else {
            f.read_to_end(&mut raw)?;
        }
        Self::parse(&raw)
    }

    pub fn parse(raw: &[u8]) -> io::Result<Self> {
        let mut c = Cursor { data: raw, pos: 0 };
        let magic = c.bytes(4)?;
        if magic == CSI_MAGIC {
            Self::parse_csi(&mut c)
        } else if magic == TBI_MAGIC {
            Self::parse_tbi(&mut c)
        } else {
            Err(bad("not a CSI or TBI index"))
        }
    }

    fn parse_csi(c: &mut Cursor<'_>) -> io::Result<Self> {
        let min_shift = c.i32()?;
        let depth = c.i32()?;
        if !(1..=30).contains(&min_shift) || !(0..=10).contains(&depth) {
            return Err(bad("invalid CSI min_shift/depth"));
        }
        let (min_shift, depth) = (min_shift as u8, depth as u8);
        let l_aux = count(c.i32()?, "aux length")?;
        let aux = c.bytes(l_aux)?;
        let header = if l_aux >= 28 {
            let mut ac = Cursor { data: aux, pos: 0 };
            TabixHeader::parse(&mut ac).ok()
        } else {
            None
        };
        let n_ref = count(c.i32()?, "reference count")?;
        let meta_id = metadata_bin(depth);
        let mut refs = Vec::with_capacity(n_ref);
        for _ in 0..n_ref {
            let n_bin = count(c.i32()?, "bin count")?;
            let mut r = RefIndex::default();
            for _ in 0..n_bin {
                let id = c.u32()?;
                let loffset = c.u64()?;
                let n_chunk = count(c.i32()?, "chunk count")?;
                if id == meta_id {
                    if n_chunk != 2 {
                        return Err(bad("metadata bin must have 2 chunks"));
                    }
                    let beg_off = c.u64()?;
                    let end_off = c.u64()?;
                    let n_mapped = c.u64()?;
                    let n_unmapped = c.u64()?;
                    r.meta = Some(RefMeta { beg_off, end_off, n_mapped, n_unmapped });
                    continue;
                }
                let mut chunks = Vec::with_capacity(n_chunk);
                for _ in 0..n_chunk {
                    let s = c.u64()?;
                    let e = c.u64()?;
                    chunks.push((s, e));
                }
                r.bins.insert(id, BinEntry { loffset, chunks });
            }
            refs.push(r);
        }
        let n_no_coor = if c.remaining() >= 8 { Some(c.u64()?) } else { None };
        Ok(Self { kind: IndexKind::Csi, min_shift, depth, header, refs, n_no_coor })
    }

    fn parse_tbi(c: &mut Cursor<'_>) -> io::Result<Self> {
        let n_ref = count(c.i32()?, "reference count")?;
        let header = TabixHeader::parse(c)?;
        let (min_shift, depth) = (14u8, 5u8);
        let meta_id = metadata_bin(depth);
        let mut refs = Vec::with_capacity(n_ref);
        for _ in 0..n_ref {
            let n_bin = count(c.i32()?, "bin count")?;
            let mut r = RefIndex::default();
            for _ in 0..n_bin {
                let id = c.u32()?;
                let n_chunk = count(c.i32()?, "chunk count")?;
                if id == meta_id {
                    if n_chunk != 2 {
                        return Err(bad("metadata bin must have 2 chunks"));
                    }
                    let beg_off = c.u64()?;
                    let end_off = c.u64()?;
                    let n_mapped = c.u64()?;
                    let n_unmapped = c.u64()?;
                    r.meta = Some(RefMeta { beg_off, end_off, n_mapped, n_unmapped });
                    continue;
                }
                let mut chunks = Vec::with_capacity(n_chunk);
                for _ in 0..n_chunk {
                    let s = c.u64()?;
                    let e = c.u64()?;
                    chunks.push((s, e));
                }
                r.bins.insert(id, BinEntry { loffset: 0, chunks });
            }
            let n_intv = count(c.i32()?, "linear index size")?;
            let mut linear = Vec::with_capacity(n_intv);
            for _ in 0..n_intv {
                linear.push(c.u64()?);
            }
            r.linear = linear;
            refs.push(r);
        }
        let n_no_coor = if c.remaining() >= 8 { Some(c.u64()?) } else { None };
        Ok(Self { kind: IndexKind::Tbi, min_shift, depth, header: Some(header), refs, n_no_coor })
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        let meta_id = metadata_bin(self.depth);
        match self.kind {
            IndexKind::Csi => {
                out.extend_from_slice(&CSI_MAGIC);
                out.extend_from_slice(&(self.min_shift as i32).to_le_bytes());
                out.extend_from_slice(&(self.depth as i32).to_le_bytes());
                let mut aux = Vec::new();
                if let Some(h) = &self.header {
                    h.write_to(&mut aux);
                }
                out.extend_from_slice(&(aux.len() as i32).to_le_bytes());
                out.extend_from_slice(&aux);
                out.extend_from_slice(&(self.refs.len() as i32).to_le_bytes());
                for r in &self.refs {
                    let n_bin = r.bins.len() + usize::from(r.meta.is_some());
                    out.extend_from_slice(&(n_bin as i32).to_le_bytes());
                    for (&id, b) in &r.bins {
                        out.extend_from_slice(&id.to_le_bytes());
                        out.extend_from_slice(&b.loffset.to_le_bytes());
                        out.extend_from_slice(&(b.chunks.len() as i32).to_le_bytes());
                        for (s, e) in &b.chunks {
                            out.extend_from_slice(&s.to_le_bytes());
                            out.extend_from_slice(&e.to_le_bytes());
                        }
                    }
                    if let Some(m) = &r.meta {
                        out.extend_from_slice(&meta_id.to_le_bytes());
                        out.extend_from_slice(&0u64.to_le_bytes());
                        out.extend_from_slice(&2i32.to_le_bytes());
                        for v in [m.beg_off, m.end_off, m.n_mapped, m.n_unmapped] {
                            out.extend_from_slice(&v.to_le_bytes());
                        }
                    }
                }
            }
            IndexKind::Tbi => {
                out.extend_from_slice(&TBI_MAGIC);
                out.extend_from_slice(&(self.refs.len() as i32).to_le_bytes());
                let header = self.header.clone().unwrap_or_else(|| TabixHeader::vcf(Vec::new()));
                header.write_to(&mut out);
                for r in &self.refs {
                    let n_bin = r.bins.len() + usize::from(r.meta.is_some());
                    out.extend_from_slice(&(n_bin as i32).to_le_bytes());
                    for (&id, b) in &r.bins {
                        out.extend_from_slice(&id.to_le_bytes());
                        out.extend_from_slice(&(b.chunks.len() as i32).to_le_bytes());
                        for (s, e) in &b.chunks {
                            out.extend_from_slice(&s.to_le_bytes());
                            out.extend_from_slice(&e.to_le_bytes());
                        }
                    }
                    if let Some(m) = &r.meta {
                        out.extend_from_slice(&meta_id.to_le_bytes());
                        out.extend_from_slice(&2i32.to_le_bytes());
                        for v in [m.beg_off, m.end_off, m.n_mapped, m.n_unmapped] {
                            out.extend_from_slice(&v.to_le_bytes());
                        }
                    }
                    out.extend_from_slice(&(r.linear.len() as i32).to_le_bytes());
                    for v in &r.linear {
                        out.extend_from_slice(&v.to_le_bytes());
                    }
                }
            }
        }
        if let Some(n) = self.n_no_coor {
            out.extend_from_slice(&n.to_le_bytes());
        }
        out
    }

    /// Write the index BGZF-compressed, as htslib does.
    pub fn save<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let f = File::create(path)?;
        let mut w = noodles_bgzf::io::Writer::new(f);
        w.write_all(&self.to_bytes())?;
        w.finish()?;
        Ok(())
    }

    pub fn names(&self) -> &[String] {
        self.header.as_ref().map(|h| h.names.as_slice()).unwrap_or(&[])
    }

    pub fn ref_id(&self, name: &str) -> Option<usize> {
        self.names().iter().position(|n| n == name)
    }

    pub fn n_refs(&self) -> usize {
        self.refs.len()
    }

    pub fn n_records(&self, ref_id: usize) -> Option<u64> {
        self.refs.get(ref_id).and_then(|r| r.meta.map(|m| m.n_mapped))
    }

    pub fn total_records(&self) -> u64 {
        self.refs.iter().filter_map(|r| r.meta.map(|m| m.n_mapped)).sum()
    }

    pub fn max_position(&self) -> u64 {
        max_pos(self.min_shift, self.depth)
    }

    /// Virtual-position chunks that may contain records overlapping
    /// `[beg0, end0)` (0-based, half-open), merged and sorted (htslib
    /// `hts_itr_query`).
    pub fn query(&self, ref_id: usize, beg0: u64, end0: u64) -> Vec<(u64, u64)> {
        let Some(r) = self.refs.get(ref_id) else { return Vec::new() };
        let end0 = end0.min(self.max_position());
        if beg0 >= end0 {
            return Vec::new();
        }
        let min_off = match self.kind {
            IndexKind::Tbi => {
                if r.linear.is_empty() {
                    0
                } else {
                    let mut i = (beg0 >> self.min_shift) as usize;
                    if i >= r.linear.len() {
                        i = r.linear.len() - 1;
                    }
                    let mut off = r.linear[i];
                    while off == 0 && i > 0 {
                        i -= 1;
                        off = r.linear[i];
                    }
                    off
                }
            }
            IndexKind::Csi => {
                let mut bin = bin_first(self.depth) + (beg0 >> self.min_shift) as u32;
                let mut found = None;
                loop {
                    if let Some(e) = r.bins.get(&bin) {
                        found = Some(e.loffset);
                        break;
                    }
                    if bin == 0 {
                        break;
                    }
                    let first = (bin_parent(bin) << 3) + 1;
                    if bin > first {
                        bin -= 1;
                    } else {
                        bin = bin_parent(bin);
                    }
                }
                if found.is_none() {
                    if let Some(e) = r.bins.get(&0) {
                        found = Some(e.loffset);
                    }
                }
                found.unwrap_or(0)
            }
        };

        let mut bins = Vec::new();
        reg2bins(beg0, end0, self.min_shift, self.depth, &mut bins);
        let mut chunks: Vec<(u64, u64)> = Vec::new();
        for b in bins {
            if let Some(e) = r.bins.get(&b) {
                for &(s, en) in &e.chunks {
                    if en > min_off {
                        chunks.push((s.max(min_off), en));
                    }
                }
            }
        }
        if chunks.is_empty() {
            return chunks;
        }
        chunks.sort_unstable();
        let mut merged: Vec<(u64, u64)> = Vec::with_capacity(chunks.len());
        for (s, e) in chunks {
            if let Some(last) = merged.last_mut() {
                if s <= last.1 || (s >> 16) == (last.1 >> 16) {
                    if e > last.1 {
                        last.1 = e;
                    }
                    continue;
                }
            }
            merged.push((s, e));
        }
        merged
    }
}

#[cfg(test)]
#[path = "../../tests/unit/csi_binning.rs"]
mod tests;
