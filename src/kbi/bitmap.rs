//! Position bitmap (`.kbm`): the set of record positions per contig as a
//! roaring-style bitmap, next to the `.kbi` point index. Each contig is
//! split into 65536-position chunks; a chunk with few positions is a sorted
//! array of low halves, a dense one a 1024-word bitset. Membership,
//! "any record in [a, b]" and set algebra (and, or, and-not) work on the
//! containers directly, so `isec`, `filter --mask` and `gtcheck` can decide
//! which records matter without decoding either file. A second, usually
//! tiny, set holds the reference bases records cover beyond their POS, so a
//! mask query sees whole REF spans.

use std::collections::BTreeMap;
use std::io::{Read, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};

const MAGIC: &[u8; 4] = b"KBM1";
/// Arrays turn into bitsets past this many positions (roaring's threshold).
const ARRAY_MAX: usize = 4096;
/// Longest REF span whose covered bases are recorded.
const MAX_SPAN: u32 = 1 << 20;

#[derive(Clone, Debug, PartialEq, Eq)]
enum Container {
    Array(Vec<u16>),
    Bitset(Box<[u64; 1024]>),
}

impl Container {
    fn contains(&self, low: u16) -> bool {
        match self {
            Self::Array(v) => v.binary_search(&low).is_ok(),
            Self::Bitset(b) => b[(low >> 6) as usize] >> (low & 63) & 1 == 1,
        }
    }

    fn insert(&mut self, low: u16) {
        match self {
            Self::Array(v) => {
                if let Err(i) = v.binary_search(&low) {
                    v.insert(i, low);
                    if v.len() > ARRAY_MAX {
                        let mut b = Box::new([0u64; 1024]);
                        for &x in v.iter() {
                            b[(x >> 6) as usize] |= 1u64 << (x & 63);
                        }
                        *self = Self::Bitset(b);
                    }
                }
            }
            Self::Bitset(b) => b[(low >> 6) as usize] |= 1u64 << (low & 63),
        }
    }

    fn len(&self) -> usize {
        match self {
            Self::Array(v) => v.len(),
            Self::Bitset(b) => b.iter().map(|w| w.count_ones() as usize).sum(),
        }
    }

    /// Any position with `lo <= low <= hi`.
    fn any_between(&self, lo: u16, hi: u16) -> bool {
        match self {
            Self::Array(v) => {
                let i = v.partition_point(|&x| x < lo);
                i < v.len() && v[i] <= hi
            }
            Self::Bitset(b) => {
                let (wl, wh) = ((lo >> 6) as usize, (hi >> 6) as usize);
                for w in wl..=wh {
                    let mut word = b[w];
                    if w == wl {
                        word &= u64::MAX << (lo & 63);
                    }
                    if w == wh {
                        word &= u64::MAX >> (63 - (hi & 63));
                    }
                    if word != 0 {
                        return true;
                    }
                }
                false
            }
        }
    }

    fn iter(&self) -> Box<dyn Iterator<Item = u16> + '_> {
        match self {
            Self::Array(v) => Box::new(v.iter().copied()),
            Self::Bitset(b) => Box::new((0..1024usize).flat_map(move |w| {
                let word = b[w];
                (0..64u32)
                    .filter(move |&i| word >> i & 1 == 1)
                    .map(move |i| (w * 64) as u16 + i as u16)
            })),
        }
    }

    fn to_bitset(&self) -> Box<[u64; 1024]> {
        match self {
            Self::Bitset(b) => b.clone(),
            Self::Array(v) => {
                let mut b = Box::new([0u64; 1024]);
                for &x in v {
                    b[(x >> 6) as usize] |= 1u64 << (x & 63);
                }
                b
            }
        }
    }

    fn from_bitset(b: Box<[u64; 1024]>) -> Option<Self> {
        let n: usize = b.iter().map(|w| w.count_ones() as usize).sum();
        if n == 0 {
            return None;
        }
        if n > ARRAY_MAX {
            return Some(Self::Bitset(b));
        }
        let mut v = Vec::with_capacity(n);
        for (w, &word) in b.iter().enumerate() {
            let mut bits = word;
            while bits != 0 {
                v.push((w * 64) as u16 + bits.trailing_zeros() as u16);
                bits &= bits - 1;
            }
        }
        Some(Self::Array(v))
    }
}

/// One contig's positions: containers keyed by the high 16 bits.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ContigBitmap {
    chunks: BTreeMap<u16, Container>,
}

impl ContigBitmap {
    pub fn insert(&mut self, pos: u32) {
        self.chunks
            .entry((pos >> 16) as u16)
            .or_insert_with(|| Container::Array(Vec::new()))
            .insert(pos as u16);
    }

    pub fn contains(&self, pos: u32) -> bool {
        self.chunks
            .get(&((pos >> 16) as u16))
            .is_some_and(|c| c.contains(pos as u16))
    }

    /// Any position in `[beg, end]` (inclusive).
    pub fn any_in_range(&self, beg: u32, end: u32) -> bool {
        if end < beg {
            return false;
        }
        let (hb, he) = ((beg >> 16) as u16, (end >> 16) as u16);
        self.chunks.range(hb..=he).any(|(&h, c)| {
            let lo = if h == hb { beg as u16 } else { 0 };
            let hi = if h == he { end as u16 } else { u16::MAX };
            c.any_between(lo, hi)
        })
    }

    pub fn len(&self) -> usize {
        self.chunks.values().map(Container::len).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.chunks.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = u32> + '_ {
        self.chunks
            .iter()
            .flat_map(|(&h, c)| c.iter().map(move |l| ((h as u32) << 16) | l as u32))
    }

    fn combine(&self, other: &Self, op: impl Fn(u64, u64) -> u64, keep_left_only: bool) -> Self {
        let mut out = Self::default();
        for (&h, a) in &self.chunks {
            match other.chunks.get(&h) {
                Some(b) => {
                    let (ba, bb) = (a.to_bitset(), b.to_bitset());
                    let mut r = Box::new([0u64; 1024]);
                    for i in 0..1024 {
                        r[i] = op(ba[i], bb[i]);
                    }
                    if let Some(c) = Container::from_bitset(r) {
                        out.chunks.insert(h, c);
                    }
                }
                None => {
                    if keep_left_only {
                        out.chunks.insert(h, a.clone());
                    }
                }
            }
        }
        out
    }

    pub fn and(&self, other: &Self) -> Self {
        self.combine(other, |a, b| a & b, false)
    }

    pub fn and_not(&self, other: &Self) -> Self {
        self.combine(other, |a, b| a & !b, true)
    }

    pub fn or(&self, other: &Self) -> Self {
        let mut out = self.combine(other, |a, b| a | b, true);
        for (&h, b) in &other.chunks {
            out.chunks.entry(h).or_insert_with(|| b.clone());
        }
        out
    }
}

/// Positions of every record, per contig name, plus the reference bases
/// their REF alleles cover beyond POS.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PosBitmap {
    contigs: BTreeMap<String, ContigBitmap>,
    tails: BTreeMap<String, ContigBitmap>,
}

fn insert_into(map: &mut BTreeMap<String, ContigBitmap>, chrom: &str, pos: u32) {
    match map.get_mut(chrom) {
        Some(c) => c.insert(pos),
        None => {
            let mut c = ContigBitmap::default();
            c.insert(pos);
            map.insert(chrom.to_string(), c);
        }
    }
}

fn union(
    a: &BTreeMap<String, ContigBitmap>,
    b: &BTreeMap<String, ContigBitmap>,
) -> BTreeMap<String, ContigBitmap> {
    let mut out = a.clone();
    for (name, y) in b {
        match out.get_mut(name) {
            Some(x) => *x = x.or(y),
            None => {
                out.insert(name.clone(), y.clone());
            }
        }
    }
    out
}

impl PosBitmap {
    pub fn insert(&mut self, chrom: &str, pos: u32) {
        insert_into(&mut self.contigs, chrom, pos);
    }

    /// A record at `pos` whose REF allele is `ref_len` bases long.
    pub fn insert_span(&mut self, chrom: &str, pos: u32, ref_len: u32) {
        self.insert(chrom, pos);
        let last = pos.saturating_add(ref_len.min(MAX_SPAN).saturating_sub(1));
        for p in pos + 1..=last {
            insert_into(&mut self.tails, chrom, p);
        }
    }

    pub fn contig(&self, chrom: &str) -> Option<&ContigBitmap> {
        self.contigs.get(chrom)
    }

    pub fn contains(&self, chrom: &str, pos: u32) -> bool {
        self.contig(chrom).is_some_and(|c| c.contains(pos))
    }

    /// Any record *starting* in `[beg, end]`.
    pub fn any_in_range(&self, chrom: &str, beg: u32, end: u32) -> bool {
        self.contig(chrom).is_some_and(|c| c.any_in_range(beg, end))
    }

    /// Any record whose REF span touches `[beg, end]`.
    pub fn any_span_in_range(&self, chrom: &str, beg: u32, end: u32) -> bool {
        self.any_in_range(chrom, beg, end)
            || self
                .tails
                .get(chrom)
                .is_some_and(|c| c.any_in_range(beg, end))
    }

    /// Number of records (positions).
    pub fn len(&self) -> usize {
        self.contigs.values().map(ContigBitmap::len).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Reference bases covered by REF alleles beyond their POS.
    pub fn tail_len(&self) -> usize {
        self.tails.values().map(ContigBitmap::len).sum()
    }

    pub fn contig_names(&self) -> impl Iterator<Item = &str> {
        self.contigs.keys().map(String::as_str)
    }

    /// Positions present in both (spans are dropped).
    pub fn and(&self, other: &Self) -> Self {
        let mut out = Self::default();
        for (name, a) in &self.contigs {
            if let Some(b) = other.contigs.get(name) {
                let c = a.and(b);
                if !c.is_empty() {
                    out.contigs.insert(name.clone(), c);
                }
            }
        }
        out
    }

    /// Positions and spans present in either.
    pub fn or(&self, other: &Self) -> Self {
        Self {
            contigs: union(&self.contigs, &other.contigs),
            tails: union(&self.tails, &other.tails),
        }
    }

    /// Positions of `self` absent from `other` (spans are dropped).
    pub fn and_not(&self, other: &Self) -> Self {
        let mut out = Self::default();
        for (name, a) in &self.contigs {
            let c = match other.contigs.get(name) {
                Some(b) => a.and_not(b),
                None => a.clone(),
            };
            if !c.is_empty() {
                out.contigs.insert(name.clone(), c);
            }
        }
        out
    }

    /// The positions and REF spans of every record of a VCF/BCF.
    pub fn from_vcf(path: &Path) -> Result<Self> {
        let mut reader = crate::vcf::UnifiedVcfReader::open(path)
            .with_context(|| format!("open {}", path.display()))?;
        let _ = reader.header()?;
        let mut out = Self::default();
        while let Some(line) = reader.read_line()? {
            if line.is_empty() || line.as_bytes()[0] == b'#' {
                continue;
            }
            let mut it = line.split('\t');
            let (Some(chrom), Some(pos)) = (it.next(), it.next()) else {
                continue;
            };
            let ref_len = it.nth(1).map_or(1, |r| r.len() as u32);
            if let Ok(p) = pos.parse::<u32>() {
                out.insert_span(chrom, p, ref_len);
            }
        }
        Ok(out)
    }

    /// `KBM1`, then the position table and the span table, each: contig
    /// count and per contig its name, chunk count and each chunk as (high
    /// half, kind, count, payload); all integers little-endian.
    pub fn save(&self, path: &Path) -> Result<()> {
        let mut w = std::io::BufWriter::new(
            std::fs::File::create(path).with_context(|| format!("create {}", path.display()))?,
        );
        w.write_all(MAGIC)?;
        for table in [&self.contigs, &self.tails] {
            w.write_all(&(table.len() as u32).to_le_bytes())?;
            for (name, c) in table {
                w.write_all(&(name.len() as u32).to_le_bytes())?;
                w.write_all(name.as_bytes())?;
                w.write_all(&(c.chunks.len() as u32).to_le_bytes())?;
                for (&h, cont) in &c.chunks {
                    w.write_all(&h.to_le_bytes())?;
                    match cont {
                        Container::Array(v) => {
                            w.write_all(&[0u8])?;
                            w.write_all(&(v.len() as u32).to_le_bytes())?;
                            for x in v {
                                w.write_all(&x.to_le_bytes())?;
                            }
                        }
                        Container::Bitset(b) => {
                            w.write_all(&[1u8])?;
                            w.write_all(&(b.len() as u32).to_le_bytes())?;
                            for x in b.iter() {
                                w.write_all(&x.to_le_bytes())?;
                            }
                        }
                    }
                }
            }
        }
        w.flush()?;
        Ok(())
    }

    pub fn load(path: &Path) -> Result<Self> {
        let mut bytes = Vec::new();
        std::fs::File::open(path)
            .with_context(|| format!("open {}", path.display()))?
            .read_to_end(&mut bytes)?;
        let mut c = std::io::Cursor::new(bytes.as_slice());
        let mut magic = [0u8; 4];
        c.read_exact(&mut magic)?;
        if &magic != MAGIC {
            bail!("{}: not a position bitmap", path.display());
        }
        let contigs = read_table(&mut c, path)?;
        let tails = read_table(&mut c, path)?;
        Ok(Self { contigs, tails })
    }

    /// `<path>.kbm`, the sidecar written next to a file's `.kbi`.
    pub fn sidecar_path(input: &Path) -> std::path::PathBuf {
        input.with_extension("kbm")
    }

    /// The sidecar of `input` if it exists.
    pub fn load_sidecar(input: &Path) -> Option<Self> {
        let p = Self::sidecar_path(input);
        if p.exists() {
            Self::load(&p).ok()
        } else {
            None
        }
    }
}

fn read_u32(c: &mut std::io::Cursor<&[u8]>) -> Result<u32> {
    let mut b = [0u8; 4];
    c.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}

fn read_table(
    c: &mut std::io::Cursor<&[u8]>,
    path: &Path,
) -> Result<BTreeMap<String, ContigBitmap>> {
    let n_contigs = read_u32(c)?;
    let mut out = BTreeMap::new();
    for _ in 0..n_contigs {
        let n = read_u32(c)? as usize;
        let mut name = vec![0u8; n];
        c.read_exact(&mut name)?;
        let name = String::from_utf8(name).context("contig name")?;
        let n_chunks = read_u32(c)?;
        let mut cb = ContigBitmap::default();
        for _ in 0..n_chunks {
            let mut h = [0u8; 2];
            c.read_exact(&mut h)?;
            let h = u16::from_le_bytes(h);
            let mut kind = [0u8; 1];
            c.read_exact(&mut kind)?;
            let count = read_u32(c)? as usize;
            let cont = match kind[0] {
                0 => {
                    let mut v = Vec::with_capacity(count.min(1 << 16));
                    for _ in 0..count {
                        let mut x = [0u8; 2];
                        c.read_exact(&mut x)?;
                        v.push(u16::from_le_bytes(x));
                    }
                    Container::Array(v)
                }
                1 => {
                    if count != 1024 {
                        bail!("{}: corrupt bitset container", path.display());
                    }
                    let mut b = Box::new([0u64; 1024]);
                    for x in b.iter_mut() {
                        let mut w = [0u8; 8];
                        c.read_exact(&mut w)?;
                        *x = u64::from_le_bytes(w);
                    }
                    Container::Bitset(b)
                }
                k => bail!("{}: unknown container kind {k}", path.display()),
            };
            cb.chunks.insert(h, cont);
        }
        out.insert(name, cb);
    }
    Ok(out)
}

#[cfg(test)]
#[path = "../../tests/unit/kbi_bitmap.rs"]
mod tests;
