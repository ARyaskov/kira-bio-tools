//! Reference FASTA access. Contigs are loaded on first use through the
//! `.fai` index (built in memory when the file has none), so a
//! chromosome-sorted pass holds one contig at a time. Gzipped input has no
//! `.gzi` and is loaded whole.

use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use anyhow::{Context, Result, bail};
use fxhash::FxHashMap;

use crate::bam::reader::FastaLike;

struct Contig {
    name: String,
    length: u64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

pub struct IndexedFasta {
    path: PathBuf,
    contigs: Vec<Contig>,
    by_name: FxHashMap<String, usize>,
    /// Uppercase sequence per contig; `None` once loading failed.
    seqs: Vec<OnceLock<Option<Box<[u8]>>>>,
}

fn fai_path(path: &Path) -> PathBuf {
    let mut p = path.as_os_str().to_os_string();
    p.push(".fai");
    PathBuf::from(p)
}

fn read_fai(path: &Path) -> Result<Vec<Contig>> {
    let mut out = Vec::new();
    for (ln, line) in BufReader::new(File::open(path)?).lines().enumerate() {
        let line = line?;
        let t = line.trim_end();
        if t.is_empty() {
            continue;
        }
        let f: Vec<&str> = t.split('\t').collect();
        if f.len() < 5 {
            bail!("{}:{}: expected 5 columns", path.display(), ln + 1);
        }
        let num = |i: usize| -> Result<u64> {
            f[i].parse::<u64>().with_context(|| format!("{}:{}: bad column {}", path.display(), ln + 1, i + 1))
        };
        out.push(Contig { name: f[0].to_string(), length: num(1)?, offset: num(2)?, line_bases: num(3)?, line_width: num(4)? });
    }
    Ok(out)
}

/// A `.fai` written for another line-ending convention (a CRLF checkout of an
/// LF file) points into the wrong bytes; check the newline geometry before
/// trusting it.
fn fai_matches_file(path: &Path, contigs: &[Contig]) -> bool {
    let Ok(mut f) = File::open(path) else { return false };
    let mut b = [0u8; 2];
    for c in contigs {
        if c.line_bases == 0 || c.length == 0 {
            continue;
        }
        if c.offset == 0 || c.line_width < c.line_bases {
            return false;
        }
        let nl_before = f.seek(SeekFrom::Start(c.offset - 1)).is_ok() && f.read_exact(&mut b[..1]).is_ok() && b[0] == b'\n';
        if !nl_before {
            return false;
        }
        if c.length > c.line_bases {
            if f.seek(SeekFrom::Start(c.offset + c.line_bases)).is_err() || f.read_exact(&mut b).is_err() {
                return false;
            }
            let ok = match c.line_width - c.line_bases {
                1 => b[0] == b'\n',
                2 => b == [b'\r', b'\n'],
                _ => false,
            };
            if !ok {
                return false;
            }
        }
    }
    true
}

/// Scan a plain FASTA once: index every contig, keeping in memory only the
/// ones whose line layout is irregular (they cannot be re-read by offset).
fn build_fai(path: &Path) -> Result<(Vec<Contig>, Vec<Option<Box<[u8]>>>)> {
    let mut r = BufReader::with_capacity(1 << 20, File::open(path)?);
    let mut contigs: Vec<Contig> = Vec::new();
    let mut preload: Vec<Option<Box<[u8]>>> = Vec::new();
    let mut raw = Vec::new();
    let mut offset = 0u64;
    let mut cur: Option<(Contig, Vec<u8>, bool)> = None;

    fn finish(cur: Option<(Contig, Vec<u8>, bool)>, contigs: &mut Vec<Contig>, preload: &mut Vec<Option<Box<[u8]>>>) {
        if let Some((c, seq, irregular)) = cur {
            preload.push(if irregular || c.line_bases == 0 { Some(seq.into_boxed_slice()) } else { None });
            contigs.push(c);
        }
    }

    loop {
        raw.clear();
        let n = r.read_until(b'\n', &mut raw)?;
        if n == 0 {
            break;
        }
        let line_start = offset;
        offset += n as u64;
        if raw.first() == Some(&b'>') {
            finish(cur.take(), &mut contigs, &mut preload);
            let rest = &raw[1..];
            let end = rest.iter().position(|&b| b == b' ' || b == b'\t' || b == b'\r' || b == b'\n').unwrap_or(rest.len());
            let name = std::str::from_utf8(&rest[..end]).context("non-UTF-8 contig name")?.to_string();
            cur = Some((Contig { name, length: 0, offset, line_bases: 0, line_width: 0 }, Vec::new(), false));
            continue;
        }
        let Some((c, seq, irregular)) = cur.as_mut() else { continue };
        let bases: usize = raw.iter().filter(|&&b| b != b'\n' && b != b'\r').count();
        if bases == 0 {
            // A blank line ends the regular layout unless it is the last one.
            *irregular = true;
            continue;
        }
        if c.line_bases == 0 {
            c.line_bases = bases as u64;
            c.line_width = n as u64;
            c.offset = line_start;
        } else if c.length % c.line_bases != 0 || (bases as u64 > c.line_bases) || (bases as u64 == c.line_bases && n as u64 != c.line_width) {
            // A short line followed by more sequence, or a changed width.
            *irregular = true;
        }
        c.length += bases as u64;
        seq.extend(raw.iter().filter(|&&b| b != b'\n' && b != b'\r').map(|b| b.to_ascii_uppercase()));
    }
    finish(cur.take(), &mut contigs, &mut preload);
    Ok((contigs, preload))
}

fn parse_fasta_bytes(data: &[u8], mut on_contig: impl FnMut(String, Vec<u8>)) -> Result<()> {
    let mut name: Option<String> = None;
    let mut cur: Vec<u8> = Vec::new();
    for line in data.split(|&b| b == b'\n') {
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' {
            if let Some(n) = name.take() {
                on_contig(n, std::mem::take(&mut cur));
            }
            let rest = &line[1..];
            let end = rest.iter().position(|&b| b == b' ' || b == b'\t' || b == b'\r').unwrap_or(rest.len());
            name = Some(std::str::from_utf8(&rest[..end]).context("non-UTF-8 contig name")?.to_string());
        } else {
            cur.extend(line.iter().filter(|&&b| b != b'\r').map(|b| b.to_ascii_uppercase()));
        }
    }
    if let Some(n) = name {
        on_contig(n, cur);
    }
    Ok(())
}

impl IndexedFasta {
    pub fn open(path: &Path) -> Result<Self> {
        let mut magic = [0u8; 2];
        let n = File::open(path).with_context(|| format!("open fasta {}", path.display()))?.read(&mut magic)?;
        if n == 2 && magic == [0x1f, 0x8b] {
            return Self::load_gzip(path);
        }
        let fai = fai_path(path);
        if fai.exists() {
            let contigs = read_fai(&fai)?;
            if fai_matches_file(path, &contigs) {
                return Ok(Self::from_contigs(path, contigs, None));
            }
            eprintln!(
                "[fasta] {}: .fai does not match the file's line layout (line endings changed?); indexing in memory",
                path.display()
            );
        }
        let (contigs, preload) = build_fai(path)?;
        Ok(Self::from_contigs(path, contigs, Some(preload)))
    }

    fn from_contigs(path: &Path, contigs: Vec<Contig>, preload: Option<Vec<Option<Box<[u8]>>>>) -> Self {
        let by_name = contigs.iter().enumerate().map(|(i, c)| (c.name.clone(), i)).collect();
        let mut seqs: Vec<OnceLock<Option<Box<[u8]>>>> = (0..contigs.len()).map(|_| OnceLock::new()).collect();
        if let Some(pre) = preload {
            for (lock, seq) in seqs.iter_mut().zip(pre) {
                if let Some(s) = seq {
                    let _ = lock.set(Some(s));
                }
            }
        }
        Self { path: path.to_path_buf(), contigs, by_name, seqs }
    }

    fn load_gzip(path: &Path) -> Result<Self> {
        let mut data = Vec::new();
        flate2::read::MultiGzDecoder::new(BufReader::new(File::open(path)?))
            .read_to_end(&mut data)
            .context("decompress fasta")?;
        let mut contigs = Vec::new();
        let mut preload = Vec::new();
        parse_fasta_bytes(&data, |name, seq| {
            contigs.push(Contig { name, length: seq.len() as u64, offset: 0, line_bases: 0, line_width: 0 });
            preload.push(Some(seq.into_boxed_slice()));
        })?;
        Ok(Self::from_contigs(path, contigs, Some(preload)))
    }

    pub fn has(&self, chr: &str) -> bool {
        self.by_name.contains_key(chr)
    }

    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.contigs.iter().map(|c| c.name.as_str())
    }

    pub fn length(&self, chr: &str) -> Option<u64> {
        self.by_name.get(chr).map(|&i| self.contigs[i].length)
    }

    /// Whole contig, uppercase, loaded on first use.
    pub fn contig(&self, chr: &str) -> Option<&[u8]> {
        let i = *self.by_name.get(chr)?;
        self.seqs[i]
            .get_or_init(|| match self.read_contig(i) {
                Ok(s) => Some(s),
                Err(e) => {
                    eprintln!("[fasta] {}: {e:#}", self.contigs[i].name);
                    None
                }
            })
            .as_deref()
    }

    pub fn base(&self, chr: &str, pos1: u32) -> Option<u8> {
        self.contig(chr)?.get((pos1 as usize).checked_sub(1)?).copied()
    }

    /// Exactly `len` bases from 1-based `pos1`.
    pub fn slice(&self, chr: &str, pos1: u32, len: usize) -> Option<&[u8]> {
        let s = self.contig(chr)?;
        let start = (pos1 as usize).saturating_sub(1);
        s.get(start..start + len)
    }

    /// Up to `len` bases from 1-based `pos1`, clamped at the contig end.
    pub fn slice_bytes(&self, chr: &str, pos1: u32, len: usize) -> Option<&[u8]> {
        let s = self.contig(chr)?;
        let start = (pos1 as usize).saturating_sub(1);
        let end = (start + len).min(s.len());
        (end > start).then(|| &s[start..end])
    }

    /// Drop every cached contig except `keep` (for chromosome-sorted passes).
    pub fn evict_except(&mut self, keep: &str) {
        for (i, c) in self.contigs.iter().enumerate() {
            if c.name != keep && c.line_bases > 0 && self.seqs[i].get().is_some_and(|s| s.is_some()) {
                self.seqs[i] = OnceLock::new();
            }
        }
    }

    fn read_contig(&self, i: usize) -> Result<Box<[u8]>> {
        let c = &self.contigs[i];
        if c.line_bases == 0 || c.length == 0 {
            return Ok(Box::default());
        }
        let n_lines = c.length.div_ceil(c.line_bases);
        let raw_len = c.length + (n_lines - 1) * (c.line_width - c.line_bases);
        let mut f = File::open(&self.path).with_context(|| format!("open {}", self.path.display()))?;
        f.seek(SeekFrom::Start(c.offset))?;
        let mut raw = vec![0u8; raw_len as usize];
        f.read_exact(&mut raw).with_context(|| format!("read contig {} at offset {}", c.name, c.offset))?;
        let mut out = Vec::with_capacity(c.length as usize);
        out.extend(raw.iter().filter(|&&b| b != b'\n' && b != b'\r').map(|b| b.to_ascii_uppercase()));
        out.truncate(c.length as usize);
        Ok(out.into_boxed_slice())
    }
}

impl FastaLike for IndexedFasta {
    fn slice(&self, chr: &str, pos: u32, len: usize) -> Option<&[u8]> {
        self.slice_bytes(chr, pos, len)
    }
}

#[cfg(test)]
#[path = "../tests/unit/fasta.rs"]
mod tests;
