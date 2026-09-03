//! Genomic region sets (`-r/-R/-t/-T`): sorted, merged intervals per contig
//! with binary-search lookups. Every command filters through this type.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result, bail};
use fxhash::FxHashMap;
use memchr::memchr;

use crate::util::Region;

/// Overlap mode of bcftools `--regions-overlap`: 0 = POS inside a region,
/// 1 = the REF span overlaps, 2 = the REF or ALT span overlaps.
pub type OverlapMode = u8;

#[derive(Clone, Default, Debug)]
pub struct RegionSet {
    /// Sorted, non-overlapping 1-based inclusive intervals per contig.
    pub by_chr: FxHashMap<String, Vec<(u32, u32)>>,
    /// Contigs in order of first appearance.
    order: Vec<String>,
}

/// Parse one `chr[:beg[-end]]` token with bcftools semantics (`chr:pos` is a
/// single position; thousands separators and k/M/G suffixes are accepted).
pub fn parse_region_spec(s: &str) -> Result<(String, u32, u32)> {
    let r = Region::parse(s).ok_or_else(|| anyhow::anyhow!("invalid region '{s}'"))?;
    let (beg, end) = r.bounds();
    if end < beg {
        bail!("invalid region '{s}': end precedes start");
    }
    Ok((r.chr, beg, end))
}

impl RegionSet {
    pub fn from_cli(spec: &str) -> Result<Self> {
        let mut f = RegionSet::default();
        for raw in spec.split(',') {
            let item = raw.trim();
            if item.is_empty() {
                continue;
            }
            let (chr, beg, end) = parse_region_spec(item)?;
            f.add(&chr, beg, end);
        }
        f.finalize();
        Ok(f)
    }

    /// BED (`.bed`/`.bed.gz`: 0-based half-open) or tab-separated
    /// `CHROM [POS [END]]` (1-based inclusive) file, optionally gzipped.
    pub fn from_file<P: AsRef<Path>>(p: P) -> Result<Self> {
        let mut f = RegionSet::default();
        let path = p.as_ref();
        let lower = path.to_string_lossy().to_ascii_lowercase();
        let is_gz = lower.ends_with(".gz");
        let is_bed = lower.trim_end_matches(".gz").ends_with(".bed");
        let file = File::open(path).with_context(|| format!("open regions {}", path.display()))?;
        let reader: Box<dyn BufRead> = if is_gz {
            Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        for (ln, line) in reader.lines().enumerate() {
            let line = line?;
            let t = line.trim();
            if t.is_empty() || t.starts_with('#') || t.starts_with("track") || t.starts_with("browser") {
                continue;
            }
            let mut parts = t.split('\t');
            let chr = parts.next().unwrap_or("");
            if chr.is_empty() {
                continue;
            }
            let num = |v: Option<&str>, what: &str| -> Result<Option<u32>> {
                match v {
                    None => Ok(None),
                    Some(s) => Ok(Some(s.trim().parse::<u32>().with_context(|| {
                        format!("{}:{}: invalid {what} '{s}'", path.display(), ln + 1)
                    })?)),
                }
            };
            let b = num(parts.next(), "start")?;
            let e = num(parts.next(), "end")?;
            let (beg, end) = match (b, e) {
                (None, _) => (1, u32::MAX),
                (Some(b), None) => {
                    if is_bed { (b + 1, u32::MAX) } else { (b, b) }
                }
                (Some(b), Some(e)) => {
                    if is_bed { (b + 1, e) } else { (b, e) }
                }
            };
            if end < beg {
                bail!("{}:{}: end precedes start", path.display(), ln + 1);
            }
            f.add(chr, beg, end);
        }
        f.finalize();
        Ok(f)
    }

    /// `-r`/`-R` (or `-t`/`-T`) pair; `None` when neither is given.
    pub fn from_args(spec: Option<&str>, file: Option<&Path>) -> Result<Option<Self>> {
        let mut set = match spec {
            Some(s) => Self::from_cli(s)?,
            None => RegionSet::default(),
        };
        if let Some(p) = file {
            let f = Self::from_file(p)?;
            for (chr, beg, end) in f.iter() {
                set.add(chr, beg, end);
            }
            set.finalize();
        }
        Ok(if spec.is_none() && file.is_none() { None } else { Some(set) })
    }

    pub fn add(&mut self, chr: &str, beg: u32, end: u32) {
        match self.by_chr.get_mut(chr) {
            Some(v) => v.push((beg.max(1), end)),
            None => {
                self.order.push(chr.to_string());
                self.by_chr.insert(chr.to_string(), vec![(beg.max(1), end)]);
            }
        }
    }

    /// Sort and merge touching intervals; call after a batch of `add`s.
    pub fn finalize(&mut self) {
        for v in self.by_chr.values_mut() {
            v.sort_unstable();
            let mut merged: Vec<(u32, u32)> = Vec::with_capacity(v.len());
            for r in v.drain(..) {
                if let Some(last) = merged.last_mut() {
                    if r.0 <= last.1.saturating_add(1) {
                        last.1 = last.1.max(r.1);
                        continue;
                    }
                }
                merged.push(r);
            }
            *v = merged;
        }
    }

    pub fn is_empty(&self) -> bool {
        self.by_chr.is_empty()
    }

    /// Contigs in order of first appearance.
    pub fn contigs(&self) -> &[String] {
        &self.order
    }

    pub fn intervals(&self, chr: &str) -> Option<&[(u32, u32)]> {
        self.by_chr.get(chr).map(Vec::as_slice)
    }

    /// All intervals, contigs in first-appearance order, positions ascending.
    pub fn iter(&self) -> impl Iterator<Item = (&str, u32, u32)> + '_ {
        self.order.iter().flat_map(move |c| {
            self.by_chr.get(c).into_iter().flat_map(move |v| v.iter().map(move |&(b, e)| (c.as_str(), b, e)))
        })
    }

    pub fn contains(&self, chr: &str, pos: u32) -> bool {
        self.overlaps_range(chr, pos, pos)
    }

    pub fn overlaps_range(&self, chr: &str, beg: u32, end: u32) -> bool {
        let Some(ranges) = self.by_chr.get(chr) else { return false };
        let idx = ranges.partition_point(|r| r.1 < beg);
        idx < ranges.len() && ranges[idx].0 <= end
    }

    /// First interval end at or after `pos` on `chr`, for skipping ahead.
    pub fn next_interval(&self, chr: &str, pos: u32) -> Option<(u32, u32)> {
        let ranges = self.by_chr.get(chr)?;
        let idx = ranges.partition_point(|r| r.1 < pos);
        ranges.get(idx).copied()
    }

    /// Whether a record overlaps the set under `mode` (see [`OverlapMode`]).
    pub fn record_passes(&self, chrom: &str, pos: u32, refa: &str, alt: &str, mode: OverlapMode) -> bool {
        if mode == 0 {
            return self.contains(chrom, pos);
        }
        let span = if mode == 2 {
            let max_alt = alt
                .split(',')
                .filter(|a| !a.starts_with('<'))
                .map(str::len)
                .max()
                .unwrap_or(refa.len());
            refa.len().max(max_alt)
        } else {
            refa.len()
        };
        let end = pos + (span.max(1) as u32) - 1;
        self.overlaps_range(chrom, pos, end)
    }

    pub fn line_passes(&self, line: &str) -> bool {
        self.line_passes_mode(line, 1)
    }

    /// [`record_passes`] on a raw VCF line.
    pub fn line_passes_mode(&self, line: &str, mode: OverlapMode) -> bool {
        let bytes = line.as_bytes();
        let Some(t1) = memchr(b'\t', bytes) else { return false };
        let chr = &line[..t1];
        let rest = &line[t1 + 1..];
        let Some(t2) = memchr(b'\t', rest.as_bytes()) else { return false };
        let Ok(pos) = rest[..t2].parse::<u32>() else { return false };
        if mode == 0 {
            return self.contains(chr, pos);
        }
        let mut cols = rest[t2 + 1..].splitn(4, '\t');
        let _id = cols.next();
        let (Some(refa), Some(alt)) = (cols.next(), cols.next()) else {
            return self.contains(chr, pos);
        };
        self.record_passes(chr, pos, refa, alt, mode)
    }

    /// Stream records overlapping the set through `cb`, reading only the
    /// BGZF blocks the `.csi`/`.tbi` index points at. Fails without an index.
    pub fn stream_with_index<P: AsRef<Path>>(
        &self,
        input: P,
        mode: OverlapMode,
        mut cb: impl FnMut(&str) -> Result<()>,
    ) -> Result<()> {
        let mut reader = crate::csi::IndexedVcfReader::open(input.as_ref())?;
        let mut chroms: Vec<&String> = self.by_chr.keys().collect();
        // Index order keeps the output in file order.
        chroms.sort_by_key(|c| reader.ref_id(c).unwrap_or(usize::MAX));
        for chrom in chroms {
            let Some(ranges) = self.by_chr.get(chrom) else { continue };
            for &(beg, end) in ranges {
                reader.query(chrom, beg, end, |line| {
                    if self.line_passes_mode(line, mode) {
                        cb(line)?;
                    }
                    Ok(true)
                })?;
            }
        }
        Ok(())
    }

    pub fn has_index_for<P: AsRef<Path>>(input: P) -> Option<std::path::PathBuf> {
        crate::csi::find_index_for(input.as_ref())
    }
}

#[cfg(test)]
#[path = "../tests/unit/regions.rs"]
mod tests;
