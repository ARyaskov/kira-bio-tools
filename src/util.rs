use memmap2::MmapMut;
use std::fs::OpenOptions;
use std::path::Path;
use std::time::{Duration, Instant};

/// Contig id: index into the source's contig dictionary (header order, then
/// order of first appearance). Never derived from the contig name.
pub type ChrId = u32;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GenomicKey(u64);

impl GenomicKey {
    #[inline]
    pub fn new(chr: ChrId, position: u32) -> Self {
        Self(((chr as u64) << 32) | (position as u64))
    }

    #[inline]
    pub fn chr(&self) -> ChrId {
        (self.0 >> 32) as ChrId
    }

    #[inline]
    pub fn position(&self) -> u32 {
        self.0 as u32
    }

    #[inline]
    pub fn as_u64(&self) -> u64 {
        self.0
    }

    #[inline]
    pub fn from_u64(v: u64) -> Self {
        Self(v)
    }
}

/// Canonical human chromosome number (1-22, X=23, Y=24, MT=25). Only used by
/// the annotation engine as a fallback key space; everything else resolves
/// contigs through the header dictionary.
pub fn chr_name_to_id(name: &str) -> Option<u8> {
    let n = name.trim().trim_start_matches("chr");

    match n {
        "1" => Some(1),
        "2" => Some(2),
        "3" => Some(3),
        "4" => Some(4),
        "5" => Some(5),
        "6" => Some(6),
        "7" => Some(7),
        "8" => Some(8),
        "9" => Some(9),
        "10" => Some(10),
        "11" => Some(11),
        "12" => Some(12),
        "13" => Some(13),
        "14" => Some(14),
        "15" => Some(15),
        "16" => Some(16),
        "17" => Some(17),
        "18" => Some(18),
        "19" => Some(19),
        "20" => Some(20),
        "21" => Some(21),
        "22" => Some(22),
        "X" => Some(23),
        "Y" => Some(24),
        "MT" | "M" => Some(25),
        _ => None,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VcfFormat {
    Plain,
    Gzip,
    Bgzf,
}

pub fn detect_format<P: AsRef<Path>>(path: P) -> std::io::Result<VcfFormat> {
    use std::fs::File;
    use std::io::Read;

    let mut file = File::open(path)?;
    let mut header = [0u8; 18];
    let n = file.read(&mut header)?;

    if n < 2 {
        return Ok(VcfFormat::Plain);
    }

    if header[0] == 0x1f && header[1] == 0x8b {
        if n >= 18 && header[12] == b'B' && header[13] == b'C' {
            return Ok(VcfFormat::Bgzf);
        }
        return Ok(VcfFormat::Gzip);
    }

    Ok(VcfFormat::Plain)
}

pub struct Timer {
    start: Instant,
    label: String,
}

impl Timer {
    pub fn new(label: &str) -> Self {
        Self {
            start: Instant::now(),
            label: label.to_string(),
        }
    }

    pub fn elapsed(&self) -> Duration {
        self.start.elapsed()
    }

    pub fn print_elapsed(&self) {
        let elapsed = self.elapsed();
        if elapsed.as_secs() > 0 {
            eprintln!("{}: {:.2}s", self.label, elapsed.as_secs_f64());
        } else {
            eprintln!("{}: {:.2}ms", self.label, elapsed.as_secs_f64() * 1000.0);
        }
    }
}

impl Drop for Timer {
    fn drop(&mut self) {
        self.print_elapsed();
    }
}

/// C `%g` with six significant digits, as htslib prints floats.
pub fn fmt_g(v: f64) -> String {
    if v == 0.0 {
        return "0".to_string();
    }
    if !v.is_finite() {
        return if v.is_nan() { "nan".into() } else if v > 0.0 { "inf".into() } else { "-inf".into() };
    }
    let trim = |s: &str| -> String {
        if s.contains('.') { s.trim_end_matches('0').trim_end_matches('.').to_string() } else { s.to_string() }
    };
    let exp = v.abs().log10().floor() as i32;
    if exp < -4 || exp >= 6 {
        let s = format!("{:.5e}", v);
        let (mant, e) = s.split_once('e').unwrap_or((s.as_str(), "0"));
        let e: i32 = e.parse().unwrap_or(0);
        format!("{}e{}{:02}", trim(mant), if e < 0 { '-' } else { '+' }, e.abs())
    } else {
        trim(&format!("{:.*}", (5 - exp).max(0) as usize, v))
    }
}

pub fn format_size(bytes: u64) -> String {
    const KB: u64 = 1024;
    const MB: u64 = KB * 1024;
    const GB: u64 = MB * 1024;

    if bytes >= GB {
        format!("{:.2} GB", bytes as f64 / GB as f64)
    } else if bytes >= MB {
        format!("{:.2} MB", bytes as f64 / MB as f64)
    } else if bytes >= KB {
        format!("{:.2} KB", bytes as f64 / KB as f64)
    } else {
        format!("{} B", bytes)
    }
}

/// A genomic region, 1-based inclusive. `start`/`end` are `None` when the
/// region names a whole contig.
#[derive(Debug, Clone)]
pub struct Region {
    pub chr: String,
    pub start: Option<u32>,
    pub end: Option<u32>,
}

/// Parse a decimal that may contain thousands separators (`1,000,000`) or an
/// SI suffix (`1k`, `2M`, `3G`), as htslib does.
pub fn parse_coordinate(s: &str) -> Option<u32> {
    let s = s.trim();
    if s.is_empty() {
        return None;
    }
    let (num, mult): (&str, u64) = match s.as_bytes()[s.len() - 1].to_ascii_lowercase() {
        b'k' => (&s[..s.len() - 1], 1_000),
        b'm' => (&s[..s.len() - 1], 1_000_000),
        b'g' => (&s[..s.len() - 1], 1_000_000_000),
        _ => (s, 1),
    };
    let cleaned: String = num.chars().filter(|c| *c != ',').collect();
    if cleaned.is_empty() || !cleaned.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let v: u64 = cleaned.parse().ok()?;
    u32::try_from(v.checked_mul(mult)?).ok()
}

impl Region {
    /// Parse `chr`, `chr:beg-end`, `chr:beg-` or `chr:pos`. With
    /// `one_coord`, `chr:pos` is a single position (bcftools); otherwise it
    /// runs to the end of the contig (tabix/samtools).
    pub fn parse_with(s: &str, one_coord: bool) -> Option<Self> {
        let s = s.trim();
        if s.is_empty() {
            return None;
        }
        // The last ':' separates the range so contig names may contain ':'.
        let Some(colon) = s.rfind(':') else {
            return Some(Region { chr: s.to_string(), start: None, end: None });
        };
        let (chr, range) = (&s[..colon], &s[colon + 1..]);
        if chr.is_empty() {
            return None;
        }
        if range.is_empty() {
            return Some(Region { chr: chr.to_string(), start: None, end: None });
        }
        match range.split_once('-') {
            Some((b, e)) => {
                let start = if b.is_empty() { 1 } else { parse_coordinate(b)? };
                let end = if e.is_empty() { u32::MAX } else { parse_coordinate(e)? };
                Some(Region { chr: chr.to_string(), start: Some(start.max(1)), end: Some(end) })
            }
            None => {
                let pos = parse_coordinate(range)?;
                let end = if one_coord { pos } else { u32::MAX };
                Some(Region { chr: chr.to_string(), start: Some(pos.max(1)), end: Some(end) })
            }
        }
    }

    /// bcftools semantics (`chr:pos` is one position).
    pub fn parse(s: &str) -> Option<Self> {
        Self::parse_with(s, true)
    }

    /// 1-based inclusive bounds, `1..=u32::MAX` for a whole contig.
    pub fn bounds(&self) -> (u32, u32) {
        (self.start.unwrap_or(1), self.end.unwrap_or(u32::MAX))
    }
}

pub fn url_encode_info_value(val: &str) -> String {
    let mut result = String::with_capacity(val.len() + 10);
    for ch in val.chars() {
        match ch {
            ' ' => result.push_str("%20"),
            ';' => result.push_str("%3B"),
            '=' => result.push_str("%3D"),
            '%' => result.push_str("%25"),
            ',' => result.push_str("%2C"),
            '\r' => result.push_str("%0D"),
            '\n' => result.push_str("%0A"),
            '\t' => result.push_str("%09"),
            _ => result.push(ch),
        }
    }
    result
}

pub fn url_decode_info_value(val: &str) -> String {
    let mut result = String::with_capacity(val.len());
    let mut chars = val.chars().peekable();

    while let Some(c) = chars.next() {
        if c == '%' {
            let hex1 = chars.next();
            let hex2 = chars.next();

            if let (Some(h1), Some(h2)) = (hex1, hex2) {
                let hex_str = format!("{}{}", h1, h2);
                if let Ok(byte) = u8::from_str_radix(&hex_str, 16) {
                    result.push(byte as char);
                    continue;
                }
            }

            result.push(c);
        } else {
            result.push(c);
        }
    }

    result
}

pub fn append_cstr(pool: &mut Vec<u8>, s: &str) -> usize {
    let ofs = pool.len();
    pool.extend_from_slice(s.as_bytes());
    pool.push(0);
    ofs
}

pub fn clean_info_values(val: &str) -> String {
    let v: Vec<&str> = val.split(',').filter(|s| *s != ".").collect();

    let v = match v.len() {
        0 => return String::new(),
        1 => {
            if v[0] == "." {
                return String::new();
            }
            v
        }
        _ => {
            let mut v = v;
            while v.last() == Some(&".") {
                v.pop();
            }
            v
        }
    };

    if v.is_empty() || v.iter().all(|s| s == &".") {
        return String::new();
    }

    v.join(",")
}

pub fn read_cstring(pool: &[u8], offset: usize) -> &str {
    if offset >= pool.len() {
        return "";
    }
    let end = pool[offset..]
        .iter()
        .position(|&b| b == 0)
        .map(|p| offset + p)
        .unwrap_or(pool.len());
    std::str::from_utf8(&pool[offset..end]).unwrap_or("")
}

use crate::annotate::structs::bundle::FieldNumber;
use fxhash::hash64;

#[inline]
pub fn fast_hash64(bytes: &[u8]) -> u64 {
    hash64(bytes)
}

pub fn extract_info_key(line: &str) -> Option<String> {
    if let Some(start) = line.find("ID=") {
        let rest = &line[start + 3..];
        if let Some(end) = rest.find(',') {
            return Some(rest[..end].to_string());
        }
    }
    None
}

pub fn extract_info_number(line: &str) -> Option<FieldNumber> {
    if let Some(start) = line.find("Number=") {
        let rest = &line[start + 7..];
        if let Some(end) = rest.find(',') {
            let num_str = &rest[..end];
            return match num_str {
                "0" => Some(FieldNumber::Zero),
                "1" => Some(FieldNumber::One),
                "." => Some(FieldNumber::Many),
                "A" => Some(FieldNumber::A),
                "R" => Some(FieldNumber::R),
                "G" => Some(FieldNumber::G),
                _ => Some(FieldNumber::Many),
            };
        }
    }
    None
}

pub fn choose_best_number(numbers: &[FieldNumber]) -> FieldNumber {
    use FieldNumber::*;

    let mut has_r = false;
    let mut has_a = false;
    let mut has_many = false;
    let mut has_one = false;

    for n in numbers {
        match n {
            R => has_r = true,
            A => has_a = true,
            Many => has_many = true,
            One => has_one = true,
            _ => {}
        }
    }

    if has_r {
        return R;
    }
    if has_a {
        return A;
    }
    if has_many {
        return Many;
    }
    if has_one {
        return One;
    }

    One
}

pub struct MmapWriter {
    file: std::fs::File,
    map: MmapMut,
    map_size: usize,
    grow_step: usize,
    offset: usize,
}

impl MmapWriter {
    pub fn create(path: &Path, initial_size: usize) -> std::io::Result<Self> {
        let size = initial_size.max(1);
        let file = OpenOptions::new()
            .create(true)
            .read(true)
            .write(true)
            .truncate(true)
            .open(path)?;
        file.set_len(size as u64)?;
        // SAFETY: the file was just created and sized by this process; no other
        // mapping or writer touches it while the map is alive.
        let map = unsafe { MmapMut::map_mut(&file)? };
        Ok(Self {
            file,
            map,
            map_size: size,
            grow_step: size,
            offset: 0,
        })
    }

    pub fn write_all(&mut self, data: &[u8]) -> std::io::Result<()> {
        self.ensure_capacity(data.len())?;
        let end = self.offset + data.len();
        self.map[self.offset..end].copy_from_slice(data);
        self.offset = end;
        Ok(())
    }

    pub fn flush(&mut self) -> std::io::Result<()> {
        self.map.flush()
    }

    pub fn finish(mut self, flush: bool) -> std::io::Result<()> {
        if flush {
            self.flush()?;
        }
        self.file.set_len(self.offset as u64)?;
        Ok(())
    }

    fn ensure_capacity(&mut self, additional: usize) -> std::io::Result<()> {
        let needed = self.offset + additional;
        if needed <= self.map_size {
            return Ok(());
        }
        let mut new_size = self.map_size;
        while new_size < needed {
            new_size += self.grow_step;
        }
        self.map.flush()?;
        self.file.set_len(new_size as u64)?;
        // SAFETY: the previous map was flushed and is replaced here; the file is
        // owned by this writer and not shared.
        self.map = unsafe { MmapMut::map_mut(&self.file)? };
        self.map_size = new_size;
        Ok(())
    }
}

#[cfg(test)]
#[path = "../tests/unit/util.rs"]
mod tests;
