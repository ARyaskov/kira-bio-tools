use std::path::Path;
use std::time::{Duration, Instant};

pub type ChrId = u8;

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

pub fn chr_id_to_name(id: ChrId) -> Option<&'static str> {
    static NAMES: [&str; 25] = [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
        "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
        "chr22", "chrX", "chrY", "chrM",
    ];
    if id >= 1 && id <= 25 {
        Some(NAMES[(id - 1) as usize])
    } else {
        None
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

#[inline]
fn parse_chromosome_fast(bytes: &[u8]) -> Option<ChrId> {
    let bytes = if bytes.len() > 3 && &bytes[..3] == b"chr" {
        &bytes[3..]
    } else {
        bytes
    };

    match bytes {
        b"X" | b"x" => Some(23),
        b"Y" | b"y" => Some(24),
        b"M" | b"MT" | b"m" | b"mt" => Some(25),
        _ if bytes.len() <= 2 && bytes.iter().all(|b| b.is_ascii_digit()) => {
            let mut num = 0u8;
            for &byte in bytes {
                num = num.wrapping_mul(10).wrapping_add(byte - b'0');
            }
            if (1..=22).contains(&num) {
                Some(num)
            } else {
                None
            }
        }
        _ => None,
    }
}

#[inline]
fn parse_u32_fast(bytes: &[u8]) -> Option<u32> {
    if bytes.is_empty() || bytes.len() > 10 {
        return None;
    }
    let mut result = 0u32;
    for &byte in bytes {
        if !byte.is_ascii_digit() {
            return None;
        }
        result = result.wrapping_mul(10).wrapping_add((byte - b'0') as u32);
    }
    Some(result)
}

pub struct Region {
    pub chr: String,
    pub start: Option<u32>,
    pub end: Option<u32>,
}

impl Region {
    pub fn parse(s: &str) -> Option<Self> {
        if let Some((chr, range)) = s.split_once(':') {
            if let Some((start, end)) = range.split_once('-') {
                Some(Region {
                    chr: chr.to_string(),
                    start: start.parse().ok(),
                    end: end.parse().ok(),
                })
            } else {
                let pos: u32 = range.parse().ok()?;
                Some(Region {
                    chr: chr.to_string(),
                    start: Some(pos),
                    end: Some(pos),
                })
            }
        } else {
            Some(Region {
                chr: s.to_string(),
                start: None,
                end: None,
            })
        }
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
