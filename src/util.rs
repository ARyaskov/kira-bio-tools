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

/// htslib's `kputd`: `%g` with six significant digits, as bcftools prints
/// every float (QUAL, Float INFO and FORMAT values).
pub fn fmt_g(v: f64) -> String {
    let mut s = String::with_capacity(16);
    push_g(&mut s, v);
    s
}

/// [`fmt_g`] appended to `out`.
pub fn push_g(out: &mut String, v: f64) {
    if v == 0.0 {
        out.push_str(if v.is_sign_negative() { "-0" } else { "0" });
        return;
    }
    if v.is_nan() {
        out.push_str("nan");
        return;
    }
    let d = if v < 0.0 {
        out.push('-');
        -v
    } else {
        v
    };
    if d.is_infinite() {
        out.push_str("inf");
        return;
    }
    if !(0.0001..=999999.0).contains(&d) {
        push_percent_g(out, d);
        return;
    }
    // Six significant digits: rint(d * 10^k), k by magnitude, as kputd does.
    const LIMIT: [f64; 9] = [0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0];
    const SCALE: [f64; 10] = [1e9, 1e8, 1e7, 1e6, 1e5, 1e4, 1e3, 1e2, 1e1, 1.0];
    let idx = LIMIT.iter().position(|&l| d < l).unwrap_or(9);
    let k = 9 - idx;
    let i = (d * SCALE[idx]).round_ties_even() as u64;
    let mut buf = itoa::Buffer::new();
    let digits = buf.format(i);
    if k == 0 {
        out.push_str(digits);
        return;
    }
    let (int_part, frac) = if digits.len() > k { digits.split_at(digits.len() - k) } else { ("", digits) };
    out.push_str(if int_part.is_empty() { "0" } else { int_part });
    let frac = frac.trim_end_matches('0');
    if frac.is_empty() {
        return;
    }
    out.push('.');
    for _ in digits.len()..k {
        out.push('0');
    }
    out.push_str(frac);
}

/// C `%g` for a positive finite value outside kputd's fixed range.
fn push_percent_g(out: &mut String, d: f64) {
    let sci = format!("{:.5e}", d);
    let (mant, exp) = sci.split_once('e').unwrap_or((sci.as_str(), "0"));
    let exp: i32 = exp.parse().unwrap_or(0);
    if (-4..6).contains(&exp) {
        let s = format!("{:.*}", (5 - exp) as usize, d);
        out.push_str(if s.contains('.') { s.trim_end_matches('0').trim_end_matches('.') } else { &s });
    } else {
        out.push_str(if mant.contains('.') { mant.trim_end_matches('0').trim_end_matches('.') } else { mant });
        out.push('e');
        out.push(if exp < 0 { '-' } else { '+' });
        let e = exp.unsigned_abs();
        if e < 10 {
            out.push('0');
        }
        let mut b = itoa::Buffer::new();
        out.push_str(b.format(e));
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

/// Current UTC time in `ctime` layout (`Tue Sep 16 12:00:00 2026`), for
/// `##…Command=…; Date=` header lines.
pub fn now_string() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    let days = secs.div_euclid(86_400);
    let rem = secs.rem_euclid(86_400);
    let (h, m, s) = (rem / 3600, rem % 3600 / 60, rem % 60);
    // Civil-from-days (Howard Hinnant), proleptic Gregorian calendar.
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let mo = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if mo <= 2 { y + 1 } else { y };
    const WD: [&str; 7] = ["Thu", "Fri", "Sat", "Sun", "Mon", "Tue", "Wed"];
    const MON: [&str; 12] = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
    format!("{} {} {:2} {:02}:{:02}:{:02} {}", WD[days.rem_euclid(7) as usize], MON[(mo - 1) as usize], d, h, m, s, y)
}

/// `erfc` as in fdlibm `s_erf.c` (the algorithm behind glibc and msvcrt).
pub fn erfc(x: f64) -> f64 {
    const TINY: f64 = 1e-300;
    const ERX: f64 = 8.45062911510467529297e-01;
    const PP: [f64; 5] = [
        1.28379167095512558561e-01,
        -3.25042107247001499370e-01,
        -2.84817495755985104766e-02,
        -5.77027029648944159157e-03,
        -2.37630166566501626084e-05,
    ];
    const QQ: [f64; 5] = [
        3.97917223959155352819e-01,
        6.50222499887672944485e-02,
        5.08130628187576562776e-03,
        1.32494738004321644526e-04,
        -3.96022827877536812320e-06,
    ];
    const PA: [f64; 7] = [
        -2.36211856075265944077e-03,
        4.14856118683748331666e-01,
        -3.72207876035701323847e-01,
        3.18346619901161753674e-01,
        -1.10894694282396677476e-01,
        3.54783043256182359371e-02,
        -2.16637559486879084300e-03,
    ];
    const QA: [f64; 6] = [
        1.06420880400844228286e-01,
        5.40397917702171048937e-01,
        7.18286544141962662868e-02,
        1.26171219808761642112e-01,
        1.36370839120290507362e-02,
        1.19844998467991074170e-02,
    ];
    const RA: [f64; 8] = [
        -9.86494403484714822705e-03,
        -6.93858572707181764372e-01,
        -1.05586262253232909814e+01,
        -6.23753324503260060396e+01,
        -1.62396669462573470355e+02,
        -1.84605092906711035994e+02,
        -8.12874355063065934246e+01,
        -9.81432934416914548592e+00,
    ];
    const SA: [f64; 8] = [
        1.96512716674392571292e+01,
        1.37657754143519042600e+02,
        4.34565877475229228821e+02,
        6.45387271733267880336e+02,
        4.29008140027567833386e+02,
        1.08635005541779435134e+02,
        6.57024977031928170135e+00,
        -6.04244152148580987438e-02,
    ];
    const RB: [f64; 7] = [
        -9.86494292470009928597e-03,
        -7.99283237680523006574e-01,
        -1.77579549177547519889e+01,
        -1.60636384855821916062e+02,
        -6.37566443368389627722e+02,
        -1.02509513161107724954e+03,
        -4.83519191608651397019e+02,
    ];
    const SB: [f64; 7] = [
        3.03380607434824582924e+01,
        3.25792512996573918826e+02,
        1.53672958608443695994e+03,
        3.19985821950859553908e+03,
        2.55305040643316442583e+03,
        4.74528541206955367215e+02,
        -2.24409524465858183362e+01,
    ];
    let hx = (x.to_bits() >> 32) as i32;
    let ix = hx & 0x7fff_ffff;
    if ix >= 0x7ff0_0000 {
        return (((hx as u32) >> 31) << 1) as f64 + 1.0 / x;
    }
    if ix < 0x3feb_0000 {
        if ix < 0x3c70_0000 {
            return 1.0 - x;
        }
        let z = x * x;
        let r = PP[0] + z * (PP[1] + z * (PP[2] + z * (PP[3] + z * PP[4])));
        let s = 1.0 + z * (QQ[0] + z * (QQ[1] + z * (QQ[2] + z * (QQ[3] + z * QQ[4]))));
        let y = r / s;
        if hx < 0x3fd0_0000 {
            return 1.0 - (x + x * y);
        }
        let mut r = x * y;
        r += x - 0.5;
        return 0.5 - r;
    }
    if ix < 0x3ff4_0000 {
        let s = x.abs() - 1.0;
        let p = PA[0] + s * (PA[1] + s * (PA[2] + s * (PA[3] + s * (PA[4] + s * (PA[5] + s * PA[6])))));
        let q = 1.0 + s * (QA[0] + s * (QA[1] + s * (QA[2] + s * (QA[3] + s * (QA[4] + s * QA[5])))));
        if hx >= 0 {
            let z = 1.0 - ERX;
            return z - p / q;
        }
        let z = ERX + p / q;
        return 1.0 + z;
    }
    if ix < 0x403c_0000 {
        let ax = x.abs();
        let s = 1.0 / (ax * ax);
        let (r, ss);
        if ix < 0x4006_db6d {
            r = RA[0] + s * (RA[1] + s * (RA[2] + s * (RA[3] + s * (RA[4] + s * (RA[5] + s * (RA[6] + s * RA[7]))))));
            ss = 1.0 + s * (SA[0] + s * (SA[1] + s * (SA[2] + s * (SA[3] + s * (SA[4] + s * (SA[5] + s * (SA[6] + s * SA[7])))))));
        } else {
            if hx < 0 && ix >= 0x4018_0000 {
                return 2.0 - TINY;
            }
            r = RB[0] + s * (RB[1] + s * (RB[2] + s * (RB[3] + s * (RB[4] + s * (RB[5] + s * RB[6])))));
            ss = 1.0 + s * (SB[0] + s * (SB[1] + s * (SB[2] + s * (SB[3] + s * (SB[4] + s * (SB[5] + s * SB[6]))))));
        }
        let z = f64::from_bits(ax.to_bits() & 0xffff_ffff_0000_0000);
        let rr = (-z * z - 0.5625).exp() * ((z - ax) * (z + ax) + r / ss).exp();
        return if hx > 0 { rr / ax } else { 2.0 - rr / ax };
    }
    if hx > 0 { TINY * TINY } else { 2.0 - TINY }
}

/// bcftools `phred_score`: `-4.3429 * ln(prob)`, capped at 99.
pub fn phred_score(prob: f64) -> f64 {
    if prob == 0.0 {
        return 99.0;
    }
    let p = -4.3429 * prob.ln();
    if p > 99.0 { 99.0 } else { p }
}

/// C `%.*e`: mantissa with `prec` decimals and a signed two-digit exponent.
pub fn fmt_e(v: f64, prec: usize) -> String {
    if !v.is_finite() {
        return if v.is_nan() { "nan".into() } else if v > 0.0 { "inf".into() } else { "-inf".into() };
    }
    let s = format!("{:.*e}", prec, v);
    let (mant, exp) = s.split_once('e').unwrap_or((&s, "0"));
    let e: i32 = exp.parse().unwrap_or(0);
    format!("{}e{}{:02}", mant, if e < 0 { '-' } else { '+' }, e.abs())
}

/// Line reader over a plain or gzip/BGZF-compressed text file.
pub fn open_text_maybe_gz(path: &Path) -> anyhow::Result<Box<dyn std::io::BufRead>> {
    use anyhow::Context as _;
    use std::io::Read as _;
    let mut f = std::fs::File::open(path).with_context(|| format!("Failed to read {}", path.display()))?;
    let mut magic = [0u8; 2];
    let n = f.read(&mut magic)?;
    let f = std::fs::File::open(path)?;
    if n == 2 && magic == [0x1f, 0x8b] {
        Ok(Box::new(std::io::BufReader::new(flate2::read::MultiGzDecoder::new(f))))
    } else {
        Ok(Box::new(std::io::BufReader::new(f)))
    }
}
