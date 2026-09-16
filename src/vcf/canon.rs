//! Record text the way bcftools writes it. htslib parses every record it
//! reads and re-renders it: floats go through `kputd` (`%g`, six significant
//! digits), so `1e-1` comes out as `0.1` and `1e6` as `1e+06`, and a sample
//! column shorter than FORMAT is padded with missing fields (`.` becomes
//! `.:.:.`). kira copies text records verbatim; [`RecordCanon`] rewrites
//! the QUAL column, the Float-typed INFO/FORMAT values and short sample
//! columns of a line only when they differ from that rendering, leaving
//! already-canonical lines untouched.

use std::borrow::Cow;

use crate::fx::FxHashSet;
use memchr::{memchr, memchr_iter};

use super::header::{FieldType, HeaderInfo};
use crate::util::fmt_g;

/// One float value as bcftools prints it: parsed like `strtod`, narrowed to
/// `float`, rendered with `kputd`. Text that is not a float is kept.
pub fn canonical_float(s: &str) -> Cow<'_, str> {
    if looks_canonical(s.as_bytes()) {
        return Cow::Borrowed(s);
    }
    match s.parse::<f64>() {
        Ok(d) => {
            let r = fmt_g(d as f32 as f64);
            if r == s {
                Cow::Borrowed(s)
            } else {
                Cow::Owned(r)
            }
        }
        Err(_) => Cow::Borrowed(s),
    }
}

/// A comma-separated float vector; `.` entries stay.
pub fn canonical_floats(s: &str) -> Cow<'_, str> {
    let b = s.as_bytes();
    let mut out: Option<String> = None;
    let mut start = 0usize;
    loop {
        let end = memchr(b',', &b[start..]).map_or(s.len(), |i| start + i);
        let part = &s[start..end];
        match canonical_float(part) {
            Cow::Owned(r) => {
                let o = out.get_or_insert_with(|| {
                    let mut o = String::with_capacity(s.len() + 8);
                    o.push_str(&s[..start]);
                    o
                });
                o.push_str(&r);
            }
            Cow::Borrowed(_) => {
                if let Some(o) = out.as_mut() {
                    o.push_str(part);
                }
            }
        }
        if end == s.len() {
            break;
        }
        if let Some(o) = out.as_mut() {
            o.push(',');
        }
        start = end + 1;
    }
    out.map_or(Cow::Borrowed(s), Cow::Owned)
}

/// True when `b` is text `kputd` could have produced (or `.`), so parsing
/// and re-rendering it would give it back unchanged.
fn looks_canonical(b: &[u8]) -> bool {
    if b == b"." {
        return true;
    }
    let b = b.strip_prefix(b"-").unwrap_or(b);
    if b.is_empty() || b == b"0" {
        return b == b"0";
    }
    let (mant, exp) = match memchr(b'e', b) {
        Some(i) => (&b[..i], Some(&b[i + 1..])),
        None => (b, None),
    };
    let (int_part, frac, has_dot) = match memchr(b'.', mant) {
        Some(i) => (&mant[..i], &mant[i + 1..], true),
        None => (mant, &b""[..], false),
    };
    if int_part.is_empty()
        || !int_part.iter().all(u8::is_ascii_digit)
        || !frac.iter().all(u8::is_ascii_digit)
    {
        return false;
    }
    if int_part.len() > 1 && int_part[0] == b'0' {
        return false;
    }
    if has_dot && (frac.is_empty() || frac[frac.len() - 1] == b'0') {
        return false;
    }
    let zero_int = int_part == b"0";
    let lead_zeros = if zero_int {
        frac.iter().take_while(|&&c| c == b'0').count()
    } else {
        0
    };
    let sig = if zero_int {
        frac.len() - lead_zeros
    } else {
        int_part.len() + frac.len()
    };
    if sig == 0 || sig > 6 {
        return false;
    }
    match exp {
        // Fixed notation covers 0.0001 ..= 999999.
        None => {
            if zero_int {
                lead_zeros <= 3
            } else {
                int_part.len() <= 6
            }
        }
        // d.ddddde[+-]XX with at least two exponent digits, outside that range.
        Some(e) => {
            if zero_int || int_part.len() != 1 || e.len() < 3 {
                return false;
            }
            let (sign, digits) = (e[0], &e[1..]);
            if (sign != b'+' && sign != b'-')
                || !digits.iter().all(u8::is_ascii_digit)
                || (digits.len() > 2 && digits[0] == b'0')
            {
                return false;
            }
            let x: i32 = std::str::from_utf8(digits)
                .ok()
                .and_then(|d| d.parse().ok())
                .unwrap_or(0);
            let x = if sign == b'-' { -x } else { x };
            !(-4..6).contains(&x)
        }
    }
}

type Edit = (usize, usize, String);

fn splice(line: &str, edits: Vec<Edit>) -> Option<String> {
    if edits.is_empty() {
        return None;
    }
    let mut out = String::with_capacity(line.len() + 16);
    let mut last = 0;
    for (s, e, v) in edits {
        out.push_str(&line[last..s]);
        out.push_str(&v);
        last = e;
    }
    out.push_str(&line[last..]);
    Some(out)
}

/// Rewrites the Float-typed values and short sample columns of text
/// records, given the header.
#[derive(Clone, Debug)]
pub struct RecordCanon {
    info: FxHashSet<Box<[u8]>>,
    info_first: [bool; 256],
    format: FxHashSet<Box<[u8]>>,
    /// Last FORMAT column seen and the indices of its Float keys.
    fmt_cache: (Vec<u8>, Vec<usize>),
}

impl Default for RecordCanon {
    fn default() -> Self {
        Self {
            info: FxHashSet::default(),
            info_first: [false; 256],
            format: FxHashSet::default(),
            fmt_cache: (Vec::new(), Vec::new()),
        }
    }
}

impl RecordCanon {
    pub fn from_headers<S: AsRef<str>>(headers: &[S]) -> Self {
        Self::from_header_info(&HeaderInfo::parse(headers))
    }

    pub fn from_header_info(h: &HeaderInfo) -> Self {
        let mut c = Self::default();
        for (k, d) in &h.info {
            if d.ty == FieldType::Float && !k.is_empty() {
                c.info_first[k.as_bytes()[0] as usize] = true;
                c.info.insert(k.as_bytes().into());
            }
        }
        for (k, d) in &h.format {
            if d.ty == FieldType::Float && !k.is_empty() {
                c.format.insert(k.as_bytes().into());
            }
        }
        c
    }

    pub fn info_is_float(&self, key: &str) -> bool {
        self.info.contains(key.as_bytes())
    }

    pub fn format_is_float(&self, key: &str) -> bool {
        self.format.contains(key.as_bytes())
    }

    /// The Float INFO values of an INFO column, re-rendered when needed.
    pub fn canon_info<'a>(&self, info: &'a str) -> Cow<'a, str> {
        let mut edits = Vec::new();
        self.info_edits(info, 0, info.len(), &mut edits);
        splice(info, edits).map_or(Cow::Borrowed(info), Cow::Owned)
    }

    /// The record with QUAL, every Float value and every short sample column
    /// re-rendered, or `None` when nothing had to change. Header lines are
    /// never touched.
    pub fn canon_line(&mut self, line: &str) -> Option<String> {
        let b = line.as_bytes();
        if b.first().is_none_or(|&c| c == b'#') {
            return None;
        }
        let mut tabs = memchr_iter(b'\t', b);
        let t5 = tabs.nth(4)?;
        let t6 = tabs.next()?;
        let mut edits: Vec<Edit> = Vec::new();
        if let Cow::Owned(q) = canonical_float(&line[t5 + 1..t6]) {
            edits.push((t5 + 1, t6, q));
        }
        let Some(t7) = tabs.next() else {
            return splice(line, edits);
        };
        let t8 = tabs.next();
        if !self.info.is_empty() {
            self.info_edits(line, t7 + 1, t8.unwrap_or(b.len()), &mut edits);
        }
        if let Some(t8) = t8 {
            let t9 = tabs.next();
            let fmt = &b[t8 + 1..t9.unwrap_or(b.len())];
            let n_fields = if fmt.is_empty() {
                0
            } else {
                memchr_iter(b':', fmt).count() + 1
            };
            let positions = self.float_positions(fmt);
            if let Some(t9) = t9 {
                let mut start = t9 + 1;
                loop {
                    let end = memchr(b'\t', &b[start..]).map_or(b.len(), |i| start + i);
                    sample_edits(line, start, end, positions, n_fields, &mut edits);
                    if end == b.len() {
                        break;
                    }
                    start = end + 1;
                }
            }
        }
        splice(line, edits)
    }

    fn info_edits(&self, line: &str, start: usize, end: usize, edits: &mut Vec<Edit>) {
        let b = line.as_bytes();
        let mut s = start;
        while s < end {
            let e = memchr(b';', &b[s..end]).map_or(end, |i| s + i);
            if let Some(eq) = memchr(b'=', &b[s..e]) {
                let key = &b[s..s + eq];
                if !key.is_empty() && self.info_first[key[0] as usize] && self.info.contains(key) {
                    let vs = s + eq + 1;
                    if let Cow::Owned(v) = canonical_floats(&line[vs..e]) {
                        edits.push((vs, e, v));
                    }
                }
            }
            s = e + 1;
        }
    }

    fn float_positions(&mut self, fmt: &[u8]) -> &[usize] {
        if self.fmt_cache.0 != fmt {
            self.fmt_cache.0.clear();
            self.fmt_cache.0.extend_from_slice(fmt);
            self.fmt_cache.1.clear();
            for (i, key) in fmt.split(|&c| c == b':').enumerate() {
                if self.format.contains(key) {
                    self.fmt_cache.1.push(i);
                }
            }
        }
        &self.fmt_cache.1
    }
}

/// Re-render the fields at `positions` (ascending) of one sample column and
/// pad it to `n_fields` fields with `.` the way htslib prints a short one.
fn sample_edits(
    line: &str,
    start: usize,
    end: usize,
    positions: &[usize],
    n_fields: usize,
    edits: &mut Vec<Edit>,
) {
    let b = line.as_bytes();
    let (mut pi, mut idx, mut s) = (0usize, 0usize, start);
    loop {
        let e = memchr(b':', &b[s..end]).map_or(end, |i| s + i);
        if pi < positions.len() && idx == positions[pi] {
            if let Cow::Owned(v) = canonical_floats(&line[s..e]) {
                edits.push((s, e, v));
            }
            pi += 1;
        }
        idx += 1;
        if e == end {
            break;
        }
        s = e + 1;
    }
    if idx < n_fields && end > start {
        edits.push((end, end, ":.".repeat(n_fields - idx)));
    }
}

#[cfg(test)]
#[path = "../../tests/unit/vcf_canon.rs"]
mod tests;
