//! Zero-copy VCF parser with minimal allocations
//!
//! Optimizations:
//! - No String allocations for fields
//! - Direct byte slice operations
//! - SIMD-friendly scanning where applicable

pub struct FastVcfParser<'a> {
    line: &'a str,
    pos: usize,
}

impl<'a> FastVcfParser<'a> {
    #[inline]
    pub fn new(line: &'a str) -> Self {
        Self { line, pos: 0 }
    }

    #[inline]
    pub fn parse_standard_fields(&mut self) -> Option<VcfFields<'a>> {
        let fields = VcfFields {
            chrom: self.next_field()?,
            pos: self.next_field()?,
            id: self.next_field()?,
            ref_allele: self.next_field()?,
            alt: self.next_field()?,
            qual: self.next_field()?,
            filter: self.next_field()?,
            info: self.next_field()?,
        };

        Some(fields)
    }

    #[inline]
    fn next_field(&mut self) -> Option<&'a str> {
        let start = self.pos;
        let bytes = self.line.as_bytes();

        while self.pos < bytes.len() && bytes[self.pos] != b'\t' {
            self.pos += 1;
        }

        let field = &self.line[start..self.pos];

        if self.pos < bytes.len() {
            self.pos += 1; // Skip tab
        }

        Some(field)
    }

    #[inline]
    pub fn rest(&self) -> &'a str {
        &self.line[self.pos..]
    }
}

pub struct VcfFields<'a> {
    pub chrom: &'a str,
    pub pos: &'a str,
    pub id: &'a str,
    pub ref_allele: &'a str,
    pub alt: &'a str,
    pub qual: &'a str,
    pub filter: &'a str,
    pub info: &'a str,
}

impl<'a> VcfFields<'a> {
    #[inline]
    pub fn position(&self) -> Option<u32> {
        fast_parse_u32(self.pos.as_bytes())
    }
}

#[inline]
fn fast_parse_u32(bytes: &[u8]) -> Option<u32> {
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

pub struct InfoParser<'a> {
    data: &'a str,
    pos: usize,
}

impl<'a> InfoParser<'a> {
    pub fn new(data: &'a str) -> Self {
        Self { data, pos: 0 }
    }

    pub fn iter(&self) -> InfoIterator<'a> {
        InfoIterator {
            data: self.data,
            pos: 0,
        }
    }
}

pub struct InfoIterator<'a> {
    data: &'a str,
    pos: usize,
}

impl<'a> Iterator for InfoIterator<'a> {
    type Item = (&'a str, Option<&'a str>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.data.len() {
            return None;
        }

        let start = self.pos;
        let bytes = self.data.as_bytes();

        // Find semicolon
        while self.pos < bytes.len() && bytes[self.pos] != b';' {
            self.pos += 1;
        }

        let segment = &self.data[start..self.pos];

        if self.pos < bytes.len() {
            self.pos += 1; // Skip semicolon
        }

        // Parse key=value
        if let Some(eq_pos) = segment.find('=') {
            Some((&segment[..eq_pos], Some(&segment[eq_pos + 1..])))
        } else {
            Some((segment, None))
        }
    }
}
