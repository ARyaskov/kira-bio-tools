//! Compact line-batch buffer for the CPU annotate pipeline.
//!
//! Single contiguous `Vec<u8>` byte pool + `Vec<(u32, u32)>` line index.
//! Sendable through channels; lines are borrowed via [`ReadBatch::line`].

/// Owned line batch. Lines are stored contiguously without separators in
/// `bytes`; `lines[i] = (offset, length)` indexes line `i`.
pub struct ReadBatch {
    bytes: Vec<u8>,
    lines: Vec<(u32, u32)>,
}

impl ReadBatch {
    pub fn with_capacity(bytes_cap: usize, lines_cap: usize) -> Self {
        Self {
            bytes: Vec::with_capacity(bytes_cap),
            lines: Vec::with_capacity(lines_cap),
        }
    }

    pub fn len(&self) -> usize {
        self.lines.len()
    }

    pub fn is_empty(&self) -> bool {
        self.lines.is_empty()
    }

    pub fn byte_len(&self) -> usize {
        self.bytes.len()
    }

    /// Appends a line. Strips trailing `\n` / `\r\n` / `\r`.
    pub fn push_line(&mut self, line: &str) {
        let trimmed = line.trim_end_matches(['\n', '\r']);
        let offset = self.bytes.len();
        self.bytes.extend_from_slice(trimmed.as_bytes());
        let len = self.bytes.len() - offset;
        debug_assert!(offset <= u32::MAX as usize, "ReadBatch offset overflow");
        debug_assert!(len <= u32::MAX as usize, "ReadBatch line length overflow");
        self.lines.push((offset as u32, len as u32));
    }

    /// Appends a line from raw bytes (zero-copy fast path). Strips trailing `\n` / `\r`.
    pub fn push_line_bytes(&mut self, bytes: &[u8]) {
        let mut end = bytes.len();
        while end > 0 && (bytes[end - 1] == b'\n' || bytes[end - 1] == b'\r') {
            end -= 1;
        }
        let trimmed = &bytes[..end];
        let offset = self.bytes.len();
        self.bytes.extend_from_slice(trimmed);
        let len = self.bytes.len() - offset;
        debug_assert!(offset <= u32::MAX as usize, "ReadBatch offset overflow");
        debug_assert!(len <= u32::MAX as usize, "ReadBatch line length overflow");
        self.lines.push((offset as u32, len as u32));
    }

    /// Removes the most recently pushed line, rolling back its bytes and index entry.
    pub fn pop_last_line(&mut self) {
        if let Some((offset, _len)) = self.lines.pop() {
            self.bytes.truncate(offset as usize);
        }
    }

    /// First byte of the most recently pushed line, or `None` if empty.
    pub fn last_line_first_byte(&self) -> Option<u8> {
        let &(offset, len) = self.lines.last()?;
        if len == 0 {
            return None;
        }
        self.bytes.get(offset as usize).copied()
    }

    /// Borrows line `idx`. Returns `""` for out-of-range; `len()` is authoritative.
    pub fn line(&self, idx: usize) -> &str {
        let Some(&(offset, len)) = self.lines.get(idx) else {
            return "";
        };
        let start = offset as usize;
        let end = start + len as usize;
        // SAFETY: bytes were pushed from a `&str` so the slice is valid UTF-8.
        debug_assert!(self.bytes.get(start..end).is_some());
        unsafe { std::str::from_utf8_unchecked(&self.bytes[start..end]) }
    }

    pub fn iter(&self) -> ReadBatchIter<'_> {
        ReadBatchIter {
            batch: self,
            cursor: 0,
        }
    }
}

/// Sequential iterator over [`ReadBatch`] lines.
pub struct ReadBatchIter<'a> {
    batch: &'a ReadBatch,
    cursor: usize,
}

impl<'a> Iterator for ReadBatchIter<'a> {
    type Item = &'a str;
    fn next(&mut self) -> Option<&'a str> {
        if self.cursor >= self.batch.len() {
            return None;
        }
        let line = self.batch.line(self.cursor);
        self.cursor += 1;
        Some(line)
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.batch.len() - self.cursor;
        (remaining, Some(remaining))
    }
}

impl<'a> ExactSizeIterator for ReadBatchIter<'a> {}

#[cfg(test)]
#[path = "../../../tests/unit/annotate_cpu_v2_read_batch.rs"]
mod tests;
