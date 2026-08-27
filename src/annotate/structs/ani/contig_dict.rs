//! Header-derived contig name ↔ id mapping for the ANI index.
//!
//! Built from the source VCF's `##contig` header lines at index-build time,
//! serialised into the `.ani` file, and loaded at lookup time.

use std::collections::HashMap;
use std::io::{self, Read, Write};

/// Contig name table. `id_to_name[i]` is the chr name with id `i`. The reverse
/// `name_to_id` is built lazily on the lookup side.
#[derive(Debug, Clone, Default)]
pub struct ContigDict {
    id_to_name: Vec<String>,
    name_to_id: HashMap<String, u32>,
}

impl ContigDict {
    /// Build a dict from an ordered list of unique chr names. The position
    /// in the list becomes the contig id (`u32`). Returns `None` on duplicates.
    pub fn from_names<I, S>(names: I) -> Option<Self>
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        let id_to_name: Vec<String> = names.into_iter().map(Into::into).collect();
        let mut name_to_id = HashMap::with_capacity(id_to_name.len());
        for (id, name) in id_to_name.iter().enumerate() {
            if name_to_id.insert(name.clone(), id as u32).is_some() {
                return None;
            }
        }
        Some(Self {
            id_to_name,
            name_to_id,
        })
    }

    /// Parse `##contig=<ID=…>` header lines from a VCF, returning a dict in
    /// the order the contigs were declared.
    pub fn from_header_lines<'a, I>(lines: I) -> Self
    where
        I: IntoIterator<Item = &'a str>,
    {
        let mut names = Vec::new();
        let mut seen = HashMap::new();
        for line in lines {
            let Some(id) = extract_contig_id(line) else {
                continue;
            };
            if seen.insert(id.clone(), names.len()).is_some() {
                continue;
            }
            names.push(id);
        }
        Self::from_names(names).unwrap_or_default()
    }

    /// O(1) lookup by chr name with bcftools-style chr-prefix normalisation.
    ///
    /// Tries exact match, then strip-`chr` / add-`chr`, then `M`/`MT`/`chrM`
    /// for the mitochondrion.
    #[inline]
    pub fn id(&self, name: &str) -> Option<u32> {
        if let Some(&id) = self.name_to_id.get(name) {
            return Some(id);
        }
        if let Some(stripped) = name.strip_prefix("chr") {
            if let Some(&id) = self.name_to_id.get(stripped) {
                return Some(id);
            }
            if stripped == "M"
                && let Some(&id) = self.name_to_id.get("MT")
            {
                return Some(id);
            }
        } else {
            let with_prefix = format!("chr{name}");
            if let Some(&id) = self.name_to_id.get(with_prefix.as_str()) {
                return Some(id);
            }
            if name == "MT"
                && let Some(&id) = self.name_to_id.get("chrM")
            {
                return Some(id);
            }
            if name == "M"
                && let Some(&id) = self.name_to_id.get("chrM")
            {
                return Some(id);
            }
        }
        None
    }

    /// O(1) reverse lookup.
    #[inline]
    pub fn name(&self, id: u32) -> Option<&str> {
        self.id_to_name.get(id as usize).map(String::as_str)
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.id_to_name.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.id_to_name.is_empty()
    }

    /// Append a new contig (idempotent: returns the existing id on duplicate).
    pub fn insert(&mut self, name: &str) -> u32 {
        if let Some(&id) = self.name_to_id.get(name) {
            return id;
        }
        let id = self.id_to_name.len() as u32;
        self.id_to_name.push(name.to_string());
        self.name_to_id.insert(name.to_string(), id);
        id
    }

    /// Wire format:
    ///   u32 LE: n_contigs
    ///   for each contig:
    ///     u16 LE: name_len
    ///     [u8; name_len]: name bytes (UTF-8, no NUL terminator)
    pub fn write_to<W: Write>(&self, mut out: W) -> io::Result<()> {
        out.write_all(&(self.id_to_name.len() as u32).to_le_bytes())?;
        for name in &self.id_to_name {
            let bytes = name.as_bytes();
            assert!(bytes.len() <= u16::MAX as usize, "contig name too long");
            out.write_all(&(bytes.len() as u16).to_le_bytes())?;
            out.write_all(bytes)?;
        }
        Ok(())
    }

    /// Serialised byte size.
    pub fn serialized_len(&self) -> usize {
        let mut total = 4;
        for name in &self.id_to_name {
            total += 2 + name.len();
        }
        total
    }

    pub fn read_from<R: Read>(mut input: R) -> io::Result<Self> {
        let mut buf4 = [0u8; 4];
        input.read_exact(&mut buf4)?;
        let n = u32::from_le_bytes(buf4) as usize;
        let mut names = Vec::with_capacity(n);
        let mut buf2 = [0u8; 2];
        for _ in 0..n {
            input.read_exact(&mut buf2)?;
            let len = u16::from_le_bytes(buf2) as usize;
            let mut name_bytes = vec![0u8; len];
            input.read_exact(&mut name_bytes)?;
            names.push(String::from_utf8(name_bytes).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "contig name not valid UTF-8")
            })?);
        }
        Self::from_names(names).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "duplicate contig name in .ani")
        })
    }

    /// Zero-copy parse from a slice (mmap-load path).
    pub fn parse_bytes(bytes: &[u8]) -> io::Result<Self> {
        if bytes.len() < 4 {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "contig dict too short"));
        }
        let n = u32::from_le_bytes(bytes[..4].try_into().unwrap()) as usize;
        let mut pos = 4usize;
        let mut names = Vec::with_capacity(n);
        for _ in 0..n {
            if pos + 2 > bytes.len() {
                return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "contig name length"));
            }
            let len = u16::from_le_bytes(bytes[pos..pos + 2].try_into().unwrap()) as usize;
            pos += 2;
            if pos + len > bytes.len() {
                return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "contig name bytes"));
            }
            let name = std::str::from_utf8(&bytes[pos..pos + len])
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "contig name UTF-8"))?
                .to_string();
            names.push(name);
            pos += len;
        }
        Self::from_names(names).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "duplicate contig name in .ani")
        })
    }
}

/// Pull the `ID=` value from a `##contig=<ID=...,...>` header line.
fn extract_contig_id(line: &str) -> Option<String> {
    let body = line.strip_prefix("##contig=<")?.strip_suffix('>')?;
    for kv in body.split(',') {
        let Some(rest) = kv.strip_prefix("ID=") else {
            continue;
        };
        let trimmed = rest.trim_matches('"');
        return Some(trimmed.to_string());
    }
    None
}

#[cfg(test)]
#[path = "../../../../tests/unit/annotate_structs_ani_contig_dict.rs"]
mod tests;
