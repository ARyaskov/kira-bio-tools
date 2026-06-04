//! Ktile freshness check.
//!
//! Compares the source-file mtime + size captured in the ktile header
//! against the current state of the source on disk.

use std::path::Path;

use super::reader::KtileReader;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KtileFreshness {
    /// Source size + mtime match the current source on disk.
    Fresh,
    /// Source size differs.
    StaleSize { ktile: u64, source: u64 },
    /// Source mtime differs.
    StaleMtime { ktile: u64, source: u64 },
    /// Metadata missing on one side; freshness undetermined.
    Unknown,
}

impl KtileFreshness {
    pub fn is_fresh(self) -> bool {
        matches!(self, Self::Fresh)
    }
}

/// Validates a ktile sidecar against its source VCF on disk.
pub fn check_ktile_freshness(
    ktile_path: &Path,
    source_path: &Path,
) -> std::io::Result<KtileFreshness> {
    let reader = KtileReader::open(ktile_path)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    let header = reader.header();

    // Sentinel 0 means writer couldn't stat the source (e.g. stdin); treat as Unknown.
    if header.source_size == 0 || header.source_mtime_unix == 0 {
        return Ok(KtileFreshness::Unknown);
    }

    let meta = std::fs::metadata(source_path)?;
    let cur_size = meta.len();
    if cur_size != header.source_size {
        return Ok(KtileFreshness::StaleSize {
            ktile: header.source_size,
            source: cur_size,
        });
    }
    let cur_mtime = meta
        .modified()?
        .duration_since(std::time::UNIX_EPOCH)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?
        .as_secs();
    if cur_mtime != header.source_mtime_unix {
        return Ok(KtileFreshness::StaleMtime {
            ktile: header.source_mtime_unix,
            source: cur_mtime,
        });
    }
    Ok(KtileFreshness::Fresh)
}

/// Canonical sidecar path next to `source_path` (`foo.vcf.gz` → `foo.vcf.gz.ktile`).
pub fn ktile_path_for(source_path: &Path) -> std::path::PathBuf {
    let mut p = source_path.to_path_buf();
    let mut name = p
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("input")
        .to_string();
    name.push_str(".ktile");
    p.set_file_name(name);
    p
}

#[cfg(test)]
#[path = "../../../tests/unit/annotate_ktile_freshness.rs"]
mod tests;
