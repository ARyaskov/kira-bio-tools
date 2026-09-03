//! Output sink shared by every VCF-writing command: plain VCF, BGZF VCF or
//! BCF, to a file or to standard output. `finish` finalizes the stream; a
//! sink dropped early still flushes what it has (BGZF finalizes on Drop) but
//! only `finish` reports errors.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::Compression;

pub use crate::annotate::postproc::{OutputKind, parse_output_type};
use crate::bcf::BcfWriter;
use crate::bgzf::{BgzfWriter, FILE_BUFFER_SIZE, STREAM_BUFFER_SIZE};

pub enum VcfSink {
    Plain(BufWriter<Box<dyn Write + Send>>),
    Bgzf(BgzfWriter),
    Bcf { writer: BcfWriter, partial: Vec<u8> },
}

fn is_stdout(path: Option<&Path>) -> bool {
    path.is_none_or(|p| p == Path::new("-"))
}

impl VcfSink {
    /// Open a sink. `path == None` or `-` means standard output. `headers` is
    /// needed for BCF (its dictionary is built from the header).
    pub fn open(path: Option<&Path>, kind: OutputKind, headers: &[String]) -> Result<Self> {
        Self::open_with_level(path, kind, None, headers)
    }

    /// Like [`open`] with an explicit compression level overriding the one
    /// embedded in `kind` (`-l`).
    pub fn open_with_level(
        path: Option<&Path>,
        kind: OutputKind,
        level: Option<u32>,
        headers: &[String],
    ) -> Result<Self> {
        let streaming = is_stdout(path);
        let raw: Box<dyn Write + Send> = if streaming {
            Box::new(io::stdout())
        } else {
            let p = path.unwrap();
            Box::new(File::create(p).with_context(|| format!("create {}", p.display()))?)
        };
        Ok(match kind {
            OutputKind::Vcf => Self::Plain(BufWriter::with_capacity(1 << 20, raw)),
            OutputKind::VcfGz(l) => {
                let lvl = level.unwrap_or(l).min(9);
                let buf = if streaming { STREAM_BUFFER_SIZE } else { FILE_BUFFER_SIZE };
                Self::Bgzf(BgzfWriter::from_writer_buffered(raw, Compression::new(lvl), buf)?)
            }
            OutputKind::Bcf(l) => {
                let lvl = level.unwrap_or(l).min(9);
                let writer = BcfWriter::from_writer(raw, lvl > 0, lvl, headers, streaming)?;
                Self::Bcf { writer, partial: Vec::new() }
            }
        })
    }

    pub fn is_bcf(&self) -> bool {
        matches!(self, Self::Bcf { .. })
    }

    /// Write header lines (no-op for BCF, whose header was written at open).
    pub fn write_header(&mut self, headers: &[String]) -> Result<()> {
        if self.is_bcf() {
            return Ok(());
        }
        for h in headers {
            self.write_line(h)?;
        }
        Ok(())
    }

    /// Write one line (without trailing newline).
    pub fn write_line(&mut self, line: &str) -> Result<()> {
        match self {
            Self::Plain(w) => {
                w.write_all(line.as_bytes())?;
                w.write_all(b"\n")?;
            }
            Self::Bgzf(w) => {
                w.write_all(line.as_bytes())?;
                w.write_all(b"\n")?;
            }
            Self::Bcf { writer, .. } => {
                if !line.is_empty() && line.as_bytes()[0] != b'#' {
                    writer.write_vcf_line(line)?;
                }
            }
        }
        Ok(())
    }

    pub fn finish(self) -> Result<()> {
        match self {
            Self::Plain(mut w) => w.flush().context("flush output"),
            Self::Bgzf(w) => w.finish().context("finalize BGZF output"),
            Self::Bcf { writer, partial } => {
                if !partial.is_empty() {
                    let line = String::from_utf8_lossy(&partial).into_owned();
                    let line = line.trim_end_matches(['\r', '\n']);
                    if !line.is_empty() {
                        // Header lines are ignored by write_vcf_line.
                        let mut writer = writer;
                        writer.write_vcf_line(line)?;
                        return writer.finish();
                    }
                }
                writer.finish()
            }
        }
    }
}

impl Write for VcfSink {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Self::Plain(w) => w.write_all(buf)?,
            Self::Bgzf(w) => w.write_all(buf)?,
            Self::Bcf { writer, partial } => {
                let mut rest = buf;
                while let Some(nl) = memchr::memchr(b'\n', rest) {
                    let (head, tail) = rest.split_at(nl + 1);
                    let line: Vec<u8> = if partial.is_empty() {
                        head.to_vec()
                    } else {
                        let mut v = std::mem::take(partial);
                        v.extend_from_slice(head);
                        v
                    };
                    let s = std::str::from_utf8(&line)
                        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "non-UTF-8 VCF line"))?;
                    let s = s.trim_end_matches(['\r', '\n']);
                    if !s.is_empty() {
                        writer.write_vcf_line(s).map_err(io::Error::other)?;
                    }
                    rest = tail;
                }
                partial.extend_from_slice(rest);
            }
        }
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Plain(w) => w.flush(),
            Self::Bgzf(w) => w.flush(),
            Self::Bcf { .. } => Ok(()),
        }
    }
}
