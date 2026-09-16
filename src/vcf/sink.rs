//! Output sink shared by every VCF-writing command: plain VCF, BGZF VCF or
//! BCF, to a file or to standard output. `finish` finalizes the stream; a
//! sink dropped early still flushes what it has (BGZF finalizes on Drop) but
//! only `finish` reports errors. Text records pass through [`FloatCanon`],
//! so floats print as bcftools prints them; BCF records are re-rendered by
//! the decoder anyway.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::Compression;

pub use crate::annotate::postproc::{OutputKind, parse_output_type};
use crate::bcf::BcfWriter;
use crate::bgzf::{BgzfWriter, FILE_BUFFER_SIZE, STREAM_BUFFER_SIZE};
use crate::vcf::floats::FloatCanon;

enum Out {
    Plain(BufWriter<Box<dyn Write + Send>>),
    Bgzf(BgzfWriter),
    Bcf(BcfWriter),
}

pub struct VcfSink {
    out: Out,
    canon: Option<FloatCanon>,
    /// Bytes of an unterminated line received through `Write`.
    partial: Vec<u8>,
}

fn is_stdout(path: Option<&Path>) -> bool {
    path.is_none_or(|p| p == Path::new("-"))
}

impl VcfSink {
    /// Open a sink. `path == None` or `-` means standard output. `headers` is
    /// needed for BCF (its dictionary is built from the header) and for the
    /// Float tags a text sink re-renders.
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
        let out = match kind {
            OutputKind::Vcf => Out::Plain(BufWriter::with_capacity(1 << 20, raw)),
            OutputKind::VcfGz(l) => {
                let lvl = level.unwrap_or(l).min(9);
                let buf = if streaming {
                    STREAM_BUFFER_SIZE
                } else {
                    FILE_BUFFER_SIZE
                };
                Out::Bgzf(BgzfWriter::from_writer_buffered(
                    raw,
                    Compression::new(lvl),
                    buf,
                )?)
            }
            OutputKind::Bcf(l) => {
                let lvl = level.unwrap_or(l).min(9);
                Out::Bcf(BcfWriter::from_writer(
                    raw,
                    lvl > 0,
                    lvl,
                    headers,
                    streaming,
                )?)
            }
        };
        let canon = if matches!(out, Out::Bcf(_)) {
            None
        } else {
            Some(FloatCanon::from_headers(headers))
        };
        Ok(Self {
            out,
            canon,
            partial: Vec::new(),
        })
    }

    pub fn is_bcf(&self) -> bool {
        matches!(self.out, Out::Bcf(_))
    }

    /// The BCF writer's dictionary, for raw passthrough checks.
    pub fn bcf_dict(&self) -> Option<&crate::bcf::BcfHeaderDict> {
        match &self.out {
            Out::Bcf(writer) => Some(&writer.dict),
            _ => None,
        }
    }

    /// Copy an encoded BCF record verbatim (BCF sinks only).
    pub fn write_raw_bcf(&mut self, shared: &[u8], indiv: &[u8]) -> Result<()> {
        match &mut self.out {
            Out::Bcf(writer) => writer.write_raw(shared, indiv),
            _ => anyhow::bail!("raw BCF records can only be written to a BCF sink"),
        }
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
        self.put_line(line, true).map_err(Into::into)
    }

    fn put_line(&mut self, line: &str, newline: bool) -> io::Result<()> {
        match &mut self.out {
            Out::Bcf(writer) => {
                if !line.is_empty() && line.as_bytes()[0] != b'#' {
                    writer.write_vcf_line(line).map_err(io::Error::other)?;
                }
                Ok(())
            }
            out => {
                let canon = self.canon.as_mut().and_then(|c| c.canon_line(line));
                let bytes = canon.as_deref().unwrap_or(line).as_bytes();
                let w: &mut dyn Write = match out {
                    Out::Plain(w) => w,
                    Out::Bgzf(w) => w,
                    Out::Bcf(_) => unreachable!(),
                };
                w.write_all(bytes)?;
                if newline {
                    w.write_all(b"\n")?;
                }
                Ok(())
            }
        }
    }

    fn put_bytes_line(&mut self, bytes: &[u8], newline: bool) -> io::Result<()> {
        let s = std::str::from_utf8(bytes)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "non-UTF-8 VCF line"))?;
        let s = s.trim_end_matches(['\r', '\n']);
        if self.is_bcf() && s.is_empty() {
            return Ok(());
        }
        self.put_line(s, newline)
    }

    pub fn finish(mut self) -> Result<()> {
        if !self.partial.is_empty() {
            let partial = std::mem::take(&mut self.partial);
            self.put_bytes_line(&partial, false)?;
        }
        match self.out {
            Out::Plain(mut w) => w.flush().context("flush output"),
            Out::Bgzf(w) => w.finish().context("finalize BGZF output"),
            Out::Bcf(writer) => writer.finish(),
        }
    }
}

impl Write for VcfSink {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let mut rest = buf;
        while let Some(nl) = memchr::memchr(b'\n', rest) {
            let (head, tail) = rest.split_at(nl + 1);
            if self.partial.is_empty() {
                self.put_bytes_line(head, true)?;
            } else {
                let mut v = std::mem::take(&mut self.partial);
                v.extend_from_slice(head);
                self.put_bytes_line(&v, true)?;
            }
            rest = tail;
        }
        self.partial.extend_from_slice(rest);
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        match &mut self.out {
            Out::Plain(w) => w.flush(),
            Out::Bgzf(w) => w.flush(),
            Out::Bcf(_) => Ok(()),
        }
    }
}
