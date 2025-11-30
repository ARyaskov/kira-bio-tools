use std::fs::File;
use std::io::{self, BufRead, BufReader, Seek, SeekFrom};
use std::path::Path;

use flate2::read::GzDecoder;
use memmap2::Mmap;
use rayon::prelude::*;
use thiserror::Error;

use crate::bgzf::{BgzfLineReader, BgzfReader, VirtualPosition};
use crate::util::{detect_format, parse_vcf_line_fast, ChrId, GenomicKey, VcfFormat};

#[derive(Debug, Error)]
pub enum VcfError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
    #[error("BGZF error: {0}")]
    Bgzf(#[from] crate::bgzf::BgzfError),
    #[error("Invalid VCF format: {0}")]
    InvalidFormat(String),
    #[error("Parse error at line {line}: {message}")]
    ParseError { line: usize, message: String },
}

pub type Result<T> = std::result::Result<T, VcfError>;

#[derive(Debug, Clone)]
pub struct VcfRecord {
    pub chr_id: ChrId,
    pub position: u32,
    pub offset: u64,
}

impl VcfRecord {
    pub fn key(&self) -> GenomicKey {
        GenomicKey::new(self.chr_id, self.position)
    }
}

pub enum VcfReader {
    Plain(PlainVcfReader),
    Gzip(GzipVcfReader),
    Bgzf(BgzfVcfReader),
}

impl VcfReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let format = detect_format(path)?;

        match format {
            VcfFormat::Plain => Ok(VcfReader::Plain(PlainVcfReader::open(path)?)),
            VcfFormat::Gzip => Ok(VcfReader::Gzip(GzipVcfReader::open(path)?)),
            VcfFormat::Bgzf => Ok(VcfReader::Bgzf(BgzfVcfReader::open(path)?)),
        }
    }

    pub fn is_bgzf(&self) -> bool {
        matches!(self, VcfReader::Bgzf(_))
    }

    pub fn header(&mut self) -> Result<Vec<String>> {
        match self {
            VcfReader::Plain(r) => r.header(),
            VcfReader::Gzip(r) => r.header(),
            VcfReader::Bgzf(r) => r.header(),
        }
    }

    pub fn records(&mut self) -> RecordIterator<'_> {
        RecordIterator { reader: self }
    }

    pub fn next_record(&mut self) -> Result<Option<VcfRecord>> {
        match self {
            VcfReader::Plain(r) => r.next_record(),
            VcfReader::Gzip(r) => r.next_record(),
            VcfReader::Bgzf(r) => r.next_record(),
        }
    }

    pub fn reference_sequences(&self) -> &[String] {
        match self {
            VcfReader::Plain(r) => &r.contigs,
            VcfReader::Gzip(r) => &r.contigs,
            VcfReader::Bgzf(r) => &r.contigs,
        }
    }
}

pub struct RecordIterator<'a> {
    reader: &'a mut VcfReader,
}

impl<'a> Iterator for RecordIterator<'a> {
    type Item = Result<VcfRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.next_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

pub struct PlainVcfReader {
    reader: BufReader<File>,
    buf: String,
    offset: u64,
    contigs: Vec<String>,
    header_parsed: bool,
}

impl PlainVcfReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        Ok(Self {
            reader: BufReader::with_capacity(8 * 1024 * 1024, file),
            buf: String::with_capacity(4096),
            offset: 0,
            contigs: Vec::new(),
            header_parsed: false,
        })
    }

    pub fn header(&mut self) -> Result<Vec<String>> {
        let mut headers = Vec::new();

        loop {
            self.buf.clear();
            let bytes = self.reader.read_line(&mut self.buf)?;
            if bytes == 0 {
                break;
            }

            if self.buf.starts_with('#') {
                if self.buf.starts_with("##contig=") {
                    if let Some(id) = extract_contig_id(&self.buf) {
                        self.contigs.push(id);
                    }
                }
                self.offset += bytes as u64;
                headers.push(self.buf.trim_end().to_string());

                if self.buf.starts_with("#CHROM") {
                    break;
                }
            } else {
                break;
            }
        }

        self.header_parsed = true;
        Ok(headers)
    }

    pub fn next_record(&mut self) -> Result<Option<VcfRecord>> {
        if !self.header_parsed {
            self.header()?;
        }

        loop {
            self.buf.clear();
            let start_offset = self.offset;
            let bytes = self.reader.read_line(&mut self.buf)?;

            if bytes == 0 {
                return Ok(None);
            }

            self.offset += bytes as u64;

            if self.buf.starts_with('#') || self.buf.trim().is_empty() {
                continue;
            }

            if let Some((chr_id, position)) = parse_vcf_line_fast(self.buf.as_bytes()) {
                return Ok(Some(VcfRecord {
                    chr_id,
                    position,
                    offset: start_offset,
                }));
            }
        }
    }
}

pub struct GzipVcfReader {
    reader: BufReader<GzDecoder<File>>,
    buf: String,
    offset: u64,
    contigs: Vec<String>,
    header_parsed: bool,
}

impl GzipVcfReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        let decoder = GzDecoder::new(file);
        Ok(Self {
            reader: BufReader::with_capacity(8 * 1024 * 1024, decoder),
            buf: String::with_capacity(4096),
            offset: 0,
            contigs: Vec::new(),
            header_parsed: false,
        })
    }

    pub fn header(&mut self) -> Result<Vec<String>> {
        let mut headers = Vec::new();

        loop {
            self.buf.clear();
            let bytes = self.reader.read_line(&mut self.buf)?;
            if bytes == 0 {
                break;
            }

            if self.buf.starts_with('#') {
                if self.buf.starts_with("##contig=") {
                    if let Some(id) = extract_contig_id(&self.buf) {
                        self.contigs.push(id);
                    }
                }
                self.offset += bytes as u64;
                headers.push(self.buf.trim_end().to_string());

                if self.buf.starts_with("#CHROM") {
                    break;
                }
            } else {
                break;
            }
        }

        self.header_parsed = true;
        Ok(headers)
    }

    pub fn next_record(&mut self) -> Result<Option<VcfRecord>> {
        if !self.header_parsed {
            self.header()?;
        }

        loop {
            self.buf.clear();
            let start_offset = self.offset;
            let bytes = self.reader.read_line(&mut self.buf)?;

            if bytes == 0 {
                return Ok(None);
            }

            self.offset += bytes as u64;

            if self.buf.starts_with('#') || self.buf.trim().is_empty() {
                continue;
            }

            if let Some((chr_id, position)) = parse_vcf_line_fast(self.buf.as_bytes()) {
                return Ok(Some(VcfRecord {
                    chr_id,
                    position,
                    offset: start_offset,
                }));
            }
        }
    }
}

pub struct BgzfVcfReader {
    reader: BgzfLineReader<File>,
    contigs: Vec<String>,
    header_parsed: bool,
}

impl BgzfVcfReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let bgzf_reader = BgzfReader::open(path)?;
        Ok(Self {
            reader: BgzfLineReader::new(bgzf_reader),
            contigs: Vec::new(),
            header_parsed: false,
        })
    }

    pub fn header(&mut self) -> Result<Vec<String>> {
        let mut headers = Vec::new();

        loop {
            match self.reader.read_line()? {
                Some((line, _)) => {
                    if line.starts_with('#') {
                        if line.starts_with("##contig=") {
                            if let Some(id) = extract_contig_id(line) {
                                self.contigs.push(id);
                            }
                        }
                        headers.push(line.to_string());

                        if line.starts_with("#CHROM") {
                            break;
                        }
                    } else {
                        break;
                    }
                }
                None => break,
            }
        }

        self.header_parsed = true;
        Ok(headers)
    }

    pub fn next_record(&mut self) -> Result<Option<VcfRecord>> {
        if !self.header_parsed {
            self.header()?;
        }

        loop {
            match self.reader.read_line()? {
                Some((line, vpos)) => {
                    if line.starts_with('#') || line.is_empty() {
                        continue;
                    }

                    if let Some((chr_id, position)) = parse_vcf_line_fast(line.as_bytes()) {
                        return Ok(Some(VcfRecord {
                            chr_id,
                            position,
                            offset: vpos.as_u64(),
                        }));
                    }
                }
                None => return Ok(None),
            }
        }
    }

    pub fn virtual_position(&self) -> VirtualPosition {
        self.reader.virtual_position()
    }
}

fn extract_contig_id(line: &str) -> Option<String> {
    let start = line.find("ID=")? + 3;
    let rest = &line[start..];
    let end = rest.find(|c| c == ',' || c == '>')?;
    Some(rest[..end].to_string())
}

pub struct MmapVcfParser<'a> {
    data: &'a [u8],
}

impl<'a> MmapVcfParser<'a> {
    pub fn new(mmap: &'a Mmap) -> Self {
        Self { data: mmap }
    }

    pub fn parse_parallel(&self, num_threads: usize) -> Vec<VcfRecord> {
        let chunk_size = self.data.len() / num_threads;
        let chunks: Vec<_> = (0..num_threads)
            .map(|i| {
                let start = i * chunk_size;
                let end = if i == num_threads - 1 {
                    self.data.len()
                } else {
                    (i + 1) * chunk_size
                };
                (start, end)
            })
            .collect();

        chunks
            .into_par_iter()
            .flat_map(|(start, end)| {
                let adjusted_start = if start == 0 {
                    0
                } else {
                    self.data[start..]
                        .iter()
                        .position(|&b| b == b'\n')
                        .map(|p| start + p + 1)
                        .unwrap_or(end)
                };

                let adjusted_end = if end >= self.data.len() {
                    self.data.len()
                } else {
                    self.data[..end]
                        .iter()
                        .rposition(|&b| b == b'\n')
                        .map(|p| p + 1)
                        .unwrap_or(end)
                };

                self.parse_chunk(adjusted_start, adjusted_end)
            })
            .collect()
    }

    fn parse_chunk(&self, start: usize, end: usize) -> Vec<VcfRecord> {
        let mut records = Vec::new();
        let mut pos = start;

        while pos < end {
            let line_end = self.data[pos..end]
                .iter()
                .position(|&b| b == b'\n')
                .map(|p| pos + p)
                .unwrap_or(end);

            let line = &self.data[pos..line_end];

            if !line.is_empty() && line[0] != b'#' {
                if let Some((chr_id, position)) = parse_vcf_line_fast(line) {
                    records.push(VcfRecord {
                        chr_id,
                        position,
                        offset: pos as u64,
                    });
                }
            }

            pos = line_end + 1;
        }

        records
    }
}

pub fn fetch_line<P: AsRef<Path>>(path: P, offset: u64) -> Result<String> {
    let path = path.as_ref();
    let format = detect_format(path)?;

    match format {
        VcfFormat::Plain => {
            let mut file = File::open(path)?;
            file.seek(SeekFrom::Start(offset))?;
            let mut reader = BufReader::new(file);
            let mut line = String::new();
            reader.read_line(&mut line)?;
            Ok(line.trim_end().to_string())
        }
        VcfFormat::Gzip => {
            Err(VcfError::InvalidFormat(
                "Cannot seek in gzip file. Use BGZF compression.".into(),
            ))
        }
        VcfFormat::Bgzf => {
            let file = File::open(path)?;
            let mut bgzf_reader = BgzfReader::new(file);
            bgzf_reader.seek(VirtualPosition::from_u64(offset))?;

            let mut line = String::new();
            bgzf_reader.read_line(&mut line)?;

            if line.ends_with('\n') {
                line.pop();
            }
            if line.ends_with('\r') {
                line.pop();
            }
            Ok(line)
        }
    }
}