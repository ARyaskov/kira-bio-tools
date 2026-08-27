use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::annotate::ktile::KtileReader;
use crate::bgzf::{
    BgzfReader as NoodlesBgzfReader, MtBgzfReader as ParallelBgzfReader, VirtualPosition,
};
use crate::util::{VcfFormat, chr_name_to_id, detect_format};
use crate::vcf::parser::extract_contig_id;
use crate::vcf::structs::{Result, VcfError, VcfRecord};

pub enum UnifiedVcfReader {
    Plain(PlainReader),
    Bgzf(BgzfReader),
    BgzfIndexing(BgzfIndexingReader),
    Ktile(KtileSourceReader),
    Bcf(BcfSourceReader),
}

pub struct BcfSourceReader {
    inner: crate::bcf::BcfReader,
    header_lines: Vec<String>,
    header_emitted: bool,
    cursor_in_header: usize,
}

impl BcfSourceReader {
    pub fn open(path: &Path) -> Result<Self> {
        let inner = crate::bcf::BcfReader::open(path)
            .map_err(|e| VcfError::Io(std::io::Error::other(format!("bcf open: {e}"))))?;
        let header_lines = inner.header_lines.clone();
        Ok(Self { inner, header_lines, header_emitted: false, cursor_in_header: 0 })
    }
    pub fn header(&self) -> &[String] { &self.header_lines }
    pub fn read_line(&mut self) -> Result<Option<String>> {
        if !self.header_emitted {
            if self.cursor_in_header < self.header_lines.len() {
                let l = self.header_lines[self.cursor_in_header].clone();
                self.cursor_in_header += 1;
                return Ok(Some(l));
            }
            self.header_emitted = true;
        }
        self.inner.read_record_line().map_err(|e| VcfError::Io(std::io::Error::other(format!("bcf read: {e}"))))
    }
}

impl UnifiedVcfReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();

        // `.ktile` sidecars are extension-dispatched.
        if path
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| e.eq_ignore_ascii_case("ktile"))
        {
            return Self::open_ktile(path);
        }
        if path
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| e.eq_ignore_ascii_case("bcf"))
        {
            return Ok(Self::Bcf(BcfSourceReader::open(path)?));
        }

        let format = detect_format(path)?;
        match format {
            VcfFormat::Plain => Ok(Self::Plain(PlainReader::open(path)?)),
            VcfFormat::Bgzf => Ok(Self::Bgzf(BgzfReader::open(path)?)),
            VcfFormat::Gzip => Err(VcfError::InvalidFormat),
        }
    }

    /// Opens an annotation source from a pre-built `.ktile` sidecar.
    pub fn open_ktile<P: AsRef<Path>>(path: P) -> Result<Self> {
        Ok(Self::Ktile(KtileSourceReader::open(path.as_ref())?))
    }

    pub fn open_for_indexing<P: AsRef<Path>>(path: P) -> Result<Self> {
        Ok(Self::BgzfIndexing(BgzfIndexingReader::open(path)?))
    }

    pub fn read_record(&mut self) -> Result<Option<VcfRecord>> {
        match self {
            Self::Plain(r) => r.read_record(),
            Self::Bgzf(r) => r.read_record(),
            Self::BgzfIndexing(r) => r.read_record(),
            Self::Ktile(r) => r.read_record(),
            Self::Bcf(_) => Err(VcfError::InvalidFormat),
        }
    }

    pub fn read_line(&mut self) -> Result<Option<String>> {
        match self {
            Self::Plain(r) => r.read_line(),
            Self::Bgzf(r) => r.read_line(),
            Self::BgzfIndexing(r) => r.read_line(),
            Self::Ktile(r) => r.read_line(),
            Self::Bcf(r) => {
                while let Some(l) = r.read_line()? {
                    if l.starts_with('#') { continue; }
                    return Ok(Some(l));
                }
                Ok(None)
            }
        }
    }

    /// Read a line and, if the source has pre-parsed metadata,
    /// the (chr_id, pos) pair. Returns `(line, None)` for sources without
    /// side-channel metadata; `Ktile` returns `(line, Some((chr_id, pos)))`.
    pub fn read_line_with_meta(&mut self) -> Result<Option<(String, Option<(u32, u32)>)>> {
        match self {
            Self::Ktile(r) => r.read_line_with_meta(),
            other => Ok(other.read_line()?.map(|line| (line, None))),
        }
    }

    /// Zero-copy fast path: appends the next line directly into `batch`.
    /// Returns `Ok(true)` on success, `Ok(false)` at EOF.
    pub fn read_line_into_batch(
        &mut self,
        batch: &mut crate::annotate::cpu_v2::ReadBatch,
    ) -> Result<bool> {
        match self {
            Self::Ktile(r) => r.read_line_into_batch(batch),
            other => match other.read_line()? {
                Some(line) => {
                    batch.push_line(&line);
                    Ok(true)
                }
                None => Ok(false),
            },
        }
    }

    pub fn header(&self) -> Result<Vec<String>> {
        match self {
            Self::Plain(r) => Ok(r.headers.clone()),
            Self::Bgzf(r) => Ok(r.headers.clone()),
            Self::BgzfIndexing(r) => Ok(r.headers.clone()),
            Self::Ktile(r) => Ok(r.headers.clone()),
            Self::Bcf(r) => Ok(r.header_lines.clone()),
        }
    }

    pub fn contigs(&self) -> Vec<String> {
        match self {
            Self::Plain(r) => r.contigs.clone(),
            Self::Bgzf(r) => r.contigs.clone(),
            Self::BgzfIndexing(r) => r.contigs.clone(),
            Self::Ktile(r) => r.contigs.clone(),
            Self::Bcf(r) => r.inner.dict.contigs.clone(),
        }
    }

    pub fn reference_sequences(&self) -> Result<Vec<String>> {
        Ok(self.contigs())
    }

    pub fn virtual_position(&self) -> Option<VirtualPosition> {
        match self {
            Self::Plain(_) => None,
            Self::Bgzf(_) => None,
            Self::BgzfIndexing(r) => Some(r.vpos),
            Self::Ktile(_) => None,
            Self::Bcf(_) => None,
        }
    }

    pub fn next_record_with_vpos(&mut self) -> Result<Option<(VcfRecord, VirtualPosition)>> {
        match self {
            Self::BgzfIndexing(r) => r.next_record_with_vpos(),
            _ => Err(VcfError::InvalidFormat),
        }
    }
}

pub struct PlainReader {
    reader: BufReader<File>,
    buffer: String,
    contigs: Vec<String>,
    headers: Vec<String>,
    offset: u64,
    first_data_line: Option<String>,
}

impl PlainReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);
        let mut buffer = String::new();
        let mut contigs = Vec::new();
        let mut headers = Vec::new();
        let mut offset = 0u64;
        let mut first_data_line = None;

        loop {
            buffer.clear();
            let n = reader.read_line(&mut buffer)?;
            if n == 0 {
                break;
            }
            offset += n as u64;
            if !buffer.starts_with('#') {
                first_data_line = Some(buffer.trim_end().to_string());
                break;
            }
            headers.push(buffer.trim_end().to_string());
            if buffer.starts_with("##contig=") {
                if let Some(id) = extract_contig_id(&buffer) {
                    contigs.push(id);
                }
            }
        }

        Ok(Self {
            reader,
            buffer: String::new(),
            contigs,
            headers,
            offset,
            first_data_line,
        })
    }

    pub fn read_record(&mut self) -> Result<Option<VcfRecord>> {
        if let Some(line) = self.first_data_line.take() {
            return parse_vcf_record(&line, self.offset);
        }

        let start_offset = self.offset;
        self.buffer.clear();
        let n = self.reader.read_line(&mut self.buffer)?;
        if n == 0 {
            return Ok(None);
        }
        self.offset += n as u64;

        if self.buffer.starts_with('#') {
            return self.read_record();
        }

        parse_vcf_record(&self.buffer, start_offset)
    }

    pub fn read_line(&mut self) -> Result<Option<String>> {
        if let Some(line) = self.first_data_line.take() {
            return Ok(Some(line));
        }

        self.buffer.clear();
        let n = self.reader.read_line(&mut self.buffer)?;
        if n == 0 {
            return Ok(None);
        }
        self.offset += n as u64;
        Ok(Some(self.buffer.trim_end().to_string()))
    }
}

/// Streaming BGZF reader using the multithreaded inflater pool. For
/// tabix-indexing use [`BgzfIndexingReader`].
pub struct BgzfReader {
    reader: ParallelBgzfReader<File>,
    buffer: String,
    contigs: Vec<String>,
    headers: Vec<String>,
    first_data_line: Option<String>,
}

impl BgzfReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let mut reader = ParallelBgzfReader::open(path)?;
        let mut buffer = String::new();
        let mut contigs = Vec::new();
        let mut headers = Vec::new();
        let mut first_data_line = None;

        loop {
            buffer.clear();
            let n = reader.read_line(&mut buffer)?;
            if n == 0 {
                break;
            }
            if !buffer.starts_with('#') {
                first_data_line = Some(buffer.trim_end().to_string());
                break;
            }
            headers.push(buffer.trim_end().to_string());
            if buffer.starts_with("##contig=") {
                if let Some(id) = extract_contig_id(&buffer) {
                    contigs.push(id);
                }
            }
        }

        Ok(Self {
            reader,
            buffer: String::new(),
            contigs,
            headers,
            first_data_line,
        })
    }

    pub fn read_record(&mut self) -> Result<Option<VcfRecord>> {
        if let Some(line) = self.first_data_line.take() {
            return parse_vcf_record(&line, 0);
        }

        self.buffer.clear();
        let n = self.reader.read_line(&mut self.buffer)?;
        if n == 0 {
            return Ok(None);
        }

        if self.buffer.starts_with('#') {
            return self.read_record();
        }

        parse_vcf_record(&self.buffer, 0)
    }

    pub fn read_line(&mut self) -> Result<Option<String>> {
        if let Some(line) = self.first_data_line.take() {
            return Ok(Some(line));
        }

        self.buffer.clear();
        let n = self.reader.read_line(&mut self.buffer)?;
        if n == 0 {
            return Ok(None);
        }
        Ok(Some(self.buffer.trim_end().to_string()))
    }
}

pub struct BgzfIndexingReader {
    reader: NoodlesBgzfReader<File>,
    buffer: String,
    contigs: Vec<String>,
    headers: Vec<String>,
    vpos: VirtualPosition,
    first_data_line: Option<String>,
}

impl BgzfIndexingReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let mut reader = NoodlesBgzfReader::open(path)?;
        let mut buffer = String::new();
        let mut contigs = Vec::new();
        let mut headers = Vec::new();
        let mut first_data_line = None;
        // Virtual position at the START of the first data record (captured before that line is read,
        // since the header scan consumes it). This is the record's chunk-start offset for indexing.
        let vpos = loop {
            let vpos_before = reader.virtual_position();
            buffer.clear();
            let n = reader.read_line(&mut buffer)?;
            if n == 0 {
                break vpos_before;
            }
            if !buffer.starts_with('#') {
                first_data_line = Some(buffer.trim_end().to_string());
                break vpos_before;
            }
            headers.push(buffer.trim_end().to_string());
            if buffer.starts_with("##contig=") {
                if let Some(id) = extract_contig_id(&buffer) {
                    contigs.push(id);
                }
            }
        };

        Ok(Self {
            reader,
            buffer: String::new(),
            contigs,
            headers,
            vpos,
            first_data_line,
        })
    }

    pub fn read_record(&mut self) -> Result<Option<VcfRecord>> {
        if let Some(line) = self.first_data_line.take() {
            return parse_vcf_record(&line, self.vpos.as_u64());
        }

        self.vpos = self.reader.virtual_position();
        self.buffer.clear();
        let n = self.reader.read_line(&mut self.buffer)?;
        if n == 0 {
            return Ok(None);
        }

        if self.buffer.starts_with('#') {
            return self.read_record();
        }

        parse_vcf_record(&self.buffer, self.vpos.as_u64())
    }

    pub fn read_line(&mut self) -> Result<Option<String>> {
        if let Some(line) = self.first_data_line.take() {
            return Ok(Some(line));
        }

        self.buffer.clear();
        let n = self.reader.read_line(&mut self.buffer)?;
        if n == 0 {
            return Ok(None);
        }
        Ok(Some(self.buffer.trim_end().to_string()))
    }

    pub fn current_vpos(&self) -> VirtualPosition {
        self.vpos
    }

    pub fn next_record_with_vpos(&mut self) -> Result<Option<(VcfRecord, VirtualPosition)>> {
        // The record carries its own start vpos (`offset`); here we return the position AFTER the
        // record so the indexer can form a proper [start, end) chunk (not a zero-length one).
        match self.read_record()? {
            Some(rec) => {
                let after = self.reader.virtual_position();
                Ok(Some((rec, after)))
            }
            None => Ok(None),
        }
    }
}

fn parse_vcf_record(line: &str, offset: u64) -> Result<Option<VcfRecord>> {
    let cols: Vec<&str> = line.trim_end().split('\t').collect();
    if cols.len() < 8 {
        return Ok(None);
    }

    let chrom = cols[0];
    let pos = cols[1]
        .parse::<u32>()
        .map_err(|_| VcfError::InvalidFormat)?;

    let chr_id = chr_name_to_id(chrom).unwrap_or(0);

    let format = if cols.len() > 8 {
        Some(cols[8].to_string())
    } else {
        None
    };

    let samples = if cols.len() > 9 {
        cols[9..].iter().map(|s| s.to_string()).collect()
    } else {
        Vec::new()
    };

    Ok(Some(VcfRecord {
        chrom: chrom.to_string(),
        pos,
        id: cols[2].to_string(),
        ref_allele: cols[3].to_string(),
        alt: cols[4].to_string(),
        qual: cols[5].to_string(),
        filter: cols[6].to_string(),
        info: cols[7].to_string(),
        format,
        samples,
        chr_id,
        position: pos,
        offset,
    }))
}

/// Streaming reader over a `.ktile` sidecar.
pub struct KtileSourceReader {
    reader: KtileReader,
    cursor: usize,
    pub(crate) headers: Vec<String>,
    pub(crate) contigs: Vec<String>,
}

impl KtileSourceReader {
    pub fn open(path: &Path) -> Result<Self> {
        let reader = KtileReader::open(path).map_err(|_| VcfError::InvalidFormat)?;
        let mut headers: Vec<String> = Vec::new();
        let mut contigs: Vec<String> = Vec::new();
        for line in reader.headers_block().split('\n') {
            if line.is_empty() {
                continue;
            }
            headers.push(line.to_string());
            if line.starts_with("##contig=") {
                if let Some(id) = extract_contig_id(line) {
                    contigs.push(id);
                }
            }
        }
        Ok(Self {
            reader,
            cursor: 0,
            headers,
            contigs,
        })
    }

    pub fn read_line(&mut self) -> Result<Option<String>> {
        if self.cursor >= self.reader.n_records() {
            return Ok(None);
        }
        let line = self.reader.line_owned(self.cursor);
        self.cursor += 1;
        Ok(Some(line))
    }

    /// Like [`read_line`] but also returns the pre-stored (chr_id, pos).
    pub fn read_line_with_meta(
        &mut self,
    ) -> Result<Option<(String, Option<(u32, u32)>)>> {
        if self.cursor >= self.reader.n_records() {
            return Ok(None);
        }
        let i = self.cursor;
        let line = self.reader.line_owned(i);
        let chr_id = self.reader.chr_id(i);
        let pos = self.reader.position(i);
        self.cursor += 1;
        Ok(Some((line, Some((chr_id, pos)))))
    }

    /// Zero-copy fast path — appends the next line directly into `batch`.
    pub fn read_line_into_batch(
        &mut self,
        batch: &mut crate::annotate::cpu_v2::ReadBatch,
    ) -> Result<bool> {
        if self.cursor >= self.reader.n_records() {
            return Ok(false);
        }
        self.reader.push_line_into_batch(self.cursor, batch);
        self.cursor += 1;
        Ok(true)
    }

    pub fn read_record(&mut self) -> Result<Option<VcfRecord>> {
        if self.cursor >= self.reader.n_records() {
            return Ok(None);
        }
        let i = self.cursor;
        self.cursor += 1;
        let line = self.reader.line_owned(i);
        let mut rec = parse_vcf_record(&line, 0)?;
        if let Some(ref mut r) = rec {
            r.chr_id = self.reader.chr_id(i) as u8;
        }
        Ok(rec)
    }
}
