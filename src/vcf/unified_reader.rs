use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use crate::annotate::ktile::KtileReader;
use crate::bcf::BCF_MAGIC;
use crate::bgzf::{
    BgzfReader as NoodlesBgzfReader, MtBgzfReader, VirtualPosition, is_bgzf_header, is_gzip_header,
};
use crate::vcf::header::ContigDict;
use crate::vcf::structs::{Result, VcfError, VcfRecord};

/// Reader over any VCF-like source: plain/gzip/BGZF VCF text, BCF (detected by
/// magic, compressed or not), `.ktile` sidecars, files or standard input.
pub enum UnifiedVcfReader {
    Text(TextReader),
    BgzfIndexing(BgzfIndexingReader),
    Ktile(KtileSourceReader),
    Bcf(BcfSourceReader),
}

pub struct BcfSourceReader {
    inner: crate::bcf::BcfReader,
    header_lines: Vec<String>,
    contigs: ContigDict,
    header_emitted: bool,
    cursor_in_header: usize,
}

impl BcfSourceReader {
    pub fn open(path: &Path) -> Result<Self> {
        let inner = crate::bcf::BcfReader::open(path)
            .map_err(|e| VcfError::Io(io::Error::other(format!("bcf open: {e}"))))?;
        Ok(Self::from_reader(inner))
    }

    pub fn from_bufread(r: Box<dyn BufRead + Send>) -> Result<Self> {
        let inner = crate::bcf::BcfReader::from_bufread(r)
            .map_err(|e| VcfError::Io(io::Error::other(format!("bcf open: {e}"))))?;
        Ok(Self::from_reader(inner))
    }

    fn from_reader(inner: crate::bcf::BcfReader) -> Self {
        let header_lines = inner.header_lines.clone();
        // Contig ids follow BCF `rid` order.
        let mut contigs = ContigDict::new();
        let max_rid = inner.dict.contig_idx.values().copied().max();
        if let Some(max) = max_rid {
            for rid in 0..=max {
                match inner.dict.contig_name(rid) {
                    Some(n) => {
                        contigs.insert(n);
                    }
                    None => {
                        contigs.insert(&format!("__unnamed_rid_{rid}"));
                    }
                }
            }
        }
        Self { inner, header_lines, contigs, header_emitted: false, cursor_in_header: 0 }
    }

    pub fn header(&self) -> &[String] { &self.header_lines }

    /// Header lines first (once), then records.
    pub fn read_line(&mut self) -> Result<Option<String>> {
        if !self.header_emitted {
            if self.cursor_in_header < self.header_lines.len() {
                let l = self.header_lines[self.cursor_in_header].clone();
                self.cursor_in_header += 1;
                return Ok(Some(l));
            }
            self.header_emitted = true;
        }
        self.read_data_line()
    }

    pub fn read_data_line(&mut self) -> Result<Option<String>> {
        self.inner
            .read_record_line()
            .map_err(|e| VcfError::Io(io::Error::other(format!("bcf read: {e}"))))
    }

    pub fn read_record(&mut self) -> Result<Option<VcfRecord>> {
        match self.read_data_line()? {
            Some(line) => parse_vcf_record(&line, 0, &mut self.contigs),
            None => Ok(None),
        }
    }
}

impl UnifiedVcfReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        if path == Path::new("-") {
            return Self::from_stream(Box::new(io::stdin()));
        }
        if path
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| e.eq_ignore_ascii_case("ktile"))
        {
            return Self::open_ktile(path);
        }
        let file = File::open(path)?;
        Self::from_stream(Box::new(file))
    }

    /// Open from any byte stream: format is sniffed from the first bytes
    /// (BGZF, plain gzip or plain text; BCF by magic after decompression).
    pub fn from_stream(r: Box<dyn Read + Send>) -> Result<Self> {
        let (head, r) = peek_prefix(r, 18)?;
        let mut buf: Box<dyn BufRead + Send> = if is_bgzf_header(&head) {
            Box::new(MtBgzfReader::new(r))
        } else if is_gzip_header(&head) {
            Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(r)))
        } else {
            Box::new(BufReader::with_capacity(1 << 20, r))
        };
        let is_bcf = {
            let b = buf.fill_buf()?;
            b.len() >= BCF_MAGIC.len() && &b[..BCF_MAGIC.len()] == BCF_MAGIC
        };
        if is_bcf {
            Ok(Self::Bcf(BcfSourceReader::from_bufread(buf)?))
        } else {
            Ok(Self::Text(TextReader::new(buf)?))
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
            Self::Text(r) => r.read_record(),
            Self::BgzfIndexing(r) => r.read_record(),
            Self::Ktile(r) => r.read_record(),
            Self::Bcf(r) => r.read_record(),
        }
    }

    /// Next data line (header lines are never returned).
    pub fn read_line(&mut self) -> Result<Option<String>> {
        match self {
            Self::Text(r) => r.read_line(),
            Self::BgzfIndexing(r) => r.read_line(),
            Self::Ktile(r) => r.read_line(),
            Self::Bcf(r) => r.read_data_line(),
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
            Self::Text(r) => Ok(r.headers.clone()),
            Self::BgzfIndexing(r) => Ok(r.headers.clone()),
            Self::Ktile(r) => Ok(r.headers.clone()),
            Self::Bcf(r) => Ok(r.header_lines.clone()),
        }
    }

    /// Contig dictionary: header contigs plus any contig met in the data so far.
    pub fn contigs(&self) -> &ContigDict {
        match self {
            Self::Text(r) => &r.contigs,
            Self::BgzfIndexing(r) => &r.contigs,
            Self::Ktile(r) => &r.contigs,
            Self::Bcf(r) => &r.contigs,
        }
    }

    pub fn reference_sequences(&self) -> Result<Vec<String>> {
        Ok(self.contigs().names().to_vec())
    }

    pub fn virtual_position(&self) -> Option<VirtualPosition> {
        match self {
            Self::BgzfIndexing(r) => Some(r.vpos),
            _ => None,
        }
    }

    pub fn next_record_with_vpos(&mut self) -> Result<Option<(VcfRecord, VirtualPosition)>> {
        match self {
            Self::BgzfIndexing(r) => r.next_record_with_vpos(),
            _ => Err(VcfError::InvalidFormat),
        }
    }

    /// Next data line with its start/end virtual positions (BGZF indexing source only).
    pub fn next_line_with_vpos(&mut self) -> Result<Option<(String, VirtualPosition, VirtualPosition)>> {
        match self {
            Self::BgzfIndexing(r) => r.next_line_with_vpos(),
            _ => Err(VcfError::InvalidFormat),
        }
    }
}

fn peek_prefix(mut r: Box<dyn Read + Send>, n: usize) -> io::Result<(Vec<u8>, Box<dyn Read + Send>)> {
    let mut head = vec![0u8; n];
    let mut got = 0usize;
    while got < n {
        let k = r.read(&mut head[got..])?;
        if k == 0 {
            break;
        }
        got += k;
    }
    head.truncate(got);
    let chained: Box<dyn Read + Send> = Box::new(Cursor::new(head.clone()).chain(r));
    Ok((head, chained))
}

fn read_header_lines<R: BufRead + ?Sized>(
    reader: &mut R,
    buffer: &mut String,
    offset: &mut u64,
) -> Result<(Vec<String>, ContigDict, Option<String>)> {
    let mut headers = Vec::new();
    let mut first_data_line = None;
    loop {
        buffer.clear();
        let n = reader.read_line(buffer)?;
        if n == 0 {
            break;
        }
        *offset += n as u64;
        if !buffer.starts_with('#') {
            first_data_line = Some(buffer.trim_end_matches(['\r', '\n']).to_string());
            break;
        }
        headers.push(buffer.trim_end_matches(['\r', '\n']).to_string());
    }
    let contigs = ContigDict::from_header_lines(headers.iter().map(String::as_str));
    Ok((headers, contigs, first_data_line))
}

/// Streaming text VCF (plain, gzip or BGZF; file or stdin).
pub struct TextReader {
    reader: Box<dyn BufRead + Send>,
    buffer: String,
    headers: Vec<String>,
    contigs: ContigDict,
    offset: u64,
    first_data_line: Option<String>,
}

impl TextReader {
    pub fn new(mut reader: Box<dyn BufRead + Send>) -> Result<Self> {
        let mut buffer = String::new();
        let mut offset = 0u64;
        let (headers, contigs, first_data_line) = read_header_lines(&mut *reader, &mut buffer, &mut offset)?;
        Ok(Self { reader, buffer: String::new(), headers, contigs, offset, first_data_line })
    }

    pub fn read_record(&mut self) -> Result<Option<VcfRecord>> {
        if let Some(line) = self.first_data_line.take() {
            let start = self.offset - (line.len() as u64 + 1);
            return parse_vcf_record(&line, start, &mut self.contigs);
        }
        loop {
            let start_offset = self.offset;
            self.buffer.clear();
            let n = self.reader.read_line(&mut self.buffer)?;
            if n == 0 {
                return Ok(None);
            }
            self.offset += n as u64;
            if self.buffer.starts_with('#') || self.buffer.trim_end().is_empty() {
                continue;
            }
            return parse_vcf_record(&self.buffer, start_offset, &mut self.contigs);
        }
    }

    pub fn read_line(&mut self) -> Result<Option<String>> {
        if let Some(line) = self.first_data_line.take() {
            return Ok(Some(line));
        }
        loop {
            self.buffer.clear();
            let n = self.reader.read_line(&mut self.buffer)?;
            if n == 0 {
                return Ok(None);
            }
            self.offset += n as u64;
            if self.buffer.starts_with('#') {
                continue;
            }
            return Ok(Some(self.buffer.trim_end_matches(['\r', '\n']).to_string()));
        }
    }
}

/// Single-threaded BGZF reader that tracks virtual positions, for index building.
pub struct BgzfIndexingReader {
    reader: NoodlesBgzfReader<File>,
    buffer: String,
    headers: Vec<String>,
    contigs: ContigDict,
    vpos: VirtualPosition,
    first_data_line: Option<String>,
}

impl BgzfIndexingReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let mut reader = NoodlesBgzfReader::open(path)?;
        let mut buffer = String::new();
        let mut headers = Vec::new();
        let mut first_data_line = None;
        // Virtual position at the START of the first data record, captured before that line is
        // read, since the header scan consumes it.
        let vpos = loop {
            let vpos_before = reader.virtual_position();
            buffer.clear();
            let n = reader.read_line(&mut buffer)?;
            if n == 0 {
                break vpos_before;
            }
            if !buffer.starts_with('#') {
                first_data_line = Some(buffer.trim_end_matches(['\r', '\n']).to_string());
                break vpos_before;
            }
            headers.push(buffer.trim_end_matches(['\r', '\n']).to_string());
        };
        let contigs = ContigDict::from_header_lines(headers.iter().map(String::as_str));
        Ok(Self { reader, buffer: String::new(), headers, contigs, vpos, first_data_line })
    }

    pub fn read_record(&mut self) -> Result<Option<VcfRecord>> {
        if let Some(line) = self.first_data_line.take() {
            return parse_vcf_record(&line, self.vpos.as_u64(), &mut self.contigs);
        }
        loop {
            self.vpos = self.reader.virtual_position();
            self.buffer.clear();
            let n = self.reader.read_line(&mut self.buffer)?;
            if n == 0 {
                return Ok(None);
            }
            if self.buffer.starts_with('#') || self.buffer.trim_end().is_empty() {
                continue;
            }
            return parse_vcf_record(&self.buffer, self.vpos.as_u64(), &mut self.contigs);
        }
    }

    pub fn read_line(&mut self) -> Result<Option<String>> {
        if let Some(line) = self.first_data_line.take() {
            return Ok(Some(line));
        }
        loop {
            self.vpos = self.reader.virtual_position();
            self.buffer.clear();
            let n = self.reader.read_line(&mut self.buffer)?;
            if n == 0 {
                return Ok(None);
            }
            if self.buffer.starts_with('#') {
                continue;
            }
            return Ok(Some(self.buffer.trim_end_matches(['\r', '\n']).to_string()));
        }
    }

    pub fn current_vpos(&self) -> VirtualPosition {
        self.vpos
    }

    pub fn next_record_with_vpos(&mut self) -> Result<Option<(VcfRecord, VirtualPosition)>> {
        // The record carries its own start vpos (`offset`); the returned position is the one
        // AFTER the record so the indexer can form a proper [start, end) chunk.
        match self.read_record()? {
            Some(rec) => {
                let after = self.reader.virtual_position();
                Ok(Some((rec, after)))
            }
            None => Ok(None),
        }
    }

    /// Next data line with `(start, end)` virtual positions.
    pub fn next_line_with_vpos(&mut self) -> Result<Option<(String, VirtualPosition, VirtualPosition)>> {
        if let Some(line) = self.first_data_line.take() {
            let start = self.vpos;
            let after = self.reader.virtual_position();
            return Ok(Some((line, start, after)));
        }
        loop {
            let start = self.reader.virtual_position();
            self.buffer.clear();
            let n = self.reader.read_line(&mut self.buffer)?;
            if n == 0 {
                return Ok(None);
            }
            self.vpos = start;
            if self.buffer.starts_with('#') {
                continue;
            }
            let after = self.reader.virtual_position();
            return Ok(Some((self.buffer.trim_end_matches(['\r', '\n']).to_string(), start, after)));
        }
    }
}

/// Parse a single VCF data line held in memory. `offset` is meaningless for a
/// record that never came from a file, so it is reported as 0; the contig id is
/// 0 because there is no dictionary to resolve it against.
pub fn parse_vcf_line(line: &str) -> Option<VcfRecord> {
    let mut dict = ContigDict::new();
    parse_vcf_record(line, 0, &mut dict).ok().flatten()
}

pub fn parse_vcf_record(line: &str, offset: u64, contigs: &mut ContigDict) -> Result<Option<VcfRecord>> {
    let cols: Vec<&str> = line.trim_end_matches(['\r', '\n']).split('\t').collect();
    if cols.len() < 8 {
        return Ok(None);
    }

    let chrom = cols[0];
    let pos = cols[1]
        .parse::<u32>()
        .map_err(|_| VcfError::ParseError(format!("invalid POS {:?} on {chrom}", cols[1])))?;

    let chr_id = contigs.insert(chrom);

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
    pub(crate) contigs: ContigDict,
}

impl KtileSourceReader {
    pub fn open(path: &Path) -> Result<Self> {
        let reader = KtileReader::open(path).map_err(|_| VcfError::InvalidFormat)?;
        let headers: Vec<String> = reader
            .headers_block()
            .split('\n')
            .filter(|l| !l.is_empty())
            .map(|l| l.to_string())
            .collect();
        let contigs = ContigDict::from_header_lines(headers.iter().map(String::as_str));
        Ok(Self { reader, cursor: 0, headers, contigs })
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
    pub fn read_line_with_meta(&mut self) -> Result<Option<(String, Option<(u32, u32)>)>> {
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

    /// Zero-copy fast path: appends the next line directly into `batch`.
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
        parse_vcf_record(&line, 0, &mut self.contigs)
    }
}
