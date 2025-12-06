use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::bgzf::{BgzfReader as NoodlesBgzfReader, VirtualPosition};
use crate::bgzf_parallel::{BatchedLineReader, ParallelBgzfReader};
use crate::util::{detect_format, parse_vcf_line_fast, VcfFormat};
use crate::vcf::parser::extract_contig_id;
use crate::vcf::structs::{Result, VcfError, VcfRecord};
use crate::vcf_parser_fast::{FastVcfParser, VcfFields};

pub enum UnifiedVcfReader {
    Plain(PlainReader),
    Bgzf(BgzfReader),
    BgzfIndexing(BgzfIndexingReader),
}

impl UnifiedVcfReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let format = detect_format(path)?;

        match format {
            VcfFormat::Plain => Ok(Self::Plain(PlainReader::open(path)?)),
            VcfFormat::Bgzf => Ok(Self::Bgzf(BgzfReader::open(path)?)),
            VcfFormat::Gzip => Err(VcfError::InvalidFormat(
                "Plain gzip not supported. Use BGZF compression (bgzip).".into(),
            )),
        }
    }

    pub fn open_for_indexing<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let format = detect_format(path)?;

        match format {
            VcfFormat::Plain => Ok(Self::Plain(PlainReader::open(path)?)),
            VcfFormat::Bgzf => Ok(Self::BgzfIndexing(BgzfIndexingReader::open(path)?)),
            VcfFormat::Gzip => Err(VcfError::InvalidFormat(
                "Plain gzip not supported. Use BGZF compression (bgzip).".into(),
            )),
        }
    }

    pub fn header(&mut self) -> Result<Vec<String>> {
        match self {
            Self::Plain(r) => r.header(),
            Self::Bgzf(r) => r.header(),
            Self::BgzfIndexing(r) => r.header(),
        }
    }

    pub fn read_line(&mut self) -> Result<Option<String>> {
        match self {
            Self::Plain(r) => r.read_line(),
            Self::Bgzf(r) => r.read_line(),
            Self::BgzfIndexing(r) => r.read_line(),
        }
    }

    pub fn read_batch(&mut self, size: usize) -> Result<Vec<String>> {
        match self {
            Self::Plain(r) => r.read_batch(size),
            Self::Bgzf(r) => r.read_batch(size),
            Self::BgzfIndexing(r) => r.read_batch(size),
        }
    }

    pub fn parse_line<'a>(&self, line: &'a str) -> Option<VcfFields<'a>> {
        let mut parser = FastVcfParser::new(line);
        parser.parse_standard_fields()
    }

    pub fn reference_sequences(&self) -> &[String] {
        match self {
            Self::Plain(r) => &r.contigs,
            Self::Bgzf(r) => &r.contigs,
            Self::BgzfIndexing(r) => &r.contigs,
        }
    }

    pub fn virtual_position(&self) -> Option<VirtualPosition> {
        match self {
            Self::Plain(_) => None,
            Self::Bgzf(_) => None,
            Self::BgzfIndexing(r) => Some(r.virtual_position()),
        }
    }

    pub fn next_record_with_vpos(&mut self) -> Result<Option<(VcfRecord, VirtualPosition)>> {
        match self {
            Self::BgzfIndexing(r) => r.next_record_with_vpos(),
            _ => Err(VcfError::InvalidFormat(
                "virtual_position only available for BGZF indexing mode".into(),
            )),
        }
    }
}

pub struct PlainReader {
    reader: BufReader<File>,
    pub contigs: Vec<String>,
    header_parsed: bool,
}

impl PlainReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        Ok(Self {
            reader: BufReader::with_capacity(8 * 1024 * 1024, file),
            contigs: Vec::new(),
            header_parsed: false,
        })
    }

    pub fn header(&mut self) -> Result<Vec<String>> {
        if self.header_parsed {
            return Ok(Vec::new());
        }

        let mut headers = Vec::new();
        let mut line = String::new();

        loop {
            line.clear();
            let bytes = self.reader.read_line(&mut line)?;
            if bytes == 0 {
                break;
            }

            if line.starts_with('#') {
                if line.starts_with("##contig=") {
                    if let Some(id) = extract_contig_id(&line) {
                        self.contigs.push(id);
                    }
                }
                headers.push(line.trim_end().to_string());

                if line.starts_with("#CHROM") {
                    break;
                }
            } else {
                break;
            }
        }

        self.header_parsed = true;
        Ok(headers)
    }

    pub fn read_line(&mut self) -> Result<Option<String>> {
        if !self.header_parsed {
            self.header()?;
        }

        let mut line = String::new();
        loop {
            line.clear();
            let bytes = self.reader.read_line(&mut line)?;
            if bytes == 0 {
                return Ok(None);
            }

            if !line.starts_with('#') && !line.trim().is_empty() {
                if line.ends_with('\n') {
                    line.pop();
                }
                if line.ends_with('\r') {
                    line.pop();
                }
                return Ok(Some(line));
            }
        }
    }

    pub fn read_batch(&mut self, size: usize) -> Result<Vec<String>> {
        if !self.header_parsed {
            self.header()?;
        }

        let mut batch = Vec::with_capacity(size);
        for _ in 0..size {
            match self.read_line()? {
                Some(line) => batch.push(line),
                None => break,
            }
        }
        Ok(batch)
    }
}

pub struct BgzfReader {
    reader: BatchedLineReader,
    pub contigs: Vec<String>,
    header_parsed: bool,
    line_buffer: Vec<(String, u64)>,
    buffer_pos: usize,
}

impl BgzfReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let parallel_reader = ParallelBgzfReader::open(path)?;
        Ok(Self {
            reader: BatchedLineReader::new(parallel_reader, 200_000),
            contigs: Vec::new(),
            header_parsed: false,
            line_buffer: Vec::new(),
            buffer_pos: 0,
        })
    }

    pub fn header(&mut self) -> Result<Vec<String>> {
        if self.header_parsed {
            return Ok(Vec::new());
        }

        let mut headers = Vec::new();

        loop {
            let line = match self.read_line_internal()? {
                Some(l) => l,
                None => break,
            };

            if line.starts_with('#') {
                if line.starts_with("##contig=") {
                    if let Some(id) = extract_contig_id(&line) {
                        self.contigs.push(id);
                    }
                }
                headers.push(line);

                if headers.last().unwrap().starts_with("#CHROM") {
                    break;
                }
            } else {
                self.line_buffer.insert(0, (line, 0));
                break;
            }
        }

        self.header_parsed = true;
        Ok(headers)
    }

    fn read_line_internal(&mut self) -> Result<Option<String>> {
        if self.buffer_pos < self.line_buffer.len() {
            let line = self.line_buffer[self.buffer_pos].0.clone();
            self.buffer_pos += 1;
            return Ok(Some(line));
        }

        self.line_buffer.clear();
        self.buffer_pos = 0;

        let batch = self.reader.read_batch()?;
        if batch.is_empty() {
            return Ok(None);
        }

        self.line_buffer = batch;
        if self.line_buffer.is_empty() {
            return Ok(None);
        }

        let line = self.line_buffer[0].0.clone();
        self.buffer_pos = 1;
        Ok(Some(line))
    }

    pub fn read_line(&mut self) -> Result<Option<String>> {
        if !self.header_parsed {
            self.header()?;
        }

        self.read_line_internal()
    }

    pub fn read_batch(&mut self, size: usize) -> Result<Vec<String>> {
        if !self.header_parsed {
            self.header()?;
        }

        let mut result = Vec::with_capacity(size);

        while result.len() < size {
            if self.buffer_pos < self.line_buffer.len() {
                result.push(self.line_buffer[self.buffer_pos].0.clone());
                self.buffer_pos += 1;
            } else {
                self.line_buffer.clear();
                self.buffer_pos = 0;

                let batch = self.reader.read_batch()?;
                if batch.is_empty() {
                    break;
                }

                self.line_buffer = batch;
            }
        }

        Ok(result)
    }
}

pub struct BgzfIndexingReader {
    reader: NoodlesBgzfReader<File>,
    pub contigs: Vec<String>,
    header_parsed: bool,
    line_buffer: String,
}

impl BgzfIndexingReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let bgzf_reader = NoodlesBgzfReader::open(path)?;
        Ok(Self {
            reader: bgzf_reader,
            contigs: Vec::new(),
            header_parsed: false,
            line_buffer: String::with_capacity(4096),
        })
    }

    pub fn header(&mut self) -> Result<Vec<String>> {
        if self.header_parsed {
            return Ok(Vec::new());
        }

        let mut headers = Vec::new();

        loop {
            self.line_buffer.clear();
            let bytes = self.reader.read_line(&mut self.line_buffer)?;
            if bytes == 0 {
                break;
            }

            if self.line_buffer.starts_with('#') {
                if self.line_buffer.starts_with("##contig=") {
                    if let Some(id) = extract_contig_id(&self.line_buffer) {
                        self.contigs.push(id);
                    }
                }
                headers.push(self.line_buffer.trim_end().to_string());

                if self.line_buffer.starts_with("#CHROM") {
                    break;
                }
            } else {
                break;
            }
        }

        self.header_parsed = true;
        Ok(headers)
    }

    pub fn read_line(&mut self) -> Result<Option<String>> {
        if !self.header_parsed {
            self.header()?;
        }

        loop {
            self.line_buffer.clear();
            let bytes = self.reader.read_line(&mut self.line_buffer)?;
            if bytes == 0 {
                return Ok(None);
            }

            if !self.line_buffer.starts_with('#') && !self.line_buffer.trim().is_empty() {
                if self.line_buffer.ends_with('\n') {
                    self.line_buffer.pop();
                }
                if self.line_buffer.ends_with('\r') {
                    self.line_buffer.pop();
                }
                return Ok(Some(self.line_buffer.clone()));
            }
        }
    }

    pub fn read_batch(&mut self, size: usize) -> Result<Vec<String>> {
        if !self.header_parsed {
            self.header()?;
        }

        let mut batch = Vec::with_capacity(size);
        for _ in 0..size {
            match self.read_line()? {
                Some(line) => batch.push(line),
                None => break,
            }
        }
        Ok(batch)
    }

    pub fn virtual_position(&self) -> VirtualPosition {
        self.reader.virtual_position()
    }

    pub fn next_record_with_vpos(&mut self) -> Result<Option<(VcfRecord, VirtualPosition)>> {
        if !self.header_parsed {
            self.header()?;
        }

        loop {
            let vpos = self.virtual_position();

            self.line_buffer.clear();
            let bytes = self.reader.read_line(&mut self.line_buffer)?;
            if bytes == 0 {
                return Ok(None);
            }

            if self.line_buffer.starts_with('#') || self.line_buffer.trim().is_empty() {
                continue;
            }

            if let Some((chr_id, position)) = parse_vcf_line_fast(self.line_buffer.as_bytes()) {
                return Ok(Some((
                    VcfRecord {
                        chr_id,
                        position,
                        offset: vpos.as_u64(),
                    },
                    vpos,
                )));
            }
        }
    }
}
