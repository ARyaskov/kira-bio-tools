use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use flate2::read::GzDecoder;

use crate::bgzf::{BgzfLineReader, BgzfReader, VirtualPosition};
use crate::util::{detect_format, parse_vcf_line_fast, VcfFormat};
use crate::vcf::parser::extract_contig_id;
use crate::vcf::structs::{Result, VcfRecord};

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

    pub fn next_raw_line(&mut self) -> Result<Option<(String, u64)>> {
        match self {
            VcfReader::Plain(r) => r.next_raw_line(),
            VcfReader::Gzip(r) => r.next_raw_line(),
            VcfReader::Bgzf(r) => r.next_raw_line(),
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
    pub contigs: Vec<String>,
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

    pub fn next_raw_line(&mut self) -> Result<Option<(String, u64)>> {
        if !self.header_parsed {
            self.header()?;
        }

        let mut line = String::new();
        let offset = self.offset;
        let bytes = self.reader.read_line(&mut line)?;
        if bytes == 0 {
            return Ok(None);
        }
        self.offset += bytes as u64;

        Ok(Some((line, offset)))
    }
}

pub struct GzipVcfReader {
    reader: BufReader<GzDecoder<File>>,
    buf: String,
    offset: u64,
    pub contigs: Vec<String>,
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

    pub fn next_raw_line(&mut self) -> Result<Option<(String, u64)>> {
        if !self.header_parsed {
            self.header()?;
        }

        let mut line = String::new();
        let offset = self.offset;
        let bytes = self.reader.read_line(&mut line)?;
        if bytes == 0 {
            return Ok(None);
        }
        self.offset += bytes as u64;

        Ok(Some((line, offset)))
    }
}

pub struct BgzfVcfReader {
    reader: BgzfLineReader<File>,
    pub contigs: Vec<String>,
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

    pub fn next_raw_line(&mut self) -> Result<Option<(String, u64)>> {
        let (line, vpos) = match self.reader.read_line()? {
            Some(v) => v,
            None => return Ok(None),
        };
        Ok(Some((line.to_string(), vpos.as_u64())))
    }
}
