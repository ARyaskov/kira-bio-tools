use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::bgzf::{BgzfReader as NoodlesBgzfReader, VirtualPosition};
use crate::util::{chr_name_to_id, detect_format, VcfFormat};
use crate::vcf::parser::extract_contig_id;
use crate::vcf::structs::{Result, VcfError, VcfRecord};

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
            VcfFormat::Gzip => Err(VcfError::InvalidFormat),
        }
    }

    pub fn open_for_indexing<P: AsRef<Path>>(path: P) -> Result<Self> {
        Ok(Self::BgzfIndexing(BgzfIndexingReader::open(path)?))
    }

    pub fn read_record(&mut self) -> Result<Option<VcfRecord>> {
        match self {
            Self::Plain(r) => r.read_record(),
            Self::Bgzf(r) => r.read_record(),
            Self::BgzfIndexing(r) => r.read_record(),
        }
    }

    pub fn read_line(&mut self) -> Result<Option<String>> {
        match self {
            Self::Plain(r) => r.read_line(),
            Self::Bgzf(r) => r.read_line(),
            Self::BgzfIndexing(r) => r.read_line(),
        }
    }

    pub fn header(&self) -> Result<Vec<String>> {
        match self {
            Self::Plain(r) => Ok(r.headers.clone()),
            Self::Bgzf(r) => Ok(r.headers.clone()),
            Self::BgzfIndexing(r) => Ok(r.headers.clone()),
        }
    }

    pub fn contigs(&self) -> Vec<String> {
        match self {
            Self::Plain(r) => r.contigs.clone(),
            Self::Bgzf(r) => r.contigs.clone(),
            Self::BgzfIndexing(r) => r.contigs.clone(),
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
}

impl PlainReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);
        let mut buffer = String::new();
        let mut contigs = Vec::new();
        let mut headers = Vec::new();
        let mut offset = 0u64;

        loop {
            buffer.clear();
            let n = reader.read_line(&mut buffer)?;
            if n == 0 {
                break;
            }
            offset += n as u64;
            if !buffer.starts_with('#') {
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
            buffer,
            contigs,
            headers,
            offset,
        })
    }

    pub fn read_record(&mut self) -> Result<Option<VcfRecord>> {
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
        self.buffer.clear();
        let n = self.reader.read_line(&mut self.buffer)?;
        if n == 0 {
            return Ok(None);
        }
        self.offset += n as u64;
        Ok(Some(self.buffer.trim_end().to_string()))
    }
}

pub struct BgzfReader {
    reader: NoodlesBgzfReader<File>,
    buffer: String,
    contigs: Vec<String>,
    headers: Vec<String>,
}

impl BgzfReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let mut reader = NoodlesBgzfReader::open(path)?;
        let mut buffer = String::new();
        let mut contigs = Vec::new();
        let mut headers = Vec::new();

        loop {
            buffer.clear();
            let n = reader.read_line(&mut buffer)?;
            if n == 0 {
                break;
            }
            if !buffer.starts_with('#') {
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
            buffer,
            contigs,
            headers,
        })
    }

    pub fn read_record(&mut self) -> Result<Option<VcfRecord>> {
        self.buffer.clear();
        let n = self.reader.read_line(&mut self.buffer)?;
        if n == 0 {
            return Ok(None);
        }

        if self.buffer.starts_with('#') {
            return self.read_record();
        }

        let vpos = self.reader.virtual_position();
        parse_vcf_record(&self.buffer, vpos.as_u64())
    }

    pub fn read_line(&mut self) -> Result<Option<String>> {
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
}

impl BgzfIndexingReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let mut reader = NoodlesBgzfReader::open(path)?;
        let mut buffer = String::new();
        let mut contigs = Vec::new();
        let mut headers = Vec::new();

        loop {
            buffer.clear();
            let n = reader.read_line(&mut buffer)?;
            if n == 0 {
                break;
            }
            if !buffer.starts_with('#') {
                break;
            }
            headers.push(buffer.trim_end().to_string());
            if buffer.starts_with("##contig=") {
                if let Some(id) = extract_contig_id(&buffer) {
                    contigs.push(id);
                }
            }
        }

        let vpos = reader.virtual_position();

        Ok(Self {
            reader,
            buffer,
            contigs,
            headers,
            vpos,
        })
    }

    pub fn read_record(&mut self) -> Result<Option<VcfRecord>> {
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
        let vpos = self.reader.virtual_position();
        match self.read_record()? {
            Some(rec) => Ok(Some((rec, vpos))),
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
