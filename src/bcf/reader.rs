use super::BCF_MAGIC;
use super::header::{BcfHeaderDict, parse_header_to_dict};
use super::record::{decode_blocks_to_vcf, read_record_raw};
use anyhow::{Context, Result, bail};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

/// Where BCF bytes come from. `File` keeps BGZF virtual positions for
/// indexing and random access; `Stream` is any already-decompressed source.
pub enum BcfInput {
    Stream(Box<dyn BufRead + Send>),
    File(noodles_bgzf::io::Reader<File>),
}

impl Read for BcfInput {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            BcfInput::Stream(r) => r.read(buf),
            BcfInput::File(r) => r.read(buf),
        }
    }
}

/// Remove a `,IDX=N` attribute from a structured header line.
pub fn strip_idx(line: &str) -> String {
    if !line.starts_with("##") || !line.contains(",IDX=") {
        return line.to_string();
    }
    let Some(s) = line.find(",IDX=") else { return line.to_string() };
    let rest = &line[s + 5..];
    let e = rest.find(|c: char| !c.is_ascii_digit()).unwrap_or(rest.len());
    format!("{}{}", &line[..s], &rest[e..])
}

pub struct BcfReader {
    inner: BcfInput,
    pub dict: BcfHeaderDict,
    pub header_lines: Vec<String>,
}

impl BcfReader {
    pub fn open<P: AsRef<Path>>(p: P) -> Result<Self> {
        let path = p.as_ref();
        let mut f = File::open(path).with_context(|| format!("open {:?}", path))?;
        let mut first2 = [0u8; 2];
        f.read_exact(&mut first2).context("probe magic")?;
        use std::io::Seek as _;
        f.seek(io::SeekFrom::Start(0))?;
        let inner = if first2 == [0x1F, 0x8B] {
            BcfInput::File(noodles_bgzf::io::Reader::new(f))
        } else {
            BcfInput::Stream(Box::new(BufReader::with_capacity(1 << 20, f)))
        };
        Self::from_input(inner)
    }

    /// Reader over an already-decompressed BCF byte stream.
    pub fn from_bufread(r: Box<dyn BufRead + Send>) -> Result<Self> {
        Self::from_input(BcfInput::Stream(r))
    }

    fn from_input(mut inner: BcfInput) -> Result<Self> {
        let mut magic = [0u8; 5];
        inner.read_exact(&mut magic).context("read BCF magic")?;
        if magic != BCF_MAGIC { bail!("not a BCF file (bad magic)"); }
        let mut l_text_buf = [0u8; 4];
        inner.read_exact(&mut l_text_buf).context("read l_text")?;
        let l_text = u32::from_le_bytes(l_text_buf);
        if l_text > (1 << 30) { bail!("corrupt BCF header length {l_text}"); }
        let mut header_bytes = vec![0u8; l_text as usize];
        inner.read_exact(&mut header_bytes).context("read header text")?;
        while header_bytes.last() == Some(&0) { header_bytes.pop(); }
        let header_text = String::from_utf8_lossy(&header_bytes);
        let raw_lines: Vec<String> = header_text.lines().map(|s| s.to_string()).collect();
        let dict = parse_header_to_dict(&raw_lines);
        // The IDX attributes only describe this file's dictionary; text output drops them.
        let header_lines: Vec<String> = raw_lines.iter().map(|l| strip_idx(l)).collect();
        Ok(Self { inner, dict, header_lines })
    }

    /// Number of header lines (for callers that want the raw dictionary size).
    pub fn header_len(&self) -> usize {
        self.header_lines.len()
    }

    pub fn read_record_line(&mut self) -> Result<Option<String>> {
        match read_record_raw(&mut self.inner)? {
            Some((shared, indiv)) => decode_blocks_to_vcf(&shared, &indiv, &self.dict).map(Some),
            None => Ok(None),
        }
    }

    pub fn read_record_raw(&mut self) -> Result<Option<(Vec<u8>, Vec<u8>)>> {
        read_record_raw(&mut self.inner)
    }

    pub fn decode(&self, shared: &[u8], indiv: &[u8]) -> Result<String> {
        decode_blocks_to_vcf(shared, indiv, &self.dict)
    }

    pub fn is_seekable(&self) -> bool {
        matches!(self.inner, BcfInput::File(_))
    }

    /// BGZF virtual position of the next record (BGZF files only).
    pub fn virtual_position(&self) -> Option<u64> {
        match &self.inner {
            BcfInput::File(r) => Some(u64::from(r.virtual_position())),
            BcfInput::Stream(_) => None,
        }
    }

    pub fn seek(&mut self, vpos: u64) -> Result<()> {
        match &mut self.inner {
            BcfInput::File(r) => {
                r.seek(noodles_bgzf::VirtualPosition::from(vpos)).context("bgzf seek")?;
                Ok(())
            }
            BcfInput::Stream(_) => bail!("BCF stream is not seekable"),
        }
    }
}
