use super::header::{BcfHeaderDict, parse_header_to_dict};
use super::record::decode_record_to_vcf;
use super::BCF_MAGIC;
use anyhow::{Context, Result, bail};
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

pub struct BcfReader {
    inner: Box<dyn Read>,
    pub dict: BcfHeaderDict,
    pub header_lines: Vec<String>,
}

impl BcfReader {
    pub fn open<P: AsRef<Path>>(p: P) -> Result<Self> {
        let path = p.as_ref();
        let mut f = File::open(path).with_context(|| format!("open {:?}", path))?;
        let mut first2 = [0u8; 2];
        use std::io::Read as _;
        f.read_exact(&mut first2).context("probe magic")?;
        use std::io::Seek as _;
        f.seek(std::io::SeekFrom::Start(0))?;
        let is_bgzf = first2 == [0x1F, 0x8B];
        let mut inner: Box<dyn Read> = if is_bgzf {
            Box::new(BufReader::with_capacity(1 << 20, noodles_bgzf::io::Reader::new(f)))
        } else {
            Box::new(BufReader::with_capacity(1 << 20, f))
        };
        let mut magic = [0u8; 5];
        inner.read_exact(&mut magic).context("read BCF magic")?;
        if &magic != BCF_MAGIC { bail!("not a BCF file (bad magic)"); }
        let mut l_text_buf = [0u8; 4];
        inner.read_exact(&mut l_text_buf).context("read l_text")?;
        let l_text = u32::from_le_bytes(l_text_buf);
        let mut header_bytes = vec![0u8; l_text as usize];
        inner.read_exact(&mut header_bytes).context("read header text")?;
        while header_bytes.last() == Some(&0) { header_bytes.pop(); }
        let header_text = String::from_utf8_lossy(&header_bytes);
        let header_lines: Vec<String> = header_text.lines().map(|s| s.to_string()).collect();
        let dict = parse_header_to_dict(&header_lines);
        Ok(Self { inner, dict, header_lines })
    }

    fn open_bgzf(path: &Path) -> Result<Box<dyn Read>> {
        let f = File::open(path)?;
        let rd = noodles_bgzf::io::Reader::new(f);
        Ok(Box::new(BufReader::with_capacity(1 << 20, rd)))
    }

    pub fn read_record_line(&mut self) -> Result<Option<String>> {
        decode_record_to_vcf(&mut self.inner, &self.dict)
    }
}
