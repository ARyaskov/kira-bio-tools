use anyhow::Result;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::bgzf::{BgzfLineReader, BgzfReader};
use crate::detect_format;
use crate::VcfFormat;

pub enum VcfAnnotationReader {
    Plain(BufReader<File>),
    Bgzf(BgzfLineReader<File>),
}

impl VcfAnnotationReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let format = detect_format(path)?;

        match format {
            VcfFormat::Plain => {
                let file = File::open(path)?;
                Ok(Self::Plain(BufReader::with_capacity(8 * 1024 * 1024, file)))
            }
            VcfFormat::Gzip => {
                anyhow::bail!("Gzip format not supported. Use BGZF compression (bgzip).")
            }
            VcfFormat::Bgzf => {
                let file = File::open(path)?;
                let bgzf = BgzfReader::new(file);
                Ok(Self::Bgzf(BgzfLineReader::new(bgzf)))
            }
        }
    }

    pub fn read_line(&mut self) -> Result<Option<String>> {
        match self {
            Self::Plain(reader) => {
                let mut line = String::new();
                match reader.read_line(&mut line)? {
                    0 => Ok(None),
                    _ => {
                        if line.ends_with('\n') {
                            line.pop();
                        }
                        if line.ends_with('\r') {
                            line.pop();
                        }
                        Ok(Some(line))
                    }
                }
            }
            Self::Bgzf(reader) => match reader.read_line()? {
                Some((line, _vpos)) => Ok(Some(line.to_string())),
                None => Ok(None),
            },
        }
    }

    pub fn read_header(&mut self) -> Result<Vec<String>> {
        let mut headers = Vec::new();

        loop {
            match self.read_line()? {
                Some(line) => {
                    if line.starts_with('#') {
                        let is_chrom_line = line.starts_with("#CHROM");
                        headers.push(line);
                        if is_chrom_line {
                            break;
                        }
                    } else {
                        break;
                    }
                }
                None => break,
            }
        }

        Ok(headers)
    }
}

pub struct BatchVcfReader {
    reader: VcfAnnotationReader,
    batch_size: usize,
    headers_read: bool,
}

impl BatchVcfReader {
    pub fn new(reader: VcfAnnotationReader, batch_size: usize) -> Self {
        Self {
            reader,
            batch_size,
            headers_read: false,
        }
    }

    pub fn read_batch(&mut self) -> Result<Vec<String>> {
        let mut batch = Vec::with_capacity(self.batch_size);

        for _ in 0..self.batch_size {
            match self.reader.read_line()? {
                Some(line) => {
                    if !line.starts_with('#') && !line.trim().is_empty() {
                        batch.push(line);
                    }
                }
                None => break,
            }
        }

        Ok(batch)
    }

    pub fn into_headers_and_self(mut self) -> Result<(Vec<String>, Self)> {
        if !self.headers_read {
            let headers = self.reader.read_header()?;
            self.headers_read = true;
            Ok((headers, self))
        } else {
            Ok((Vec::new(), self))
        }
    }
}
