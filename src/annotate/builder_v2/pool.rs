use anyhow::Result;
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};

pub struct StringPool {
    mem: Vec<u8>,
    file: Option<BufWriter<File>>,
    path: Option<PathBuf>,
    size: u64,
    limit: Option<usize>,
    spilled: bool,
}

impl StringPool {
    pub fn new() -> Self {
        Self {
            mem: Vec::new(),
            file: None,
            path: None,
            size: 0,
            limit: None,
            spilled: false,
        }
    }

    pub fn with_limit(limit: Option<usize>, spill_path: Option<PathBuf>) -> Self {
        Self {
            mem: Vec::new(),
            file: None,
            path: spill_path,
            size: 0,
            limit,
            spilled: false,
        }
    }

    pub fn len(&self) -> usize {
        self.size as usize
    }

    pub fn is_in_memory(&self) -> bool {
        self.file.is_none()
    }

    pub fn append_cstr(&mut self, s: &str) -> usize {
        let ofs = self.size as usize;
        if self.file.is_none() {
            if let Some(limit) = self.limit {
                let need = self.size as usize + s.len() + 1;
                if need > limit {
                    if let Some(ref path) = self.path {
                        let path = path.clone();
                        let _ = self.spill_to_file(&path);
                        self.spilled = true;
                    }
                }
            }
        }
        if let Some(ref mut file) = self.file {
            let _ = file.write_all(s.as_bytes());
            let _ = file.write_all(&[0]);
            self.size += (s.len() + 1) as u64;
        } else {
            self.mem.extend_from_slice(s.as_bytes());
            self.mem.push(0);
            self.size += (s.len() + 1) as u64;
        }
        ofs
    }

    pub fn spilled(&self) -> bool {
        self.spilled
    }

    pub fn spill_to_file(&mut self, path: &Path) -> Result<()> {
        if self.file.is_some() {
            return Ok(());
        }

        let file = OpenOptions::new()
            .create(true)
            .write(true)
            .truncate(true)
            .open(path)?;
        let mut writer = BufWriter::new(file);
        writer.write_all(&self.mem)?;
        writer.flush()?;

        self.size = self.mem.len() as u64;
        self.mem.clear();
        self.mem.shrink_to_fit();

        self.path = Some(path.to_path_buf());
        self.file = Some(writer);
        Ok(())
    }

    pub fn write_to(&mut self, out: &mut File) -> Result<()> {
        if let Some(ref mut file) = self.file {
            file.flush()?;
            let path = self
                .path
                .as_ref()
                .ok_or_else(|| anyhow::anyhow!("String pool file path missing"))?;
            let mut input = File::open(path)?;
            let mut buf = vec![0u8; 8 * 1024 * 1024];
            loop {
                let n = input.read(&mut buf)?;
                if n == 0 {
                    break;
                }
                out.write_all(&buf[..n])?;
            }
        } else {
            out.write_all(&self.mem)?;
        }
        Ok(())
    }

    pub fn cleanup(&mut self) {
        if let Some(ref path) = self.path {
            let _ = std::fs::remove_file(path);
        }
    }
}
