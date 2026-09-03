//! Region queries over an indexed BGZF VCF or BCF.

use std::fs::File;
use std::io::{BufRead, Read};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};

use super::binning::BinIndex;
use super::builder::{find_index_for, vcf_line_interval};
use crate::bcf::record::{record_end0, record_meta};
use crate::bcf::{BCF_MAGIC, BcfReader};

enum Source {
    Vcf(noodles_bgzf::io::Reader<File>),
    Bcf(BcfReader),
}

pub struct IndexedVcfReader {
    src: Source,
    index: BinIndex,
    index_path: PathBuf,
    headers: Vec<String>,
    names: Vec<String>,
}

impl IndexedVcfReader {
    pub fn open(path: &Path) -> Result<Self> {
        let Some(idx) = find_index_for(path) else {
            bail!("no .csi or .tbi index found for {}", path.display());
        };
        Self::open_with_index(path, &idx)
    }

    pub fn open_with_index(path: &Path, index_path: &Path) -> Result<Self> {
        let index = BinIndex::load(index_path).with_context(|| format!("read index {}", index_path.display()))?;
        let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
        let mut probe = noodles_bgzf::io::Reader::new(f);
        let mut magic = [0u8; 5];
        let n = probe.read(&mut magic)?;
        if n == 5 && magic == BCF_MAGIC {
            let r = BcfReader::open(path)?;
            let headers = r.header_lines.clone();
            let max_rid = r.dict.contig_idx.values().copied().max();
            let mut names = Vec::new();
            if let Some(max) = max_rid {
                for rid in 0..=max {
                    names.push(r.dict.contig_name(rid).unwrap_or("").to_string());
                }
            }
            return Ok(Self { src: Source::Bcf(r), index, index_path: index_path.to_path_buf(), headers, names });
        }
        let f = File::open(path)?;
        let mut r = noodles_bgzf::io::Reader::new(f);
        let mut headers = Vec::new();
        let mut line = String::new();
        loop {
            line.clear();
            let n = r.read_line(&mut line)?;
            if n == 0 || !line.starts_with('#') {
                break;
            }
            headers.push(line.trim_end_matches(['\r', '\n']).to_string());
        }
        let names = index.names().to_vec();
        Ok(Self { src: Source::Vcf(r), index, index_path: index_path.to_path_buf(), headers, names })
    }

    pub fn headers(&self) -> &[String] {
        &self.headers
    }

    pub fn index(&self) -> &BinIndex {
        &self.index
    }

    pub fn index_path(&self) -> &Path {
        &self.index_path
    }

    /// Contig names in index order.
    pub fn names(&self) -> &[String] {
        &self.names
    }

    pub fn ref_id(&self, chrom: &str) -> Option<usize> {
        self.names.iter().position(|n| n == chrom)
    }

    /// Visit records overlapping `[beg1, end1]` (1-based, inclusive) on
    /// `chrom`, in file order. The callback returns `false` to stop early.
    pub fn query<F>(&mut self, chrom: &str, beg1: u32, end1: u32, mut f: F) -> Result<()>
    where
        F: FnMut(&str) -> Result<bool>,
    {
        let Some(rid) = self.ref_id(chrom) else { return Ok(()) };
        let beg0 = beg1.saturating_sub(1) as u64;
        let end0 = end1 as u64;
        let chunks = self.index.query(rid, beg0, end0);
        match &mut self.src {
            Source::Vcf(r) => {
                let mut line = String::new();
                for (s, e) in chunks {
                    r.seek(noodles_bgzf::VirtualPosition::from(s)).context("bgzf seek")?;
                    loop {
                        if u64::from(r.virtual_position()) >= e {
                            break;
                        }
                        line.clear();
                        let n = r.read_line(&mut line)?;
                        if n == 0 {
                            break;
                        }
                        let l = line.trim_end_matches(['\r', '\n']);
                        if l.is_empty() || l.starts_with('#') {
                            continue;
                        }
                        let Some((c, rb, re)) = vcf_line_interval(l) else { continue };
                        if c != chrom {
                            continue;
                        }
                        if rb >= end0 {
                            return Ok(());
                        }
                        if re > beg0 && !f(l)? {
                            return Ok(());
                        }
                    }
                }
            }
            Source::Bcf(r) => {
                for (s, e) in chunks {
                    r.seek(s)?;
                    loop {
                        if r.virtual_position().unwrap_or(u64::MAX) >= e {
                            break;
                        }
                        let Some((shared, indiv)) = r.read_record_raw()? else { break };
                        let Some(m) = record_meta(&shared) else { continue };
                        if m.rid < 0 || m.rid as usize != rid {
                            continue;
                        }
                        let rb = m.pos.max(0) as u64;
                        if rb >= end0 {
                            return Ok(());
                        }
                        if record_end0(&m) > beg0 {
                            let line = r.decode(&shared, &indiv)?;
                            if !f(&line)? {
                                return Ok(());
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }
}
