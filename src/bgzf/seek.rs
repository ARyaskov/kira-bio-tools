use anyhow::{Context, Result};
use noodles_bgzf::io::Reader as NoodlesBgzf;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub struct BgzfSeekReader {
    inner: NoodlesBgzf<File>,
}

impl BgzfSeekReader {
    pub fn open<P: AsRef<Path>>(p: P) -> Result<Self> {
        let f = File::open(p.as_ref()).with_context(|| format!("open {:?}", p.as_ref()))?;
        Ok(Self { inner: NoodlesBgzf::new(f) })
    }

    pub fn seek_to(&mut self, vpos: u64) -> Result<()> {
        let bgzf_vpos = noodles_bgzf::VirtualPosition::from(vpos);
        self.inner.seek(bgzf_vpos)?;
        Ok(())
    }

    pub fn read_line(&mut self) -> Result<Option<String>> {
        let mut br = BufReader::with_capacity(1 << 16, &mut self.inner);
        let mut buf = String::new();
        let n = br.read_line(&mut buf)?;
        if n == 0 { return Ok(None); }
        if buf.ends_with('\n') { buf.pop(); }
        if buf.ends_with('\r') { buf.pop(); }
        Ok(Some(buf))
    }
}
