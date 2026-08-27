//! gVCF block collapse: emit consecutive REF-only sites as single &lt;NON_REF&gt;
//! record per block, partitioned by DP bins (matches bcftools `--gvcf MIN_DP[,DP_BIN...]`).

use anyhow::Result;
use std::io::Write;

#[derive(Clone, Debug)]
pub struct GvcfBlocker {
    bins: Vec<u32>,
    cur: Option<GvcfBlock>,
}

#[derive(Clone, Debug)]
struct GvcfBlock {
    chrom: String,
    start: u64,
    end: u64,
    bin_lo: u32,
    bin_hi: u32,
    ref_base: String,
    min_dp: u32,
    min_qual: f64,
    n_sites: u32,
}

impl GvcfBlocker {
    pub fn new(bins: Vec<u32>) -> Self {
        let mut b = bins;
        b.sort();
        if b.is_empty() { b.push(0); }
        Self { bins: b, cur: None }
    }

    /// Try to extend current block with a REF-only site. Returns Some(flush)
    /// if current block must be flushed first.
    pub fn add_ref_site<W: Write>(
        &mut self,
        chrom: &str,
        pos: u64,
        ref_base: &str,
        dp: u32,
        qual: f64,
        out: &mut W,
    ) -> Result<()> {
        let bin = self.dp_bin(dp);
        let (bin_lo, bin_hi) = (self.bins[bin], self.bins.get(bin + 1).copied().unwrap_or(u32::MAX));
        let should_extend = self.cur.as_ref()
            .map(|c| c.chrom == chrom && c.end + 1 == pos && c.bin_lo == bin_lo && c.bin_hi == bin_hi)
            .unwrap_or(false);
        if should_extend {
            let c = self.cur.as_mut().unwrap();
            c.end = pos;
            c.min_dp = c.min_dp.min(dp);
            c.min_qual = c.min_qual.min(qual);
            c.n_sites += 1;
        } else {
            self.flush(out)?;
            self.cur = Some(GvcfBlock {
                chrom: chrom.to_string(),
                start: pos, end: pos,
                bin_lo, bin_hi,
                ref_base: ref_base.to_string(),
                min_dp: dp, min_qual: qual,
                n_sites: 1,
            });
        }
        Ok(())
    }

    /// Force flush of current block (called on chromosome change or end-of-stream).
    pub fn flush<W: Write>(&mut self, out: &mut W) -> Result<()> {
        if let Some(b) = self.cur.take() {
            writeln!(out, "{}\t{}\t.\t{}\t<NON_REF>\t{}\t.\tEND={};MIN_DP={};BIN={}",
                b.chrom, b.start, b.ref_base, format_qual(b.min_qual), b.end, b.min_dp, b.bin_lo)?;
        }
        Ok(())
    }

    /// Reset on chromosome change without writing (used internally by add_ref_site).
    pub fn reset_chrom<W: Write>(&mut self, out: &mut W) -> Result<()> {
        self.flush(out)
    }

    fn dp_bin(&self, dp: u32) -> usize {
        let mut idx = 0;
        for (i, &b) in self.bins.iter().enumerate() {
            if dp >= b { idx = i; } else { break; }
        }
        idx
    }
}

fn format_qual(q: f64) -> String {
    if q.is_finite() { format!("{:.2}", q) } else { ".".into() }
}

#[cfg(test)]
#[path = "../../tests/unit/call_gvcf.rs"]
mod tests;
