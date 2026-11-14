//! gVCF block collapse: emit consecutive REF-only sites as single &lt;NON_REF&gt;
//! record per block, partitioned by DP bins (matches bcftools `--gvcf MIN_DP[,DP_BIN...]`).

use anyhow::Result;
use std::io::Write;

#[derive(Clone, Debug)]
pub struct GvcfBlocker {
    bins: Vec<u32>,
    cur: Option<GvcfBlock>,
    likelihoods: Option<LikelihoodBlock>,
}

#[derive(Clone, Debug)]
struct LikelihoodBlock {
    chrom: String,
    start: u64,
    end: u64,
    range: usize,
    ref_base: String,
    alt: String,
    min_dp: u32,
    depths: Vec<u32>,
    pls: Vec<Vec<i32>>,
    qs: Option<String>,
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
    pub fn from_spec(spec: &str) -> Result<Self> {
        let bins = spec
            .split(',')
            .map(|value| value.parse::<u32>().map_err(anyhow::Error::from))
            .collect::<Result<Vec<_>>>()?;
        anyhow::ensure!(
            !bins.is_empty() && bins.windows(2).all(|pair| pair[0] < pair[1]),
            "gVCF depth thresholds must be strictly increasing"
        );
        Ok(Self::new(bins))
    }
    pub fn new(bins: Vec<u32>) -> Self {
        let mut b = bins;
        b.sort();
        if b.is_empty() {
            b.push(0);
        }
        Self {
            bins: b,
            cur: None,
            likelihoods: None,
        }
    }

    pub fn write_likelihoods<W: Write>(&mut self, text: &str, out: &mut W) -> Result<()> {
        for line in text.lines().filter(|line| !line.is_empty()) {
            let cols: Vec<&str> = line.split('\t').collect();
            anyhow::ensure!(cols.len() >= 9, "invalid mpileup record");
            let pos = cols[1].parse::<u64>()?;
            let is_ref = matches!(cols[4], "<*>" | "<NON_REF>" | ".");
            let keys: Vec<&str> = cols[8].split(':').collect();
            let dp_index = keys.iter().position(|key| *key == "DP");
            let pl_index = keys.iter().position(|key| *key == "PL");
            let samples: Vec<Vec<&str>> = cols[9..]
                .iter()
                .map(|sample| sample.split(':').collect())
                .collect();
            let depths = dp_index.and_then(|index| {
                samples
                    .iter()
                    .map(|sample| sample.get(index)?.parse::<u32>().ok())
                    .collect::<Option<Vec<_>>>()
            });
            let min_dp = depths
                .as_ref()
                .and_then(|values| values.iter().min().copied())
                .unwrap_or(0);
            let range = self.bins.partition_point(|threshold| *threshold <= min_dp);
            let can_collapse = is_ref && depths.is_some() && range > 0;
            if self.likelihoods.as_ref().is_some_and(|block| {
                !can_collapse
                    || block.chrom != cols[0]
                    || pos > block.end + 1
                    || range != block.range
            }) {
                if let Some(block) = &mut self.likelihoods
                    && block.chrom == cols[0]
                    && block.end == pos
                {
                    block.end = block.end.saturating_sub(1);
                }
                self.flush_likelihoods(out)?;
            }
            if can_collapse {
                let pls = pl_index
                    .and_then(|index| {
                        samples
                            .iter()
                            .map(|sample| {
                                sample
                                    .get(index)?
                                    .split(',')
                                    .map(|value| value.parse::<i32>().ok())
                                    .collect::<Option<Vec<_>>>()
                            })
                            .collect::<Option<Vec<_>>>()
                    })
                    .unwrap_or_default();
                let end = cols[7]
                    .split(';')
                    .find_map(|field| {
                        field
                            .strip_prefix("END=")
                            .and_then(|value| value.parse().ok())
                    })
                    .unwrap_or(pos);
                match &mut self.likelihoods {
                    Some(block) => {
                        block.end = end;
                        block.min_dp = block.min_dp.min(min_dp);
                        for (current, next) in block.depths.iter_mut().zip(depths.unwrap()) {
                            *current = (*current).min(next);
                        }
                        if pls.is_empty() {
                            block.pls.clear();
                        } else {
                            for (current, next) in block.pls.iter_mut().zip(pls) {
                                anyhow::ensure!(
                                    current.len() == 3 && next.len() == 3,
                                    "unexpected number of reference PL fields"
                                );
                                if (next[1], next[2]) < (current[1], current[2]) {
                                    current[1] = next[1];
                                    current[2] = next[2];
                                }
                            }
                        }
                    }
                    None => {
                        self.likelihoods = Some(LikelihoodBlock {
                            chrom: cols[0].to_string(),
                            start: pos,
                            end,
                            range,
                            ref_base: cols[3].to_string(),
                            alt: cols[4].to_string(),
                            min_dp,
                            depths: depths.unwrap(),
                            pls,
                            qs: cols[7]
                                .split(';')
                                .find_map(|field| field.strip_prefix("QS=").map(String::from)),
                        });
                    }
                }
            } else if is_ref && min_dp > 0 {
                let mut cols: Vec<String> = cols.iter().map(|value| value.to_string()).collect();
                if cols[7] == "." {
                    cols[7].clear();
                } else {
                    cols[7].push(';');
                }
                cols[7].push_str(&format!("MIN_DP={min_dp}"));
                writeln!(out, "{}", cols.join("\t"))?;
            } else {
                writeln!(out, "{line}")?;
            }
        }
        Ok(())
    }

    fn flush_likelihoods<W: Write>(&mut self, out: &mut W) -> Result<()> {
        if let Some(block) = self.likelihoods.take() {
            write!(
                out,
                "{}\t{}\t.\t{}\t{}\t.\t.\t",
                block.chrom, block.start, block.ref_base, block.alt
            )?;
            if block.start < block.end {
                write!(out, "END={};", block.end)?;
            }
            write!(out, "MIN_DP={}", block.min_dp)?;
            if let Some(qs) = block.qs {
                write!(out, ";QS={qs}")?;
            }
            write!(
                out,
                "\t{}",
                if block.pls.is_empty() { "DP" } else { "PL:DP" }
            )?;
            for (i, depth) in block.depths.iter().enumerate() {
                write!(out, "\t")?;
                if let Some(pls) = block.pls.get(i) {
                    for (j, value) in pls.iter().enumerate() {
                        if j > 0 {
                            write!(out, ",")?;
                        }
                        write!(out, "{value}")?;
                    }
                    write!(out, ":")?;
                }
                write!(out, "{depth}")?;
            }
            writeln!(out)?;
        }
        Ok(())
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
        let (bin_lo, bin_hi) = (
            self.bins[bin],
            self.bins.get(bin + 1).copied().unwrap_or(u32::MAX),
        );
        let should_extend = self
            .cur
            .as_ref()
            .map(|c| {
                c.chrom == chrom && c.end + 1 == pos && c.bin_lo == bin_lo && c.bin_hi == bin_hi
            })
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
                start: pos,
                end: pos,
                bin_lo,
                bin_hi,
                ref_base: ref_base.to_string(),
                min_dp: dp,
                min_qual: qual,
                n_sites: 1,
            });
        }
        Ok(())
    }

    /// Force flush of current block (called on chromosome change or end-of-stream).
    pub fn flush<W: Write>(&mut self, out: &mut W) -> Result<()> {
        self.flush_likelihoods(out)?;
        if let Some(b) = self.cur.take() {
            writeln!(
                out,
                "{}\t{}\t.\t{}\t<NON_REF>\t{}\t.\tEND={};MIN_DP={};BIN={}",
                b.chrom,
                b.start,
                b.ref_base,
                format_qual(b.min_qual),
                b.end,
                b.min_dp,
                b.bin_lo
            )?;
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
            if dp >= b {
                idx = i;
            } else {
                break;
            }
        }
        idx
    }
}

fn format_qual(q: f64) -> String {
    if q.is_finite() {
        format!("{:.2}", q)
    } else {
        ".".into()
    }
}

#[cfg(test)]
#[path = "../../tests/unit/call_gvcf.rs"]
mod tests;
