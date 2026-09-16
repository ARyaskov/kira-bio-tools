//! SimHash sketches of per-sample dosage vectors (`gtcheck --sketch`).
//!
//! All-vs-all sample comparison scores every pair at every site, O(n²·m).
//! The sketch pass instead projects each sample's centred dosage vector onto
//! `bits` Rademacher (±1) random directions, one draw per site, and keeps the
//! sign of every projection. The Hamming distance between two sketches is a
//! binomial estimate of the angle between the dosage vectors (P(bit differs)
//! = θ/π), so a duplicate or close relative sits far nearer to a sample than
//! anyone else and the k nearest sketches of each sample are the candidate
//! pairs worth scoring exactly. One pass costs O(n·m·bits) adds and the
//! candidate search O(n²·bits/64) popcounts, no genotype decoding.

/// Streaming projection accumulator.
pub struct DosageSketch {
    n: usize,
    bits: usize,
    acc: Vec<f32>,
    signs: Vec<f32>,
    sites: usize,
    seed: u64,
}

fn splitmix64(x: &mut u64) -> u64 {
    *x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

impl DosageSketch {
    /// `bits` is rounded up to a multiple of 64 (at least 64).
    pub fn new(n: usize, bits: usize) -> Self {
        let bits = bits.max(64).div_ceil(64) * 64;
        Self {
            n,
            bits,
            acc: vec![0.0; n * bits],
            signs: vec![0.0; bits],
            sites: 0,
            seed: 0x005E_ED0F_DA7A,
        }
    }

    pub fn bits(&self) -> usize {
        self.bits
    }

    /// Sites that contributed (called in at least two samples and not
    /// monomorphic among them).
    pub fn sites(&self) -> usize {
        self.sites
    }

    /// One site: dosage 0/1/2 per sample, negative = missing.
    pub fn add_site(&mut self, dosage: &[i8]) {
        debug_assert_eq!(dosage.len(), self.n);
        let (mut sum, mut called) = (0i64, 0usize);
        let mut lo = i8::MAX;
        let mut hi = i8::MIN;
        for &d in dosage {
            if d >= 0 {
                sum += d as i64;
                called += 1;
                lo = lo.min(d);
                hi = hi.max(d);
            }
        }
        if called < 2 || lo == hi {
            return;
        }
        let mean = sum as f32 / called as f32;
        // Fresh ±1 directions for this site.
        let mut state = self.seed ^ (self.sites as u64).wrapping_mul(0xD6E8_FEB8_6659_FD93);
        for w in 0..self.bits / 64 {
            let word = splitmix64(&mut state);
            for b in 0..64 {
                self.signs[w * 64 + b] = if (word >> b) & 1 == 1 { 1.0 } else { -1.0 };
            }
        }
        let bits = self.bits;
        for (s, &d) in dosage.iter().enumerate() {
            if d < 0 {
                continue;
            }
            let v = d as f32 - mean;
            let row = &mut self.acc[s * bits..(s + 1) * bits];
            for (a, sg) in row.iter_mut().zip(&self.signs) {
                *a += v * sg;
            }
        }
        self.sites += 1;
    }

    /// The sign bits of every sample, `bits / 64` words each.
    pub fn finish(&self) -> Vec<Vec<u64>> {
        let words = self.bits / 64;
        (0..self.n)
            .map(|s| {
                let row = &self.acc[s * self.bits..(s + 1) * self.bits];
                (0..words)
                    .map(|w| {
                        let mut word = 0u64;
                        for b in 0..64 {
                            if row[w * 64 + b] > 0.0 {
                                word |= 1u64 << b;
                            }
                        }
                        word
                    })
                    .collect()
            })
            .collect()
    }
}

pub fn hamming(a: &[u64], b: &[u64]) -> u32 {
    a.iter().zip(b).map(|(x, y)| (x ^ y).count_ones()).sum()
}

/// For every sample its `k` nearest others by Hamming distance (ties by
/// index), as `(hi, lo)` pairs with `lo < hi`, sorted and deduplicated.
pub fn nearest_pairs(sketches: &[Vec<u64>], k: usize) -> Vec<(usize, usize)> {
    let n = sketches.len();
    let mut out: Vec<(usize, usize)> = Vec::new();
    let mut dist: Vec<(u32, usize)> = Vec::with_capacity(n);
    for i in 0..n {
        dist.clear();
        dist.extend(
            (0..n)
                .filter(|&j| j != i)
                .map(|j| (hamming(&sketches[i], &sketches[j]), j)),
        );
        let take = k.min(dist.len());
        if take < dist.len() {
            dist.select_nth_unstable(take);
            dist.truncate(take);
        }
        dist.sort_unstable();
        for &(_, j) in &dist {
            out.push((i.max(j), i.min(j)));
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

#[cfg(test)]
#[path = "../tests/unit/sketch.rs"]
mod tests;
