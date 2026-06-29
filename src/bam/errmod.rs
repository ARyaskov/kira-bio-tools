//! Genotype likelihoods + Hardy-Weinberg prior + QUAL.
//!
//! Likelihoods are kept as natural-log values (`log_p`), normalised so the most
//! likely genotype is 0. Integer PLs are derived for output; QUAL is computed
//! from the (unsaturated) log posterior so it is a true float, not a `u8`.

const MAX_QUAL: usize = 64;
const PL_CAP: u32 = 255;

pub struct ErrorModel {
    log_factorial: Vec<f64>,
}

#[derive(Clone)]
pub struct GenotypeLikelihoods {
    /// Per-genotype log-likelihood (natural log), normalised so max == 0.
    pub log_p: Vec<f64>,
    /// Phred-scaled likelihoods for output (min == 0, capped at 255).
    pub pl: Vec<u32>,
}

impl GenotypeLikelihoods {
    /// Build from raw natural-log likelihoods (any scale); normalises so the
    /// best genotype is 0 and derives capped integer PLs.
    fn from_log_p(mut log_p: Vec<f64>) -> Self {
        let max = log_p.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        if max.is_finite() {
            for l in log_p.iter_mut() {
                *l -= max;
            }
        } else {
            for l in log_p.iter_mut() {
                *l = 0.0;
            }
        }
        let pl = log_p
            .iter()
            .map(|&l| {
                let phred = (-10.0 * l / std::f64::consts::LN_10).round();
                phred.clamp(0.0, PL_CAP as f64) as u32
            })
            .collect();
        Self { log_p, pl }
    }
}

impl ErrorModel {
    pub fn new() -> Self {
        let mut lf = vec![0.0; 4096];
        for i in 2..lf.len() {
            lf[i] = lf[i - 1] + (i as f64).ln();
        }
        Self { log_factorial: lf }
    }

    /// Count-based likelihoods: one averaged quality per allele. Used by the
    /// indel / synthetic-quality paths where per-read bases are not available.
    pub fn likelihoods(
        &self,
        n_alleles: usize,
        counts: &[u32],
        quals: &[u32],
    ) -> GenotypeLikelihoods {
        let n_gt = n_alleles * (n_alleles + 1) / 2;
        let total: u32 = counts.iter().sum();
        if total == 0 {
            return GenotypeLikelihoods::from_log_p(vec![0.0; n_gt]);
        }

        let avg_q = |i: usize| -> f64 {
            if counts[i] == 0 {
                30.0
            } else {
                (quals[i] as f64 / counts[i] as f64).min(MAX_QUAL as f64)
            }
        };

        let mut log_p: Vec<f64> = vec![f64::NEG_INFINITY; n_gt];
        let mut gt_idx = 0usize;
        for j in 0..n_alleles {
            for i in 0..=j {
                let mut ll = 0.0f64;
                for k in 0..n_alleles {
                    let c = counts[k] as f64;
                    if c == 0.0 {
                        continue;
                    }
                    let e = 10f64.powf(-avg_q(k) / 10.0);
                    ll += c * allele_prob(i, j, k, e).max(1e-30).ln();
                }
                log_p[gt_idx] = ll;
                gt_idx += 1;
            }
        }
        let _ = &self.log_factorial;
        GenotypeLikelihoods::from_log_p(log_p)
    }

    /// Exact per-read genotype likelihoods. `reads` is `(allele_index, qual)`
    /// per observed base, with quality already capped by mapping quality
    /// (samtools-style `min(BQ, MAPQ)`). This is the correct diploid model —
    /// a product over reads, not a single averaged quality per allele.
    pub fn likelihoods_per_read(
        &self,
        n_alleles: usize,
        reads: &[(usize, u8)],
    ) -> GenotypeLikelihoods {
        let n_gt = n_alleles * (n_alleles + 1) / 2;
        if reads.is_empty() {
            return GenotypeLikelihoods::from_log_p(vec![0.0; n_gt]);
        }
        // Pre-compute ln-probabilities per read once.
        let err: Vec<f64> = reads
            .iter()
            .map(|&(_, q)| 10f64.powf(-(q as f64) / 10.0))
            .collect();

        let mut log_p: Vec<f64> = vec![0.0; n_gt];
        let mut gt_idx = 0usize;
        for j in 0..n_alleles {
            for i in 0..=j {
                let mut ll = 0.0f64;
                for (r, &(a, _)) in reads.iter().enumerate() {
                    ll += allele_prob(i, j, a, err[r]).max(1e-30).ln();
                }
                log_p[gt_idx] = ll;
                gt_idx += 1;
            }
        }
        GenotypeLikelihoods::from_log_p(log_p)
    }
}

/// P(observed allele `a` | diploid genotype `i/j`, per-base error `e`).
#[inline]
fn allele_prob(i: usize, j: usize, a: usize, e: f64) -> f64 {
    if i == a && j == a {
        1.0 - e
    } else if i == a || j == a {
        0.5 * (1.0 - e) + 0.5 * (e / 3.0)
    } else {
        e / 3.0
    }
}

impl GenotypeLikelihoods {
    pub fn to_pl_string(&self) -> String {
        self.pl
            .iter()
            .map(u32::to_string)
            .collect::<Vec<_>>()
            .join(",")
    }

    pub fn most_likely_gt(&self, n_alleles: usize) -> (usize, usize) {
        let mut best_idx = 0;
        let mut best = f64::NEG_INFINITY;
        for (k, &v) in self.log_p.iter().enumerate() {
            if v > best {
                best = v;
                best_idx = k;
            }
        }
        for j in 0..n_alleles {
            for i in 0..=j {
                if best_idx == 0 {
                    return (i, j);
                }
                best_idx -= 1;
            }
        }
        (0, 0)
    }

    /// Apply Hardy-Weinberg prior with given per-allele variant rate.
    /// Returns posterior likelihoods (data + prior), renormalised.
    pub fn with_prior(&self, n_alleles: usize, alt_rate: f64) -> GenotypeLikelihoods {
        let ref_freq = (1.0 - alt_rate).max(1e-12);
        let alt_freq_each = if n_alleles > 1 {
            alt_rate / (n_alleles as f64 - 1.0)
        } else {
            0.0
        };
        let freq = |i: usize| if i == 0 { ref_freq } else { alt_freq_each };

        let mut post: Vec<f64> = Vec::with_capacity(self.log_p.len());
        let mut gt_idx = 0;
        for j in 0..n_alleles {
            for i in 0..=j {
                let p_prior = freq(i) * freq(j) * if i == j { 1.0 } else { 2.0 };
                post.push(self.log_p[gt_idx] + p_prior.max(1e-30).ln());
                gt_idx += 1;
            }
        }
        GenotypeLikelihoods::from_log_p(post)
    }

    /// Variant QUAL = -10·log10 P(genotype 0/0 | data), as a true float.
    /// 0 ⇒ reference is the best call; large ⇒ strong variant signal.
    /// Call after [`with_prior`] for a posterior-based QUAL.
    pub fn qual(&self) -> f64 {
        if self.log_p.is_empty() {
            return 0.0;
        }
        // log_p is normalised so max == 0; convert to posterior over genotypes.
        let z: f64 = self.log_p.iter().map(|&l| l.exp()).sum();
        if z <= 0.0 {
            return 0.0;
        }
        let p_ref = self.log_p[0].exp() / z;
        (-10.0 * p_ref.max(1e-30).log10()).max(0.0)
    }

    /// True if the best posterior GT is not 0/0.
    pub fn is_variant(&self, n_alleles: usize) -> bool {
        let (gi, gj) = self.most_likely_gt(n_alleles);
        gi != 0 || gj != 0
    }
}

#[cfg(test)]
#[path = "../../tests/unit/bam_errmod.rs"]
mod tests;
