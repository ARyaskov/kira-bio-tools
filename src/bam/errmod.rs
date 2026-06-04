//! Genotype likelihoods + Hardy-Weinberg prior + QUAL.

const MAX_QUAL: usize = 64;

pub struct ErrorModel {
    log_factorial: Vec<f64>,
}

#[derive(Clone)]
pub struct GenotypeLikelihoods {
    pub pl: Vec<u8>,
}

impl ErrorModel {
    pub fn new() -> Self {
        let mut lf = vec![0.0; 4096];
        for i in 2..lf.len() {
            lf[i] = lf[i - 1] + (i as f64).ln();
        }
        Self { log_factorial: lf }
    }

    pub fn likelihoods(
        &self,
        n_alleles: usize,
        counts: &[u32],
        quals: &[u32],
    ) -> GenotypeLikelihoods {
        let n_gt = n_alleles * (n_alleles + 1) / 2;
        let mut log_p: Vec<f64> = vec![f64::NEG_INFINITY; n_gt];
        let total: u32 = counts.iter().sum();
        if total == 0 {
            return GenotypeLikelihoods { pl: vec![0; n_gt] };
        }

        let avg_q = |i: usize| -> f64 {
            if counts[i] == 0 {
                30.0
            } else {
                (quals[i] as f64 / counts[i] as f64).min(MAX_QUAL as f64)
            }
        };

        let mut gt_idx = 0usize;
        for j in 0..n_alleles {
            for i in 0..=j {
                let mut ll = 0.0f64;
                for k in 0..n_alleles {
                    let c = counts[k] as f64;
                    if c == 0.0 {
                        continue;
                    }
                    let qk = avg_q(k);
                    let e = 10f64.powf(-qk / 10.0);
                    let p = if i == k && j == k {
                        1.0 - e
                    } else if i == k || j == k {
                        0.5 * (1.0 - e) + 0.5 * (e / 3.0)
                    } else {
                        e / 3.0
                    };
                    ll += c * p.max(1e-30).ln();
                }
                log_p[gt_idx] = ll;
                gt_idx += 1;
            }
        }

        let max_ll = log_p.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let pl: Vec<u8> = log_p
            .iter()
            .map(|&l| {
                let phred = (-10.0 * (l - max_ll) / std::f64::consts::LN_10).round();
                phred.max(0.0).min(255.0) as u8
            })
            .collect();
        let _ = &self.log_factorial;
        GenotypeLikelihoods { pl }
    }
}

impl GenotypeLikelihoods {
    pub fn to_pl_string(&self) -> String {
        self.pl
            .iter()
            .map(u8::to_string)
            .collect::<Vec<_>>()
            .join(",")
    }

    pub fn most_likely_gt(&self, n_alleles: usize) -> (usize, usize) {
        let mut best_idx = 0;
        let mut best_pl: u8 = u8::MAX;
        for (k, &v) in self.pl.iter().enumerate() {
            if v < best_pl {
                best_pl = v;
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
    /// Returns posterior PL = data PL + prior PL, normalised so min = 0.
    pub fn with_prior(&self, n_alleles: usize, alt_rate: f64) -> GenotypeLikelihoods {
        let ref_freq = (1.0 - alt_rate).max(1e-12);
        let alt_freq_each = if n_alleles > 1 {
            alt_rate / (n_alleles as f64 - 1.0)
        } else {
            0.0
        };
        let freq = |i: usize| if i == 0 { ref_freq } else { alt_freq_each };

        let n_gt = n_alleles * (n_alleles + 1) / 2;
        let mut post: Vec<f64> = Vec::with_capacity(n_gt);
        let mut gt_idx = 0;
        for j in 0..n_alleles {
            for i in 0..=j {
                let p_prior = freq(i) * freq(j) * if i == j { 1.0 } else { 2.0 };
                let prior_pl = -10.0 * p_prior.max(1e-30).log10();
                let pl_post = self.pl[gt_idx] as f64 + prior_pl;
                post.push(pl_post);
                gt_idx += 1;
            }
        }
        let min_pl = post.iter().cloned().fold(f64::INFINITY, f64::min);
        let pl: Vec<u8> = post
            .iter()
            .map(|&p| (p - min_pl).round().clamp(0.0, 255.0) as u8)
            .collect();
        GenotypeLikelihoods { pl }
    }

    /// Variant QUAL = posterior PL of GT 0/0 (after `with_prior`).
    /// 0 = ref is the best call; large = strong variant signal.
    pub fn qual(&self) -> u8 {
        self.pl.first().copied().unwrap_or(0)
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
