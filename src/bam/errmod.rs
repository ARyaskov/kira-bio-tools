//! Genotype likelihoods with dependent errors: port of htslib `errmod.c`
//! (Li 2011, the model behind `bcftools mpileup` PLs). Bases are packed as
//! `q << 5 | strand << 4 | base`; repeated observations of the same base on
//! the same strand are down-weighted by `fk[k] = (1 - depcorr)^k (1 - eps) + eps`.

use std::sync::OnceLock;

/// `bcf_call_init(theta = 0.83)`: dependency coefficient `1 - theta`.
const DEPCORR: f64 = 1.0 - 0.83;
const EPS: f64 = 0.03;
/// Sites deeper than this are subsampled, as htslib does.
const MAX_BASES: usize = 255;

struct Coef {
    fk: Vec<f64>,
    /// `beta[q << 16 | n << 8 | k]`.
    beta: Vec<f64>,
    /// `lhet[n << 8 | k]`.
    lhet: Vec<f64>,
}

fn lgamma(x: f64) -> f64 {
    // Lanczos approximation (g = 7), accurate to ~1e-13 for x > 0.
    const G: f64 = 7.0;
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        let pi = std::f64::consts::PI;
        return (pi / (pi * x).sin()).ln() - lgamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = C[0];
    let t = x + G + 0.5;
    for (i, c) in C.iter().enumerate().skip(1) {
        a += c / (x + i as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

fn cal_coef(depcorr: f64, eps: f64) -> Coef {
    let mut fk = vec![0.0; 256];
    fk[0] = 1.0;
    for (k, v) in fk.iter_mut().enumerate().skip(1) {
        *v = (1.0 - depcorr).powi(k as i32) * (1.0 - eps) + eps;
    }
    let mut lc = vec![0.0f64; 256 * 256];
    for n in 1..256usize {
        let lgn = lgamma(n as f64 + 1.0);
        for k in 1..=n {
            lc[n << 8 | k] = lgn - lgamma(k as f64 + 1.0) - lgamma((n - k) as f64 + 1.0);
        }
    }
    // Binomial tails in log space, as errmod.c does, so deep sites with
    // high qualities do not underflow.
    let mut beta = vec![0.0f64; 256 * 256 * 64];
    for q in 1..64usize {
        let e = 10f64.powf(-(q as f64) / 10.0);
        let le = e.ln();
        let le1 = (1.0 - e).ln();
        for n in 1..=255usize {
            let base = q << 16 | n << 8;
            let mut sum1 = lc[n << 8 | n] + n as f64 * le;
            beta[base + n] = f64::INFINITY;
            let mut k = n as i64 - 1;
            while k >= 0 {
                let ku = k as usize;
                let sum = sum1 + (lc[n << 8 | ku] + ku as f64 * le + (n - ku) as f64 * le1 - sum1).exp().ln_1p();
                beta[base + ku] = -10.0 / std::f64::consts::LN_10 * (sum1 - sum);
                sum1 = sum;
                k -= 1;
            }
        }
    }
    let mut lhet = vec![0.0f64; 256 * 256];
    for n in 0..256usize {
        for k in 0..256usize {
            lhet[n << 8 | k] = lc[n << 8 | k] - std::f64::consts::LN_2 * n as f64;
        }
    }
    Coef { fk, beta, lhet }
}

static COEF: OnceLock<Coef> = OnceLock::new();

/// `drand48` (48-bit LCG), seeded like `srand48(0)`, for the deep-site subsample.
struct Drand48(u64);

impl Drand48 {
    fn new(seed: u32) -> Self {
        Self(((seed as u64) << 16) | 0x330E)
    }
    fn next(&mut self) -> f64 {
        self.0 = (self.0.wrapping_mul(0x5DEE_CE66D).wrapping_add(0xB)) & ((1u64 << 48) - 1);
        self.0 as f64 / (1u64 << 48) as f64
    }
}

thread_local! {
    static RNG: std::cell::RefCell<Drand48> = std::cell::RefCell::new(Drand48::new(0));
}

/// Pack one observation: quality (clamped to 4..=63 like `bcf_call_glfgen`),
/// strand (1 = reverse) and base code.
#[inline]
pub fn pack_base(q: u8, reverse: bool, base: u8) -> u16 {
    let q = q.clamp(4, 63) as u16;
    (q << 5) | ((reverse as u16) << 4) | (base as u16 & 0xf)
}

#[derive(Clone, Copy)]
pub struct ErrorModel {
    coef: &'static Coef,
}

impl Default for ErrorModel {
    fn default() -> Self {
        Self::new()
    }
}

impl ErrorModel {
    pub fn new() -> Self {
        Self { coef: COEF.get_or_init(|| cal_coef(DEPCORR, EPS)) }
    }

    /// `errmod_cal`: phred-scaled genotype likelihoods `q[i * m + j]` over
    /// `m` base codes. `bases` is reordered (sorted, subsampled when deep).
    pub fn cal(&self, bases: &mut Vec<u16>, m: usize) -> Vec<f32> {
        let mut q = vec![0.0f32; m * m];
        if bases.is_empty() || m == 0 || m > 16 {
            return q;
        }
        let mut n = bases.len();
        if n > MAX_BASES {
            RNG.with(|r| {
                let mut r = r.borrow_mut();
                let mut i = n;
                while i > 1 {
                    let j = (r.next() * i as f64) as usize;
                    bases.swap(j, i - 1);
                    i -= 1;
                }
            });
            n = MAX_BASES;
        }
        bases[..n].sort_unstable();
        let coef = self.coef;
        let mut w = [0u32; 32];
        let mut fsum = [0.0f64; 16];
        let mut bsum = [0.0f64; 16];
        let mut c = [0u32; 16];
        for j in (0..n).rev() {
            let b = bases[j];
            let mut qv = (b >> 5) as usize;
            if qv < 4 {
                qv = 4;
            }
            if qv > 63 {
                qv = 63;
            }
            let k = (b & 0x1f) as usize;
            let kb = k & 0xf;
            let fkw = coef.fk[(w[k] as usize).min(255)];
            fsum[kb] += fkw;
            bsum[kb] += fkw * coef.beta[qv << 16 | n << 8 | (c[kb] as usize).min(255)];
            c[kb] += 1;
            w[k] += 1;
        }
        for j in 0..m {
            let (mut tmp1, mut tmp3, mut tmp2) = (0.0f64, 0.0f64, 0u32);
            for k in 0..m {
                if k == j {
                    continue;
                }
                tmp1 += bsum[k];
                tmp2 += c[k];
                tmp3 += fsum[k];
            }
            let _ = tmp3;
            if tmp2 > 0 {
                q[j * m + j] = tmp1 as f32;
            }
            for k in j + 1..m {
                let cjk = (c[j] + c[k]) as usize;
                let (mut tmp1, mut tmp2) = (0.0f64, 0u32);
                for i in 0..m {
                    if i == j || i == k {
                        continue;
                    }
                    tmp1 += bsum[i];
                    tmp2 += c[i];
                }
                let lhet = coef.lhet[cjk.min(255) << 8 | (c[k] as usize).min(255)];
                let v = if tmp2 > 0 { -4.343 * lhet + tmp1 } else { -4.343 * lhet };
                q[j * m + k] = v as f32;
                q[k * m + j] = v as f32;
            }
            for k in 0..m {
                if q[j * m + k] < 0.0 {
                    q[j * m + k] = 0.0;
                }
            }
        }
        q
    }

    /// PLs for the genotypes over `alleles` (base codes in output allele
    /// order, VCF genotype order), min-normalised and capped at 255 like
    /// `bcf_call_combine`.
    pub fn pls(&self, bases: &mut Vec<u16>, alleles: &[u8]) -> Vec<u32> {
        let m = 1 + alleles.iter().copied().max().unwrap_or(0) as usize;
        let m = m.max(alleles.len()).min(16);
        let q = self.cal(bases, m);
        let n = alleles.len();
        let mut raw: Vec<f32> = Vec::with_capacity(n * (n + 1) / 2);
        for j in 0..n {
            for i in 0..=j {
                let (a, b) = (alleles[i] as usize, alleles[j] as usize);
                raw.push(q[a * m + b]);
            }
        }
        let min = raw.iter().copied().fold(f32::INFINITY, f32::min);
        raw.iter().map(|&x| ((x - min + 0.499) as i64).clamp(0, 255) as u32).collect()
    }
}

#[cfg(test)]
#[path = "../../tests/unit/bam_errmod.rs"]
mod tests;
