pub const PL2P_SIZE: usize = 256;
pub const LN10_OVER_10: f64 = 0.23025850929940458;
pub const M10_OVER_LN10: f64 = -4.342944819032518;

pub fn init_pl2p() -> [f64; PL2P_SIZE] {
    let mut t = [0.0f64; PL2P_SIZE];
    for i in 0..PL2P_SIZE { t[i] = 10f64.powf(-(i as f64) / 10.0); }
    t
}

#[inline]
pub fn pl_to_prob(pl: i32, table: &[f64; PL2P_SIZE]) -> f64 {
    if pl < 0 { 1.0 }
    else if (pl as usize) < PL2P_SIZE { table[pl as usize] }
    else { 10f64.powf(-(pl as f64) / 10.0) }
}

#[inline]
pub fn log10_sum_exp(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY { return b; }
    if b == f64::NEG_INFINITY { return a; }
    let m = a.max(b);
    m + (10f64.powf(a - m) + 10f64.powf(b - m)).log10()
}

/// Watterson factor: 1 + 1/2 + 1/3 + ... + 1/(2N-1) for 2N alleles.
/// Used to scale theta prior in mcall_init().
pub fn watterson_factor(n_alleles_total: usize) -> f64 {
    if n_alleles_total < 2 { return 1.0; }
    let mut s = 1.0;
    for i in 2..n_alleles_total { s += 1.0 / (i as f64); }
    s
}

/// Index for diploid genotype (i,j) where i<=j in canonical VCF PL order:
/// 0/0, 0/1, 1/1, 0/2, 1/2, 2/2, 0/3, 1/3, 2/3, 3/3, ...
#[inline]
pub fn gt_index(i: usize, j: usize) -> usize {
    let (lo, hi) = if i <= j { (i, j) } else { (j, i) };
    hi * (hi + 1) / 2 + lo
}

/// Number of diploid genotypes for n alleles: n*(n+1)/2.
#[inline]
pub fn n_genotypes(n_alleles: usize) -> usize { n_alleles * (n_alleles + 1) / 2 }

#[cfg(test)]
#[path = "../../tests/unit/call_math.rs"]
mod tests;
