//! logsumexp hot path for BAQ HMM. Uses std::f32::{exp,ln_1p}; LLVM lowers these
//! to libmvec/SVML where available, which beats hand-rolled polynomial wrappers
//! once you account for inlining and register pressure inside the HMM kernels.

/// log(exp(a) + exp(b)). Stable for |a-b| up to ~80.
#[inline(always)]
pub fn logsumexp2(a: f32, b: f32) -> f32 {
    if !a.is_finite() {
        return b;
    }
    if !b.is_finite() {
        return a;
    }
    let (m, d) = if a >= b { (a, b - a) } else { (b, a - b) };
    if d < -80.0 {
        m
    } else {
        m + (d.exp()).ln_1p()
    }
}

#[inline(always)]
pub fn logsumexp3(a: f32, b: f32, c: f32) -> f32 {
    logsumexp2(logsumexp2(a, b), c)
}
