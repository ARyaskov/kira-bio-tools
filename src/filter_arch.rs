// Unified SIMD module for:
//  - FAST split key=value
//  - FAST integer parsing
//  - FAST float parsing
//  - FAST compare
//  - FAST info=...; parser (SIMD / fallback)
//

use std::collections::HashMap;

//
// ---------------------------------------------------
// PUBLIC API
// ---------------------------------------------------
//

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CmpOp {
    Eq,
    Ne,
    Lt,
    Gt,
    Le,
    Ge,
}

pub trait FilterArch {
    fn split_key_value<'a>(&self, s: &'a str) -> Option<(&'a str, &'a str)>;

    fn parse_int(&self, s: &str) -> Option<i64>;
    fn parse_float(&self, s: &str) -> Option<f64>;

    fn cmp_i64(&self, a: i64, op: CmpOp, b: i64) -> bool;
    fn cmp_f64(&self, a: f64, op: CmpOp, b: f64) -> bool;

    fn parse_info_simd(&self, s: &str) -> HashMap<String, String>;
}

//
// ---------------------------------------------------
// SELECT RUNTIME ARCH
// ---------------------------------------------------
//

// Inline module — fallback
mod arch_fallback {
    use super::*;

    pub struct FallbackImpl;

    impl FilterArch for FallbackImpl {
        #[inline]
        fn split_key_value<'a>(&self, s: &'a str) -> Option<(&'a str, &'a str)> {
            s.split_once('=')
        }

        #[inline]
        fn parse_int(&self, s: &str) -> Option<i64> {
            s.parse().ok()
        }

        #[inline]
        fn parse_float(&self, s: &str) -> Option<f64> {
            s.parse().ok()
        }

        #[inline]
        fn cmp_i64(&self, a: i64, op: CmpOp, b: i64) -> bool {
            match op {
                CmpOp::Eq => a == b,
                CmpOp::Ne => a != b,
                CmpOp::Lt => a < b,
                CmpOp::Gt => a > b,
                CmpOp::Le => a <= b,
                CmpOp::Ge => a >= b,
            }
        }

        #[inline]
        fn cmp_f64(&self, a: f64, op: CmpOp, b: f64) -> bool {
            match op {
                CmpOp::Eq => a == b,
                CmpOp::Ne => a != b,
                CmpOp::Lt => a < b,
                CmpOp::Gt => a > b,
                CmpOp::Le => a <= b,
                CmpOp::Ge => a >= b,
            }
        }

        #[inline]
        fn parse_info_simd(&self, s: &str) -> HashMap<String, String> {
            let mut map = HashMap::new();
            for item in s.split(';') {
                if let Some((k, v)) = item.split_once('=') {
                    map.insert(k.to_string(), v.to_string());
                }
            }
            map
        }
    }
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
mod arch_avx2 {
    use super::*;
    use core::arch::x86_64::*;

    pub struct Avx2Impl;

    impl FilterArch for Avx2Impl {
        #[inline]
        fn split_key_value<'a>(&self, s: &'a str) -> Option<(&'a str, &'a str)> {
            unsafe { split_eq_avx2(s) }
        }

        #[inline]
        fn parse_int(&self, s: &str) -> Option<i64> {
            fast_atoi(s)
        }

        #[inline]
        fn parse_float(&self, s: &str) -> Option<f64> {
            fast_atof(s)
        }

        #[inline]
        fn cmp_i64(&self, a: i64, op: CmpOp, b: i64) -> bool {
            arch_fallback::FallbackImpl.cmp_i64(a, op, b)
        }

        #[inline]
        fn cmp_f64(&self, a: f64, op: CmpOp, b: f64) -> bool {
            arch_fallback::FallbackImpl.cmp_f64(a, op, b)
        }

        #[inline]
        fn parse_info_simd(&self, s: &str) -> HashMap<String, String> {
            arch_fallback::FallbackImpl.parse_info_simd(s)
        }
    }

    #[target_feature(enable = "avx2")]
    pub unsafe fn split_eq_avx2<'a>(s: &'a str) -> Option<(&'a str, &'a str)> {
        use core::arch::x86_64::*;

        let bytes = s.as_bytes();
        let len = bytes.len();

        let needle = _mm256_set1_epi8(b'=' as i8);

        let mut i = 0;

        while i + 32 <= len {
            // Load 32 bytes
            let ptr = bytes.as_ptr().add(i) as *const __m256i;
            let chunk = _mm256_loadu_si256(ptr);

            // Compare each byte with '='
            let mask = _mm256_movemask_epi8(_mm256_cmpeq_epi8(chunk, needle));

            if mask != 0 {
                // Position of the first '=' inside this 32-byte window
                let pos_in_chunk = mask.trailing_zeros() as usize;
                let pos = i + pos_in_chunk;

                return Some((&s[..pos], &s[pos + 1..]));
            }

            i += 32;
        }

        // Tail fallback
        if let Some(p) = bytes[i..].iter().position(|&c| c == b'=') {
            let pos = i + p;
            return Some((&s[..pos], &s[pos + 1..]));
        }

        None
    }

    #[inline]
    fn fast_atoi(s: &str) -> Option<i64> {
        let mut neg = false;
        let bytes = s.as_bytes();
        let mut i = 0;

        if bytes.get(0) == Some(&b'-') {
            neg = true;
            i += 1;
        }

        let mut n: i64 = 0;
        while i < bytes.len() {
            let c = bytes[i];
            if c < b'0' || c > b'9' {
                return None;
            }
            n = n * 10 + (c - b'0') as i64;
            i += 1;
        }

        Some(if neg { -n } else { n })
    }

    #[inline]
    fn fast_atof(s: &str) -> Option<f64> {
        let mut neg = false;
        let bytes = s.as_bytes();
        let mut i = 0;

        if bytes.get(0) == Some(&b'-') {
            neg = true;
            i += 1;
        }

        let mut int: i64 = 0;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            int = int * 10 + (bytes[i] - b'0') as i64;
            i += 1;
        }

        let mut frac = 0.0;
        let mut base = 0.1;

        if i < bytes.len() && bytes[i] == b'.' {
            i += 1;
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                frac += (bytes[i] - b'0') as f64 * base;
                base *= 0.1;
                i += 1;
            }
        }

        let val = int as f64 + frac;
        Some(if neg { -val } else { val })
    }
}

#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
mod arch_neon {
    use super::*;
    use core::arch::aarch64::*;

    pub struct NeonImpl;

    impl FilterArch for NeonImpl {
        fn split_key_value(&self, s: &str) -> Option<(&str, &str)> {
            unsafe { split_eq_neon(s) }
        }

        fn parse_int(&self, s: &str) -> Option<i64> {
            fast_atoi(s)
        }

        fn parse_float(&self, s: &str) -> Option<f64> {
            fast_atof(s)
        }

        fn cmp_i64(&self, a: i64, op: CmpOp, b: i64) -> bool {
            crate::filter_arch::arch_fallback::FallbackImpl.cmp_i64(a, op, b)
        }

        fn cmp_f64(&self, a: f64, op: CmpOp, b: f64) -> bool {
            crate::filter_arch::arch_fallback::FallbackImpl.cmp_f64(a, op, b)
        }

        fn parse_info_simd(&self, s: &str) -> HashMap<String, String> {
            crate::filter_arch::arch_fallback::FallbackImpl.parse_info_simd(s)
        }
    }

    #[target_feature(enable = "neon")]
    unsafe fn split_eq_neon(s: &str) -> Option<(&str, &str)> {
        let bytes = s.as_bytes();
        let len = bytes.len();
        let needle = vdupq_n_u8(b'=');

        let mut i = 0;
        while i + 16 <= len {
            let chunk = vld1q_u8(bytes.as_ptr().add(i));
            let mask = vceqq_u8(chunk, needle);
            if vmaxvq_u8(mask) != 0 {
                for j in 0..16 {
                    if bytes[i + j] == b'=' {
                        return Some((&s[..i + j], &s[i + j + 1..]));
                    }
                }
            }
            i += 16;
        }
        s.split_once('=')
    }

    #[inline]
    fn fast_atoi(s: &str) -> Option<i64> {
        let mut neg = false;
        let bytes = s.as_bytes();
        let mut i = 0;

        if bytes.get(0) == Some(&b'-') {
            neg = true;
            i += 1;
        }

        let mut n = 0_i64;
        while i < bytes.len() {
            let c = bytes[i];
            if c < b'0' || c > b'9' {
                return None;
            }
            n = n * 10 + (c - b'0') as i64;
            i += 1;
        }

        Some(if neg { -n } else { n })
    }

    #[inline]
    fn fast_atof(s: &str) -> Option<f64> {
        let mut neg = false;
        let bytes = s.as_bytes();
        let mut i = 0;

        if bytes.get(0) == Some(&b'-') {
            neg = true;
            i += 1;
        }

        let mut int = 0_i64;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            int = int * 10 + (bytes[i] - b'0') as i64;
            i += 1;
        }

        let mut frac = 0.0_f64;
        let mut base = 0.1;

        if i < bytes.len() && bytes[i] == b'.' {
            i += 1;
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                frac += (bytes[i] - b'0') as f64 * base;
                base *= 0.1;
                i += 1;
            }
        }

        let val = int as f64 + frac;
        Some(if neg { -val } else { val })
    }
}


#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
static ARCH: arch_avx2::Avx2Impl = arch_avx2::Avx2Impl;

#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
static ARCH: arch_neon::NeonImpl = arch_neon::NeonImpl;

#[cfg(not(any(
    all(target_arch = "x86_64", target_feature = "avx2"),
    all(target_arch = "aarch64", target_feature = "neon")
)))]
static ARCH: arch_fallback::FallbackImpl = arch_fallback::FallbackImpl;

pub fn get_arch() -> &'static dyn FilterArch {
    &ARCH
}
