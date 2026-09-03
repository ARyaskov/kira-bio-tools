use crate::vcf::structs::{VcfFields, VcfFieldsFull};

#[cfg(target_arch = "x86_64")]
pub mod x86_avx2;

#[cfg(target_arch = "aarch64")]
pub mod arm_neon;

pub mod fallback;

#[cfg(target_arch = "x86_64")]
static SIMD_AVAILABLE: std::sync::OnceLock<bool> = std::sync::OnceLock::new();

#[inline(always)]
#[cfg(target_arch = "x86_64")]
fn is_simd_available() -> bool {
    *SIMD_AVAILABLE.get_or_init(|| is_x86_feature_detected!("avx2"))
}

pub struct SimdVcfParser;

impl SimdVcfParser {
    #[inline(always)]
    pub fn parse_line(line: &[u8]) -> Option<VcfFieldsFull<'_>> {
        #[cfg(target_arch = "x86_64")]
        {
            if is_simd_available() {
                // SAFETY: AVX2 support was verified at runtime just above.
                return unsafe { x86_avx2::parse_vcf_line_simd(line) };
            }
        }

        #[cfg(target_arch = "aarch64")]
        {
            // SAFETY: NEON is a baseline feature of every aarch64 target.
            return unsafe { arm_neon::parse_vcf_line_simd(line) };
        }

        #[cfg(not(target_arch = "aarch64"))]
        fallback::parse_vcf_line_scalar(line)
    }

    #[inline(always)]
    pub fn parse_fields(line: &[u8]) -> Option<VcfFields<'_>> {
        #[cfg(target_arch = "x86_64")]
        {
            if is_simd_available() {
                // SAFETY: AVX2 support was verified at runtime just above.
                return unsafe { x86_avx2::parse_vcf_fields_avx2(line) };
            }
        }

        #[cfg(target_arch = "aarch64")]
        {
            // SAFETY: NEON is a baseline feature of every aarch64 target.
            return unsafe { arm_neon::parse_vcf_fields_neon(line) };
        }

        #[cfg(not(target_arch = "aarch64"))]
        fallback::parse_vcf_fields_scalar(line)
    }

    #[inline(always)]
    pub fn parse_chr_pos(line: &[u8]) -> Option<(u8, u32)> {
        #[cfg(target_arch = "x86_64")]
        {
            if is_simd_available() {
                // SAFETY: AVX2 support was verified at runtime just above.
                return unsafe { x86_avx2::parse_chr_pos_avx2(line) };
            }
        }

        #[cfg(target_arch = "aarch64")]
        {
            // SAFETY: NEON is a baseline feature of every aarch64 target.
            return unsafe { arm_neon::parse_chr_pos_neon(line) };
        }

        #[cfg(not(target_arch = "aarch64"))]
        fallback::parse_chr_pos_scalar(line)
    }
}

/// Strict decimal parse: `None` for empty input, a non-digit or overflow
/// (`12abc` is rejected rather than truncated to 12).
#[inline]
pub fn parse_u32_strict(bytes: &[u8]) -> Option<u32> {
    if bytes.is_empty() || bytes.len() > 10 {
        return None;
    }
    let mut x: u32 = 0;
    for &c in bytes {
        if !c.is_ascii_digit() {
            return None;
        }
        x = x.checked_mul(10)?.checked_add((c - b'0') as u32)?;
    }
    Some(x)
}

#[inline(always)]
pub fn parse_u32_bytes(bytes: &[u8]) -> Option<u32> {
    parse_u32_strict(bytes)
}
