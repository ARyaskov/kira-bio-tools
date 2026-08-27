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
                return unsafe { x86_avx2::parse_vcf_line_simd(line) };
            }
        }

        #[cfg(target_arch = "aarch64")]
        {
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
                return unsafe { x86_avx2::parse_vcf_fields_avx2(line) };
            }
        }

        #[cfg(target_arch = "aarch64")]
        {
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
                return unsafe { x86_avx2::parse_chr_pos_avx2(line) };
            }
        }

        #[cfg(target_arch = "aarch64")]
        {
            return unsafe { arm_neon::parse_chr_pos_neon(line) };
        }

        #[cfg(not(target_arch = "aarch64"))]
        fallback::parse_chr_pos_scalar(line)
    }
}

#[inline(always)]
pub fn parse_u32_bytes(bytes: &[u8]) -> Option<u32> {
    #[cfg(target_arch = "x86_64")]
    {
        if is_simd_available() {
            return Some(unsafe { x86_avx2::parse_u32_fast(bytes.as_ptr(), bytes.len()) });
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        return Some(unsafe { arm_neon::parse_u32_fast(bytes.as_ptr(), bytes.len()) });
    }

    #[cfg(not(target_arch = "aarch64"))]
    fallback::parse_u32_scalar(bytes)
}
