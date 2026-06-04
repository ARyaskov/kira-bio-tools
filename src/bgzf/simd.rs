#![allow(unsafe_op_in_unsafe_fn)]

use super::structs::BGZF_HEADER;

#[cfg(target_arch = "x86_64")]
pub mod x86_simd {
    use super::*;
    use core::arch::x86_64::*;

    #[target_feature(enable = "avx2")]
    pub unsafe fn memcpy_avx(dst: *mut u8, src: *const u8, len: usize) {
        let mut i = 0;

        while i + 32 <= len {
            let v = _mm256_loadu_si256(src.add(i) as *const _);
            _mm256_storeu_si256(dst.add(i) as *mut _, v);
            i += 32;
        }

        while i < len {
            *dst.add(i) = *src.add(i);
            i += 1;
        }
    }

    #[target_feature(enable = "avx2")]
    pub unsafe fn copy_bgzf_header(dst: *mut u8) {
        let mut tmp = [0u8; 32];
        tmp[..18].copy_from_slice(&BGZF_HEADER);
        let v = _mm256_loadu_si256(tmp.as_ptr() as *const _);
        _mm256_storeu_si256(dst as *mut _, v);
    }
}

#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
pub mod arm_simd {
    use super::*;
    use core::arch::aarch64::*;

    #[target_feature(enable = "neon")]
    pub unsafe fn memcpy_neon(dst: *mut u8, src: *const u8, len: usize) {
        let mut i = 0;

        while i + 16 <= len {
            let v = vld1q_u8(src.add(i));
            vst1q_u8(dst.add(i), v);
            i += 16;
        }

        while i < len {
            *dst.add(i) = *src.add(i);
            i += 1;
        }
    }

    #[target_feature(enable = "neon")]
    pub unsafe fn copy_bgzf_header(dst: *mut u8) {
        let v1 = vld1q_u8(BGZF_HEADER.as_ptr());
        vst1q_u8(dst, v1);
        let v2 = vld1_u8(BGZF_HEADER.as_ptr().add(16));
        vst1_u8(dst.add(16), v2);
    }
}

#[inline(always)]
pub fn crc32_fallback(crc: u32, data: &[u8]) -> u32 {
    let mut hasher = crc32fast::Hasher::new_with_initial(crc);
    hasher.update(data);
    hasher.finalize()
}

#[inline(always)]
pub fn compute_crc32(data: &[u8]) -> u32 {
    crc32_fallback(0, data)
}

#[inline(always)]
pub unsafe fn fast_memcpy(dst: *mut u8, src: *const u8, len: usize) {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            return x86_simd::memcpy_avx(dst, src, len);
        }
    }

    #[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
    {
        return arm_simd::memcpy_neon(dst, src, len);
    }

    #[cfg(not(all(target_arch = "aarch64", target_feature = "neon")))]
    std::ptr::copy_nonoverlapping(src, dst, len);
}

#[inline(always)]
pub unsafe fn fast_copy_bgzf_header(dst: *mut u8) {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            return x86_simd::copy_bgzf_header(dst);
        }
    }

    #[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
    {
        return arm_simd::copy_bgzf_header(dst);
    }

    #[cfg(not(all(target_arch = "aarch64", target_feature = "neon")))]
    std::ptr::copy_nonoverlapping(BGZF_HEADER.as_ptr(), dst, 18);
}
