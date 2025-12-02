use core::arch::global_asm;

/// Runtime CPU feature detection
pub struct NormContext {
    use_avx2: bool,
    use_neon: bool,
}

impl NormContext {
    pub fn detect() -> Self {
        Self {
            #[cfg(target_arch = "x86_64")]
            use_avx2: is_x86_feature_detected!("avx2"),
            #[cfg(not(target_arch = "x86_64"))]
            use_avx2: false,

            #[cfg(target_arch = "aarch64")]
            use_neon: true, // NEON is always available on aarch64
            #[cfg(not(target_arch = "aarch64"))]
            use_neon: false,
        }
    }

    /// Check if AVX2 is available
    #[inline]
    pub fn has_avx2(&self) -> bool {
        self.use_avx2
    }

    /// Check if NEON is available
    #[inline]
    pub fn has_neon(&self) -> bool {
        self.use_neon
    }

    /// Get feature description string
    pub fn features(&self) -> String {
        let mut features = Vec::new();
        if self.use_avx2 {
            features.push("AVX2");
        }
        if self.use_neon {
            features.push("NEON");
        }
        if features.is_empty() {
            features.push("Scalar");
        }
        features.join(", ")
    }

    #[inline]
    pub fn normalize(&self, ref_allele: &str, alt_allele: &str) -> (String, String, usize, usize) {
        #[cfg(target_arch = "x86_64")]
        {
            if self.use_avx2 {
                return normalize_avx2(ref_allele, alt_allele);
            }
        }

        #[cfg(target_arch = "aarch64")]
        {
            if self.use_neon {
                return normalize_neon(ref_allele, alt_allele);
            }
        }

        normalize_scalar(ref_allele, alt_allele)
    }
}

/// Scalar fallback implementation
#[inline]
pub fn normalize_scalar(ref_allele: &str, alt_allele: &str) -> (String, String, usize, usize) {
    let prefix = common_prefix_scalar(ref_allele.as_bytes(), alt_allele.as_bytes());
    let (r_tail, a_tail) = (
        &ref_allele.as_bytes()[prefix..],
        &alt_allele.as_bytes()[prefix..],
    );
    let suffix = common_suffix_scalar(r_tail, a_tail);
    let new_ref = String::from_utf8(r_tail[..r_tail.len() - suffix].to_vec()).unwrap();
    let new_alt = String::from_utf8(a_tail[..a_tail.len() - suffix].to_vec()).unwrap();
    (new_ref, new_alt, prefix, suffix)
}

#[inline]
fn common_prefix_scalar(a: &[u8], b: &[u8]) -> usize {
    let len = a.len().min(b.len());
    let mut i = 0;
    while i < len && a[i] == b[i] {
        i += 1;
    }
    i
}

#[inline]
fn common_suffix_scalar(a: &[u8], b: &[u8]) -> usize {
    let mut i = 0;
    let len = a.len().min(b.len());
    while i < len {
        if a[a.len() - 1 - i] != b[b.len() - 1 - i] {
            break;
        }
        i += 1;
    }
    i
}

/// AVX2 optimized implementation
#[cfg(target_arch = "x86_64")]
#[inline]
pub fn normalize_avx2(ref_allele: &str, alt_allele: &str) -> (String, String, usize, usize) {
    use core::arch::x86_64::*;

    unsafe {
        let rb = ref_allele.as_bytes();
        let ab = alt_allele.as_bytes();

        let prefix = common_prefix_avx2(rb, ab);
        let r_tail = &rb[prefix..];
        let a_tail = &ab[prefix..];

        let suffix = common_suffix_avx2(r_tail, a_tail);
        let new_ref = String::from_utf8_unchecked(r_tail[..r_tail.len() - suffix].to_vec());
        let new_alt = String::from_utf8_unchecked(a_tail[..a_tail.len() - suffix].to_vec());

        (new_ref, new_alt, prefix, suffix)
    }
}

#[cfg(target_arch = "x86_64")]
#[inline]
unsafe fn common_prefix_avx2(a: &[u8], b: &[u8]) -> usize {
    use core::arch::x86_64::*;

    let len = a.len().min(b.len());
    let mut i = 0;

    // Process 32 bytes at a time
    while i + 32 <= len {
        let va = _mm256_loadu_si256(a.as_ptr().add(i) as *const __m256i);
        let vb = _mm256_loadu_si256(b.as_ptr().add(i) as *const __m256i);
        let cmp = _mm256_cmpeq_epi8(va, vb);
        let mask = _mm256_movemask_epi8(cmp);

        if mask != -1 {
            return i + (!mask as u32).trailing_zeros() as usize;
        }
        i += 32;
    }

    // Scalar tail
    while i < len && a[i] == b[i] {
        i += 1;
    }
    i
}

#[cfg(target_arch = "x86_64")]
#[inline]
unsafe fn common_suffix_avx2(a: &[u8], b: &[u8]) -> usize {
    use core::arch::x86_64::*;

    let la = a.len();
    let lb = b.len();
    let len = la.min(lb);
    let mut i = 0;

    // Process 32 bytes at a time from the end
    while i + 32 <= len {
        let pa = la - i - 32;
        let pb = lb - i - 32;

        let va = _mm256_loadu_si256(a.as_ptr().add(pa) as *const __m256i);
        let vb = _mm256_loadu_si256(b.as_ptr().add(pb) as *const __m256i);
        let cmp = _mm256_cmpeq_epi8(va, vb);
        let mask = _mm256_movemask_epi8(cmp);

        if mask != -1 {
            let lz = (!mask).leading_zeros() as usize;
            return i + (32 - lz);
        }
        i += 32;
    }

    // Scalar tail
    while i < len {
        if a[la - i - 1] != b[lb - i - 1] {
            return i;
        }
        i += 1;
    }
    i
}

/// NEON optimized implementation for ARM
#[cfg(target_arch = "aarch64")]
#[inline]
pub fn normalize_neon(ref_allele: &str, alt_allele: &str) -> (String, String, usize, usize) {
    use core::arch::aarch64::*;

    unsafe {
        let rb = ref_allele.as_bytes();
        let ab = alt_allele.as_bytes();

        let prefix = common_prefix_neon(rb, ab);
        let r_tail = &rb[prefix..];
        let a_tail = &ab[prefix..];

        let suffix = common_suffix_neon(r_tail, a_tail);
        let new_ref = String::from_utf8_unchecked(r_tail[..r_tail.len() - suffix].to_vec());
        let new_alt = String::from_utf8_unchecked(a_tail[..a_tail.len() - suffix].to_vec());

        (new_ref, new_alt, prefix, suffix)
    }
}

#[cfg(target_arch = "aarch64")]
#[inline]
unsafe fn common_prefix_neon(a: &[u8], b: &[u8]) -> usize {
    use core::arch::aarch64::*;

    let len = a.len().min(b.len());
    let mut i = 0;

    // Process 16 bytes at a time
    while i + 16 <= len {
        let av = vld1q_u8(a.as_ptr().add(i));
        let bv = vld1q_u8(b.as_ptr().add(i));
        let cmp = vceqq_u8(av, bv);

        // Check if all bytes match
        let min_val = vminvq_u8(cmp);
        if min_val != 0xFF {
            // Find first mismatch
            for j in 0..16 {
                if *a.get_unchecked(i + j) != *b.get_unchecked(i + j) {
                    return i + j;
                }
            }
        }
        i += 16;
    }

    // Scalar tail
    while i < len && a[i] == b[i] {
        i += 1;
    }
    i
}

#[cfg(target_arch = "aarch64")]
#[inline]
unsafe fn common_suffix_neon(a: &[u8], b: &[u8]) -> usize {
    use core::arch::aarch64::*;

    let la = a.len();
    let lb = b.len();
    let len = la.min(lb);
    let mut i = 0;

    // Process 16 bytes at a time from the end
    while i + 16 <= len {
        let pa = la - i - 16;
        let pb = lb - i - 16;

        let av = vld1q_u8(a.as_ptr().add(pa));
        let bv = vld1q_u8(b.as_ptr().add(pb));
        let cmp = vceqq_u8(av, bv);

        let min_val = vminvq_u8(cmp);
        if min_val != 0xFF {
            // Find first mismatch from the end
            for j in (0..16).rev() {
                if *a.get_unchecked(pa + j) != *b.get_unchecked(pb + j) {
                    return i + (15 - j);
                }
            }
        }
        i += 16;
    }

    // Scalar tail
    while i < len {
        if a[la - i - 1] != b[lb - i - 1] {
            return i;
        }
        i += 1;
    }
    i
}

/// Legacy API for backward compatibility
#[inline]
pub fn normalize(ref_allele: &str, alt_allele: &str) -> (String, String, usize, usize) {
    // Use runtime detection for best performance
    lazy_static::lazy_static! {
        static ref CONTEXT: NormContext = NormContext::detect();
    }
    CONTEXT.normalize(ref_allele, alt_allele)
}
