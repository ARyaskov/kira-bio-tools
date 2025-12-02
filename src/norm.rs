use core::arch::global_asm;

pub fn normalize(ref_allele: &str, alt_allele: &str) -> (String, String, usize, usize) {
    let prefix = common_prefix_simd(ref_allele.as_bytes(), alt_allele.as_bytes());
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

#[cfg(target_arch = "x86_64")]
#[inline]
fn common_prefix_simd(a: &[u8], b: &[u8]) -> usize {
    use core::arch::x86_64::*;
    let len = a.len().min(b.len());
    let mut i = 0;

    unsafe {
        while i + 32 <= len {
            let av = _mm256_loadu_si256(a.as_ptr().add(i) as *const _);
            let bv = _mm256_loadu_si256(b.as_ptr().add(i) as *const _);
            let cmp = _mm256_cmpeq_epi8(av, bv);
            let mask = _mm256_movemask_epi8(cmp);
            if mask != -1 {
                let tz = (!mask).trailing_zeros() as usize;
                return i + tz;
            }
            i += 32;
        }
    }

    while i < len && a[i] == b[i] {
        i += 1;
    }
    i
}

#[cfg(all(target_arch = "aarch64"))]
#[inline]
fn common_prefix_simd(a: &[u8], b: &[u8]) -> usize {
    use core::arch::aarch64::*;
    let len = a.len().min(b.len());
    let mut i = 0;

    unsafe {
        while i + 16 <= len {
            let av = vld1q_u8(a.as_ptr().add(i));
            let bv = vld1q_u8(b.as_ptr().add(i));
            let cmp = vceqq_u8(av, bv);
            let mask = vmaxvq_u8(vshrq_n_u8(cmp, 7));
            if mask != 0b1111_1111 {
                for j in 0..16 {
                    if *a.get_unchecked(i + j) != *b.get_unchecked(i + j) {
                        return i + j;
                    }
                }
            }
            i += 16;
        }
    }

    while i < len && a[i] == b[i] {
        i += 1;
    }
    i
}

#[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
#[inline]
fn common_prefix_simd(a: &[u8], b: &[u8]) -> usize {
    let mut i = 0;
    let len = a.len().min(b.len());
    while i < len && a[i] == b[i] {
        i += 1;
    }
    i
}
