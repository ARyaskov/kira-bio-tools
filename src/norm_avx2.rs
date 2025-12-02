use std::arch::x86_64::*;

unsafe fn common_prefix_avx2(a: &[u8], b: &[u8]) -> usize {
    let len = a.len().min(b.len());
    let mut i = 0;

    while i + 32 <= len {
        let va = _mm256_loadu_si256(a.as_ptr().add(i) as *const _);
        let vb = _mm256_loadu_si256(b.as_ptr().add(i) as *const _);
        let cmp = _mm256_cmpeq_epi8(va, vb);
        let mask = _mm256_movemask_epi8(cmp);
        if mask != -1 {
            return i + (!mask as u32).trailing_zeros() as usize;
        }
        i += 32;
    }

    while i < len {
        if a[i] != b[i] {
            return i;
        }
        i += 1;
    }
    i
}

unsafe fn common_suffix_avx2(a: &[u8], b: &[u8]) -> usize {
    let la = a.len();
    let lb = b.len();
    let len = la.min(lb);
    let mut i = 0;

    while i + 32 <= len {
        let pa = la - i - 32;
        let pb = lb - i - 32;

        let va = _mm256_loadu_si256(a.as_ptr().add(pa) as *const _);
        let vb = _mm256_loadu_si256(b.as_ptr().add(pb) as *const _);
        let cmp = _mm256_cmpeq_epi8(va, vb);
        let mask = _mm256_movemask_epi8(cmp);

        if mask != -1 {
            let lz = (!mask).leading_zeros() as usize;
            return i + (32 - lz);
        }
        i += 32;
    }

    while i < len {
        if a[la - i - 1] != b[lb - i - 1] {
            return i;
        }
        i += 1;
    }
    i
}

pub fn normalize_avx2(r: &str, a: &str) -> (String, String, usize, usize) {
    unsafe {
        let rb = r.as_bytes();
        let ab = a.as_bytes();

        let cp = common_prefix_avx2(rb, ab);
        let rb2 = &rb[cp..];
        let ab2 = &ab[cp..];

        let cs = common_suffix_avx2(rb2, ab2);
        let r2 = &rb2[..rb2.len().saturating_sub(cs)];
        let a2 = &ab2[..ab2.len().saturating_sub(cs)];

        (
            std::str::from_utf8_unchecked(r2).to_string(),
            std::str::from_utf8_unchecked(a2).to_string(),
            cp,
            cs,
        )
    }
}
