//! Checksums for BGZF blocks (`crc32fast` picks the SIMD path at runtime).

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
