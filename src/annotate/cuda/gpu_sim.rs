//! Pure-Rust simulator of the GPU MPH lookup kernel.
//!
//! Replicates [`kira_kv_engine::Index::lookup_u64`] starting from a
//! `GpuExport` POD snapshot. The CUDA and OpenCL kernels mirror this same
//! arithmetic — any divergence is a kernel bug.

#![cfg(feature = "gpu")]

use kira_kv_engine::{BloomExport, GpuExport};

// ---------- canonical hash (engine `prehash_seed` step) --------------------

/// `simd_hash::mix64` from kira_kv_engine. Three xorshift+mul rounds.
#[inline]
const fn mix64(mut x: u64) -> u64 {
    x ^= x >> 32;
    x = x.wrapping_mul(0xd6e8feb86659fd93);
    x ^= x >> 32;
    x = x.wrapping_mul(0xd6e8feb86659fd93);
    x ^= x >> 32;
    x
}

/// `canonical_hash_bytes` specialised to the 8-byte u64 path used by
/// `Index::lookup_u64`. Equivalent to `simd_hash::hash_u64_one(key, seed)`.
#[inline]
pub const fn canonical_u64(key: u64, prehash_seed: u64) -> u64 {
    mix64(key ^ prehash_seed)
}

// ---------- block-Bloom prefilter -----------------------------------------

const BLOOM_BLOCK_WORDS: usize = 8;

/// Per-block bit pattern derived from a hash (Impala/Pibiri split-block layout).
/// Must match `block_bloom::block_and_mask` exactly.
#[inline]
fn bloom_block_and_mask(hash: u64, blocks: usize) -> (usize, [u64; BLOOM_BLOCK_WORDS]) {
    let block = (((hash >> 32) as u128 * blocks as u128) >> 32) as usize;
    let seed = hash & 0xFFFF_FFFF;
    const SALT: [u32; BLOOM_BLOCK_WORDS] = [
        0x47b6_137b, 0x4476_8924, 0x1820_5237, 0x2384_8965, 0x8e6e_2354, 0x0f7c_c9b6, 0xe43d_5fa5,
        0xa4d5_2dc1,
    ];
    let mut mask = [0u64; BLOOM_BLOCK_WORDS];
    let mut i = 0;
    while i < BLOOM_BLOCK_WORDS {
        let bit = ((seed as u32).wrapping_mul(SALT[i]) >> 27) & 0x3F;
        mask[i] = 1u64 << bit;
        i += 1;
    }
    (block, mask)
}

/// Returns `true` if all 8 lane bits for `canonical` are set in the
/// corresponding block of the exported bloom filter.
pub fn bloom_contains(bloom: &BloomExport, canonical: u64) -> bool {
    let (block, mask) = bloom_block_and_mask(canonical, bloom.blocks);
    let base = block * BLOOM_BLOCK_WORDS;
    for w in 0..BLOOM_BLOCK_WORDS {
        if bloom.words[base + w] & mask[w] != mask[w] {
            return false;
        }
    }
    true
}

// ---------- PtrHash25 MPH lookup ------------------------------------------

/// Same constants the engine uses for the 2-zone bucket skew.
const ALPHA_BUCKETS: f64 = 0.30;
const BETA_KEYS: f64 = 0.60;

/// Multiply-high reduction for hash → [0, n). 64-bit mulhi.
#[inline]
const fn fast_reduce(hash: u64, n: usize) -> usize {
    ((hash as u128 * n as u128) >> 64) as usize
}

/// Bucket-for: 2-level skewed mapping from h1 → bucket id. Mirrors
/// `ptrhash25::bucket_for` byte-for-byte.
#[inline]
fn ptrhash_bucket_for(h1: u64, num_buckets: u32) -> usize {
    let num_buckets = num_buckets as usize;
    let zone_decider = (h1 >> 48) as u32;
    let beta_threshold = (BETA_KEYS * 65536.0) as u32;
    let large_buckets = ((num_buckets as f64) * ALPHA_BUCKETS) as usize;
    let small_buckets = num_buckets - large_buckets;
    let h_low = h1 & 0x0000_FFFF_FFFF_FFFF;
    if zone_decider < beta_threshold {
        fast_reduce(h_low << 16, large_buckets.max(1))
    } else {
        large_buckets + fast_reduce(h_low << 16, small_buckets.max(1))
    }
}

/// Slot-for: combine h2 with the bucket's pilot byte and reduce into [0, n).
/// Mirrors `ptrhash25::slot_for`.
#[inline]
fn ptrhash_slot_for(h2: u64, pilot: u8, n: u64) -> usize {
    let pilot_mix = (pilot as u64).wrapping_mul(0xA24B_1F6F_DA39_2B31);
    let mixed = (h2 ^ pilot_mix)
        .rotate_left(31)
        .wrapping_mul(0xD6E8_FEB8_6659_FD93);
    fast_reduce(mixed, n as usize)
}

/// PtrHash25 keyed hash: pre-rotate the key, mix with `mph_salt`, derive two
/// independent halves h1/h2. Mirrors `ptrhash25::hash_key` for the non-AES
/// path (the AES branch isn't exported via `GpuExport`).
#[inline]
fn ptrhash_hash_key(canonical: u64, mph_salt: u64, prerotate: u8) -> (u64, u64) {
    let rotated = canonical.rotate_left(prerotate as u32);
    let base = mix64(rotated ^ mph_salt);
    let h1 = base;
    let h2 = base.rotate_left(23) ^ 0xA24B_1F6F_DA39_2B31;
    (h1, h2)
}

/// u16 fingerprint check (low 16 bits of canonical hash) — matches
/// `index::fingerprint16_mph`.
#[inline]
const fn fingerprint16(canonical: u64) -> u16 {
    (canonical & 0xFFFF) as u16
}

// ---------- top-level kernel -----------------------------------------------

/// Result of a single key lookup, mirroring `Index::lookup_u64` semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GpuLookup {
    /// Key passed all filters; this is the MPH slot index.
    Found(u32),
    /// Bloom filter rejected. Foreign-key fast path.
    BloomMiss,
    /// Fingerprint check rejected. Foreign key that collided with a real one.
    FingerprintMiss,
}

impl GpuLookup {
    /// Strip the rejection-reason for the common `Option<usize>` API parity.
    #[inline]
    pub const fn as_option(self) -> Option<u32> {
        match self {
            Self::Found(idx) => Some(idx),
            _ => None,
        }
    }
}

/// Reference implementation of the GPU MPH lookup that the CUDA/OpenCL
/// kernels must match.
pub fn lookup_u64(export: &GpuExport, key: u64) -> GpuLookup {
    let canonical = canonical_u64(key, export.prehash_seed);

    if let Some(bloom) = &export.bloom
        && !bloom_contains(bloom, canonical)
    {
        return GpuLookup::BloomMiss;
    }

    let (h1, h2) = ptrhash_hash_key(canonical, export.mph_salt, export.prerotate);
    let bucket = ptrhash_bucket_for(h1, export.num_buckets);
    let pilot = export.pilots[bucket];
    let slot = ptrhash_slot_for(h2, pilot, export.num_slots);

    if let Some(fps) = &export.fingerprints {
        let expected = fingerprint16(canonical);
        if fps.get(slot).copied() != Some(expected) {
            return GpuLookup::FingerprintMiss;
        }
    }

    GpuLookup::Found(slot as u32)
}

/// Batched form: iterates the simulator on a key slice.
pub fn lookup_batch(export: &GpuExport, keys: &[u64], out: &mut [Option<u32>]) {
    assert!(out.len() >= keys.len(), "out slice too small");
    for (i, &k) in keys.iter().enumerate() {
        out[i] = lookup_u64(export, k).as_option();
    }
}

#[cfg(test)]
#[path = "../../../tests/unit/annotate_cuda_gpu_sim.rs"]
mod tests;
