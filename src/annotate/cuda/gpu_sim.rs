//! Pure-Rust simulator of the GPU MPH lookup kernel.
//!
//! Replicates [`kira_kv_engine::Index::lookup_u64`] starting from a
//! `GpuExport` POD snapshot. The CUDA and OpenCL kernels mirror this same
//! arithmetic — any divergence is a kernel bug.
//!
//! kira_kv_engine 0.6.3 splits an index into parts of ~32 768 keys, so every
//! real index is multi-part: the part is picked first, and its own salt and
//! bucket/slot ranges drive the rest. A single-part export takes the same
//! path with the identity offsets, which is the pre-0.6.3 formula.

#![cfg(feature = "gpu")]

use kira_kv_engine::{BloomExport, GpuExport, GpuPart};

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
/// Must match `block_bloom::block_and_mask` exactly. `bit_shift` is the
/// export's own lane shift: 26 (6-bit index) for filters built by 0.6.3,
/// 27 (5-bit) for ones deserialized from a 0.6 file.
#[inline]
fn bloom_block_and_mask(hash: u64, blocks: usize, bit_shift: u32) -> (usize, [u64; BLOOM_BLOCK_WORDS]) {
    let block = (((hash >> 32) as u128 * blocks as u128) >> 32) as usize;
    let seed = hash as u32;
    const SALT: [u32; BLOOM_BLOCK_WORDS] = [
        0x47b6_137b, 0x4476_8924, 0x1820_5237, 0x2384_8965, 0x8e6e_2354, 0x0f7c_c9b6, 0xe43d_5fa5,
        0xa4d5_2dc1,
    ];
    let mut mask = [0u64; BLOOM_BLOCK_WORDS];
    let mut i = 0;
    while i < BLOOM_BLOCK_WORDS {
        let bit = (seed.wrapping_mul(SALT[i]) >> bit_shift) & 0x3F;
        mask[i] = 1u64 << bit;
        i += 1;
    }
    (block, mask)
}

/// Returns `true` if all 8 lane bits for `canonical` are set in the
/// corresponding block of the exported bloom filter.
pub fn bloom_contains(bloom: &BloomExport, canonical: u64) -> bool {
    let (block, mask) = bloom_block_and_mask(canonical, bloom.blocks, bloom.bit_shift);
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
/// `ptrhash25::PART_MUL`: Fibonacci multiplier of the part selector.
const PART_MUL: u64 = 0x9E37_79B9_7F4A_7C15;
/// `ptrhash25::REMIX_MUL`: per-part rehash of the base hash.
const REMIX_MUL: u64 = 0xBF58_476D_1CE4_E5B9;
/// `ptrhash25::H2_XOR`.
const H2_XOR: u64 = 0xA24B_1F6F_DA39_2B31;

/// Multiply-high reduction for hash → [0, n). 64-bit mulhi.
#[inline]
const fn fast_reduce(hash: u64, n: usize) -> usize {
    ((hash as u128 * n as u128) >> 64) as usize
}

/// `ptrhash25::large_buckets_of`, for exports that carry no part table.
#[inline]
fn large_buckets_of(num_buckets: usize) -> usize {
    ((num_buckets as f64) * ALPHA_BUCKETS) as usize
}

/// Bucket within one part: 2-level skewed mapping from h1 → bucket id.
/// Mirrors `ptrhash25::bucket_in_part` byte-for-byte.
#[inline]
fn ptrhash_bucket_in_part(h1: u64, large_buckets: usize, small_buckets: usize) -> usize {
    let zone_decider = (h1 >> 48) as u32;
    let beta_threshold = (BETA_KEYS * 65536.0) as u32;
    let h_low = (h1 & 0x0000_FFFF_FFFF_FFFF) << 16;
    if zone_decider < beta_threshold {
        fast_reduce(h_low, large_buckets.max(1))
    } else {
        large_buckets + fast_reduce(h_low, small_buckets.max(1))
    }
}

/// Slot-for: combine h2 with the bucket's pilot byte and reduce into [0, n).
/// Mirrors `ptrhash25::slot_for`.
#[inline]
fn ptrhash_slot_for(h2: u64, pilot: u8, n: usize) -> usize {
    let pilot_mix = (pilot as u64).wrapping_mul(0xA24B_1F6F_DA39_2B31);
    let mixed = (h2 ^ pilot_mix)
        .rotate_left(31)
        .wrapping_mul(0xD6E8_FEB8_6659_FD93);
    fast_reduce(mixed, n)
}

/// `ptrhash25::part_of`: which part owns an already-rotated key.
#[inline]
fn ptrhash_part_of(rotated: u64, part_salt: u64, parts: usize) -> usize {
    fast_reduce((rotated ^ part_salt).wrapping_mul(PART_MUL), parts)
}

/// Part geometry, falling back to the whole index when the export carries no
/// part table (single-partition builds from kira_kv_engine 0.6).
#[inline]
fn part_geometry(export: &GpuExport, part: usize) -> GpuPart {
    match export.parts.get(part) {
        Some(p) => *p,
        None => GpuPart {
            slot_off: 0,
            salt: export.mph_salt,
            bucket_off: 0,
            num_slots: export.num_slots as u32,
            num_buckets: export.num_buckets,
            large_buckets: large_buckets_of(export.num_buckets as usize) as u32,
        },
    }
}

/// PtrHash25 keyed hash: pre-rotate the key, mix with `mph_salt`, pick the
/// part, then derive the part-local halves h1/h2. Mirrors
/// `ptrhash25::PtrHash25Mphf::index_u64` for the non-AES path (the AES branch
/// isn't exported via `GpuExport`).
#[inline]
fn ptrhash_hash_key(export: &GpuExport, canonical: u64) -> (usize, u64, u64) {
    let rotated = canonical.rotate_left(export.prerotate as u32);
    let base = mix64(rotated ^ export.mph_salt);
    let (part, h1) = if export.parts.len() > 1 {
        let p = ptrhash_part_of(rotated, export.part_salt, export.parts.len());
        (p, (base ^ export.parts[p].salt).wrapping_mul(REMIX_MUL))
    } else {
        (0, base)
    };
    (part, h1, h1.rotate_left(23) ^ H2_XOR)
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

    let (part, h1, h2) = ptrhash_hash_key(export, canonical);
    let geom = part_geometry(export, part);
    let large = geom.large_buckets as usize;
    let small = (geom.num_buckets as usize).saturating_sub(large);
    let bucket = geom.bucket_off as usize + ptrhash_bucket_in_part(h1, large, small);
    let pilot = export.pilots[bucket];
    let slot = geom.slot_off + ptrhash_slot_for(h2, pilot, geom.num_slots as usize) as u64;

    if let Some(fps) = &export.fingerprints {
        let expected = fingerprint16(canonical);
        if fps.get(slot as usize).copied() != Some(expected) {
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
