// PtrHash25 MPH lookup kernel (OpenCL parity of `ani_kernel.cu`).
//
// Same algorithm, same constants. Mirrors `gpu_sim::lookup_u64` and
// `ani_kernel.cu::ptrhash25_lookup` byte-for-byte. See gpu_sim.rs for the
// canonical reference + correctness test against `Index::lookup_u64`.

#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable

// ---------- PtrHash25 lookup primitives ----------------------------------

inline ulong mix64(ulong x) {
    x ^= x >> 32;
    x *= 0xd6e8feb86659fd93UL;
    x ^= x >> 32;
    x *= 0xd6e8feb86659fd93UL;
    x ^= x >> 32;
    return x;
}

inline ulong canonical_u64(ulong key, ulong prehash_seed) {
    return mix64(key ^ prehash_seed);
}

inline ulong rotl64(ulong v, uint n) {
    return (v << n) | (v >> (64u - n));
}

// OpenCL has `mul_hi(ulong, ulong)` for the upper 64 bits of a 128-bit
// multiply — this is what the Rust `fast_reduce` does via `as u128`.
inline ulong fast_reduce(ulong hash, ulong n) {
    return mul_hi(hash, n);
}

// 2-zone bucket-for. Constants:
//   ALPHA_BUCKETS = 0.30
//   BETA_KEYS = 0.60 → beta_threshold = 0.60 * 65536 = 39321
inline ulong ptrhash_bucket_for(ulong h1, uint num_buckets) {
    const uint beta_threshold = 39321u;
    ulong large_buckets = (ulong)((double)num_buckets * 0.30);
    if (large_buckets == 0UL) large_buckets = 1UL;
    ulong small_buckets = (ulong)num_buckets - large_buckets;
    if (small_buckets == 0UL) small_buckets = 1UL;

    uint zone_decider = (uint)(h1 >> 48);
    ulong h_low = h1 & 0x0000FFFFFFFFFFFFUL;
    ulong shifted = h_low << 16;
    if (zone_decider < beta_threshold) {
        return fast_reduce(shifted, large_buckets);
    } else {
        return large_buckets + fast_reduce(shifted, small_buckets);
    }
}

inline ulong ptrhash_slot_for(ulong h2, uchar pilot, ulong n) {
    ulong pilot_mix = (ulong)pilot * 0xA24B1F6FDA392B31UL;
    ulong mixed = rotl64(h2 ^ pilot_mix, 31u) * 0xD6E8FEB86659FD93UL;
    return fast_reduce(mixed, n);
}

inline void ptrhash_hash_key(ulong canonical, ulong mph_salt, uint prerotate,
                             ulong *h1_out, ulong *h2_out) {
    ulong rotated = rotl64(canonical, prerotate);
    ulong base = mix64(rotated ^ mph_salt);
    *h1_out = base;
    *h2_out = rotl64(base, 23u) ^ 0xA24B1F6FDA392B31UL;
}

// ---------- Block-Bloom prefilter ----------------------------------------

inline bool bloom_contains(__global const ulong *bloom_words,
                           uint bloom_blocks,
                           ulong canonical) {
    if (bloom_words == 0 || bloom_blocks == 0u) {
        return true; // no filter → pass
    }
    const uint SALT[8] = {
        0x47b6137bu, 0x44768924u, 0x18205237u, 0x23848965u,
        0x8e6e2354u, 0x0f7cc9b6u, 0xe43d5fa5u, 0xa4d52dc1u,
    };
    ulong h_hi = canonical >> 32;
    uint block = (uint)((h_hi * (ulong)bloom_blocks) >> 32);
    uint seed = (uint)(canonical & 0xFFFFFFFFu);
    uint base = block * 8u;
    for (int i = 0; i < 8; ++i) {
        uint bit = (seed * SALT[i] >> 27) & 0x3Fu;
        ulong mask = 1UL << bit;
        if ((bloom_words[base + i] & mask) != mask) {
            return false;
        }
    }
    return true;
}

inline ushort fingerprint16(ulong canonical) {
    return (ushort)(canonical & 0xFFFFu);
}

inline uint ptrhash25_lookup(
    ulong  key,
    ulong  prehash_seed,
    ulong  mph_salt,
    uint   num_buckets,
    ulong  num_slots,
    uint   prerotate,
    __global const uchar  *pilots,
    __global const ulong  *bloom_words,
    uint   bloom_blocks,
    __global const ushort *fingerprints  // null when lean_mph
) {
    ulong canonical = canonical_u64(key, prehash_seed);
    if (!bloom_contains(bloom_words, bloom_blocks, canonical)) {
        return 0xFFFFFFFFu;
    }
    ulong h1, h2;
    ptrhash_hash_key(canonical, mph_salt, prerotate, &h1, &h2);
    ulong bucket = ptrhash_bucket_for(h1, num_buckets);
    uchar pilot = pilots[(uint)bucket];
    ulong slot = ptrhash_slot_for(h2, pilot, num_slots);
    if (fingerprints != 0) {
        if (fingerprints[(uint)slot] != fingerprint16(canonical)) {
            return 0xFFFFFFFFu;
        }
    }
    return (uint)slot;
}

// Batched: one work-item per key.
__kernel void ani_ptrhash25_lookup_kernel(
    __global const ulong  *keys,
    ulong  prehash_seed,
    ulong  mph_salt,
    uint   num_buckets,
    ulong  num_slots,
    uint   prerotate,
    __global const uchar  *pilots,
    __global const ulong  *bloom_words,
    uint   bloom_blocks,
    __global const ushort *fingerprints,
    __global uint         *out_idx,
    int    nkeys
) {
    int gid = get_global_id(0);
    if (gid >= nkeys) return;
    out_idx[gid] = ptrhash25_lookup(
        keys[gid], prehash_seed, mph_salt, num_buckets, num_slots, prerotate,
        pilots, bloom_words, bloom_blocks, fingerprints);
}
