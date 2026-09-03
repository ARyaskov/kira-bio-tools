// PtrHash25 MPH lookup kernels for the ANI annotate path.
//
// **Algorithm.** Implements `kira_kv_engine::Index::lookup_u64` directly from
// the `GpuExport` POD snapshot. The Rust simulator at `gpu_sim.rs` runs the
// same arithmetic byte-for-byte; if their outputs ever diverge on the same
// input the kernel is wrong, and the unit tests in `gpu_sim::tests` will
// fail first.
//
// **Why this replaced the previous kernel.** The legacy `ani_lookup_kernel_v2`
// used the old CHD-style 3-pilot vector hash `(g[v0] + g[v1] + g[v2]) % n`
// with an XOR-commutative variant key — a layout that's no longer produced
// by kira_kv_engine 0.6.0. The new PtrHash25 backend has different on-disk
// constants and a different lookup formula; consuming the old layout on a
// new build would silently return random slot indices.
//
// **Partitioning (kira_kv_engine 0.6.3).** The engine now splits the key set
// into parts of ~32768 keys, so every index of real size is multi-part. The
// part selector runs first and each part carries its own hash salt plus its
// own bucket and slot ranges; a one-part export reduces to the pre-0.6.3
// formula. Getting this wrong returns plausible but wrong slots, which the
// `entry_keys` re-check downgrades to a silent miss.
//
// **Constants must match `ptrhash25.rs` exactly:**
//   - mix64: triple xorshift+mul fmix
//   - ALPHA_BUCKETS = 0.30, BETA_KEYS = 0.60 → beta_threshold = 39321
//   - fast_reduce: hi-64 of the 128-bit product
//   - part selector: mulhi((rotated ^ part_salt) * 0x9E3779B97F4A7C15, nparts)
//   - per-part remix: h1 = (base ^ part.salt) * 0xBF58476D1CE4E5B9
//   - h2 derivation: h1.rotate_left(23) ^ 0xA24B1F6FDA392B31
//   - slot mix: pilot * 0xA24B1F6FDA392B31, then rotl(31), then * 0xD6E8FEB86659FD93
//   - fingerprint: low 16 bits of canonical hash
//   - bloom: 8-word split-block, ((h>>32) * blocks) >> 32 selects block, lane
//     bit is `(seed * SALT[w]) >> bit_shift`, with bit_shift 26 for filters
//     built by 0.6.3 and 27 for ones loaded from a 0.6 file
//
// **Memory layout pushed to the device per index (allocated once at GpuAni::load):**
//   - prehash_seed (u64), mph_salt (u64), part_salt (u64), num_buckets (u32),
//     num_slots (u64), num_parts (u32), prerotate (u8), bloom_bit_shift (u32)
//     — uploaded as scalars via cudaMemcpyToSymbol or kernel args.
//   - parts: AniMphPart array of length `num_parts` (32 B each, so a few MB
//     even for 100M keys).
//   - pilots: u8 array of length `num_buckets`. ~30-50 MB for a 100M-key index.
//   - bloom_words: optional u64 array of length `bloom_blocks * 8`, a multiple
//     of 256 blocks, ~16-20 MB for 100M keys at 11 bits/key.
//   - fingerprints: optional u16 array of length `num_slots`. ~200 MB for 100M.

#include <cuda_runtime.h>
#include <stdint.h>

typedef struct {
    unsigned char  chr_id;
    unsigned int   pos;
    unsigned int   ref_ofs;
    unsigned int   alt_ofs;
    unsigned int   id_ofs;
    unsigned int   qual_ofs;
    unsigned int   filter_ofs;
    unsigned int   info_ofs;
    unsigned int   info_len;
    unsigned int   format_ofs;
    unsigned int   samples_ofs;
} AniEntry;

typedef struct {
    unsigned int chr_id;
    unsigned int min_pos;
    unsigned int max_pos;
    unsigned int block_start;
    unsigned int block_count;
} AniPosContig;

// Mirrors `kira_kv_engine::GpuPart` and the Rust `AniMphPart` upload struct.
typedef struct {
    unsigned long long slot_off;
    unsigned long long salt;
    unsigned int       bucket_off;
    unsigned int       num_slots;
    unsigned int       num_buckets;
    unsigned int       large_buckets;
} AniMphPart;

typedef struct {
    unsigned int base_pos;
    unsigned int _pad;
    unsigned long long masks[8];
    unsigned int offsets_start;
    unsigned int _pad2;
} AniPosBlock;

// ---------- PtrHash25 lookup primitives ----------------------------------

__device__ __forceinline__ unsigned long long mix64(unsigned long long x) {
    x ^= x >> 32;
    x *= 0xd6e8feb86659fd93ULL;
    x ^= x >> 32;
    x *= 0xd6e8feb86659fd93ULL;
    x ^= x >> 32;
    return x;
}

__device__ __forceinline__ unsigned long long canonical_u64(unsigned long long key,
                                                              unsigned long long prehash_seed) {
    return mix64(key ^ prehash_seed);
}

__device__ __forceinline__ unsigned long long rotl64(unsigned long long v, unsigned int n) {
    return (v << n) | (v >> (64u - n));
}

__device__ __forceinline__ unsigned long long fast_reduce(unsigned long long hash,
                                                            unsigned long long n) {
    // mulhi(hash, n) via 128-bit multiply. nvcc emits MUL + MULHI for this.
    return __umul64hi(hash, n);
}

// ---------- 2-zone bucket-for (matches ptrhash25::bucket_in_part) --------

// (BETA_KEYS * 65536) = (0.60 * 65536) = 39321.6 → truncated to 39321.
// `large_buckets` comes from the export (floor(num_buckets * 0.30)) so the
// device never re-derives it in a different float rounding mode.
__device__ __forceinline__ unsigned long long ptrhash_bucket_in_part(unsigned long long h1,
                                                                       unsigned int large,
                                                                       unsigned int small) {
    const unsigned int beta_threshold = 39321u; // BETA_KEYS * 65536
    unsigned long long large_buckets = (unsigned long long)large;
    unsigned long long small_buckets = (unsigned long long)small;

    unsigned int zone_decider = (unsigned int)(h1 >> 48);
    unsigned long long h_low = h1 & 0x0000FFFFFFFFFFFFULL;
    unsigned long long shifted = h_low << 16;
    if (zone_decider < beta_threshold) {
        return fast_reduce(shifted, large_buckets == 0ULL ? 1ULL : large_buckets);
    } else {
        return large_buckets + fast_reduce(shifted, small_buckets == 0ULL ? 1ULL : small_buckets);
    }
}

__device__ __forceinline__ unsigned long long ptrhash_slot_for(unsigned long long h2,
                                                                 unsigned char pilot,
                                                                 unsigned long long n) {
    unsigned long long pilot_mix = (unsigned long long)pilot * 0xA24B1F6FDA392B31ULL;
    unsigned long long mixed = rotl64(h2 ^ pilot_mix, 31u) * 0xD6E8FEB86659FD93ULL;
    return fast_reduce(mixed, n);
}

// Picks the owning part and derives its local h1/h2. With a single part this
// is the pre-0.6.3 formula (`h1 = base`).
__device__ __forceinline__ void ptrhash_hash_key(unsigned long long canonical,
                                                   unsigned long long mph_salt,
                                                   unsigned long long part_salt,
                                                   unsigned int prerotate,
                                                   const AniMphPart *parts,
                                                   unsigned int num_parts,
                                                   unsigned int *part_out,
                                                   unsigned long long *h1_out,
                                                   unsigned long long *h2_out) {
    unsigned long long rotated = rotl64(canonical, prerotate);
    unsigned long long base = mix64(rotated ^ mph_salt);
    unsigned int part = 0u;
    unsigned long long h1 = base;
    if (parts != 0 && num_parts > 1u) {
        part = (unsigned int)fast_reduce((rotated ^ part_salt) * 0x9E3779B97F4A7C15ULL,
                                         (unsigned long long)num_parts);
        h1 = (base ^ parts[part].salt) * 0xBF58476D1CE4E5B9ULL;
    }
    *part_out = part;
    *h1_out = h1;
    *h2_out = rotl64(h1, 23u) ^ 0xA24B1F6FDA392B31ULL;
}

// ---------- Block-Bloom prefilter ----------------------------------------

__device__ __forceinline__ bool bloom_contains(const unsigned long long *bloom_words,
                                                unsigned int bloom_blocks,
                                                unsigned int bit_shift,
                                                unsigned long long canonical) {
    if (bloom_words == 0 || bloom_blocks == 0u) {
        return true; // no filter present → pass through
    }
    // Same SALT constants as block_bloom.rs.
    const unsigned int SALT[8] = {
        0x47b6137bu, 0x44768924u, 0x18205237u, 0x23848965u,
        0x8e6e2354u, 0x0f7cc9b6u, 0xe43d5fa5u, 0xa4d52dc1u,
    };
    unsigned long long h_hi = canonical >> 32;
    unsigned long long block_u64 = __umul64hi(h_hi << 32, (unsigned long long)bloom_blocks << 32);
    // The Rust reference is `((hash >> 32) as u128 * blocks as u128) >> 32`.
    // We translate as: ((h_hi) * blocks) >> 32. Use a fused multiply on 64-bit.
    block_u64 = (h_hi * (unsigned long long)bloom_blocks) >> 32;
    unsigned int block = (unsigned int)block_u64;
    unsigned int seed = (unsigned int)(canonical & 0xFFFFFFFFu);
    unsigned int base = block * 8u;
    #pragma unroll 8
    for (int i = 0; i < 8; ++i) {
        unsigned int bit = (seed * SALT[i] >> bit_shift) & 0x3Fu;
        unsigned long long mask = 1ULL << bit;
        if ((bloom_words[base + i] & mask) != mask) {
            return false;
        }
    }
    return true;
}

__device__ __forceinline__ unsigned short fingerprint16(unsigned long long canonical) {
    return (unsigned short)(canonical & 0xFFFFu);
}

// ---------- Top-level lookup ---------------------------------------------

// Returns 0xFFFFFFFF on miss (bloom or fingerprint reject), otherwise the
// MPH slot index in [0, num_slots).
__device__ __forceinline__ unsigned int ptrhash25_lookup(
    unsigned long long      key,
    unsigned long long      prehash_seed,
    unsigned long long      mph_salt,
    unsigned long long      part_salt,
    unsigned int            num_buckets,
    unsigned long long      num_slots,
    unsigned int            prerotate,
    const AniMphPart       *parts,
    unsigned int            num_parts,
    const unsigned char    *pilots,
    const unsigned long long *bloom_words,
    unsigned int            bloom_blocks,
    unsigned int            bloom_bit_shift,
    const unsigned short   *fingerprints   // null if lean_mph
) {
    unsigned long long canonical = canonical_u64(key, prehash_seed);
    if (!bloom_contains(bloom_words, bloom_blocks, bloom_bit_shift, canonical)) {
        return 0xFFFFFFFFu;
    }

    unsigned int part;
    unsigned long long h1, h2;
    ptrhash_hash_key(canonical, mph_salt, part_salt, prerotate, parts, num_parts,
                     &part, &h1, &h2);

    // Part geometry, or the whole index when the export carries no part table.
    unsigned long long slot_off = 0ULL;
    unsigned int bucket_off = 0u;
    unsigned int part_buckets = num_buckets;
    unsigned long long part_slots = num_slots;
    unsigned int large = (unsigned int)((double)num_buckets * 0.30);
    if (parts != 0 && num_parts > 0u) {
        slot_off = parts[part].slot_off;
        bucket_off = parts[part].bucket_off;
        part_buckets = parts[part].num_buckets;
        part_slots = (unsigned long long)parts[part].num_slots;
        large = parts[part].large_buckets;
    }
    unsigned int small = part_buckets > large ? part_buckets - large : 0u;

    unsigned long long bucket = (unsigned long long)bucket_off
                              + ptrhash_bucket_in_part(h1, large, small);
    unsigned char pilot = pilots[(unsigned int)bucket];
    unsigned long long slot = slot_off + ptrhash_slot_for(h2, pilot, part_slots);

    if (fingerprints != 0) {
        if (fingerprints[(unsigned int)slot] != fingerprint16(canonical)) {
            return 0xFFFFFFFFu;
        }
    }
    return (unsigned int)slot;
}

// Batched kernel: one thread per key. Writes `0xFFFFFFFF` on miss.
extern "C" __global__ void ani_ptrhash25_lookup_kernel(
    const unsigned long long *__restrict__ keys,
    unsigned long long       prehash_seed,
    unsigned long long       mph_salt,
    unsigned long long       part_salt,
    unsigned int             num_buckets,
    unsigned long long       num_slots,
    unsigned int             prerotate,
    const AniMphPart         *__restrict__ parts,        // null when not exported
    unsigned int             num_parts,
    const unsigned char      *__restrict__ pilots,
    const unsigned long long *__restrict__ bloom_words,  // null when filter absent
    unsigned int             bloom_blocks,
    unsigned int             bloom_bit_shift,
    const unsigned short     *__restrict__ fingerprints, // null when lean_mph
    unsigned int             *__restrict__ out_idx,
    int                      nkeys
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= nkeys) return;
    out_idx[tid] = ptrhash25_lookup(
        keys[tid], prehash_seed, mph_salt, part_salt, num_buckets, num_slots,
        prerotate, parts, num_parts, pilots, bloom_words, bloom_blocks,
        bloom_bit_shift, fingerprints);
}

// ---------- pos-index kernel (unchanged) ---------------------------------

__global__ void ani_pos_lookup_kernel(
    const unsigned int* __restrict__ chr_ids,
    const unsigned int* __restrict__ positions,
    const AniPosContig* __restrict__ contigs,
    unsigned int contig_count,
    const AniPosBlock* __restrict__ blocks,
    const unsigned int* __restrict__ pos_offsets,
    const unsigned short* __restrict__ pos_counts,
    unsigned int* __restrict__ out_offsets,
    unsigned short* __restrict__ out_counts,
    int nkeys
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= nkeys) return;

    unsigned int chr_id = chr_ids[tid];
    unsigned int pos = positions[tid];

    // Linear scan by chr_id (typically ≤ 25 contigs; on GPU branching is fine
    // because all threads in a warp scan the same array in lockstep).
    const AniPosContig *contig = 0;
    for (unsigned int i = 0; i < contig_count; ++i) {
        if (contigs[i].chr_id == chr_id) {
            contig = &contigs[i];
            break;
        }
    }
    if (contig == 0 || pos < contig->min_pos || pos > contig->max_pos) {
        out_offsets[tid] = 0u;
        out_counts[tid] = 0u;
        return;
    }

    unsigned int base = (pos / 512u) * 512u;
    unsigned int lo = contig->block_start;
    unsigned int hi = lo + contig->block_count;
    // Binary search on blocks by base_pos. blocks are in-place; warp diverges
    // only on early-exit. Acceptable for the modest contig sizes.
    while (lo < hi) {
        unsigned int mid = (lo + hi) >> 1;
        unsigned int bp = blocks[mid].base_pos;
        if (bp < base) lo = mid + 1;
        else if (bp > base) hi = mid;
        else { lo = mid; hi = mid; }
    }
    // After the loop lo points to the first block with base_pos >= base.
    // Need an exact match.
    unsigned int block_idx = lo;
    if (block_idx >= contig->block_start + contig->block_count
        || blocks[block_idx].base_pos != base) {
        out_offsets[tid] = 0u;
        out_counts[tid] = 0u;
        return;
    }
    const AniPosBlock *block = &blocks[block_idx];
    unsigned int bit = pos - base;
    unsigned int word = bit >> 6;       // bit / 64
    unsigned int bit_in_word = bit & 63;
    unsigned long long mask = block->masks[word];
    if (((mask >> bit_in_word) & 1ULL) == 0ULL) {
        out_offsets[tid] = 0u;
        out_counts[tid] = 0u;
        return;
    }
    unsigned int rank = 0u;
    for (unsigned int w = 0; w < word; ++w) {
        rank += __popcll(block->masks[w]);
    }
    unsigned long long lower_mask = (bit_in_word == 0u) ? 0ULL : ((1ULL << bit_in_word) - 1ULL);
    rank += __popcll(mask & lower_mask);
    unsigned int pos_idx = block->offsets_start + rank;
    out_offsets[tid] = pos_offsets[pos_idx];
    out_counts[tid] = pos_counts[pos_idx];
}
