#include <limits.h>

extern "C" {

struct AniPosContig {
    unsigned short chr_id;
    unsigned short pad;
    unsigned int min_pos;
    unsigned int max_pos;
    unsigned int block_start;
    unsigned int block_count;
};

struct AniPosBlock {
    unsigned int base_pos;
    unsigned int pad;
    unsigned long long masks[8];
    unsigned int offsets_start;
    unsigned int pad2;
};

__device__ __forceinline__ unsigned long long fxhash_rotl(unsigned long long v) {
    return (v << 5) | (v >> (64 - 5));
}

__device__ __forceinline__ unsigned long long fxhash_word(unsigned long long hash, unsigned long long word) {
    const unsigned long long SEED64 = 0x517cc1b727220a95ULL;
    unsigned long long v = fxhash_rotl(hash) ^ word;
    return v * SEED64;
}

__device__ __forceinline__ unsigned long long fxhash_bytes(const unsigned char* data, unsigned int len) {
    unsigned long long hash = 0ULL;
    hash = fxhash_word(hash, (unsigned long long)len);
    for (unsigned int i = 0; i < len; ++i) {
        hash = fxhash_word(hash, (unsigned long long)data[i]);
    }
    return hash;
}

__device__ __forceinline__ unsigned long long wymum(unsigned long long a, unsigned long long b) {
    unsigned long long hi = __umul64hi(a, b);
    unsigned long long lo = a * b;
    return hi ^ lo;
}

__device__ __forceinline__ unsigned long long read64_swapped(unsigned long long key) {
    unsigned int lo = (unsigned int)(key & 0xffffffffULL);
    unsigned int hi = (unsigned int)(key >> 32);
    return ((unsigned long long)lo << 32) | (unsigned long long)hi;
}

__device__ __forceinline__ unsigned long long wyhash8(unsigned long long key, unsigned long long seed) {
    const unsigned long long P0 = 0xa0761d6478bd642fULL;
    const unsigned long long P1 = 0xe7037ed1a0b428dbULL;
    const unsigned long long P5 = 0xeb44accab455d165ULL;
    unsigned long long s = seed ^ P0;
    unsigned long long r = read64_swapped(key) ^ P1;
    s = wymum(s, r);
    return wymum(s, (unsigned long long)8 ^ P5);
}

__device__ __forceinline__ unsigned long long splitmix64(unsigned long long x) {
    x = x + 0x9E3779B97F4A7C15ULL;
    unsigned long long z = x;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

__global__ void ani_lookup_kernel_v2(
    const unsigned long long* __restrict__ keys,
    const unsigned int* __restrict__ g,
    unsigned int m,
    unsigned int n,
    unsigned long long salt,
    const unsigned long long* __restrict__ entry_keys,
    unsigned int* __restrict__ out_idx,
    int nkeys
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= nkeys) return;

    unsigned long long key = keys[tid];
    unsigned long long base = wyhash8(key, salt);

    unsigned int v0 = (unsigned int)(splitmix64(base ^ 0x9E3779B97F4A7C15ULL) % (unsigned long long)m);
    unsigned int v1 = (unsigned int)(splitmix64(base + 0xA24B1F6FULL) % (unsigned long long)m);
    unsigned int v2 = (unsigned int)(splitmix64(base ^ 0x853C49E60A6C9D39ULL) % (unsigned long long)m);

    unsigned int idx = (g[v0] + g[v1] + g[v2]) % n;
    if (entry_keys[idx] == key) {
        out_idx[tid] = idx;
    } else {
        out_idx[tid] = 0xffffffffU;
    }
}

__global__ void ani_lookup_from_strings_kernel(
    const unsigned char* __restrict__ ref_pool,
    const unsigned int* __restrict__ ref_offsets,
    const unsigned int* __restrict__ ref_lens,
    const unsigned char* __restrict__ alt_pool,
    const unsigned int* __restrict__ alt_offsets,
    const unsigned int* __restrict__ alt_lens,
    const unsigned int* __restrict__ key_ref_idx,
    const unsigned char* __restrict__ key_chr,
    const unsigned int* __restrict__ key_pos,
    const unsigned int* __restrict__ g,
    unsigned int m,
    unsigned int n,
    unsigned long long salt,
    const unsigned long long* __restrict__ entry_keys,
    unsigned int* __restrict__ out_idx,
    int nkeys
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= nkeys) return;

    unsigned int ref_idx = key_ref_idx[tid];
    unsigned int ref_off = ref_offsets[ref_idx];
    unsigned int ref_len = ref_lens[ref_idx];
    unsigned int alt_off = alt_offsets[tid];
    unsigned int alt_len = alt_lens[tid];

    unsigned long long ref_hash = fxhash_bytes(ref_pool + ref_off, ref_len);
    unsigned long long alt_hash = fxhash_bytes(alt_pool + alt_off, alt_len);

    unsigned long long key = ((unsigned long long)key_chr[tid] << 32) | (unsigned long long)key_pos[tid];
    key ^= ref_hash;
    key ^= alt_hash;

    unsigned long long base = wyhash8(key, salt);

    unsigned int v0 = (unsigned int)(splitmix64(base ^ 0x9E3779B97F4A7C15ULL) % (unsigned long long)m);
    unsigned int v1 = (unsigned int)(splitmix64(base + 0xA24B1F6FULL) % (unsigned long long)m);
    unsigned int v2 = (unsigned int)(splitmix64(base ^ 0x853C49E60A6C9D39ULL) % (unsigned long long)m);

    unsigned int idx = (g[v0] + g[v1] + g[v2]) % n;
    if (entry_keys[idx] == key) {
        out_idx[tid] = idx;
    } else {
        out_idx[tid] = 0xffffffffU;
    }
}

__global__ void ani_pos_lookup_kernel(
    const unsigned char* __restrict__ chr_ids,
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

    unsigned char chr = chr_ids[tid];
    unsigned int pos = positions[tid];

    int contig_idx = -1;
    for (unsigned int i = 0; i < contig_count; ++i) {
        if (contigs[i].chr_id == chr) {
            contig_idx = (int)i;
            break;
        }
    }
    if (contig_idx < 0) {
        out_counts[tid] = 0;
        out_offsets[tid] = 0;
        return;
    }

    AniPosContig contig = contigs[contig_idx];
    if (pos < contig.min_pos || pos > contig.max_pos) {
        out_counts[tid] = 0;
        out_offsets[tid] = 0;
        return;
    }

    unsigned int base = (pos / 512) * 512;
    unsigned int lo = 0;
    unsigned int hi = contig.block_count;
    unsigned int block_idx = 0xffffffffU;
    while (lo < hi) {
        unsigned int mid = (lo + hi) >> 1;
        unsigned int idx = contig.block_start + mid;
        unsigned int bpos = blocks[idx].base_pos;
        if (bpos < base) {
            lo = mid + 1;
        } else if (bpos > base) {
            hi = mid;
        } else {
            block_idx = idx;
            break;
        }
    }
    if (block_idx == 0xffffffffU) {
        out_counts[tid] = 0;
        out_offsets[tid] = 0;
        return;
    }

    AniPosBlock block = blocks[block_idx];
    unsigned int bit = pos - base;
    unsigned int word = bit >> 6;
    unsigned int bit_in_word = bit & 63U;
    unsigned long long mask = block.masks[word];
    if (((mask >> bit_in_word) & 1ULL) == 0ULL) {
        out_counts[tid] = 0;
        out_offsets[tid] = 0;
        return;
    }

    unsigned int rank = 0;
    for (unsigned int w = 0; w < word; ++w) {
        rank += __popcll(block.masks[w]);
    }
    unsigned long long lower_mask = bit_in_word == 0 ? 0ULL : ((1ULL << bit_in_word) - 1ULL);
    rank += __popcll(mask & lower_mask);

    unsigned int pos_idx = block.offsets_start + rank;
    out_offsets[tid] = pos_offsets[pos_idx];
    out_counts[tid] = pos_counts[pos_idx];
}

__global__ void ani_info_merge_kernel(
    const unsigned int* __restrict__ alt_entry_idx,
    const unsigned int* __restrict__ alt_offsets,
    const unsigned short* __restrict__ alt_counts,
    const unsigned int* __restrict__ entry_offsets,
    const unsigned short* __restrict__ entry_counts,
    const unsigned int* __restrict__ pair_tag_ids,
    const unsigned int* __restrict__ pair_value_off,
    const unsigned int* __restrict__ pair_value_len,
    const unsigned char* __restrict__ raw_values,
    const unsigned char* __restrict__ tag_types,
    const unsigned int* __restrict__ tag_ids,
    unsigned int* __restrict__ out_pair_idx,
    int n_records,
    int n_tags
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int total = n_records * n_tags;
    if (tid >= total) return;

    int rec = tid / n_tags;
    int tag_i = tid - rec * n_tags;
    unsigned int tag_id = tag_ids[tag_i];
    if (tag_id == 0xffffffffU) {
        out_pair_idx[tid] = 0xffffffffU;
        return;
    }
    unsigned char tag_type = tag_types[tag_id];

    unsigned int base = alt_offsets[rec];
    unsigned int count = alt_counts[rec];
    unsigned int out = 0xffffffffU;

    for (unsigned int i = 0; i < count; ++i) {
        unsigned int entry = alt_entry_idx[base + i];
        if (entry == 0xffffffffU) continue;
        unsigned int eoff = entry_offsets[entry];
        unsigned int ecnt = entry_counts[entry];
        for (unsigned int j = 0; j < ecnt; ++j) {
            unsigned int pidx = eoff + j;
            if (pair_tag_ids[pidx] != tag_id) continue;
            if (tag_type == 3U) {
                out = pidx;
                break;
            }
            unsigned int vlen = pair_value_len[pidx];
            if (vlen == 0U) continue;
            unsigned int voff = pair_value_off[pidx];
            if (vlen == 1U && raw_values[voff] == '.') {
                continue;
            }
            if (vlen > 0U) {
                out = pidx;
                break;
            }
        }
        if (out != 0xffffffffU) break;
    }
    out_pair_idx[tid] = out;
}

}
