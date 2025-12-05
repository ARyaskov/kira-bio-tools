#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable

// SIMD-friendly xxhash-based GPU hash
inline ulong mix64(ulong x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdUL;
    x ^= x >> 29;
    x *= 0xc4ceb9fe1a85ec53UL;
    x ^= x >> 32;
    return x;
}

typedef struct {
    uint chr_pos;   // chr_id<<24 | pos
    uint ref_ofs;
    uint alt_ofs;
    uint info_ofs;
} AniEntryCL;

__kernel void ani_lookup_kernel_v2(
    __global const uint *chr_pos_in,   // N (chr_id<<24 | pos)
    __global const ulong *ref_hash_in, // N
    __global const ulong *alt_hash_in, // N
    __global const uint *g,            // m
    uint m,
    __global long *out_idx,            // N
    int n
){
    int gid = get_global_id(0);
    if (gid >= n) return;

    ulong h = (ulong)chr_pos_in[gid];

    h ^= ref_hash_in[gid];
    h = mix64(h);

    h ^= alt_hash_in[gid];
    h = mix64(h);

    uint slot = (uint)(h % m);
    uint idx  = g[slot];

    out_idx[gid] = (long)idx;
}
