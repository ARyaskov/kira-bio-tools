#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable

// wyhash v1 constants (must match wyhash crate v0.6)
constant ulong P0 = 0xa0761d6478bd642fUL;
constant ulong P1 = 0xe7037ed1a0b428dbUL;
constant ulong P2 = 0x8ebc6af09c88c6e3UL;
constant ulong P3 = 0x589965cc75374cc3UL;
constant ulong P4 = 0x1d8e4e27c47d124fUL;
constant ulong P5 = 0xeb44accab455d165UL;
constant ulong FX_SEED64 = 0x517cc1b727220a95UL;

inline ulong wymum(ulong a, ulong b) {
    ulong hi = mul_hi(a, b);
    ulong lo = a * b;
    return hi ^ lo;
}

inline ulong read64_swapped(ulong key) {
    uint lo = (uint)(key & 0xffffffffUL);
    uint hi = (uint)(key >> 32);
    return ((ulong)lo << 32) | (ulong)hi;
}

inline ulong wyhash8(ulong key, ulong seed) {
    ulong s = seed ^ P0;
    ulong r = read64_swapped(key) ^ P1;
    s = wymum(s, r);
    return wymum(s, (ulong)8 ^ P5);
}

inline ulong splitmix64(ulong x) {
    x = x + 0x9E3779B97F4A7C15UL;
    ulong z = x;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBUL;
    return z ^ (z >> 31);
}

inline ulong fxhash_rotl(ulong v) {
    return (v << 5) | (v >> (64 - 5));
}

inline ulong fxhash_word(ulong hash, ulong word) {
    ulong v = fxhash_rotl(hash) ^ word;
    return v * FX_SEED64;
}

inline ulong fxhash_bytes(__global const uchar *data, uint len) {
    ulong hash = 0UL;
    hash = fxhash_word(hash, (ulong)len);
    for (uint i = 0; i < len; ++i) {
        hash = fxhash_word(hash, (ulong)data[i]);
    }
    return hash;
}

__kernel void ani_lookup_kernel_v2(
    __global const ulong *keys,   // N (u64 key bytes)
    __global const uint *g,        // m
    uint m,
    uint n,
    ulong salt,
    __global const ulong *entry_keys, // n entries
    __global uint *out_idx,        // N
    int nkeys
){
    int gid = get_global_id(0);
    if (gid >= nkeys) return;

    ulong key = keys[gid];
    ulong base = wyhash8(key, salt);

    uint v0 = (uint)(splitmix64(base ^ 0x9E3779B97F4A7C15UL) % (ulong)m);
    uint v1 = (uint)(splitmix64(base + 0xA24B1F6FUL) % (ulong)m);
    uint v2 = (uint)(splitmix64(base ^ 0x853C49E60A6C9D39UL) % (ulong)m);

    uint idx = (g[v0] + g[v1] + g[v2]) % n;
    if (entry_keys[idx] == key) {
        out_idx[gid] = idx;
    } else {
        out_idx[gid] = 0xffffffffU;
    }
}

__kernel void ani_lookup_from_strings_kernel(
    __global const uchar *ref_pool,
    __global const uint *ref_offsets,
    __global const uint *ref_lens,
    __global const uchar *alt_pool,
    __global const uint *alt_offsets,
    __global const uint *alt_lens,
    __global const uint *key_ref_idx,
    __global const uchar *key_chr,
    __global const uint *key_pos,
    __global const uint *g,
    uint m,
    uint n,
    ulong salt,
    __global const ulong *entry_keys,
    __global uint *out_idx,
    int nkeys
){
    int gid = get_global_id(0);
    if (gid >= nkeys) return;

    uint ref_idx = key_ref_idx[gid];
    uint ref_off = ref_offsets[ref_idx];
    uint ref_len = ref_lens[ref_idx];
    uint alt_off = alt_offsets[gid];
    uint alt_len = alt_lens[gid];

    ulong ref_hash = fxhash_bytes(ref_pool + ref_off, ref_len);
    ulong alt_hash = fxhash_bytes(alt_pool + alt_off, alt_len);

    ulong key = ((ulong)key_chr[gid] << 32) | (ulong)key_pos[gid];
    key ^= ref_hash;
    key ^= alt_hash;

    ulong base = wyhash8(key, salt);

    uint v0 = (uint)(splitmix64(base ^ 0x9E3779B97F4A7C15UL) % (ulong)m);
    uint v1 = (uint)(splitmix64(base + 0xA24B1F6FUL) % (ulong)m);
    uint v2 = (uint)(splitmix64(base ^ 0x853C49E60A6C9D39UL) % (ulong)m);

    uint idx = (g[v0] + g[v1] + g[v2]) % n;
    if (entry_keys[idx] == key) {
        out_idx[gid] = idx;
    } else {
        out_idx[gid] = 0xffffffffU;
    }
}
