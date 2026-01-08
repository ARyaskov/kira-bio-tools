extern "C" {

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

}
