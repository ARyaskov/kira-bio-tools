extern "C" {

struct AniEntry {
    unsigned char chr_id;
    unsigned int pos;
    unsigned int ref_ofs;
    unsigned int alt_ofs;
    unsigned int info_ofs;
};

__device__ __forceinline__
unsigned long long warp_reduce_xor(unsigned long long val)
{
    // FULL_MASK = 0xFFFFFFFFu
    for (int offset = 16; offset > 0; offset >>= 1)
        val ^= __shfl_down_sync(0xFFFFFFFFu, val, offset);
    return val;
}

__global__
void ani_lookup_kernel(
    const unsigned long long* __restrict__ keys,
    const unsigned int* __restrict__ g,
    unsigned int m,
    const AniEntry* __restrict__ entries,
    long long* __restrict__ out_idx,
    int n)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n) return;

    unsigned long long key = keys[tid];

    // --- Warp reduction (each warp has 32 threads) ---
    unsigned long long warp_hash = warp_reduce_xor(key);

    // Broadcast value from lane 0 to whole warp
    warp_hash = __shfl_sync(0xFFFFFFFFu, warp_hash, 0);

    // Mix key with warp-aggregated hash
    key ^= warp_hash;

    // --- MPH mapping ---
    unsigned int slot = (unsigned int)(key % m);
    unsigned int idx  = g[slot];

    out_idx[tid] = (long long)idx;
}

} // extern "C"