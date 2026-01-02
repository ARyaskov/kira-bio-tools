// Compute MPH index: same logic as CPU AniIndex::lookup()

typedef struct {
    uchar chr_id;
    uint  pos;
    uint  ref_ofs;
    uint  alt_ofs;
    uint  id_ofs;
    uint  qual_ofs;
    uint  filter_ofs;
    uint  info_ofs;
    uint  info_len;
    uint  format_ofs;
    uint  samples_ofs;
} AniEntry;

__kernel void ani_lookup(
    __global const ulong* keys,
    __global const uint* g,
    uint m,
    __global const AniEntry* entries,
    __global long* out_idx,
    int n
) {
    int tid = get_global_id(0);
    if (tid >= n) return;

    ulong key = keys[tid];

    // MPH index
    uint slot = key % m;
    uint idx = g[slot];

    out_idx[tid] = (long)idx;
}
