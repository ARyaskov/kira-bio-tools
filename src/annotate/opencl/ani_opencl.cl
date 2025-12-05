typedef struct {
    uint chr_pos;
    uint ref_ofs;
    uint alt_ofs;
    uint info_ofs;
} AniEntryCL;

__kernel void ani_lookup_kernel(
    __global const ulong *keys,
    __global const uint  *g,
    uint m,
    __global const AniEntryCL *entries,
    __global long *out_idx,
    int n)
{
    int gid = get_global_id(0);
    if (gid >= n) return;

    ulong key = keys[gid];

    uint slot = (uint)(key % m);
    uint idx  = g[slot];

    out_idx[gid] = (long)idx;
}