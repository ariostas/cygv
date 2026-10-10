/*
 * GPU (CUDA) version of cgv's extraction stage.
 *
 * Given the union U of alpha supports (keys, degrees, alpha values) and the
 * initial instanton polynomial I, compute A_C = I[C]/deg(C) for every point,
 * layer by layer in degree, subtracting s_C z^C exp(C.alpha) for each curve C
 * with nonzero residual s_C = I[C]. Same arithmetic as the CPU path (NL lanes
 * of 62-bit Montgomery residues), so results are bit-identical.
 *
 * Within a layer every curve's exp is advanced level-synchronously: kernel e
 * processes, for all curves at once, the exp coefficients of internal degree
 * e (final by then), scattering deg(u) L[u] f[b] into the entries b+u of
 * higher internal degree. Entries live in a device hash table keyed by
 * (curve, point) and accumulate with u64 atomics; an add that wraps past 2^64
 * is corrected by adding 2^64 mod p, so the accumulator stays congruent to
 * the true sum mod p.
 */
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <strings.h>
#include <vector>
#include <algorithm>
#include "gpu_compat.h"   /* CUDA or HIP */
#include <time.h>

typedef unsigned __int128 u128;
typedef __int128 i128;
typedef uint64_t u64;
typedef int64_t i64;
typedef uint32_t u32;

#ifndef NL
#define NL 2
#endif
#include "gpu.h"

/* everything below lives in a per-lane-count namespace, so builds for several
 * NL can be linked into one executable */
#define CGV_CAT2(a, b) a##b
#define CGV_CAT(a, b) CGV_CAT2(a, b)
#define GPUNS CGV_CAT(cgv_gpu_nl, NL)
namespace GPUNS {
struct Co { u64 v[NL]; };

#define CK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { fprintf(stderr, "cgv gpu: %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); exit(1); } } while (0)
/* kernel launch failures (e.g. no memory left for local stacks) are not reported by a later
 * synchronize, only by cudaGetLastError: check it at every sync point and host read */
/* streams/events/pinned memory for the concurrent emit (c028) */
#ifdef CGV_HIP
#define cudaStream_t hipStream_t
#define cudaEvent_t hipEvent_t
#define cudaStreamCreateWithFlags hipStreamCreateWithFlags
#define cudaStreamNonBlocking hipStreamNonBlocking
#define cudaEventCreateWithFlags hipEventCreateWithFlags
#define cudaEventDisableTiming hipEventDisableTiming
#define cudaEventRecord hipEventRecord
#define cudaStreamWaitEvent hipStreamWaitEvent
#define cudaStreamSynchronize hipStreamSynchronize
#define cudaStreamDestroy hipStreamDestroy
#define cudaEventDestroy hipEventDestroy
#define HOST_ALLOC(pp, n) hipHostMalloc(pp, n, hipHostMallocDefault)
#define HOST_FREE(p) hipHostFree(p)
#else
#define HOST_ALLOC(pp, n) cudaMallocHost(pp, n)
#define HOST_FREE(p) cudaFreeHost(p)
#endif
#define SYNC() do { CK(cudaGetLastError()); CK(cudaDeviceSynchronize()); } while (0)

__constant__ u64 c_p[NL], c_pinv[NL], c_r2[NL], c_r64[NL], c_one[NL];
__constant__ int c_h11, c_keybits, c_maxdeg;

/* ---------------- field ---------------- */
__device__ __forceinline__ u64 mmul(u64 a, u64 b, int l) {
    u64 lo = a * b, hi = __umul64hi(a, b);
    u64 m = lo * c_pinv[l];
    u64 mhi = __umul64hi(m, c_p[l]);
    u64 u = hi + mhi + (lo != 0);
    return u >= c_p[l] ? u - c_p[l] : u;
}
__device__ __forceinline__ u64 madd(u64 a, u64 b, int l) { u64 c = a + b; return c >= c_p[l] ? c - c_p[l] : c; }
__device__ __forceinline__ u64 mneg(u64 a, int l) { return a ? c_p[l] - a : 0; }
__device__ __forceinline__ u64 mint(i64 x, int l) { /* small integer to Montgomery form */
    u64 p = c_p[l];
    u64 r = x >= 0 ? (u64)x % p : p - ((u64)(-x) % p);
    if (r == p) r = 0;
    return mmul(r, c_r2[l], l);
}
/* acc += x (x < p), keeping acc congruent mod p across 2^64 wraparound */
__device__ __forceinline__ void acc_add(u64 *acc, u64 x, int l) {
    for (;;) {
        u64 old = atomicAdd((unsigned long long *)acc, (unsigned long long)x);
        if (old + x >= old) return;
        x = c_r64[l];
    }
}

/* ---------------- keys ---------------- */
__device__ __forceinline__ u64 hash128(u128 k) {
    u64 x = (u64)k ^ ((u64)(k >> 64) * 0x9E3779B97F4A7C15ull);
    x ^= x >> 31; x *= 0xD6E8FEB86659FD93ull; x ^= x >> 32;
    return x;
}
__device__ void unpack_d(u128 k, int *v) {
    u128 mask = ((u128)1 << c_keybits) - 1;
    for (int t = 0; t < c_h11; t++) {
        i64 f = (i64)(u64)(k & mask);
        if (f >= ((i64)1 << (c_keybits - 1))) f -= (i64)1 << c_keybits;
        v[t] = (int)f;
        k = (k - (u128)(i128)f) >> c_keybits;
    }
}

/* ---------------- tables ---------------- */
/* state: 0 empty, 1 being written, 2 ready */
/* acc accumulates I; once degree deg has been extracted it holds A instead
 * (nothing adds to a point after its degree is extracted) */
struct IEnt { u128 key; int deg; u32 state; u64 acc[NL]; };
struct XEnt { u128 key; u32 curve; int deg; u32 state; u32 pad; u64 acc[NL]; };

/* load-acquire (gpu scope): later loads of the entry see what the writer
 * published before its release of the state flag */
__device__ __forceinline__ u32 ld_state(const u32 *s) { return gpu_ld_acquire(s); }
__device__ __forceinline__ void st_release(u32 *s, u32 v) { gpu_st_release(s, v); }

/* find or insert point k in the I table; returns position or -1 if full */
__device__ i64 i_find(IEnt *tab, u64 mask, u128 k, int deg, u32 *count, u32 cap) {
    u64 h = hash128(k) & mask;
    for (u64 probes = 0; probes <= mask; probes++) {
        IEnt *e = &tab[h];
        u32 st = ld_state(&e->state);
        if (st == 0) {
            if (atomicAdd(count, 0) >= cap) return -1;
            if (atomicCAS(&e->state, 0u, 1u) == 0u) {
                e->key = k; e->deg = deg;
                for (int l = 0; l < NL; l++) e->acc[l] = 0;
                st_release(&e->state, 2u);
                atomicAdd(count, 1u);
                return (i64)h;
            }
            st = ld_state(&e->state);
        }
#ifdef CGV_HIP
        if (st == 1) { probes--; continue; } /* being written: look again next iteration. AMD waves give diverged
                                                lanes no forward-progress guarantee, so no spin inside a branch */
#else
        while (st == 1) st = ld_state(&e->state);
#endif
        if (st == 0) { probes--; continue; } /* slot was released (overflow): look at it again */
        if (*(volatile u128 *)&e->key == k) return (i64)h;
        h = (h + 1) & mask;
    }
    return -1;
}
/* find or insert (curve, k) in the exp table; new entries are appended to order */
__device__ i64 x_find(XEnt *tab, u64 mask, u32 c, u128 k, int deg, u32 *norder, u32 *order, u32 cap,
                      int direct, u32 *lvlcnt, u32 *lvlpool, u32 stride) {
    u64 h = (hash128(k) ^ ((u64)c * 0xA24BAED4963EE407ull)) & mask;
    for (u64 probes = 0; probes <= mask; probes++) {
        XEnt *e = &tab[h];
        u32 st = ld_state(&e->state);
        if (__builtin_expect(st == 2, 1)) {   /* common case: occupied, compare and move on (found by ShinkaEvolve, +3-9%) */
            if (e->curve == c && *(u128 *)&e->key == k) return (i64)h;
            h = (h + 1) & mask;
            continue;
        }
        if (st == 0) {
            if (atomicCAS(&e->state, 0u, 1u) == 0u) {
                u32 o = atomicAdd(norder, 1u);
                if (o >= cap) { atomicExch(&e->state, 0u); return -1; } /* overflow: caller retries with a smaller batch */
                e->key = k; e->curve = c; e->deg = deg;
                for (int l = 0; l < NL; l++) e->acc[l] = 0;
                st_release(&e->state, 2u);
                order[o] = (u32)h;
                if (direct) {
                    /* warp-aggregated append: lanes inserting at the same level share one atomic */
                    unsigned act = warp_active();
                    unsigned grp = warp_match_any(act, deg);
                    int lane = threadIdx.x & 31, leader = __ffs(grp) - 1;
                    u32 q0 = 0;
                    if (lane == leader) q0 = atomicAdd(&lvlcnt[deg], (u32)__popc(grp));
                    q0 = warp_shfl(grp, q0, leader);
                    u32 q = q0 + __popc(grp & ((1u << lane) - 1));
                    if (q >= stride) return -1; /* level list full: treated as overflow */
                    lvlpool[(u64)deg * stride + q] = (u32)h;
                }
                return (i64)h;
            }
            st = ld_state(&e->state);
        }
#ifdef CGV_HIP
        if (st == 1) { probes--; continue; } /* being written: look again next iteration. AMD waves give diverged
                                                lanes no forward-progress guarantee, so no spin inside a branch */
#else
        while (st == 1) st = ld_state(&e->state);
#endif
        if (st == 0) { probes--; continue; } /* slot was released (overflow): look at it again */
        if (e->curve == c && *(u128 *)&e->key == k) return (i64)h;
        h = (h + 1) & mask;
    }
    return -1;
}

/* ---------------- device state ---------------- */
struct Dev {
    /* union U */
    int NU; u128 *ukey; int *udeg; Co *ualpha; /* [NU][h11] */
    Co *INV, *DEGC; int tabn;
    /* I table */
    IEnt *itab; u64 imask; u32 *icount; u32 icap;
    /* curves of the current layer */
    u128 *ckey; int *cdeg; Co *cs; int *cnu; u64 *coff; u32 *ncurves; u32 ccap;
    /* L pool */
    Co *L; u64 Lcap;
    /* exp table */
    XEnt *xtab; u64 xmask; u32 *norder; u32 *order; u32 xcap;
    /* level buckets */
    u32 *lvlcount; u32 *lvlpool; int *overflow;
    u64 *work, *woff; void *scan_tmp; size_t scan_tmp_bytes;
    int *Lidx;                 /* union index of each compacted L term */
    int *cnl;                  /* number of compacted L terms per curve (cnu = uncompacted prefix) */
    u64 *umask;                /* which alpha_t are nonzero at each union point */
    u32 *lvlcnt; u32 lvlstride; int direct; /* direct per-level lists: lvlpool[deg * lvlstride + i] */
};

/* initial I: insert all points */
__global__ void k_init_I(Dev d, int n, const u128 *key, const int *deg, const Co *val) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    i64 h = i_find(d.itab, d.imask, key[i], deg[i], d.icount, d.icap);
    if (h < 0) { *d.overflow = 1; return; }
    for (int l = 0; l < NL; l++) acc_add(&d.itab[h].acc[l], val[i].v[l], l);
}

/* smallest degree > dd among the I points (into *out, preset to INT_MAX) */
__global__ void k_next_deg(Dev d, int dd, int *out) {
    u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (i > d.imask) return;
    const IEnt *e = &d.itab[i];
    if (e->state == 2 && e->deg > dd) atomicMin(out, e->deg);
}

/* layer [dd, dhi): compute A for all I points with degree in the range; collect curves to apply */
/* how many curves k_extract will take from layer [dd, dhi) (same conditions, no writes), so the
 * curve buffers can grow first: k_extract converts values in place and cannot be re-run */
__global__ void k_count_extract(Dev d, int dd, int dhi, u32 *cnt) {
    u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (i > d.imask) return;
    const IEnt *e = &d.itab[i];
    if (e->state != 2 || e->deg < dd || e->deg >= dhi || e->deg >= c_maxdeg) return;
    int nz = 0;
    for (int l = 0; l < NL; l++) nz |= (e->acc[l] % c_p[l]) != 0;
    if (nz) atomicAdd(cnt, 1u);
}
__global__ void k_extract(Dev d, int dd, int dhi) {
    u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (i > d.imask) return;
    IEnt *e = &d.itab[i];
    if (e->state != 2 || e->deg < dd || e->deg >= dhi) return;
    int de = e->deg;
    Co v; int nz = 0;
    for (int l = 0; l < NL; l++) { v.v[l] = e->acc[l] % c_p[l]; nz |= v.v[l] != 0; }
    for (int l = 0; l < NL; l++) e->acc[l] = mmul(v.v[l], d.INV[de].v[l], l); /* now A */
    if (!nz || de >= c_maxdeg) return;
    u32 c = atomicAdd(d.ncurves, 1u);
    if (c >= d.ccap) { *d.overflow = 2; return; }
    d.ckey[c] = e->key; d.cdeg[c] = de; d.cs[c] = v;
    int T = c_maxdeg - de, lo = 0, hi = d.NU;
    while (lo < hi) { int mid = (lo + hi) / 2; if (d.udeg[mid] <= T) lo = mid + 1; else hi = mid; }
    d.cnu[c] = lo;
}

/* L[c] = nonzero terms of deg(u) sum_t C_t alpha_t[u] for u < cnu[c], in
 * order (compacted with a block scan); their count goes to cnl[c]. Block per curve. */
__global__ void k_buildL(Dev d, u32 c0) {
    typedef cub::BlockScan<int, 256> BS;
    __shared__ typename BS::TempStorage tmp;
    __shared__ int base;
    u32 c = c0 + blockIdx.x;
    int C[64];
    unpack_d(d.ckey[c], C);
    Co cm[64]; int nzt[64], nnz = 0;
    for (int t = 0; t < c_h11; t++) if (C[t]) { for (int l = 0; l < NL; l++) cm[nnz].v[l] = mint(C[t], l); nzt[nnz++] = t; }
    Co *Lc = d.L + d.coff[c];
    int *Li = d.Lidx + d.coff[c];
    int nu = d.cnu[c];
    u64 pat = 0;
    for (int q = 0; q < nnz; q++) pat |= (u64)1 << nzt[q];
    if (threadIdx.x == 0) base = 0;
    __syncthreads();
    for (int u0 = 0; u0 < nu; u0 += blockDim.x) {
        int u = u0 + threadIdx.x;
        Co v; int nz = 0;
        for (int l = 0; l < NL; l++) v.v[l] = 0;
        if (u < nu && (d.umask[u] & pat)) {
            /* ualpha already carries the factor deg(u) */
            const Co *al = d.ualpha + (size_t)u * c_h11;
            for (int q = 0; q < nnz; q++)
                for (int l = 0; l < NL; l++) v.v[l] = madd(v.v[l], mmul(cm[q].v[l], al[nzt[q]].v[l], l), l);
            for (int l = 0; l < NL; l++) nz |= v.v[l] != 0;
        }
        int pos, tot;
        BS(tmp).ExclusiveSum(nz ? 1 : 0, pos, tot);
        if (nz) { Lc[base + pos] = v; Li[base + pos] = u; }
        __syncthreads();
        if (threadIdx.x == 0) base += tot;
        __syncthreads();
    }
    if (threadIdx.x == 0) d.cnl[c] = base;
}

/* seed: entry (c, 0) with value 1 for each curve of the batch */
__global__ void k_seed(Dev d, u32 c0, u32 nc) {
    u32 i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nc) return;
    u32 c = c0 + i;
    i64 h = x_find(d.xtab, d.xmask, c, 0, 0, d.norder, d.order, d.xcap, d.direct, d.lvlcnt, d.lvlpool, d.lvlstride);
    if (h < 0) { *d.overflow = 3; return; }
    for (int l = 0; l < NL; l++) d.xtab[h].acc[l] = c_one[l];
}

/* histogram and bucket newly inserted entries by level */
__global__ void k_hist(Dev d, u32 lo, u32 hi) {
    u32 i = lo + blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= hi) return;
    atomicAdd(&d.lvlcount[d.xtab[d.order[i]].deg], 1u);
}
__global__ void k_bucket(Dev d, u32 lo, u32 hi, u32 *cursor) {
    u32 i = lo + blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= hi) return;
    u32 pos = d.order[i];
    u32 o = atomicAdd(&cursor[d.xtab[pos].deg], 1u);
    d.lvlpool[o] = pos;
}

/* level e, step 1: finalize f = acc / e for the entries of a segment (stored in
 * acc), and count the L terms each will scatter */
__global__ void k_final(Dev d, int e, const u32 *ent, u32 n) {
    u32 i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i > n) return;
    if (i == n) { d.work[n] = 0; return; }
    XEnt *x = &d.xtab[ent[i]];
    int nz = 0;
    for (int l = 0; l < NL; l++) {
        u64 a = x->acc[l] % c_p[l];
        u64 f = e == 0 ? a : mmul(a, d.INV[e].v[l], l);
        x->acc[l] = f; nz |= f != 0;
    }
    if (!nz) { d.work[i] = 0; return; }
    u32 c = x->curve;
    int lim = c_maxdeg - d.cdeg[c] - e;
    const int *Li = d.Lidx + d.coff[c];
    int lo = 0, hi = d.cnl[c];
    while (lo < hi) { int mid = (lo + hi) / 2; if (d.udeg[Li[mid]] <= lim) lo = mid + 1; else hi = mid; }
    d.work[i] = (u64)lo;
}
/* level e, step 2: one thread per (entry, L term) */
__global__ void k_scatter(Dev d, int e, const u32 *ent, u32 n, u64 W) {
    u64 t = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    /* entry of thread t = last i with woff[i] <= t. Consecutive t share or
     * neighbor entries, so two lanes search the whole range for the warp's
     * first and last t, and every lane then searches only between those. */
    int lane = threadIdx.x & 31;
    u64 tw0 = t - lane, tw1 = tw0 + 31 < W ? tw0 + 31 : W - 1;
    u32 wlo = 0;
    if (lane == 0 || lane == 31) {
        u64 tt = lane == 0 ? tw0 : tw1;
        u32 lo = 0, hi = n;
        while (hi - lo > 1) { u32 mid = (lo + hi) / 2; if (d.woff[mid] <= tt) lo = mid; else hi = mid; }
        wlo = lo;
    }
    u32 a = warp_shfl(0xffffffffu, wlo, 0), b = warp_shfl(0xffffffffu, wlo, 31);
    if (t >= W) return;
    u32 lo = a, hi = b + 1;
    while (hi - lo > 1) { u32 mid = (lo + hi) / 2; if (d.woff[mid] <= t) lo = mid; else hi = mid; }
    u32 i = lo;

    u64 j = t - d.woff[i];
    XEnt *x = &d.xtab[ent[i]];
    u32 c = x->curve;
    Co la = d.L[d.coff[c] + j];
    int u = d.Lidx[d.coff[c] + j];
    Co f;
    for (int l = 0; l < NL; l++) f.v[l] = x->acc[l];
    i64 h = x_find(d.xtab, d.xmask, c, x->key + d.ukey[u], e + d.udeg[u], d.norder, d.order, d.xcap, d.direct, d.lvlcnt, d.lvlpool, d.lvlstride);
    if (h < 0) { *d.overflow = 4; return; }
    for (int l = 0; l < NL; l++) acc_add(&d.xtab[h].acc[l], mmul(la.v[l], f.v[l], l), l);

}

/* Runs of small levels, all on one block (no host round trips): for each level
 * e from e0 on, finalize its entries, prefix-sum their work, and scatter. Stops
 * (writing e and its work total to out[0..1]) at a level with more than nmax
 * entries (not finalized: out[2] = 0) or more than wmax work items (finalized,
 * woff ready: out[2] = 1), out[3] = its entry count; out[0] = T + 1 when all levels are done,
 * and then out[1..3] = the overflow flag, the number of exp entries, the I table count, so the
 * host learns the batch's bookkeeping from this one copy (no further round trips). */
__global__ void __launch_bounds__(1024) k_small_levels(Dev d, int e0, int T, u32 nmax, u64 wmax, u64 *out) {
    typedef cub::BlockScan<u64, 1024> BS;
    __shared__ typename BS::TempStorage tmp;
    __shared__ u64 carry, Wtot;
    __shared__ u32 n_sh;
    __shared__ int ovf_sh;
    int tid = threadIdx.x;
    for (int e = e0; e <= T; e++) {
        __syncthreads();
        if (tid == 0) { n_sh = *(volatile u32 *)&d.lvlcnt[e]; ovf_sh = *(volatile int *)d.overflow; }
        __syncthreads();
        /* early abort: once the (sticky) overflow flag is set the batch is discarded and retried,
         * so the remaining levels are wasted work; go straight to the terminal report */
        if (ovf_sh) break;
        u32 n = n_sh;
        if (!n) continue;
        if (n > nmax || n > d.lvlstride) { if (tid == 0) { out[0] = e; out[1] = 0; out[2] = 0; out[3] = n; } return; }
        const u32 *ent = d.lvlpool + (u64)e * d.lvlstride;
        /* finalize + work counts */
        for (u32 i = tid; i < n; i += blockDim.x) {
            XEnt *x = &d.xtab[ent[i]];
            int nz = 0;
            for (int l = 0; l < NL; l++) {
                u64 a = x->acc[l] % c_p[l];
                u64 f = e == 0 ? a : mmul(a, d.INV[e].v[l], l);
                x->acc[l] = f; nz |= f != 0;
            }
            u64 w = 0;
            if (nz) {
                u32 c = x->curve;
                int lim = c_maxdeg - d.cdeg[c] - e;
                const int *Li = d.Lidx + d.coff[c];
                int lo = 0, hi = d.cnl[c];
                while (lo < hi) { int mid = (lo + hi) / 2; if (d.udeg[Li[mid]] <= lim) lo = mid + 1; else hi = mid; }
                w = lo;
            }
            d.work[i] = w;
        }
        __syncthreads();
        /* block-wide exclusive scan of work[0..n) into woff[0..n] */
        if (tid == 0) carry = 0;
        __syncthreads();
        for (u32 base = 0; base < n; base += blockDim.x) {
            u32 i = base + tid;
            u64 v = i < n ? d.work[i] : 0, pre, agg;
            BS(tmp).ExclusiveSum(v, pre, agg);
            if (i < n) d.woff[i] = carry + pre;
            __syncthreads();
            if (tid == 0) carry += agg;
            __syncthreads();
        }
        if (tid == 0) { d.woff[n] = carry; Wtot = carry; }
        __syncthreads();
        u64 W = Wtot;
        if (W > wmax) { if (tid == 0) { out[0] = e; out[1] = W; out[2] = 1; out[3] = n; } return; }
        /* scatter */
        for (u64 t = tid; t < W; t += blockDim.x) {
            u32 lo = 0, hi = n;
            while (hi - lo > 1) { u32 mid = (lo + hi) / 2; if (d.woff[mid] <= t) lo = mid; else hi = mid; }
            u64 j = t - d.woff[lo];
            XEnt *x = &d.xtab[ent[lo]];
            u32 c = x->curve;
            Co la = d.L[d.coff[c] + j];
            int u = d.Lidx[d.coff[c] + j];
            i64 h = x_find(d.xtab, d.xmask, c, x->key + d.ukey[u], e + d.udeg[u], d.norder, d.order, d.xcap, d.direct, d.lvlcnt, d.lvlpool, d.lvlstride);
            if (h < 0) { *d.overflow = 4; continue; }
            for (int l = 0; l < NL; l++) acc_add(&d.xtab[h].acc[l], mmul(la.v[l], x->acc[l], l), l);
        }
        __threadfence();
    }
    __syncthreads(); /* every thread's scatter (inserts, overflow flag) is done and fenced */
    if (tid == 0) {
        out[0] = T + 1;
        out[1] = (u64)(u32)atomicAdd((unsigned *)d.overflow, 0u);
        out[2] = atomicAdd(d.norder, 0u);
        out[3] = atomicAdd(d.icount, 0u);
    }
}

/* One cooperative kernel per batch: every level of the batch, with grid-wide
 * barriers between finalize / prefix sum / scatter, so no host round trips.
 * Needs direct level lists; bsum has one slot per block (+1). */
#define COOP_BS 256
__global__ void __launch_bounds__(COOP_BS) k_levels_coop(Dev d, int T, u64 *bsum) {
    cg::grid_group grid = cg::this_grid();
    typedef cub::BlockScan<u64, COOP_BS> BS;
    __shared__ typename BS::TempStorage tmp;
    __shared__ u64 carry;
    const u32 nb = gridDim.x, b = blockIdx.x, tid = threadIdx.x;
    const u64 gtid = (u64)b * blockDim.x + tid, gsz = (u64)nb * blockDim.x;
    for (int e = 0; e <= T; e++) {
        u32 n = *(volatile u32 *)&d.lvlcnt[e];
        if (!n) continue; /* same value in every thread: all inserts into e happened before the last barrier */
        if (n > d.lvlstride) { if (gtid == 0) *d.overflow = 6; return; }
        const u32 *ent = d.lvlpool + (u64)e * d.lvlstride;
        /* finalize + work counts */
        for (u64 i = gtid; i < n; i += gsz) {
            XEnt *x = &d.xtab[ent[i]];
            int nz = 0;
            for (int l = 0; l < NL; l++) {
                u64 a = x->acc[l] % c_p[l];
                u64 f = e == 0 ? a : mmul(a, d.INV[e].v[l], l);
                x->acc[l] = f; nz |= f != 0;
            }
            u64 w = 0;
            if (nz) {
                u32 c = x->curve;
                int lim = c_maxdeg - d.cdeg[c] - e;
                const int *Li = d.Lidx + d.coff[c];
                int lo = 0, hi = d.cnl[c];
                while (lo < hi) { int mid = (lo + hi) / 2; if (d.udeg[Li[mid]] <= lim) lo = mid + 1; else hi = mid; }
                w = lo;
            }
            d.work[i] = w;
        }
        grid.sync();
        /* prefix sum: each block scans its chunk, block 0 scans the block totals, blocks add offsets */
        u64 chunk = (n + nb - 1) / nb, c0 = (u64)b * chunk, c1 = c0 + chunk < n ? c0 + chunk : n;
        if (tid == 0) carry = 0;
        __syncthreads();
        for (u64 base = c0; base < c1; base += blockDim.x) {
            u64 i = base + tid;
            u64 v = i < c1 ? d.work[i] : 0, pre, agg;
            BS(tmp).ExclusiveSum(v, pre, agg);
            if (i < c1) d.woff[i] = carry + pre;
            __syncthreads();
            if (tid == 0) carry += agg;
            __syncthreads();
        }
        if (tid == 0) bsum[b] = carry;
        grid.sync();
        if (b == 0) {
            if (tid == 0) carry = 0;
            __syncthreads();
            for (u32 base = 0; base < nb; base += blockDim.x) {
                u32 i = base + tid;
                u64 v = i < nb ? bsum[i] : 0, pre, agg;
                BS(tmp).ExclusiveSum(v, pre, agg);
                __syncthreads();
                if (i < nb) bsum[i] = carry + pre;
                __syncthreads();
                if (tid == 0) carry += agg;
                __syncthreads();
            }
            if (tid == 0) bsum[nb] = carry;
        }
        grid.sync();
        {
            u64 off = bsum[b];
            for (u64 i = c0 + tid; i < c1; i += blockDim.x) d.woff[i] += off;
            if (gtid == 0) d.woff[n] = bsum[nb];
        }
        grid.sync();
        u64 W = bsum[nb];
        /* scatter */
        for (u64 t = gtid; t < W; t += gsz) {
            u32 lo = 0, hi = n;
            while (hi - lo > 1) { u32 mid = (lo + hi) / 2; if (d.woff[mid] <= t) lo = mid; else hi = mid; }
            u64 j = t - d.woff[lo];
            XEnt *x = &d.xtab[ent[lo]];
            u32 c = x->curve;
            Co la = d.L[d.coff[c] + j];
            int u = d.Lidx[d.coff[c] + j];
            i64 h = x_find(d.xtab, d.xmask, c, x->key + d.ukey[u], e + d.udeg[u], d.norder, d.order, d.xcap, d.direct, d.lvlcnt, d.lvlpool, d.lvlstride);
            if (h < 0) { *d.overflow = 4; continue; }
            for (int l = 0; l < NL; l++) acc_add(&d.xtab[h].acc[l], mmul(la.v[l], x->acc[l], l), l);
        }
        if (gtid == 0) bsum[nb + 1] += W; /* work counter for the profile */
        grid.sync();
    }
}

/* emit: I[C + b] -= s_C f_C[b] for every finished entry with b != 0 */
__global__ void k_emit(Dev d, u32 n) {
    u32 i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    XEnt *x = &d.xtab[d.order[i]];
    if (x->deg == 0) return;
    int nz = 0;
    for (int l = 0; l < NL; l++) nz |= x->acc[l] != 0;
    if (!nz) return;
    u32 c = x->curve;
    i64 h = i_find(d.itab, d.imask, d.ckey[c] + x->key, d.cdeg[c] + x->deg, d.icount, d.icap);
    if (h < 0) { *d.overflow = 5; return; }
    for (int l = 0; l < NL; l++) acc_add(&d.itab[h].acc[l], mneg(mmul(d.cs[c].v[l], x->acc[l], l), l), l);
}

__global__ void k_emit_range(Dev d, u32 lo, u32 hi) {
    u32 i = lo + blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= hi) return;
    XEnt *x = &d.xtab[d.order[i]];
    int nz = 0;
    for (int l = 0; l < NL; l++) nz |= x->acc[l] != 0;
    if (x->deg == 0 || !nz) return;
    u32 c = x->curve;
    i64 h = i_find(d.itab, d.imask, d.ckey[c] + x->key, d.cdeg[c] + x->deg, d.icount, d.icap);
    if (h < 0) { *d.overflow = 5; return; }
    for (int l = 0; l < NL; l++) acc_add(&d.itab[h].acc[l], mneg(mmul(d.cs[c].v[l], x->acc[l], l), l), l);
}
/* emit entries [lo, hi) of the order and clear their slots in the same pass (fused k_emit_range +
 * k_clear): an entry's emit reads only its own exp slot, and every order index names a distinct slot,
 * so clearing a slot right after its emit cannot affect any other entry's emit */
__global__ void k_emit_clear_range(Dev d, u32 lo, u32 hi, int *eflag) {
    u32 i = lo + blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= hi) return;
    XEnt *x = &d.xtab[d.order[i]];
    int nz = 0;
    for (int l = 0; l < NL; l++) nz |= x->acc[l] != 0;
    if (x->deg != 0 && nz) {
        u32 c = x->curve;
        i64 h = i_find(d.itab, d.imask, d.ckey[c] + x->key, d.cdeg[c] + x->deg, d.icount, d.icap);
        if (h < 0) { *eflag = 5; } else
        for (int l = 0; l < NL; l++) acc_add(&d.itab[h].acc[l], mneg(mmul(d.cs[c].v[l], x->acc[l], l), l), l);
    }
    x->state = 0;
    for (int l = 0; l < NL; l++) x->acc[l] = 0;
}
__global__ void k_emit_clear(Dev d, u32 n) {
    u32 i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    XEnt *x = &d.xtab[d.order[i]];
    int nz = 0;
    for (int l = 0; l < NL; l++) nz |= x->acc[l] != 0;
    if (x->deg != 0 && nz) {
        u32 c = x->curve;
        i64 h = i_find(d.itab, d.imask, d.ckey[c] + x->key, d.cdeg[c] + x->deg, d.icount, d.icap);
        if (h < 0) { *d.overflow = 5; } else
        for (int l = 0; l < NL; l++) acc_add(&d.itab[h].acc[l], mneg(mmul(d.cs[c].v[l], x->acc[l], l), l), l);
    }
    x->state = 0;
    for (int l = 0; l < NL; l++) x->acc[l] = 0;
}
__global__ void k_clear(Dev d, u32 n) {
    u32 i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    XEnt *x = &d.xtab[d.order[i]];
    x->state = 0;
    for (int l = 0; l < NL; l++) x->acc[l] = 0;
}

/* move every entry of an old I table into the (empty) current one */
__global__ void k_rehash_I(Dev d, const IEnt *old, u64 oldcap) {
    u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= oldcap) return;
    const IEnt *e = &old[i];
    if (e->state != 2) return;
    u64 h = hash128(e->key) & d.imask;
    for (;;) {
        if (atomicCAS(&d.itab[h].state, 0u, 1u) == 0u) {
            IEnt *t = &d.itab[h];
            t->key = e->key; t->deg = e->deg;
            for (int l = 0; l < NL; l++) t->acc[l] = e->acc[l];
            t->state = 2;
            return;
        }
        h = (h + 1) & d.imask;
    }
}

/* collect all I points with their A */
__global__ void k_collect(Dev d, u128 *okey, int *odeg, Co *oA, u32 *n) {
    u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (i > d.imask) return;
    IEnt *e = &d.itab[i];
    if (e->state != 2) return;
    u32 o = atomicAdd(n, 1u);
    okey[o] = e->key; odeg[o] = e->deg;
    for (int l = 0; l < NL; l++) oA[o].v[l] = e->acc[l];
}

static double wall(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }
enum { T_SETUP, T_EXTRACT, T_BUILDL, T_BUCKET, T_LEVEL, T_EMIT, T_CLEAR, T_WASTE, T_N };
static const char *TNAME[T_N] = {"setup", "extract", "buildL+seed", "bucket", "level", "emit", "clear", "retried(wasted)"};
static inline unsigned nblk(u64 n, int bs) { return (unsigned)((n + bs - 1) / bs); }
template <class T> static T *dalloc(size_t n) {
    T *p;
    cudaError_t e = cudaMalloc(&p, (n ? n : 1) * sizeof(T));
    if (e != cudaSuccess) {
        size_t f = 0, t = 0; cudaMemGetInfo(&f, &t);
        fprintf(stderr, "cgv gpu: cannot allocate %.2f GB (%.2f GB free): %s\n", (n ? n : 1) * sizeof(T) / 1e9, f / 1e9, cudaGetErrorString(e));
        exit(1);
    }
    return p;
}
template <class T> static void h2d(T *d, const T *h, size_t n) { if (n) CK(cudaMemcpy(d, h, n * sizeof(T), cudaMemcpyHostToDevice)); }
template <class T> static T d2h1(const T *d) { T v; CK(cudaGetLastError()); CK(cudaMemcpy(&v, d, sizeof(T), cudaMemcpyDeviceToHost)); return v; }

/* ======================= bundled path (layers with many curves) =======================
 * Up to 32 curves with the same nonzero-coordinate pattern (hence the same L and exp
 * supports) share one hash entry per point, with a value per curve: one probe serves
 * all of them and the per-curve atomics are coalesced. Lane k of a warp = curve k. */
struct XEntB { u128 key; u32 curve; int deg; u32 state; u32 pad; u64 acc[NL][32]; };

__device__ i64 x_findB(XEntB *tab, u64 mask, u32 c, u128 k, int deg, u32 *norder, u32 *order, u32 cap,
                       int direct, u32 *lvlcnt, u32 *lvlpool, u32 stride) {
    u64 h = (hash128(k) ^ ((u64)c * 0xA24BAED4963EE407ull)) & mask;
    for (u64 probes = 0; probes <= mask; probes++) {
        XEntB *e = &tab[h];
        u32 st = ld_state(&e->state);
        if (st == 0) {
            if (atomicCAS(&e->state, 0u, 1u) == 0u) {
                u32 o = atomicAdd(norder, 1u);
                if (o >= cap) { atomicExch(&e->state, 0u); return -1; }
                e->key = k; e->curve = c; e->deg = deg;
                st_release(&e->state, 2u);
                order[o] = (u32)h;
                if (direct) {
                    u32 q = atomicAdd(&lvlcnt[deg], 1u); /* callers are single lanes: no aggregation */
                    if (q >= stride) return -1;
                    lvlpool[(u64)deg * stride + q] = (u32)h;
                }
                return (i64)h;
            }
            st = ld_state(&e->state);
        }
#ifdef CGV_HIP
        if (st == 1) { probes--; continue; } /* being written: look again next iteration. AMD waves give diverged
                                                lanes no forward-progress guarantee, so no spin inside a branch */
#else
        while (st == 1) st = ld_state(&e->state);
#endif
        if (st == 0) { probes--; continue; }
        if (*(volatile u32 *)&e->curve == c && *(volatile u128 *)&e->key == k) return (i64)h;
        h = (h + 1) & mask;
    }
    return -1;
}
/* acc of new entries must be zero: the table is cleared entry by entry after each batch */

__global__ void k_pattern(Dev d, u32 nc, u64 *pat) {
    u32 c = blockIdx.x * blockDim.x + threadIdx.x;
    if (c >= nc) return;
    int C[64];
    unpack_d(d.ckey[c], C);
    u64 m = 0;
    for (int t = 0; t < c_h11; t++) if (C[t]) m |= (u64)1 << t;
    pat[c] = m;
}
/* one warp per bundle: L terms where any of its curves is nonzero, in union order */
__global__ void k_buildL_b(Dev d, const int *bcur, u32 b0, const u64 *bofs, int nuT, int *bcnl) {
    u32 beta = b0 + blockIdx.x;
    int k = threadIdx.x;
    int c = bcur[(u64)beta * 32 + k];
    Co cm[64]; int nzt[64], nnz = 0;
    if (c >= 0) {
        int C[64];
        unpack_d(d.ckey[c], C);
        for (int t = 0; t < c_h11; t++) if (C[t]) { for (int l = 0; l < NL; l++) cm[nnz].v[l] = mint(C[t], l); nzt[nnz++] = t; }
    }
    u64 *Lv = (u64 *)d.L;
    u64 base = bofs[beta - b0];
    int cnt = 0;
    for (int u = 0; u < nuT; u++) {
        const Co *al = d.ualpha + (size_t)u * c_h11;
        u64 v[NL]; int nz = 0;
        for (int l = 0; l < NL; l++) v[l] = 0;
        for (int q = 0; q < nnz; q++)
            for (int l = 0; l < NL; l++) v[l] = madd(v[l], mmul(cm[q].v[l], al[nzt[q]].v[l], l), l);
        for (int l = 0; l < NL; l++) nz |= v[l] != 0; /* ualpha carries deg(u) */
        if (!warp_any(0xffffffffu, nz)) continue;
        u64 j = base + cnt;
        for (int l = 0; l < NL; l++) Lv[(j * NL + l) * 32 + k] = v[l];
        if (k == 0) d.Lidx[j] = u;
        cnt++;
    }
    if (k == 0) bcnl[beta - b0] = cnt;
}
__global__ void k_seed_b(Dev d, XEntB *xb, const int *bcur, u32 b0, u32 nbb) {
    u32 i = blockIdx.x;
    if (i >= nbb) return;
    __shared__ i64 hs;
    if (threadIdx.x == 0) hs = x_findB(xb, d.xmask, b0 + i, 0, 0, d.norder, d.order, d.xcap, d.direct, d.lvlcnt, d.lvlpool, d.lvlstride);
    __syncthreads();
    if (hs < 0) { if (threadIdx.x == 0) *d.overflow = 3; return; }
    int k = threadIdx.x;
    int c = bcur[(u64)(b0 + i) * 32 + k];
    for (int l = 0; l < NL; l++) xb[hs].acc[l][k] = c >= 0 ? c_one[l] : 0;
}
/* warp per entry: finalize the 32 values, count L terms (per bundle) below the limit */
__global__ void k_final_b(Dev d, XEntB *xb, int e, const u32 *ent, u32 n, int T, u32 b0, const u64 *bofs, const int *bcnl) {
    u64 w = ((u64)blockIdx.x * blockDim.x + threadIdx.x) / 32;
    int k = threadIdx.x & 31;
    if (w > n) return;
    if (w == n) { if (k == 0) d.work[n] = 0; return; }
    XEntB *x = &xb[ent[w]];
    int nz = 0;
    for (int l = 0; l < NL; l++) {
        u64 a = x->acc[l][k] % c_p[l];
        u64 f = e == 0 ? a : mmul(a, d.INV[e].v[l], l);
        x->acc[l][k] = f; nz |= f != 0;
    }
    int any = warp_any(0xffffffffu, nz);
    if (k) return;
    if (!any) { d.work[w] = 0; return; }
    u32 beta = x->curve;
    int lim = T - e;
    const int *Li = d.Lidx + bofs[beta - b0];
    int lo = 0, hi = bcnl[beta - b0];
    while (lo < hi) { int mid = (lo + hi) / 2; if (d.udeg[Li[mid]] <= lim) lo = mid + 1; else hi = mid; }
    d.work[w] = (u64)lo;
}
/* warp per (entry, term): lane 0 probes, every lane adds its curve's product */
__global__ void k_scatter_b(Dev d, XEntB *xb, int e, const u32 *ent, u32 n, u64 W, u32 b0, const u64 *bofs) {
    u64 t = ((u64)blockIdx.x * blockDim.x + threadIdx.x) / 32;
    int k = threadIdx.x & 31;
    if (t >= W) return;
    u32 lo = 0, hi = n;
    while (hi - lo > 1) { u32 mid = (lo + hi) / 2; if (d.woff[mid] <= t) lo = mid; else hi = mid; }
    u64 j = t - d.woff[lo];
    XEntB *x = &xb[ent[lo]];
    u32 beta = x->curve;
    u64 jj = bofs[beta - b0] + j;
    int u = d.Lidx[jj];
    i64 h = 0;
    if (k == 0) {
        h = x_findB(xb, d.xmask, beta, x->key + d.ukey[u], e + d.udeg[u], d.norder, d.order, d.xcap, d.direct, d.lvlcnt, d.lvlpool, d.lvlstride);
        if (h < 0) *d.overflow = 4;
    }
    h = warp_shfl(0xffffffffu, h, 0);
    if (h < 0) return;
    const u64 *Lv = (const u64 *)d.L;
    for (int l = 0; l < NL; l++) {
        u64 la = Lv[(jj * NL + l) * 32 + k], fk = x->acc[l][k];
        if (la && fk) acc_add(&xb[h].acc[l][k], mmul(la, fk, l), l);
    }
}
/* warp per entry: each lane emits its curve's contribution, then the entry is cleared */
__global__ void k_emit_clear_b(Dev d, XEntB *xb, u32 n, const int *bcur) {
    u64 w = ((u64)blockIdx.x * blockDim.x + threadIdx.x) / 32;
    int k = threadIdx.x & 31;
    if (w >= n) return;
    XEntB *x = &xb[d.order[w]];
    int c = bcur[(u64)x->curve * 32 + k];
    int nz = 0;
    for (int l = 0; l < NL; l++) nz |= x->acc[l][k] != 0;
    if (c >= 0 && x->deg != 0 && nz) {
        i64 h = i_find(d.itab, d.imask, d.ckey[c] + x->key, d.cdeg[c] + x->deg, d.icount, d.icap);
        if (h < 0) *d.overflow = 5;
        else for (int l = 0; l < NL; l++) acc_add(&d.itab[h].acc[l], mneg(mmul(d.cs[c].v[l], x->acc[l][k], l), l), l);
    }
    for (int l = 0; l < NL; l++) x->acc[l][k] = 0;
    warp_sync();
    if (k == 0) x->state = 0;
}
__global__ void k_clear_b(Dev d, XEntB *xb, u32 n) {
    u64 w = ((u64)blockIdx.x * blockDim.x + threadIdx.x) / 32;
    int k = threadIdx.x & 31;
    if (w >= n) return;
    XEntB *x = &xb[d.order[w]];
    for (int l = 0; l < NL; l++) x->acc[l][k] = 0;
    warp_sync();
    if (k == 0) x->state = 0;
}

} /* namespace */
using namespace GPUNS;

extern "C" double gpu_free_gb(int device) {
    /* a GPU driving a display is off limits unless CGV_ALLOW_DISPLAY_GPU=1: filling its
     * memory and running long kernels on it can black out the desktop */
    char bus[32];
    if (cudaDeviceGetPCIBusId(bus, sizeof(bus), device) != cudaSuccess) { cudaGetLastError(); return -1; }
    int display = 0;
#ifndef CGV_HIP  /* AMD: no display check (compute servers); nvidia-smi is NVIDIA-only */
    {   /* nvidia-smi knows which GPU drives a display; match it by PCI bus id */
        FILE *p = popen("nvidia-smi --query-gpu=pci.bus_id,display_active --format=csv,noheader 2>/dev/null", "r");
        char line[256];
        while (p && fgets(line, sizeof(line), p)) {
            /* nvidia-smi prints an 8-digit domain (00000000:01:00.0), CUDA a 4-digit one (0000:01:00.0) */
            const char *a = strchr(line, ':'), *b = strchr(bus, ':');
            if (a && b && !strncasecmp(a, b, strlen(b)) && strstr(line, "Enabled")) display = 1;
        }
        if (p) pclose(p);
    }
#endif
    if (display && !getenv("CGV_ALLOW_DISPLAY_GPU")) {
        fprintf(stderr, "  (GPU %d drives a display; set CGV_ALLOW_DISPLAY_GPU=1 to use it anyway)\n", device);
        return -1;
    }
    if (cudaSetDevice(device) != cudaSuccess) { cudaGetLastError(); return -1; }
    size_t f = 0, t = 0;
    if (cudaMemGetInfo(&f, &t) != cudaSuccess) { cudaGetLastError(); return -1; }
    return f / 1e9;
}

extern "C" int gpu_extract(const GpuIn *in, GpuOut *out) {
    double tt[T_N] = {0}, tmark = wall(), tb;
    int dev = in->device;
    CK(cudaSetDevice(dev));
    cudaDeviceProp prop; CK(cudaGetDeviceProperties(&prop, dev));
    if (prop.warpSize != 32) { fprintf(stderr, "cgv gpu: needs 32-lane warps (this device has %d); use the CPU\n", prop.warpSize); return 1; }
    size_t freemem, totmem; CK(cudaMemGetInfo(&freemem, &totmem));
    /* memory budget for the exp table and L pool: all free memory on a discrete GPU; on an
     * integrated GPU (shared with the CPU) 24 GB unless CGV_GPU_MEM_GB says otherwise. The I table
     * still grows into whatever is actually free. */
    size_t budget = getenv("CGV_GPU_MEM_GB") ? (size_t)(atof(getenv("CGV_GPU_MEM_GB")) * 1e9) : prop.integrated ? (size_t)24e9 : (size_t)-1;
    size_t free0 = freemem;
    auto budget_free = [&](size_t f) { size_t used = free0 > f ? free0 - f : 0; return used >= budget ? (size_t)0 : std::min(f, budget - used); };
    if (in->verbose) fprintf(stderr, "  gpu %d: %s, %.1f GB free\n", dev, prop.name, freemem / 1e9);

    u64 one[NL], r64[NL];
    for (int l = 0; l < NL; l++) {
        u64 p = in->p[l];
        u128 r = ((u128)1 << 64) % p;
        r64[l] = (u64)r;
        one[l] = (u64)r; /* Montgomery form of 1 is 2^64 mod p */
    }
    CK(cudaMemcpyToSymbol(c_p, in->p, sizeof(u64) * NL));
    CK(cudaMemcpyToSymbol(c_pinv, in->pinv, sizeof(u64) * NL));
    CK(cudaMemcpyToSymbol(c_r2, in->r2, sizeof(u64) * NL));
    CK(cudaMemcpyToSymbol(c_r64, r64, sizeof(u64) * NL));
    CK(cudaMemcpyToSymbol(c_one, one, sizeof(u64) * NL));
    CK(cudaMemcpyToSymbol(c_h11, &in->h11, sizeof(int)));
    CK(cudaMemcpyToSymbol(c_keybits, &in->keybits, sizeof(int)));
    CK(cudaMemcpyToSymbol(c_maxdeg, &in->maxdeg, sizeof(int)));
    int H = in->h11, D = in->maxdeg;

    Dev d; memset(&d, 0, sizeof(d));
    d.NU = in->NU;
    d.ukey = dalloc<u128>(in->NU); h2d(d.ukey, (const u128 *)in->ukey, in->NU);
    d.udeg = dalloc<int>(in->NU); h2d(d.udeg, in->udeg, in->NU);
    {
        /* alpha_t[u] * deg(u) (the Euler-recurrence weight), and each point's mask of nonzero alpha_t */
        std::vector<Co> ua((size_t)in->NU * H);
        std::vector<u64> um(in->NU);
        const Co *src = (const Co *)in->ualpha, *degc = (const Co *)in->DEGC;
        for (int u = 0; u < in->NU; u++) {
            const Co &dg = degc[in->udeg[u]];
            u64 m = 0;
            for (int t = 0; t < H; t++) {
                Co x = src[(size_t)u * H + t], y;
                int nz = 0;
                for (int l = 0; l < NL; l++) {
                    nz |= x.v[l] != 0;
                    u128 prod = (u128)x.v[l] * dg.v[l]; /* Montgomery product on the host */
                    u64 lo = (u64)prod, mm = lo * in->pinv[l];
                    u128 t2 = prod + (u128)mm * in->p[l];
                    u64 r = (u64)(t2 >> 64);
                    y.v[l] = r >= in->p[l] ? r - in->p[l] : r;
                }
                ua[(size_t)u * H + t] = y;
                if (nz) m |= (u64)1 << t;
            }
            um[u] = m;
        }
        d.ualpha = dalloc<Co>((size_t)in->NU * H); h2d(d.ualpha, ua.data(), (size_t)in->NU * H);
        d.umask = dalloc<u64>(in->NU); h2d(d.umask, um.data(), in->NU);
    }
    /* for batch sizing: which alpha_t are nonzero at each union point, and per mask
     * the number of union points of degree <= g (so a curve's L size is known from
     * its nonzero pattern without building L) */
    std::vector<u64> mvals;
    std::vector<std::vector<u32>> mcnt; /* mcnt[id][g] = # union points with that mask and degree <= g */
    {
        const Co *ua = (const Co *)in->ualpha;
        std::vector<int> mid(in->NU);
        for (int u = 0; u < in->NU; u++) {
            u64 m = 0;
            for (int t = 0; t < in->h11; t++) { int nz = 0; for (int l = 0; l < NL; l++) nz |= ua[(size_t)u * in->h11 + t].v[l] != 0; if (nz) m |= (u64)1 << t; }
            int id = -1;
            for (size_t q = 0; q < mvals.size(); q++) if (mvals[q] == m) { id = (int)q; break; }
            if (id < 0) { id = (int)mvals.size(); mvals.push_back(m); mcnt.emplace_back(in->maxdeg + 1, 0); }
            mid[u] = id;
        }
        for (int u = 0; u < in->NU; u++) if (in->udeg[u] <= in->maxdeg) mcnt[mid[u]][in->udeg[u]]++;
        for (auto &v : mcnt) for (int g = 1; g <= in->maxdeg; g++) v[g] += v[g - 1];
    }
    d.tabn = in->tabn;
    d.INV = dalloc<Co>(in->tabn); h2d(d.INV, (const Co *)in->INV, in->tabn);
    d.DEGC = dalloc<Co>(in->tabn); h2d(d.DEGC, (const Co *)in->DEGC, in->tabn);
    d.overflow = dalloc<int>(1); CK(cudaMemset(d.overflow, 0, sizeof(int)));

    /* I table: sized from an estimate, 4x the initial points (points grow ~2x) */
    /* I-table load: grown between layers once above i_tgt full; inserts fail (run stops, never a wrong
     * result) above i_capf. 0.7/0.85 (was 0.4/0.5): same speed, one fewer doubling on deep runs
     * (e.g. D32 peak 13.7-16.8 -> 10.6 GB). CGV_I_TGT / CGV_I_CAPF override. */
    const double i_tgt = getenv("CGV_I_TGT") ? atof(getenv("CGV_I_TGT")) : 0.7, i_capf = getenv("CGV_I_CAPF") ? atof(getenv("CGV_I_CAPF")) : 0.85;
    u64 icap = 1; while (icap < (u64)in->nI * 8 + 1024) icap <<= 1;
    if (in->icap_log2) { icap = (u64)1 << in->icap_log2; while (icap < (u64)in->nI * 2) icap <<= 1; }
    d.imask = icap - 1; d.icap = (u32)std::min<u64>((u64)(icap * i_capf), 0xffffffffu);
    d.itab = dalloc<IEnt>(icap); CK(cudaMemset(d.itab, 0, icap * sizeof(IEnt)));
    d.icount = dalloc<u32>(1); CK(cudaMemset(d.icount, 0, sizeof(u32)));
    {
        u128 *k = dalloc<u128>(in->nI); int *g = dalloc<int>(in->nI); Co *v = dalloc<Co>(in->nI);
        h2d(k, (const u128 *)in->ikey, in->nI); h2d(g, in->ideg, in->nI); h2d(v, (const Co *)in->ival, in->nI);
        k_init_I<<<nblk(in->nI, 256), 256>>>(d, in->nI, k, g, v);
        SYNC();
        cudaFree(k); cudaFree(g); cudaFree(v);
    }

    /* curve buffers */
    d.ccap = getenv("CGV_CCAP_LOG2") ? 1u << atoi(getenv("CGV_CCAP_LOG2")) : 1u << 22;  /* test knob: start small to exercise growth */
    auto alloc_curves = [&](u32 n) {
        d.ccap = n;
        d.ckey = dalloc<u128>(n); d.cdeg = dalloc<int>(n); d.cs = dalloc<Co>(n);
        d.cnu = dalloc<int>(n); d.coff = dalloc<u64>(n); d.cnl = dalloc<int>(n);
    };
    auto free_curves = [&]() { cudaFree(d.ckey); cudaFree(d.cdeg); cudaFree(d.cs); cudaFree(d.cnu); cudaFree(d.coff); cudaFree(d.cnl); };
    d.ckey = dalloc<u128>(d.ccap); d.cdeg = dalloc<int>(d.ccap); d.cs = dalloc<Co>(d.ccap);
    d.cnu = dalloc<int>(d.ccap); d.coff = dalloc<u64>(d.ccap); d.ncurves = dalloc<u32>(1);
    d.cnl = dalloc<int>(d.ccap);

    /* exp table (+ order, level pool, scan buffers) gets half the free memory, the L pool most of the rest */
    CK(cudaMemGetInfo(&freemem, &totmem)); freemem = budget_free(freemem);
    /* per table slot: the entry, plus order, lvlpool, work, woff at half the slot count */
    size_t per = sizeof(XEnt) + (2 * sizeof(u32) + 2 * sizeof(u64)) / 2;
    u64 xcapn = 1; while ((xcapn * 2) * per <= (size_t)(freemem * 0.7)) xcapn <<= 1;
    /* start the exp table at what batches need (it grows on demand, see the retry path): zeroing a
     * table sized to all of memory costs ~1 s per run on an integrated GPU */
    u64 xcap_max = xcapn, xinit = (u64)1 << (getenv("CGV_EXP_LOG2") ? atoi(getenv("CGV_EXP_LOG2")) : 24);
    if (xcapn > xinit) xcapn = xinit;
    if (in->xcap_log2) xcapn = xcap_max = (u64)1 << in->xcap_log2;
    d.xmask = xcapn - 1; d.xcap = (u32)std::min<u64>(xcapn / 2, 0xffffffffu);
    /* the exp table and its per-slot arrays; reallocated smaller if the I table needs the room */
    /* double-buffered exp table (c028): batch k's emit + clear runs on its own stream over one buffer
     * while batch k+1's levels run over the other; ev_emit[b] marks the end of the last emit on buffer b */
    XEnt *xtb[2] = {NULL, NULL}; u32 *ordb[2] = {NULL, NULL};
    int dbl = 0, dbl_off = getenv("CGV_NO_DBL") != NULL, cur = 0, ev_used[2] = {0, 0};
    const u64 fill_free8 = getenv("CGV_FILL_FREE8") ? atoi(getenv("CGV_FILL_FREE8")) : 2; /* batch sub-tables may fill to m - m*k/8 slots (k = 2: 3/4) before a batch counts as overflowing */
    u64 dbl_slots = 0; /* capacity of buffer 1 (slots; its order array has as many entries) */
    auto alloc_exp = [&](u64 n) {
        d.xmask = n - 1; d.xcap = (u32)std::min<u64>(n / 2, 0xffffffffu);
        d.xtab = dalloc<XEnt>(n); CK(cudaMemset(d.xtab, 0, n * sizeof(XEnt)));
        d.order = dalloc<u32>(d.xcap);
        xtb[0] = d.xtab; ordb[0] = d.order; xtb[1] = NULL; ordb[1] = NULL; cur = 0; dbl = 0; ev_used[0] = ev_used[1] = 0;
        /* the second buffer only takes batch-sized tables (batches target ~1e6 entries, i.e. a few M
         * slots): capped at 2^CGV_DBL_LOG2 slots (default 2^22, 256 MB) instead of a second full table,
         * and only allocated when it leaves most of the free memory alone. Bigger batches use buffer 0. */
        dbl_slots = std::min<u64>(n, (u64)1 << (getenv("CGV_DBL_LOG2") ? atoi(getenv("CGV_DBL_LOG2")) : 22));
        size_t f_, t_; CK(cudaMemGetInfo(&f_, &t_));
        if (!dbl_off && dbl_slots * (sizeof(XEnt) + sizeof(u32)) * 8 <= budget_free(f_)) {
            XEnt *x2 = NULL; u32 *o2 = NULL;
            if (cudaMalloc(&x2, dbl_slots * sizeof(XEnt)) == cudaSuccess && cudaMalloc(&o2, dbl_slots * sizeof(u32)) == cudaSuccess) {
                CK(cudaMemset(x2, 0, dbl_slots * sizeof(XEnt)));
                xtb[1] = x2; ordb[1] = o2; dbl = 1;
            } else {
                if (x2) cudaFree(x2);
                (void)cudaGetLastError(); /* a failed allocation is not an error here: run single-buffered */
            }
        }
        d.work = dalloc<u64>((u64)d.xcap + 1); d.woff = dalloc<u64>((u64)d.xcap + 1);
        d.scan_tmp = NULL; d.scan_tmp_bytes = 0;
        cub::DeviceScan::ExclusiveSum(NULL, d.scan_tmp_bytes, d.work, d.woff, (int)std::min<u64>((u64)d.xcap + 1, 0x7fffffff));
        CK(cudaMalloc(&d.scan_tmp, d.scan_tmp_bytes));
    };
    auto free_exp = [&]() {
        CK(cudaDeviceSynchronize()); /* no emit may still use either buffer */
        cudaFree(xtb[0]); cudaFree(ordb[0]); if (xtb[1]) cudaFree(xtb[1]); if (ordb[1]) cudaFree(ordb[1]);
        xtb[0] = xtb[1] = NULL; ordb[0] = ordb[1] = NULL; dbl = 0;
        cudaFree(d.work); cudaFree(d.woff); cudaFree(d.scan_tmp);
    };
    alloc_exp(xcapn);
    /* level lists: 4 bytes per entry, room for direct lists of every level at the table's cap */
    u64 lvl_cap = 0;
    auto alloc_lvl = [&](u64 cnt) { if (d.lvlpool) cudaFree(d.lvlpool); lvl_cap = cnt; d.lvlpool = dalloc<u32>(cnt); };
    d.lvlpool = NULL; alloc_lvl((u64)(D + 2) * d.xcap);
    d.norder = dalloc<u32>(1);
    d.lvlcount = dalloc<u32>(D + 2);
    u32 *cursor = dalloc<u32>(D + 2);
    d.lvlcnt = dalloc<u32>(D + 2);
    CK(cudaMemGetInfo(&freemem, &totmem)); freemem = budget_free(freemem);
    void *test_hold = NULL; /* test knob: occupy the device so only CGV_TEST_LEAVE_MB stays free */
    if (getenv("CGV_TEST_LEAVE_MB")) {
        size_t leave = (size_t)atoll(getenv("CGV_TEST_LEAVE_MB")) << 20;
        if (freemem > leave) { CK(cudaMalloc(&test_hold, freemem - leave)); CK(cudaMemGetInfo(&freemem, &totmem)); }
    }
    /* leave room for the I table to grow (a growth step holds the old and a 2x-4x new table) */
    size_t reserve = std::min<size_t>(icap * sizeof(IEnt) * 4, freemem / 2);
    size_t launch_room = (size_t)256 << 20; /* kernel launches need device memory for local stacks */
    size_t Lbytes = freemem > reserve + launch_room ? (size_t)((freemem - reserve - launch_room) * 0.8) : 0;
    size_t Lmax = (size_t)((getenv("CGV_LPOOL_GB") ? atof(getenv("CGV_LPOOL_GB")) : 2.0) * 1e9);  /* grows on demand */
    if (Lbytes > Lmax) Lbytes = Lmax;
    d.Lcap = std::max<u64>(Lbytes / (sizeof(Co) + sizeof(int)), 1 << 16);
    d.L = dalloc<Co>(d.Lcap);
    d.Lidx = dalloc<int>(d.Lcap);
    /* the L pool is only used while a batch builds its exp entries, so an I-table growth may
     * release it; it is reallocated (to what is then free) before the next batch */
    u64 L_want = 0;  /* entries a growth asked for */
    auto ensure_L = [&]() {
        if (d.L) return;
        size_t f = 0, t = 0; CK(cudaMemGetInfo(&f, &t));
        size_t keep = icap * sizeof(IEnt) + ((size_t)256 << 20); /* room for the next (staged) I growth */
        u64 n = f > keep ? (u64)((f - keep) * 0.8) / (sizeof(Co) + sizeof(int)) : 0;
        size_t room = (size_t)256 << 20; /* for kernel launches */
        n = std::max<u64>(n, f > room ? (u64)((f - room) * 0.5) / (sizeof(Co) + sizeof(int)) : 0);
        n = std::min<u64>(n, Lmax / (sizeof(Co) + sizeof(int)));   /* same cap as at the start */
        n = std::max<u64>(std::max<u64>(n, L_want), 1 << 16);
        d.Lcap = n; d.L = dalloc<Co>(n); d.Lidx = dalloc<int>(n);
        if (in->verbose) fprintf(stderr, "  gpu: L pool reallocated, %.2f GB\n", n * sizeof(Co) / 1e9);
    };
    if (in->verbose) fprintf(stderr, "  gpu tables: I %llu, exp %llu entries, L pool %.2f GB\n",
                             (unsigned long long)icap, (unsigned long long)xcapn, d.Lcap * sizeof(Co) / 1e9);

    tt[T_SETUP] = wall() - tmark;
    std::vector<u32> hcnt(D + 2);
    std::vector<int> hnu; std::vector<u64> hoff;
    i64 napplied = 0; int nbatches = 0, nretry = 0;
    u64 totW = 0;
    u64 full_mask = d.xmask; u32 full_cap = d.xcap;
    u64 xfloor = 0; /* largest exp table a single curve has needed: shrinking below it only causes regrowth */
    u64 force_mask = 0;
    double est_nu = 1; /* exp entries per L term */
    /* bundled layers */
    int bundle_min = getenv("CGV_BUNDLE") ? atoi(getenv("CGV_BUNDLE")) : 0;
    int xmode = 0;
    u64 *d_pat = NULL; u32 pat_cap = 0;
    int *d_bcur = NULL; size_t bc_cap = 0;
    u64 *d_bofs = NULL; int *d_bcnl = NULL; u32 bb_cap = 0;
    double est_b = 1, xtarget_b = getenv("CGV_XTARGET_B") ? atof(getenv("CGV_XTARGET_B")) : 4e6;
    double xtarget = getenv("CGV_XTARGET") ? atof(getenv("CGV_XTARGET")) : 1e6;
    /* grow the I table so that it can take `extra` more points (below the hard cap), or (extra = 0) to below i_tgt */
    /* between layers: keep below i_tgt full; before an emit: just make sure the hard cap (i_capf) holds */
    auto grow_I = [&](u64 extra, i64 known = -1) { /* known: the I count, if the host has it already */
        u32 cnt = known >= 0 ? (u32)known : d2h1(d.icount);
        if (extra ? (u64)cnt + extra <= d.icap : (double)cnt <= icap * i_tgt) return;
        u64 ncap = icap * 2;
        if (!extra) while ((double)cnt > ncap * i_tgt) ncap *= 2;
        else while ((double)(cnt + extra) > ncap * i_capf) ncap *= 2;
        IEnt *old = d.itab; u64 ocap = icap;
        size_t nbytes = ncap * sizeof(IEnt), obytes = ocap * sizeof(IEnt), slack = (size_t)64 << 20; /* left free for kernel launches */
        size_t f = 0, t = 0; CK(cudaMemGetInfo(&f, &t));
        if (nbytes + slack <= f) {
            /* room for both tables: rehash on the device */
            d.itab = dalloc<IEnt>(ncap); CK(cudaMemset(d.itab, 0, ncap * sizeof(IEnt)));
            icap = ncap; d.imask = icap - 1; d.icap = (u32)std::min<u64>((u64)(icap * i_capf), 0xffffffffu);
            k_rehash_I<<<nblk(ocap, 256), 256>>>(d, old, ocap);
            SYNC();
            cudaFree(old);
            if (in->verbose) fprintf(stderr, "  gpu: I table grown to %llu entries\n", (unsigned long long)icap);
            return;
        }
        /* staged: park the old table in host memory, so only the new one is on the device.
         * If that still does not fit, release the L pool, then (between layers only, when the
         * exp table is empty) halve the exp table until it does. */
        const char *how = "staged via host";
        if (nbytes + slack > f + obytes && dbl && cur == 0) {
            /* memory is short: go single-buffered (buffer 1 is idle: the host synced the emit stream) */
            SYNC(); cudaFree(xtb[1]); cudaFree(ordb[1]); xtb[1] = NULL; ordb[1] = NULL; dbl = 0; dbl_off = 1;
            CK(cudaMemGetInfo(&f, &t));
            if (nbytes + slack <= f) {
                d.itab = dalloc<IEnt>(ncap); CK(cudaMemset(d.itab, 0, ncap * sizeof(IEnt)));
                icap = ncap; d.imask = icap - 1; d.icap = (u32)std::min<u64>((u64)(icap * i_capf), 0xffffffffu);
                k_rehash_I<<<nblk(ocap, 256), 256>>>(d, old, ocap);
                SYNC();
                cudaFree(old);
                if (in->verbose) fprintf(stderr, "  gpu: I table grown to %llu entries (exp table single-buffered)\n", (unsigned long long)icap);
                return;
            }
        }
        if (nbytes + slack > f + obytes && d.L) {
            cudaFree(d.L); cudaFree(d.Lidx); d.L = NULL; d.Lidx = NULL; d.Lcap = 0;
            CK(cudaMemGetInfo(&f, &t)); how = "staged via host, L pool released";
        }
        while (nbytes + slack > f + obytes && !extra && full_mask + 1 > std::max<u64>((u64)1 << 20, xfloor)) {
            u64 n = (full_mask + 1) / 2;
            free_exp(); alloc_exp(n);
            full_mask = d.xmask; full_cap = d.xcap; force_mask = 0;
            CK(cudaMemGetInfo(&f, &t)); how = "staged via host, exp table halved";
            if (in->verbose) fprintf(stderr, "  gpu: exp table reduced to %llu entries\n", (unsigned long long)n);
        }
        if (nbytes + slack > f + obytes) {
            fprintf(stderr, "cgv gpu: I table needs %.2f GB but only %.2f GB can be freed on the device "
                            "(run on the CPU, or a GPU with more memory)\n", nbytes / 1e9, (f + obytes) / 1e9);
            exit(1);
        }
        double tg = wall();
        std::vector<IEnt> host(ocap);
        CK(cudaMemcpy(host.data(), old, obytes, cudaMemcpyDeviceToHost));
        cudaFree(old);
        d.itab = dalloc<IEnt>(ncap); CK(cudaMemset(d.itab, 0, ncap * sizeof(IEnt)));
        icap = ncap; d.imask = icap - 1; d.icap = (u32)std::min<u64>((u64)(icap * i_capf), 0xffffffffu);
        CK(cudaMemGetInfo(&f, &t));
        u64 chunk = std::max<u64>(1 << 16, std::min<u64>(ocap, (u64)((f - std::min(f, slack / 2)) / sizeof(IEnt))));
        chunk = std::min<u64>(chunk, (u64)1 << 24);
        IEnt *buf = dalloc<IEnt>(chunk);
        for (u64 lo = 0; lo < ocap; lo += chunk) {
            u64 n = std::min(chunk, ocap - lo);
            h2d(buf, host.data() + lo, n);
            k_rehash_I<<<nblk(n, 256), 256>>>(d, buf, n);
            SYNC();
        }
        cudaFree(buf);
        if (in->verbose) fprintf(stderr, "  gpu: I table grown to %llu entries (%s, %.2fs)\n", (unsigned long long)icap, how, wall() - tg);
    };
    int *d_next = dalloc<int>(1);
    /* layers span super_delta degrees: see gv.c (curves only change points >= deg + delta) */
    int super_delta = in->NU ? in->udeg[0] : D + 1;
    if (super_delta < 1 || getenv("CGV_NO_SUPER")) super_delta = 1;
    u64 *d_sl = dalloc<u64>(4);
    cudaStream_t s_emit; CK(cudaStreamCreateWithFlags(&s_emit, cudaStreamNonBlocking));
    cudaEvent_t ev_emit[2], ev_lev;
    CK(cudaEventCreateWithFlags(&ev_emit[0], cudaEventDisableTiming)); CK(cudaEventCreateWithFlags(&ev_emit[1], cudaEventDisableTiming));
    CK(cudaEventCreateWithFlags(&ev_lev, cudaEventDisableTiming));
    int *d_eovf = dalloc<int>(1); CK(cudaMemset(d_eovf, 0, sizeof(int)));
    u32 *h_icnt = NULL; if (HOST_ALLOC((void **)&h_icnt, sizeof(u32)) != cudaSuccess) { h_icnt = NULL; (void)cudaGetLastError(); }
    int icnt_ok = 0; /* *h_icnt holds the I count after the last queued emit (valid once the emit stream is synced) */
    /* cooperative kernel: as many blocks as can be resident at once */
    int coop_blocks = 0;
    if (getenv("CGV_COOP")) { /* off by default: slower than per-level launches on standard gradings */
        int per_sm = 0;
        CK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&per_sm, k_levels_coop, COOP_BS, 0));
        coop_blocks = per_sm * prop.multiProcessorCount;
        if (in->verbose) fprintf(stderr, "  gpu: cooperative level kernel, %d blocks of %d\n", coop_blocks, COOP_BS);
    }
    u64 *d_bsum = dalloc<u64>(coop_blocks + 2);
    int small_levels = getenv("CGV_SMALL_LEVELS") ? atoi(getenv("CGV_SMALL_LEVELS")) : 1;
    u32 sl_nmax = getenv("CGV_SL_N") ? atoi(getenv("CGV_SL_N")) : 4096;
    u64 sl_wmax = getenv("CGV_SL_W") ? atoll(getenv("CGV_SL_W")) : 8192;
    for (int dd = 1; dd <= D; ) {
        /* jump to the next degree that has points (gradings can have sparse degrees) */
        {
            int big = 0x7fffffff;
            CK(cudaMemcpy(d_next, &big, sizeof(int), cudaMemcpyHostToDevice));
            k_next_deg<<<nblk(icap, 256), 256>>>(d, dd - 1, d_next);
            int nd = d2h1(d_next);
            if (nd > D) break;
            dd = nd;
        }
        /* keep the I table below i_tgt full; points roughly double over the run */
        grow_I(0);
        /* headroom for one more (staged) doubling during this layer's emits, which cannot shrink
         * the exp table: free memory plus the releasable L pool must cover the current table */
        {
            size_t f = 0, t = 0; CK(cudaMemGetInfo(&f, &t));
            size_t need = icap * sizeof(IEnt) + ((size_t)256 << 20), Lb = d.L ? d.Lcap * (sizeof(Co) + sizeof(int)) : 0;
            while (f + Lb < need && full_mask + 1 > std::max<u64>((u64)1 << 20, xfloor)) {
                u64 n = (full_mask + 1) / 2;
                free_exp(); alloc_exp(n);
                full_mask = d.xmask; full_cap = d.xcap; force_mask = 0;
                CK(cudaMemGetInfo(&f, &t));
                if (in->verbose) fprintf(stderr, "  gpu: exp table reduced to %llu entries (headroom for I growth)\n", (unsigned long long)n);
            }
        }
        tb = wall();
        int dhi = dd + super_delta;
        {   /* grow the curve buffers first if this layer has more curves than they hold */
            CK(cudaMemset(d.ncurves, 0, sizeof(u32)));
            k_count_extract<<<nblk(icap, 256), 256>>>(d, dd, dhi, d.ncurves);
            u32 want = d2h1(d.ncurves);
            if (want > d.ccap) {
                u32 n = d.ccap; while (n < want + want / 4) n *= 2;
                free_curves(); alloc_curves(n);
                if (in->verbose) fprintf(stderr, "  gpu: curve buffers grown to %u (layer %d has %u curves)\n", n, dd, want);
            }
        }
        CK(cudaMemset(d.ncurves, 0, sizeof(u32)));
        k_extract<<<nblk(icap, 256), 256>>>(d, dd, dhi);
        SYNC();
        tt[T_EXTRACT] += wall() - tb;
        if (d2h1(d.overflow)) { fprintf(stderr, "cgv gpu: overflow %d\n", d2h1(d.overflow)); return 1; }
        u32 nc = d2h1(d.ncurves);
        if (!nc) { dd = dhi; continue; }
        napplied += nc;
        hnu.resize(nc); hoff.resize(nc + 1);
        u64 layW = 0; double layT0 = wall(), layLev0 = tt[T_LEVEL]; u64 layEnt = 0; int layB = 0;
        CK(cudaMemcpy(hnu.data(), d.cnu, nc * sizeof(int), cudaMemcpyDeviceToHost));
        /* L size of each curve from its nonzero pattern (for batch sizing) */
        std::vector<double> hnl(nc);
        {
            if (pat_cap < nc) { cudaFree(d_pat); pat_cap = nc; d_pat = dalloc<u64>(pat_cap); }
            k_pattern<<<nblk(nc, 256), 256>>>(d, nc, d_pat);
            std::vector<u64> hp(nc);
            std::vector<int> hdeg(nc);
            CK(cudaMemcpy(hp.data(), d_pat, nc * sizeof(u64), cudaMemcpyDeviceToHost));
            CK(cudaMemcpy(hdeg.data(), d.cdeg, nc * sizeof(int), cudaMemcpyDeviceToHost));
            std::vector<std::pair<std::pair<u64, int>, double>> cache;
            for (u32 c = 0; c < nc; c++) {
                int Tl = D - hdeg[c];
                double v = -1;
                for (auto &pr : cache) if (pr.first.first == hp[c] && pr.first.second == Tl) { v = pr.second; break; }
                if (v < 0) {
                    v = 0;
                    for (size_t q = 0; q < mvals.size(); q++) if (mvals[q] & hp[c]) v += mcnt[q][Tl];
                    cache.push_back({{hp[c], Tl}, v});
                }
                hnl[c] = v;
            }
        }
        if (bundle_min && nc >= (u32)bundle_min && super_delta == 1) { /* bundles assume one degree per layer */
            /* ---------------- bundled layer ---------------- */
            if (xmode != 1) { CK(cudaMemset(d.xtab, 0, (full_mask + 1) * sizeof(XEnt))); xmode = 1; }
            XEntB *xb = (XEntB *)d.xtab;
            u64 capB = 1; while (capB * 2 * sizeof(XEntB) <= (full_mask + 1) * sizeof(XEnt)) capB *= 2;
            int T = D - dd, nuT = hnu[0];
            /* bundles: curves sorted by nonzero pattern, up to 32 of one pattern */
            if (pat_cap < nc) { cudaFree(d_pat); pat_cap = nc; d_pat = dalloc<u64>(pat_cap); }
            k_pattern<<<nblk(nc, 256), 256>>>(d, nc, d_pat);
            std::vector<u64> hpat(nc);
            CK(cudaMemcpy(hpat.data(), d_pat, nc * sizeof(u64), cudaMemcpyDeviceToHost));
            std::vector<int> idx(nc);
            for (u32 i = 0; i < nc; i++) idx[i] = (int)i;
            std::stable_sort(idx.begin(), idx.end(), [&](int a, int b) { return hpat[a] < hpat[b]; });
            std::vector<int> hbc;
            for (u32 i = 0; i < nc;) {
                u32 j = i;
                while (j < nc && hpat[idx[j]] == hpat[idx[i]] && j - i < 32) j++;
                for (u32 q = i; q < i + 32; q++) hbc.push_back(q < j ? idx[q] : -1);
                i = j;
            }
            u32 nbund = (u32)(hbc.size() / 32);
            if (bc_cap < hbc.size()) { cudaFree(d_bcur); bc_cap = hbc.size(); d_bcur = dalloc<int>(bc_cap); }
            h2d(d_bcur, hbc.data(), hbc.size());
            if (bb_cap < nbund) { cudaFree(d_bofs); cudaFree(d_bcnl); bb_cap = nbund; d_bofs = dalloc<u64>(bb_cap); d_bcnl = dalloc<int>(bb_cap); }
            ensure_L();
            u64 Lterms = d.Lcap / 32; /* L pool as bundle terms (32*NL values + index each) */
            u32 b0 = 0, bmaxb = nbund;
            u64 force = 0;
            while (b0 < nbund) {
                u32 nb = std::min(bmaxb, nbund - b0);
                if (!force) nb = std::min<u64>(nb, std::max<u64>(1, (u64)(xtarget_b / (est_b * nuT + 1))));
                nb = std::min<u64>(nb, std::max<u64>(1, Lterms / (u64)(nuT ? nuT : 1)));
                if ((u64)nb * nuT > Lterms) { fprintf(stderr, "cgv gpu: L pool too small for a bundle\n"); return 1; }
                u64 need = (u64)(est_b * nuT * nb * 3) + 1024, m = 1 << 16;
                while (m < need && m < capB) m <<= 1;
                if (force) m = force;
                if (m > capB) m = capB;
                d.xmask = m - 1;
                d.xcap = (u32)std::min<u64>(m / 2, lvl_cap / (u64)(T + 1)); /* direct level lists always fit */
                d.direct = 1; d.lvlstride = d.xcap;
                CK(cudaMemset(d.lvlcnt, 0, (D + 2) * sizeof(u32)));
                CK(cudaMemset(d.norder, 0, sizeof(u32)));
                std::vector<u64> hofs(nb);
                for (u32 i = 0; i < nb; i++) hofs[i] = (u64)i * nuT;
                h2d(d_bofs, hofs.data(), nb);
                double tbatch = wall(); tb = tbatch;
                k_buildL_b<<<nb, 32>>>(d, d_bcur, b0, d_bofs, nuT, d_bcnl);
                k_seed_b<<<nb, 32>>>(d, xb, d_bcur, b0, nb);
                SYNC();
                double tbl = wall() - tb, tlev = 0;
                int ovf = d2h1(d.overflow) != 0;
                std::vector<u32> hlc(T + 2);
                for (int e = 0; e <= T && !ovf; e++) {
                    tb = wall();
                    if (e > 0) {
                        CK(cudaMemcpy(hlc.data() + e, d.lvlcnt + e, (T + 1 - e) * sizeof(u32), cudaMemcpyDeviceToHost));
                        while (e <= T && !hlc[e]) e++;
                        if (e > T) break;
                    }
                    u32 n = d2h1(d.lvlcnt + e);
                    if (n > d.lvlstride) { ovf = 1; break; }
                    const u32 *ent = d.lvlpool + (u64)e * d.lvlstride;
                    k_final_b<<<nblk(((u64)n + 1) * 32, 256), 256>>>(d, xb, e, ent, n, T, b0, d_bofs, d_bcnl);
                    CK(cub::DeviceScan::ExclusiveSum(d.scan_tmp, d.scan_tmp_bytes, d.work, d.woff, (int)(n + 1)));
                    u64 W = d2h1(d.woff + n);
                    totW += W * 32; layW += W * 32;
                    if (W) k_scatter_b<<<nblk(W * 32, 256), 256>>>(d, xb, e, ent, n, W, b0, d_bofs);
                    tlev += wall() - tb;
                }
                SYNC();
                if (d2h1(d.overflow)) ovf = 1;
                u32 nent = std::min(d2h1(d.norder), d.xcap);
                if (ovf) {
                    CK(cudaMemset(d.overflow, 0, sizeof(int)));
                    k_clear_b<<<nblk((u64)nent * 32, 256), 256>>>(d, xb, nent);
                    SYNC();
                    est_b *= 2;
                    if (d.xmask + 1 < capB) { force = std::min<u64>((d.xmask + 1) * 4, capB); bmaxb = nb; }
                    else { if (nb == 1) { fprintf(stderr, "cgv gpu: a single bundle overflows the exp table\n"); return 1; } bmaxb = std::max<u32>(1, nb / 2); force = capB; }
                    nretry++;
                    tt[T_WASTE] += wall() - tbatch;
                    continue;
                }
                tt[T_BUILDL] += tbl; tt[T_LEVEL] += tlev;
                tb = wall();
                grow_I((u64)nent * 32);
                k_emit_clear_b<<<nblk((u64)nent * 32, 256), 256>>>(d, xb, nent, d_bcur);
                SYNC();
                tt[T_EMIT] += wall() - tb;
                { int ov = d2h1(d.overflow); if (ov) { fprintf(stderr, "cgv gpu: overflow flag %d after emit (bundled layer %d)\n", ov, dd); return 1; } }
                nbatches++; layB++; layEnt += nent;
                est_b = std::max((double)nent / ((double)nb * (nuT ? nuT : 1)), 0.7 * est_b);
                b0 += nb; force = 0; bmaxb = nbund;
            }
            d.xmask = full_mask; d.xcap = full_cap;
            if (getenv("CGV_PROF") && in->verbose)
                fprintf(stderr, "    gpu layer %d (bundled, %u bundles): curves %u batches %d entries %.3g W %.3g  level %.3fs layer %.3fs\n", dd, nbund, nc, layB,
                        (double)layEnt, (double)layW, tt[T_LEVEL] - layLev0, wall() - layT0);
            dd = dhi;
            continue;
        }
        if (xmode != 0) { CK(cudaMemset(d.xtab, 0, (full_mask + 1) * sizeof(XEnt))); xmode = 0; }
        /* process the layer in batches of curves whose L fits and whose exp entries fit */
        ensure_L();
        u32 c0 = 0;
        u32 bmax = nc;
        while (c0 < nc) {
            ensure_L(); /* an I growth during the last emit may have released it */
            u32 nb = std::min(bmax, nc - c0);
            /* size the batch so its exp entries (estimated) stay near xtarget,
             * i.e. a table that mostly lives in L2 */
            if (!force_mask) {
                /* entries scale with the curves' L sizes: take curves until the estimate reaches xtarget */
                double acc_e = 0; u32 k = 0;
                for (; k < nb; k++) { acc_e += est_nu * hnl[c0 + k] + 1; if (acc_e > xtarget && k > 0) break; }
                nb = std::max<u32>(1, k);
            }
            /* L offsets within the pool */
            u64 tot = 0; u32 k = 0;
            for (; k < nb; k++) { if (tot + hnu[c0 + k] > d.Lcap) break; hoff[c0 + k] = tot; tot += hnu[c0 + k]; }
            nb = k;
            if (!nb) {
                /* one curve needs more L terms than the pool holds: grow the pool and retry */
                if ((u64)hnu[c0] <= d.Lcap || L_want >= (u64)hnu[c0]) { fprintf(stderr, "cgv gpu: L pool too small\n"); return 1; }
                L_want = (u64)(hnu[c0] * 1.25) + 1024;
                cudaFree(d.L); cudaFree(d.Lidx); d.L = NULL; d.Lidx = NULL;
                ensure_L();
                if (in->verbose) fprintf(stderr, "  gpu: L pool grown to %.2f GB (one curve needs it)\n", d.Lcap * sizeof(Co) / 1e9);
                continue;
            }
            {
                double sumnu = 0; for (u32 k = 0; k < nb; k++) sumnu += hnl[c0 + k];
                u64 need = (u64)((est_nu * sumnu + nb) * 3) + 1024, m = 1 << 16;
                while (m < need && m <= full_mask) m <<= 1;
                if (force_mask) m = force_mask;
                if (m > full_mask + 1) m = full_mask + 1;
                /* overflow tolerance (c025): the table is still sized at 3x the estimate (same load in the normal
                 * case), but a batch may fill it to 3/4 of its slots before it counts as overflowing (was 1/2).
                 * A batch that is 1.5-2.25x over its estimate now finishes at a higher load instead of being
                 * discarded and redone in a 4x table. The full-size table keeps cap = half (full_cap), so order,
                 * work, woff and the direct level lists (lvl_cap = (D+2)*full_cap) are sized as before. */
                d.xmask = m - 1; d.xcap = (u32)std::min<u64>(m <= full_mask ? m - m * fill_free8 / 8 : m / 2, full_cap);
            }
            int T = D - dd;
            d.direct = (u64)(T + 1) * d.xcap <= lvl_cap;
            d.lvlstride = d.xcap;
            if (d.direct) CK(cudaMemset(d.lvlcnt, 0, (D + 2) * sizeof(u32)));
            double tbatch = wall(); tb = tbatch;
            double tlev = 0, tbuck = 0;
            if (cur == 1 && (u64)d.xmask + 1 > dbl_slots) cur = 0; /* too big for buffer 1 */
            d.xtab = xtb[cur]; d.order = ordb[cur];
            CK(cudaMemcpy(d.coff + c0, hoff.data() + c0, nb * sizeof(u64), cudaMemcpyHostToDevice));
            /* the last emit over this buffer must have cleared it before the seed inserts */
            if (ev_used[cur]) CK(cudaStreamWaitEvent(0, ev_emit[cur], 0));
            CK(cudaMemset(d.norder, 0, sizeof(u32)));
            k_buildL<<<nb, 256>>>(d, c0);
            k_seed<<<nblk(nb, 256), 256>>>(d, c0, nb);
            /* no sync here: the levels are queued behind buildL/seed on the default stream, and
             * the (sticky) overflow flag is read after the levels as before */
            CK(cudaGetLastError());
            double tbl = wall() - tb;
            /* level buckets: segments[e] = list of (offset, count) in lvlpool */
            std::vector<std::vector<std::pair<u32, u32>>> seg(T + 1);
            u32 done = 0, poolused = 0;
            int ovf = 0;
            int fin = 0; u64 fin_so[4] = {0, 0, 0, 0}; /* terminal k_small_levels report (overflow, norder, icount) */
            std::vector<u32> hlc(T + 2);
            for (int e = 0; e <= T && !ovf; e++) {
                if (d.direct && coop_blocks) {
                    tb = wall();
                    CK(cudaMemset(d_bsum, 0, (coop_blocks + 2) * sizeof(u64)));
                    int Tc = T;
                    void *args[] = {(void *)&d, (void *)&Tc, (void *)&d_bsum};
                    CK(cudaLaunchCooperativeKernel((void *)k_levels_coop, coop_blocks, COOP_BS, args, 0, 0));
                    SYNC();
                    u64 Wb = d2h1(d_bsum + coop_blocks + 1);
                    totW += Wb; layW += Wb;
                    tlev += wall() - tb;
                    if (d2h1(d.overflow)) ovf = 1;
                    break;
                }
                if (d.direct && small_levels) {
                    /* runs of small levels on one block; big levels come back to the host */
                    tb = wall();
                    k_small_levels<<<1, 1024>>>(d, e, T, sl_nmax, sl_wmax, d_sl);
                    u64 so[4];
                    CK(cudaGetLastError());
                    CK(cudaMemcpy(so, d_sl, sizeof(so), cudaMemcpyDeviceToHost));
                    if ((int)so[0] > T) {
                        tlev += wall() - tb;
                        fin = 1; for (int q = 0; q < 4; q++) fin_so[q] = so[q];
                        if (so[1]) ovf = 1;
                        break;
                    }
                    e = (int)so[0];
                    u32 n = (u32)so[3];
                    if (n > d.lvlstride) { ovf = 1; break; }
                    const u32 *ent = d.lvlpool + (u64)e * d.lvlstride;
                    u64 W = so[1];
                    if (!so[2]) {
                        k_final<<<nblk((u64)n + 1, 256), 256>>>(d, e, ent, n);
                        CK(cub::DeviceScan::ExclusiveSum(d.scan_tmp, d.scan_tmp_bytes, d.work, d.woff, (int)(n + 1)));
                        W = d2h1(d.woff + n);
                    }
                    totW += W; layW += W;
                    if (W) k_scatter<<<nblk(W, 256), 256>>>(d, e, ent, n, W);
                    tlev += wall() - tb;
                    continue;
                }
                if (d.direct && e > 0) {
                    /* next nonempty level: entries only go to higher levels, so all counts are final up to it */
                    CK(cudaMemcpy(hlc.data() + e, d.lvlcnt + e, (T + 1 - e) * sizeof(u32), cudaMemcpyDeviceToHost));
                    while (e <= T && !hlc[e]) e++;
                    if (e > T) { SYNC(); if (d2h1(d.overflow)) ovf = 1; break; }
                }
                if (d.direct) {
                    /* two host round trips per level: the level size, and the work total;
                     * the (sticky) overflow flag is checked once after the last level */
                    tb = wall();
                    u32 n = d2h1(d.lvlcnt + e);
                    if (n > d.lvlstride) { ovf = 1; break; }
                    if (n) {
                        const u32 *ent = d.lvlpool + (u64)e * d.lvlstride;
                        k_final<<<nblk((u64)n + 1, 256), 256>>>(d, e, ent, n);
                        CK(cub::DeviceScan::ExclusiveSum(d.scan_tmp, d.scan_tmp_bytes, d.work, d.woff, (int)(n + 1)));
                        u64 W = d2h1(d.woff + n);
                        totW += W; layW += W;
                        if (W) k_scatter<<<nblk(W, 256), 256>>>(d, e, ent, n, W);
                    }
                    tlev += wall() - tb;
                    if (e == T) { SYNC(); if (d2h1(d.overflow)) ovf = 1; }
                    continue;
                }
                /* bucket entries inserted since the last round */
                tb = wall();
                u32 now = d2h1(d.norder);
                if (now > d.xcap) { ovf = 1; break; }
                if (now > done) {
                    CK(cudaMemset(d.lvlcount, 0, (D + 2) * sizeof(u32)));
                    k_hist<<<nblk(now - done, 256), 256>>>(d, done, now);
                    CK(cudaMemcpy(hcnt.data(), d.lvlcount, (D + 2) * sizeof(u32), cudaMemcpyDeviceToHost));
                    std::vector<u32> start(D + 2);
                    for (int g = 0; g <= D; g++) { start[g] = poolused; if (hcnt[g]) { if (g <= T) seg[g].push_back({poolused, hcnt[g]}); poolused += hcnt[g]; } }
                    CK(cudaMemcpy(cursor, start.data(), (D + 2) * sizeof(u32), cudaMemcpyHostToDevice));
                    k_bucket<<<nblk(now - done, 256), 256>>>(d, done, now, cursor);
                    done = now;
                }
                SYNC();
                tbuck += wall() - tb; tb = wall();
                for (auto &sg : seg[e]) {
                    const u32 *ent = d.lvlpool + sg.first;
                    u32 n = sg.second;
                    k_final<<<nblk((u64)n + 1, 256), 256>>>(d, e, ent, n);
                    CK(cub::DeviceScan::ExclusiveSum(d.scan_tmp, d.scan_tmp_bytes, d.work, d.woff, (int)(n + 1)));
                    u64 W = d2h1(d.woff + n);
                    totW += W; layW += W;
                    if (W) k_scatter<<<nblk(W, 256), 256>>>(d, e, ent, n, W);
                }
                SYNC();
                tlev += wall() - tb;
                if (d2h1(d.overflow)) { ovf = 1; }
            }
            if (!ovf && !fin) { SYNC(); if (d2h1(d.overflow)) ovf = 1; }
            u32 nent = std::min(fin ? (u32)fin_so[2] : d2h1(d.norder), d.xcap);
            if (ovf) {
                /* table full: clear and retry with half the curves */
                CK(cudaMemset(d.overflow, 0, sizeof(int)));
                k_clear<<<nblk(nent, 256), 256>>>(d, nent);
                SYNC();
                if (getenv("CGV_PROF") && in->verbose) fprintf(stderr, "    retry: layer %d batch at %u, %u curves, mask 2^%d, entries %u\n", dd, c0, nb, __builtin_ctzll(d.xmask + 1), nent);
                est_nu *= 2; /* remember that curves here are bigger than estimated */
                if (d.xmask < full_mask) {
                    /* estimate was short: same curves, 4x the table */
                    force_mask = std::min<u64>((d.xmask + 1) * 4, full_mask + 1);
                    bmax = nb;
                } else {
                    if (nb == 1) {
                        size_t f = 0, t = 0; CK(cudaMemGetInfo(&f, &t));
                        u64 n = (full_mask + 1) * 2;
                        if (n > xcap_max || (n - (full_mask + 1)) * (per + (u64)(D + 2) * 2) > f) { fprintf(stderr, "cgv gpu: a single curve overflows the exp table (%llu slots, no memory to grow)\n", (unsigned long long)(full_mask + 1)); return 1; }
                        free_exp(); alloc_exp(n); alloc_lvl((u64)(D + 2) * d.xcap);
                        full_mask = d.xmask; full_cap = d.xcap; xfloor = n;
                        if (in->verbose) fprintf(stderr, "  gpu: exp table grown to %llu entries (one curve needs it)\n", (unsigned long long)n);
                        bmax = 1; force_mask = full_mask + 1;
                    } else {
                        bmax = std::max<u32>(1, nb / 2);
                        force_mask = full_mask + 1;
                    }
                }
                nretry++;
                tt[T_WASTE] += wall() - tbatch;
                continue;
            }
            tt[T_BUILDL] += tbl; tt[T_LEVEL] += tlev; tt[T_BUCKET] += tbuck;
            tb = wall();
            /* emit in chunks; each entry adds at most one new point, so a chunk of size K
             * cannot overflow once the hard cap has room for K more */
            {
                const u32 K = 1u << 20;
                for (u32 lo = 0; lo < nent; lo += K) {
                    u32 hi = std::min(nent, lo + K);
                    /* earlier emits (the previous batch's, overlapped with this batch's levels, and this
                     * batch's earlier chunks) must be done so that the I count is final */
                    CK(cudaStreamSynchronize(s_emit));
                    grow_I(hi - lo, lo == 0 && icnt_ok ? (i64)*h_icnt : -1);
                    if (lo == 0) { CK(cudaEventRecord(ev_lev, 0)); CK(cudaStreamWaitEvent(s_emit, ev_lev, 0)); }
                    k_emit_clear_range<<<nblk(hi - lo, 256), 256, 0, s_emit>>>(d, lo, hi, d_eovf); /* emit + clear in one pass */
                }
                if (h_icnt) { CK(cudaMemcpyAsync(h_icnt, d.icount, sizeof(u32), cudaMemcpyDeviceToHost, s_emit)); icnt_ok = 1; }
                CK(cudaEventRecord(ev_emit[cur], s_emit)); ev_used[cur] = 1;
                CK(cudaGetLastError());
            }
            /* not waited for: the next batch runs on the other buffer while this emit runs */
            if (dbl) cur ^= 1;
            tt[T_EMIT] += wall() - tb;
            nbatches++;
            { double sumnu = 0; for (u32 k = 0; k < nb; k++) sumnu += hnl[c0 + k]; /* this batch's curves */
              est_nu = std::max((double)nent / (sumnu + nb), 0.7 * est_nu); } /* slowly decaying max */
            c0 += nb;
            layEnt += nent; layB++;
            force_mask = 0; bmax = nc;
        }
        /* layer end: wait for the last emit, back to buffer 0, check the emit failure flag */
        SYNC();
        cur = 0; d.xtab = xtb[0]; d.order = ordb[0]; ev_used[0] = ev_used[1] = 0; icnt_ok = 0;
        { int ov = d2h1(d_eovf); if (ov) { fprintf(stderr, "cgv gpu: overflow flag %d after emit (layer %d)\n", ov, dd); return 1; } }
        if (getenv("CGV_PROF") && in->verbose)
            fprintf(stderr, "    gpu layer %d: curves %u batches %d entries %.3g W %.3g  level %.3fs layer %.3fs  (%.2g W/s)\n", dd, nc, layB,
                    (double)layEnt, (double)layW, tt[T_LEVEL] - layLev0, wall() - layT0, layW / (tt[T_LEVEL] - layLev0 + 1e-9));
        dd = dhi;
    }
    SYNC();
    cudaFree(d_next); cudaFree(d_sl); cudaFree(d_bsum); cudaFree(d_eovf);
    if (xtb[1]) cudaFree(xtb[1]);
    if (ordb[1]) cudaFree(ordb[1]);
    if (h_icnt) HOST_FREE(h_icnt);
    cudaStreamDestroy(s_emit); cudaEventDestroy(ev_emit[0]); cudaEventDestroy(ev_emit[1]); cudaEventDestroy(ev_lev);
    if (test_hold) cudaFree(test_hold);
    cudaFree(d_pat); cudaFree(d_bcur); cudaFree(d_bofs); cudaFree(d_bcnl);
    u32 npts = d2h1(d.icount);
    /* only the I table is needed from here on: release everything else before the output arrays
     * (52-72 B per point) are allocated, or a deep run that fits can still fail at the very end */
    SYNC();
    { void *p[] = {d.xtab, d.order, d.lvlpool, d.L, d.Lidx, d.work, d.woff, d.scan_tmp, d.ckey, d.cdeg, d.cs, d.cnu,
                   d.coff, d.cnl, d.ukey, d.udeg, d.ualpha, d.INV, d.DEGC};
      for (void *q : p) cudaFree(q); }
    d.xtab = NULL; d.order = NULL; d.lvlpool = NULL; d.L = NULL; d.Lidx = NULL; d.work = NULL; d.woff = NULL;
    d.scan_tmp = NULL; d.ckey = NULL; d.cdeg = NULL; d.cs = NULL; d.cnu = NULL; d.coff = NULL; d.cnl = NULL;
    d.ukey = NULL; d.udeg = NULL; d.ualpha = NULL; d.INV = NULL; d.DEGC = NULL;
    u128 *okey = dalloc<u128>(npts); int *odeg = dalloc<int>(npts); Co *oA = dalloc<Co>(npts); u32 *on = dalloc<u32>(1);
    CK(cudaMemset(on, 0, sizeof(u32)));
    k_collect<<<nblk(icap, 256), 256>>>(d, okey, odeg, oA, on);
    SYNC();
    out->n = (int)d2h1(on);
    out->key = malloc((size_t)out->n * sizeof(u128));
    out->deg = (int *)malloc((size_t)out->n * sizeof(int));
    out->A = malloc((size_t)out->n * sizeof(Co));
    CK(cudaMemcpy(out->key, okey, (size_t)out->n * sizeof(u128), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(out->deg, odeg, (size_t)out->n * sizeof(int), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(out->A, oA, (size_t)out->n * sizeof(Co), cudaMemcpyDeviceToHost));
    out->napplied = napplied;
    if (in->verbose) {
        fprintf(stderr, "  gpu: %lld curves applied in %d batches (%d retries), %d points\n", (long long)napplied, nbatches, nretry, out->n);
        fprintf(stderr, "  gpu: %.3g scatter threads (%.3g/s in level)\n", (double)totW, totW / (tt[T_LEVEL] + 1e-9));
        fprintf(stderr, "  gpu time:");
        for (int i = 0; i < T_N; i++) fprintf(stderr, " %s %.2fs", TNAME[i], tt[i]);
        fprintf(stderr, "\n");
    }
    cudaFree(okey); cudaFree(odeg); cudaFree(oA); cudaFree(on);
    cudaFree(d.ukey); cudaFree(d.udeg); cudaFree(d.ualpha); cudaFree(d.INV); cudaFree(d.DEGC); cudaFree(d.overflow);
    cudaFree(d.itab); cudaFree(d.icount); cudaFree(d.ckey); cudaFree(d.cdeg); cudaFree(d.cs); cudaFree(d.cnu); cudaFree(d.coff);
    cudaFree(d.ncurves); cudaFree(d.xtab); cudaFree(d.order); cudaFree(d.lvlpool); cudaFree(d.norder); cudaFree(d.lvlcount);
    cudaFree(cursor); cudaFree(d.L); cudaFree(d.work); cudaFree(d.woff); cudaFree(d.scan_tmp); cudaFree(d.Lidx); cudaFree(d.lvlcnt); cudaFree(d.cnl); cudaFree(d.umask);
    return 0;
}
