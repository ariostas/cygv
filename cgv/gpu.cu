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
#include <cstddef>
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
#ifdef CGV_WIDE
#define GPUNS CGV_CAT(CGV_CAT(cgv_gpu_nl, NL), CGV_CAT(_w, CGV_WIDE))
#else
#define GPUNS CGV_CAT(cgv_gpu_nl, NL)
#endif
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
__constant__ int c_kbv, c_kbw[64];   /* per-coordinate key field widths (c_kbv = 1), else c_keybits each */
/* Global keys (I table, curves) are GKey: u128, or KW words in the CGV_WIDE = 256/512 builds. The exp tables always key
 * by u128 offset codes: the packed offset itself in the 128-bit build; in the wide builds the host's level code (the
 * 64/128-bit packed form lt_setup chose, linear, injective within a level), decoded to a key only at emit. */
#ifdef CGV_WIDE
#define GKW (CGV_WIDE / 64)
struct GKey { u64 w[GKW]; };
__device__ __forceinline__ GKey gk_add(const GKey &a, const GKey &b) {
    GKey r; u64 c = 0;
    for (int i = 0; i < GKW; i++) { u64 s = a.w[i] + c; u64 c1 = s < c; r.w[i] = s + b.w[i]; c = c1 + (r.w[i] < s); }
    return r;
}
/* k + K v (neg: k - K v) */
__device__ __forceinline__ GKey gk_madd(const GKey &k, const GKey &K, u64 v, int neg) {
    GKey p; u64 c = 0;
    for (int i = 0; i < GKW; i++) { u64 lo = K.w[i] * v, hi = __umul64hi(K.w[i], v); lo += c; hi += lo < c; p.w[i] = lo; c = hi; }
    GKey r; u64 b = 0;
    if (!neg) return gk_add(k, p);
    for (int i = 0; i < GKW; i++) { u64 x = k.w[i] - p.w[i], b1 = k.w[i] < p.w[i]; r.w[i] = x - b; b = b1 | (x < b); }
    return r;
}
__device__ __forceinline__ int gk_eq(const GKey &a, const GKey &b) { u64 x = 0; for (int i = 0; i < GKW; i++) x |= a.w[i] ^ b.w[i]; return !x; }
__device__ __forceinline__ int gk_eq_vol(const GKey *a, const GKey &b) { const volatile u64 *p = a->w; u64 x = 0; for (int i = 0; i < GKW; i++) x |= p[i] ^ b.w[i]; return !x; }
/* level-code decoder (gv.c lt_setup): key = base + sum_f v_f K_f (+ r K_t0), v_f = bits [off, off + b) of the code minus
 * plo (mod 2^64 or 2^128), r = (e - wlo - sum_f w_f v_f) / w0 */
__constant__ int c_d128, c_dn, c_dt0, c_dw0, c_doff[64], c_db[64], c_dw[64]; __constant__ i64 c_dwlo; __constant__ u128 c_dplo;
__constant__ GKey c_dK[64], c_dbase, c_dKt0;
__device__ __forceinline__ u128 ok_add(u128 a, u128 b) { u128 s = a + b; return c_d128 ? s : (u128)(u64)s; }
__device__ GKey dec_off(u128 o, int e) {
    u128 y = c_d128 ? o - c_dplo : (u128)((u64)o - (u64)c_dplo);
    GKey k = c_dbase; u64 ur = (u64)(e - c_dwlo);
    for (int f = 0; f < c_dn; f++) {
        u64 v = (u64)(y >> c_doff[f]) & (((u64)1 << c_db[f]) - 1);
        k = gk_madd(k, c_dK[f], v, 0); ur -= (u64)(i64)c_dw[f] * v;
    }
    if (c_dt0) {
        i64 r = (i64)ur; r = c_dw0 == 1 ? r : c_dw0 == -1 ? -r : r / c_dw0;
        k = gk_madd(k, c_dKt0, r >= 0 ? (u64)r : -(u64)r, r < 0);
    }
    return k;
}
#define OK_ADD(a, b) ok_add(a, b)
#define EMIT_KEY(kc, o, e) gk_add(kc, dec_off(o, e))
#define GK_EQ(a, b) gk_eq(a, b)
#define GK_EQ_VOL(p, b) gk_eq_vol(p, b)
#define XDEG(e) && (e)->deg == deg   /* level codes may coincide across levels (the degree fixes one coordinate) */
#define XHDEG ^ ((u64)(u32)deg * 0x9E3779B97F4A7C15ull) /* so the hash separates them too */
#else
typedef u128 GKey;
#define OK_ADD(a, b) ((a) + (b))
#define EMIT_KEY(kc, o, e) ((kc) + (o))
#define GK_EQ(a, b) ((a) == (b))
#define GK_EQ_VOL(p, b) (*(volatile u128 *)(p) == (b))
#define XDEG(e)
#define XHDEG
#endif

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
/* acc[l] += x[l] for all lanes: the NL atomics are issued back to back (one round trip, not NL),
 * the rare wrap corrections afterwards */
/* the wrap corrections (frequent: acc passes 2^64 every few adds of ~62-bit values, so most warps have a lane
 * that needs one) are likewise issued together, so a warp waits one extra round trip, not one per lane */
__device__ __forceinline__ void acc_addv(u64 *acc, const u64 *x) {
    u64 old[NL], cx[NL];
    for (int l = 0; l < NL; l++) old[l] = atomicAdd((unsigned long long *)&acc[l], (unsigned long long)x[l]);
    int any = 0;
    for (int l = 0; l < NL; l++) { cx[l] = old[l] + x[l] < old[l] ? c_r64[l] : 0; any |= cx[l] != 0; }
#ifdef CGV_HIP
    (void)any; /* AMD: one correction at a time (batching them was slower there) */
    for (int l = 0; l < NL; l++) if (cx[l]) acc_add(&acc[l], cx[l], l);
#else
    while (any) {
#ifndef __CUDA_ARCH__
        for (int l = 0; l < NL; l++) old[l] = cx[l] ? atomicAdd((unsigned long long *)&acc[l], (unsigned long long)cx[l]) : 0;
#else
        /* predicated atomics, so that nvcc issues them back to back instead of one branch (and wait) per lane */
        for (int l = 0; l < NL; l++) {
            old[l] = 0;
            asm volatile("{ .reg .pred p; setp.ne.u64 p, %2, 0; @p atom.global.add.u64 %0, [%1], %2; }"
                         : "+l"(old[l]) : "l"(&acc[l]), "l"(cx[l]) : "memory");
        }
#endif
        any = 0;
        for (int l = 0; l < NL; l++) { cx[l] = cx[l] && old[l] + cx[l] < old[l] ? c_r64[l] : 0; any |= cx[l] != 0; }
    }
#endif
}

/* ---------------- keys ---------------- */
__device__ __forceinline__ u64 hash128(u128 k) {
    u64 x = (u64)k ^ ((u64)(k >> 64) * 0x9E3779B97F4A7C15ull);
    x ^= x >> 31; x *= 0xD6E8FEB86659FD93ull; x ^= x >> 32;
    return x;
}
#ifdef CGV_WIDE
__device__ __forceinline__ u64 gk_hash(const GKey &k) {
    const u64 C[8] = {1, 0x9E3779B97F4A7C15ull, 0xC2B2AE3D27D4EB4Full, 0x165667B19E3779F9ull,
                      0x27D4EB2F165667C5ull, 0x94D049BB133111EBull, 0xBF58476D1CE4E5B9ull, 0xFF51AFD7ED558CCDull};
    u64 x = 0;
    for (int i = 0; i < GKW; i++) x ^= k.w[i] * C[i];
    x ^= x >> 31; x *= 0xD6E8FEB86659FD93ull; x ^= x >> 32;
    return x;
}
__device__ void unpack_d(GKey k, int *v) {   /* (as gv.c's wide unpack) */
    for (int t = 0; t < c_h11; t++) {
        int w = c_kbv ? c_kbw[t] : c_keybits;
        i64 f = (i64)(k.w[0] & (((u64)1 << w) - 1));
        if (f >= ((i64)1 << (w - 1))) f -= (i64)1 << w;
        v[t] = (int)f;
        u64 b = 0;   /* k = (k - f) >> w */
        GKey s;
        for (int i = 0; i < GKW; i++) { u64 x = (u64)(i ? (f < 0 ? -1 : 0) : f), y = k.w[i] - x, b1 = k.w[i] < x; s.w[i] = y - b; b = b1 | (y < b); }
        for (int i = 0; i < GKW; i++) k.w[i] = s.w[i] >> w | (i + 1 < GKW ? s.w[i + 1] << (64 - w) : 0);
    }
}
#else
#define gk_hash hash128
__device__ void unpack_d(u128 k, int *v) {
    if (c_kbv) {
        for (int t = 0; t < c_h11; t++) {
            int w = c_kbw[t];
            i64 f = (i64)(u64)(k & (((u128)1 << w) - 1));
            if (f >= ((i64)1 << (w - 1))) f -= (i64)1 << w;
            v[t] = (int)f;
            k = (k - (u128)(i128)f) >> w;
        }
        return;
    }
    u128 mask = ((u128)1 << c_keybits) - 1;
    for (int t = 0; t < c_h11; t++) {
        i64 f = (i64)(u64)(k & mask);
        if (f >= ((i64)1 << (c_keybits - 1))) f -= (i64)1 << c_keybits;
        v[t] = (int)f;
        k = (k - (u128)(i128)f) >> c_keybits;
    }
}
#endif

/* ---------------- tables ---------------- */
/* state: 0 empty, 1 being written, 2 ready */
/* acc accumulates I; once degree deg has been extracted it holds A instead
 * (nothing adds to a point after its degree is extracted) */
struct IEnt { GKey key; int deg; u32 state; u64 acc[NL]; };
/* loff: the curve's L pool offset, set at finalize (saves the scatter a dependent load of coff) */
struct XEnt { u128 key; u32 curve; int deg; u32 state; u32 loff; u64 acc[NL]; };
/* ld_rlx_head reads (curve, deg) and (state, loff) as the entry's 3rd and 4th u64; its 32-byte load also needs
 * 32-byte entries (XF_SCATTER) in a 256-byte aligned table (cudaMalloc) */
static_assert(offsetof(XEnt, curve) == 16 && offsetof(XEnt, state) == 24, "XEnt layout");

/* load-acquire (gpu scope): later loads of the entry see what the writer
 * published before its release of the state flag */
__device__ __forceinline__ u32 ld_state(const u32 *s) { return gpu_ld_acquire(s); }
__device__ __forceinline__ void st_release(u32 *s, u32 v) { gpu_st_release(s, v); }
/* exp-table probe modes (see x_find) */
enum { XF_GPU, XF_EPOCH, XF_CTA };
#if !defined(CGV_HIP) && defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 1000
/* the relaxed head load needs 32-byte entries (NL 3, 4); elsewhere (AMD, NL 2, older NVIDIA) it did
 * not pay off, and the acquire path is kept */
#define XF_SCATTER (sizeof(XEnt) % 32 == 0 ? XF_EPOCH : XF_GPU)
#else
#define XF_SCATTER XF_GPU
#endif
#ifdef CGV_HIP
#define XF_SMALL XF_GPU
#else
#define XF_SMALL XF_CTA
#endif

/* one relaxed (strong, L1-bypassing) 32-byte load of an exp entry's key, (curve, deg), (state, loff) */
__device__ __forceinline__ void ld_rlx_head(const XEnt *e, u128 *k, u32 *c, u32 *st) {
#if !defined(CGV_HIP) && defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 1000
    u64 k0, k1, w2, w3;
    asm volatile("ld.relaxed.gpu.global.v4.u64 {%0, %1, %2, %3}, [%4];" : "=l"(k0), "=l"(k1), "=l"(w2), "=l"(w3) : "l"(&e->key) : "memory");
    *k = (u128)k1 << 64 | k0; *c = (u32)w2; *st = (u32)w3;
#elif defined(CGV_HIP)
    __builtin_trap(); /* not used here (XF_SCATTER) */
#else
    __trap();
#endif
}
__device__ __forceinline__ u32 ld_acquire_cta(const u32 *s) {
#ifdef CGV_HIP
    return __hip_atomic_load(s, __ATOMIC_ACQUIRE, __HIP_MEMORY_SCOPE_WORKGROUP);
#else
    u32 v;
    asm volatile("ld.acquire.cta.global.u32 %0, [%1];" : "=r"(v) : "l"(s) : "memory");
    return v;
#endif
}


/* find or insert point k in the I table; returns position or -1 if full */
__device__ i64 i_find(IEnt *tab, u64 mask, GKey k, int deg, u32 *count, u32 cap) {
    u64 h = gk_hash(k) & mask;
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
        if (GK_EQ_VOL(&e->key, k)) return (i64)h;
        h = (h + 1) & mask;
    }
    return -1;
}
/* find or insert (curve, k) in the exp table; new entries are appended to order.
 * Slot states carry the batch epoch B (1 <= B < 2^24): B<<8|1 being written, B<<8|w ready (w >= 2),
 * any other value is empty, so a batch's table needs no clearing afterwards (an insert resets acc).
 * w names the writer: e + 3 for the k_scatter launch of level e in XF_EPOCH mode, else 2.
 * XF_GPU: acquire every probe. XF_CTA: the calling block is the only writer during the launch
 * (k_small_levels), so a block-scope acquire suffices (no L1 invalidation). XF_EPOCH (k_scatter of
 * level e <= 252, w = e + 3): an entry that is ready but not B<<8|w was completed before this launch
 * and is immutable during it, so its key and header, from one relaxed load (one L2 request instead
 * of three, no acquire), can be compared directly; only entries inserted in this launch take the
 * acquire. (Several launches of one level, as on the bucketed path, only make the acquire path more common.) */
/* y0 (optional): a new entry starts with acc = y0 instead of 0, and *ins is set, so the inserting thread
 * needs no atomic add of its own. The insert issues its stores and both counter atomics before the release,
 * so the release fence overlaps the atomics' round trip instead of adding its own. NVIDIA only (XF_Y0):
 * on AMD it was slower, and the insert there keeps the earlier order and ignores y0. */
#ifdef CGV_HIP
#define XF_Y0 0
#else
#define XF_Y0 1
#endif
template <int M>
__device__ i64 x_find(XEnt *tab, u64 mask, u32 c, u128 k, int deg, u32 *norder, u32 *order, u32 cap,
                      int direct, u32 *lvlcnt, u32 *lvlpool, u32 stride, u32 B, u32 w,
                      const u64 *y0 = NULL, int *ins = NULL) {
    const u32 S1 = B << 8 | 1, SW = B << 8 | w;
#define XF_READY(st) (((st) ^ (B << 8)) < 256u && ((st) & 255u) >= 2u)
    u64 h = (hash128(k) ^ ((u64)c * 0xA24BAED4963EE407ull) XHDEG) & mask;
    for (u64 probes = 0; probes <= mask; probes++) {
        XEnt *e = &tab[h];
        u32 st;
        if (M == XF_EPOCH) {
            u128 ek; u32 ec;
            ld_rlx_head(e, &ek, &ec, &st);
            if (__builtin_expect(XF_READY(st) && st != SW, 1)) {
                if (ec == c && ek == k XDEG(e)) return (i64)h;
                h = (h + 1) & mask;
                continue;
            }
            if (st == SW) st = ld_state(&e->state);
        } else st = M == XF_CTA ? ld_acquire_cta(&e->state) : ld_state(&e->state);
        if (__builtin_expect(XF_READY(st), 1)) {   /* common case: occupied, compare and move on (found by ShinkaEvolve, +3-9%) */
            if (e->curve == c && *(u128 *)&e->key == k XDEG(e)) return (i64)h;
            h = (h + 1) & mask;
            continue;
        }
        if (st != S1) {
            if (atomicCAS(&e->state, st, S1) == st) {
#ifdef CGV_HIP
                /* AMD: the earlier insert (counter first, acc = 0, no y0): the reordered one was slower there */
                (void)y0; (void)ins;
                u32 o = atomicAdd(norder, 1u);
                if (o >= cap) { atomicExch(&e->state, 0u); return -1; } /* overflow: caller retries with a smaller batch */
                e->key = k; e->curve = c; e->deg = deg;
                for (int l = 0; l < NL; l++) e->acc[l] = 0;
                st_release(&e->state, SW);
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
#else
                e->key = k; e->curve = c; e->deg = deg;
                for (int l = 0; l < NL; l++) e->acc[l] = y0 ? y0[l] : 0;
                u32 o = atomicAdd(norder, 1u), q = 0;
                int ok = o < cap;
                if (direct) {
                    /* warp-aggregated append: lanes inserting at the same level share one atomic;
                       overflowing lanes match on -1 and reserve nothing */
                    unsigned act = warp_active();
                    unsigned grp = warp_match_any(act, ok ? deg : -1);
                    int lane = threadIdx.x & 31, leader = __ffs(grp) - 1;
                    u32 q0 = 0;
                    if (ok && lane == leader) q0 = atomicAdd(&lvlcnt[deg], (u32)__popc(grp));
                    q0 = warp_shfl(grp, q0, leader);
                    q = q0 + __popc(grp & ((1u << lane) - 1));
                }
                if (!ok) { atomicExch(&e->state, 0u); return -1; } /* overflow: caller retries with a smaller batch */
                st_release(&e->state, SW);
                if (ins) *ins = 1;
                order[o] = (u32)h;
                if (direct) {
                    if (q >= stride) return -1; /* level list full: treated as overflow */
                    lvlpool[(u64)deg * stride + q] = (u32)h;
                }
                return (i64)h;
#endif
            }
            st = M == XF_CTA ? ld_acquire_cta(&e->state) : ld_state(&e->state);
        }
#ifdef CGV_HIP
        if (st == S1) { probes--; continue; } /* being written: look again next iteration. AMD waves give diverged
                                                lanes no forward-progress guarantee, so no spin inside a branch */
#else
        while (st == S1) st = M == XF_CTA ? ld_acquire_cta(&e->state) : ld_state(&e->state);
#endif
        if (!XF_READY(st)) { probes--; continue; } /* slot was released (overflow) or taken by a CAS since: look at it again */
        if (e->curve == c && *(u128 *)&e->key == k XDEG(e)) return (i64)h;
        h = (h + 1) & mask;
    }
#undef XF_READY
    return -1;
}

/* ---------------- device state ---------------- */
struct Dev {
    /* union U */
    int NU; u128 *ukey; int *udeg; Co *ualpha; /* [NU][h11]; ukey: offset codes (see GKey) */
    Co *INV, *DEGC; int tabn;
    /* I table */
    IEnt *itab; u64 imask; u32 *icount; u32 icap;
    /* curves of the current layer */
    GKey *ckey; int *cdeg; Co *cs; int *cnu; u64 *coff; u32 *ncurves; u32 ccap;
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
    u32 xep;                   /* exp table epoch of the current batch (see x_find) */
    /* class-pattern build (patmode): L = every support term, no zero tests, no values; each scatter
     * pair is recorded as (target slot, source slot << 32 | L index) */
    int patmode; u32 *npair; u32 pcap; u32 *pk; u64 *pv;
};
/* record one pattern pair (warp-aggregated append) */
__device__ __forceinline__ void pat_rec(Dev &d, i64 h, u32 src, u64 j) {
    unsigned act = warp_active();
    int lane = threadIdx.x & 31, leader = __ffs(act) - 1;
    u32 b0 = 0;
    if (lane == leader) b0 = atomicAdd(d.npair, (u32)__popc(act));
    b0 = warp_shfl(act, b0, leader);
    u32 x = b0 + __popc(act & ((1u << lane) - 1));
    if (x < d.pcap) { d.pk[x] = (u32)h; d.pv[x] = ((u64)src << 32) | j; } else *d.overflow = 7;
}

/* initial I: insert all points */
__global__ void k_init_I(Dev d, int n, const GKey *key, const int *deg, const Co *val) {
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
        if (u < nu && (d.umask[u] & pat) && d.patmode) { nz = 1; v.v[0] = 1; }
        else if (u < nu && (d.umask[u] & pat)) {
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
    /* breakpoints after the curve's L terms: Li[nu + g] = # terms with udeg <= g, g <= maxdeg - deg(c),
     * so a level's work count is one load instead of a search (terms are in union = degree order) */
    for (int g = threadIdx.x; g <= c_maxdeg - d.cdeg[c]; g += blockDim.x) {
        int lo = 0, hi = base;
        while (lo < hi) { int mid = (lo + hi) / 2; if (d.udeg[Li[mid]] <= g) lo = mid + 1; else hi = mid; }
        Li[nu + g] = lo;
    }
}

/* seed: entry (c, 0) with value 1 for each curve of the batch */
__global__ void k_seed(Dev d, u32 c0, u32 nc) {
    u32 i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nc) return;
    u32 c = c0 + i;
    i64 h = x_find<XF_GPU>(d.xtab, d.xmask, c, 0, 0, d.norder, d.order, d.xcap, d.direct, d.lvlcnt, d.lvlpool, d.lvlstride, d.xep, 2);
    if (h < 0) { *d.overflow = 3; return; }
    for (int l = 0; l < NL; l++) d.xtab[h].acc[l] = c_one[l];
}

/* the same L build spread over a 2D grid (curve, chunk of BL_CH u) so that a batch of few curves still fills
 * the GPU: k_buildL_cnt counts each chunk's nonzero terms, k_buildL_put recomputes them and writes them after
 * the earlier chunks' counts (cnt[x * gridDim.y + y]), with the breakpoints; the (curve, 0) block also seeds */
#define BL_CH 1024
__device__ __forceinline__ void bl_setup(Dev d, u32 c, Co *cm, int *nzt, int *nnz, u64 *pat) {
    if (threadIdx.x == 0) {
        int C[64]; unpack_d(d.ckey[c], C);
        int n = 0; u64 m = 0;
        for (int t = 0; t < c_h11; t++) if (C[t]) { for (int l = 0; l < NL; l++) cm[n].v[l] = mint(C[t], l); nzt[n++] = t; m |= (u64)1 << t; }
        *nnz = n; *pat = m;
    }
    __syncthreads();
}
__device__ __forceinline__ int bl_term(Dev d, int u, int nu, const Co *cm, const int *nzt, int nnz, u64 pat, Co &v) {
    int nz = 0;
    for (int l = 0; l < NL; l++) v.v[l] = 0;
    if (u < nu && (d.umask[u] & pat)) {
        const Co *al = d.ualpha + (size_t)u * c_h11; /* ualpha already carries the factor deg(u) */
        for (int q = 0; q < nnz; q++) {
            Co a = al[nzt[q]];
            for (int l = 0; l < NL; l++) v.v[l] = madd(v.v[l], mmul(cm[q].v[l], a.v[l], l), l);
        }
        for (int l = 0; l < NL; l++) nz |= v.v[l] != 0;
    }
    return nz;
}
__global__ void __launch_bounds__(256) k_buildL_cnt(Dev d, u32 c0, u32 *cnt) {
    __shared__ Co cm[64]; __shared__ int nzt[64], nnz; __shared__ u64 pat;
    u32 c = c0 + blockIdx.x;
    int nu = d.cnu[c], u0 = blockIdx.y * BL_CH;
    if (u0 >= nu) return;
    bl_setup(d, c, cm, nzt, &nnz, &pat);
    int n = 0;
    for (int i = 0; i < BL_CH; i += 256) { Co v; n += __syncthreads_count(bl_term(d, u0 + i + threadIdx.x, nu, cm, nzt, nnz, pat, v)); }
    if (threadIdx.x == 0) cnt[(u64)blockIdx.x * gridDim.y + blockIdx.y] = (u32)n;
}
__global__ void __launch_bounds__(256) k_buildL_put(Dev d, u32 c0, const u32 *cnt) {
    typedef cub::BlockScan<int, 256> BS;
    __shared__ typename BS::TempStorage tmp;
    __shared__ Co cm[64]; __shared__ int nzt[64], nnz; __shared__ u64 pat; __shared__ int base;
    u32 c = c0 + blockIdx.x;
    int nu = d.cnu[c], u0 = blockIdx.y * BL_CH, T = c_maxdeg - d.cdeg[c];
    Co *Lc = d.L + d.coff[c];
    int *Li = d.Lidx + d.coff[c];
    if (blockIdx.y == 0 && threadIdx.x == 0) {
        i64 h = x_find<XF_GPU>(d.xtab, d.xmask, c, 0, 0, d.norder, d.order, d.xcap, d.direct, d.lvlcnt, d.lvlpool, d.lvlstride, d.xep, 2);
        if (h < 0) *d.overflow = 3;
        else for (int l = 0; l < NL; l++) d.xtab[h].acc[l] = c_one[l];
        for (int g = 0, g1 = nu ? d.udeg[0] : T + 1; g < g1 && g <= T; g++) Li[nu + g] = 0; /* below the first term's degree */
        if (!nu) d.cnl[c] = 0;
    }
    if (u0 >= nu) return;
    bl_setup(d, c, cm, nzt, &nnz, &pat);
    const u32 *cc = cnt + (u64)blockIdx.x * gridDim.y;
    int s = 0, pre;
    for (u32 j = threadIdx.x; j < blockIdx.y; j += 256) s += (int)cc[j];
    BS(tmp).ExclusiveSum(s, pre, s); /* s = sum of the earlier chunks' counts */
    if (threadIdx.x == 0) base = s;
    __syncthreads();
    for (int i = 0; i < BL_CH; i += 256) {
        Co v;
        int u = u0 + i + threadIdx.x;
        int nz = bl_term(d, u, nu, cm, nzt, nnz, pat, v);
        int pos, tot;
        BS(tmp).ExclusiveSum(nz, pos, tot);
        int b = base;
        if (nz) { Lc[b + pos] = v; Li[b + pos] = u; }
        /* breakpoints Li[nu + g] = # terms with udeg <= g (g <= T): set by the last u of each degree */
        if (u < nu) {
            int g0 = d.udeg[u], g1 = u + 1 < nu ? d.udeg[u + 1] : T + 1;
            for (int g = g0; g < g1; g++) Li[nu + g] = b + pos + nz;
        }
        __syncthreads();
        if (threadIdx.x == 0) base += tot;
        __syncthreads();
    }
    if (threadIdx.x == 0 && u0 + BL_CH >= nu) d.cnl[c] = base; /* last chunk */
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
    if (!nz && !d.patmode) { d.work[i] = 0; return; }
    u32 c = x->curve;
    int lim = c_maxdeg - d.cdeg[c] - e;
    u64 co = d.coff[c];
    x->loff = (u32)co;
    d.work[i] = lim < 0 ? 0 : (u64)d.Lidx[co + d.cnu[c] + lim];
}
/* level e, step 2: one thread per (entry, L term); called by whole warps (t - lane is warp-uniform) */
__device__ __forceinline__ void scatter_t(Dev &d, int e, const u32 *ent, u32 n, u64 W, u64 t) {
    /* entry of thread t = last i with woff[i] <= t. Consecutive t share or
     * neighbor entries, so the warp searches the whole range for its first
     * and last t, and every lane then searches only between those. */
    int lane = threadIdx.x & 31;
    u64 tw0 = t - lane, tw1 = tw0 + 31 < W ? tw0 + 31 : W - 1;
    if (tw0 >= W) return; /* whole warp */
    /* 16-ary searches, lanes 0-15 for tw0 and 16-31 for tw1: a step loads 15 pivots at once, so
     * ~log16(n) dependent loads instead of log2(n). Invariant woff[lo] <= tt < woff[hi]. */
    u32 lo0, lo1;
    {
        int q = lane & 15, sh = lane & 16;
        u64 tt = sh ? tw1 : tw0;
        u32 lo = 0, hi = n;
        while (warp_any(0xffffffffu, hi - lo > 1)) {
            u32 len = hi - lo;
            u32 p = lo + (u32)((u64)len * q / 16);
            int pr = q && len > 1 && d.woff[p] <= tt;
            int k = __popc((warp_ballot(pr) >> sh) & 0xffffu);
            if (len > 1) {
                u32 nlo = lo + (u32)((u64)len * k / 16);
                hi = k < 15 ? lo + (u32)((u64)len * (k + 1) / 16) : hi;
                lo = nlo;
            }
        }
        lo0 = warp_shfl(0xffffffffu, lo, 0); lo1 = warp_shfl(0xffffffffu, lo, 16);
    }
    u32 a = lo0, b = lo1;
    if (t >= W) return;
    u32 lo = a, hi = b + 1;
    while (hi - lo > 1) { u32 mid = (lo + hi) / 2; if (d.woff[mid] <= t) lo = mid; else hi = mid; }
    u32 i = lo;

    u64 j = t - d.woff[i];
    XEnt *x = &d.xtab[ent[i]];
    u32 c = x->curve;
    Co la = d.L[x->loff + j];
    int u = d.Lidx[x->loff + j];
    u64 y[NL];
    for (int l = 0; l < NL; l++) y[l] = mmul(la.v[l], x->acc[l], l);
    int ins = 0;
    i64 h = (XF_SCATTER == XF_EPOCH && e <= 252 && !d.patmode)   /* (pattern builds: the generic writer) */
        ? x_find<XF_SCATTER>(d.xtab, d.xmask, c, OK_ADD(x->key, d.ukey[u]), e + d.udeg[u], d.norder, d.order, d.xcap, d.direct, d.lvlcnt, d.lvlpool, d.lvlstride, d.xep, e + 3, XF_Y0 ? y : NULL, &ins)
        : x_find<XF_GPU>(d.xtab, d.xmask, c, OK_ADD(x->key, d.ukey[u]), e + d.udeg[u], d.norder, d.order, d.xcap, d.direct, d.lvlcnt, d.lvlpool, d.lvlstride, d.xep, 2, XF_Y0 && !d.patmode ? y : NULL, &ins);
    if (h < 0) { *d.overflow = 4; return; }
    if (d.patmode) { pat_rec(d, h, ent[i], j); return; }
    if (!ins) acc_addv(d.xtab[h].acc, y);
}
__global__ void k_scatter(Dev d, int e, const u32 *ent, u32 n, u64 W) {
    scatter_t(d, e, ent, n, W, (u64)blockIdx.x * blockDim.x + threadIdx.x);
}

/* Device-driven levels: once a batch reaches a big level, the host queues every remaining level as
 * k_dfin, k_dscan, k_doff, k_dscatter with no round trips; level sizes and work totals stay on the
 * device (fixed grids, grid-stride loops). dw: [0] W of the level, [1] sum of W, [2..] block sums. */
#define DL_G 512
__global__ void __launch_bounds__(256) k_dfin(Dev d, int e, u64 *dw) {
    typedef cub::BlockScan<u64, 256> BS;
    __shared__ typename BS::TempStorage tmp;
    __shared__ u64 carry;
    u32 n = *(volatile u32 *)&d.lvlcnt[e];
    if (*(volatile int *)d.overflow || !n) { if (blockIdx.x == 0 && threadIdx.x == 0) dw[0] = 0; return; }
    if (n > d.lvlstride) { if (blockIdx.x == 0 && threadIdx.x == 0) { *d.overflow = 6; dw[0] = 0; } return; }
    const u32 *ent = d.lvlpool + (u64)e * d.lvlstride;
    u32 ch = (n + gridDim.x - 1) / gridDim.x, c0 = blockIdx.x * ch, c1 = min(n, c0 + ch);
    if (threadIdx.x == 0) carry = 0;
    __syncthreads();
    for (u32 base = c0; base < c1; base += blockDim.x) {
        u32 i = base + threadIdx.x;
        u64 w = 0, pre, agg;
        if (i < c1) {
            XEnt *x = &d.xtab[ent[i]];
            int nz = 0;
            for (int l = 0; l < NL; l++) {
                u64 a = x->acc[l] % c_p[l];
                u64 f = e == 0 ? a : mmul(a, d.INV[e].v[l], l);
                x->acc[l] = f; nz |= f != 0;
            }
            if (nz || d.patmode) {
                u32 c = x->curve;
                int lim = c_maxdeg - d.cdeg[c] - e;
                u64 co = d.coff[c];
                x->loff = (u32)co;
                w = lim < 0 ? 0 : (u64)d.Lidx[co + d.cnu[c] + lim];
            }
        }
        BS(tmp).ExclusiveSum(w, pre, agg);
        if (i < c1) d.woff[i] = carry + pre;
        __syncthreads();
        if (threadIdx.x == 0) carry += agg;
        __syncthreads();
    }
    if (threadIdx.x == 0) dw[2 + blockIdx.x] = carry;
}
__global__ void __launch_bounds__(DL_G) k_dscan(Dev d, int e, u64 *dw) {
    typedef cub::BlockScan<u64, DL_G> BS;
    __shared__ typename BS::TempStorage tmp;
    u32 n = *(volatile u32 *)&d.lvlcnt[e];
    if (*(volatile int *)d.overflow || !n) return;   /* dw[0] = 0 from k_dfin */
    u64 v = dw[2 + threadIdx.x], pre, agg;
    BS(tmp).ExclusiveSum(v, pre, agg);
    dw[2 + threadIdx.x] = pre;
    if (threadIdx.x == 0) { dw[0] = agg; dw[1] += agg; d.woff[n] = agg; dw[DL_G + 2] = 0; }
}
__global__ void __launch_bounds__(256) k_doff(Dev d, int e, const u64 *dw) {
    u32 n = *(volatile u32 *)&d.lvlcnt[e];
    if (!dw[0] || blockIdx.x == 0) return;   /* block 0's chunk has offset 0 */
    u32 ch = (n + gridDim.x - 1) / gridDim.x, c0 = blockIdx.x * ch, c1 = min(n, c0 + ch);
    u64 off = dw[2 + blockIdx.x];
    for (u32 i = c0 + threadIdx.x; i < c1; i += blockDim.x) d.woff[i] += off;
}
/* resident blocks take 256-thread chunks in order from a counter (dw[DL_G + 2], zeroed by k_dscan), so the
 * active threads stay a narrow window of t as with one block per chunk */
__global__ void __launch_bounds__(256) k_dscatter(Dev d, int e, u64 *dw) {
    __shared__ u64 b_sh;
    u64 W = *(volatile u64 *)&dw[0];
    if (!W) return;
    u32 n = *(volatile u32 *)&d.lvlcnt[e];
    const u32 *ent = d.lvlpool + (u64)e * d.lvlstride;
    for (;;) {
        __syncthreads();
        if (threadIdx.x == 0) b_sh = (u64)atomicAdd((unsigned long long *)&dw[DL_G + 2], 1ull) * blockDim.x;
        __syncthreads();
        u64 b = b_sh;
        if (b >= W) return;
        scatter_t(d, e, ent, n, W, b + threadIdx.x);
    }
}
/* the batch's bookkeeping after device-driven levels, as k_small_levels reports it when done */
__global__ void k_dreport(Dev d, int T, u64 *dw, u64 *out) {
    out[0] = T + 1; out[4] = dw[1]; dw[1] = 0;
    out[1] = (u64)(u32)atomicAdd((unsigned *)d.overflow, 0u);
    out[2] = atomicAdd(d.norder, 0u);
    out[3] = atomicAdd(d.icount, 0u);
}

/* Runs of small levels, all on one block (no host round trips): for each level
 * e from e0 on, finalize its entries, prefix-sum their work, and scatter. Stops
 * (writing e and its work total to out[0..1]) at a level with more than nmax
 * entries (not finalized: out[2] = 0) or more than wmax work items (finalized,
 * woff ready: out[2] = 1), out[3] = its entry count; out[0] = T + 1 when all levels are done,
 * and then out[1..3] = the overflow flag, the number of exp entries, the I table count, so the
 * host learns the batch's bookkeeping from this one copy (no further round trips). */
#define SL_SH 4096
__global__ void __launch_bounds__(1024) k_small_levels(Dev d, int e0, int T, u32 nmax, u64 wmax, u64 *out) {
    typedef cub::BlockScan<u64, 1024> BS;
    __shared__ typename BS::TempStorage tmp;
    __shared__ u64 carry, Wtot;
    __shared__ u32 n_sh;
    __shared__ int ovf_sh;
    __shared__ u64 s_woff[SL_SH + 1]; /* copy of woff for the scatter's searches (levels with n <= SL_SH) */
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
            if (nz || d.patmode) {
                u32 c = x->curve;
                int lim = c_maxdeg - d.cdeg[c] - e;
                u64 co = d.coff[c];
                x->loff = (u32)co;
                w = lim < 0 ? 0 : (u64)d.Lidx[co + d.cnu[c] + lim];
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
            if (i < n) { d.woff[i] = carry + pre; if (i < SL_SH) s_woff[i] = carry + pre; }
            __syncthreads();
            if (tid == 0) carry += agg;
            __syncthreads();
        }
        if (tid == 0) { d.woff[n] = carry; Wtot = carry; }
        __syncthreads();
        u64 W = Wtot;
        if (W > wmax) { if (tid == 0) { out[0] = e; out[1] = W; out[2] = 1; out[3] = n; } return; }
        const u64 *wo = n <= SL_SH ? s_woff : d.woff;
        /* scatter */
        for (u64 t = tid; t < W; t += blockDim.x) {
            u32 lo = 0, hi = n;
            while (hi - lo > 1) { u32 mid = (lo + hi) / 2; if (wo[mid] <= t) lo = mid; else hi = mid; }
            u64 j = t - wo[lo];
            XEnt *x = &d.xtab[ent[lo]];
            u32 c = x->curve;
            Co la = d.L[x->loff + j];
            int u = d.Lidx[x->loff + j];
            i64 h = x_find<XF_SMALL>(d.xtab, d.xmask, c, OK_ADD(x->key, d.ukey[u]), e + d.udeg[u], d.norder, d.order, d.xcap, d.direct, d.lvlcnt, d.lvlpool, d.lvlstride, d.xep, 2);
            if (h < 0) { *d.overflow = 4; continue; }
            if (d.patmode) { pat_rec(d, h, ent[lo], j); continue; }
            u64 y[NL];
            for (int l = 0; l < NL; l++) y[l] = mmul(la.v[l], x->acc[l], l);
            acc_addv(d.xtab[h].acc, y);
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
        if (*(volatile int *)d.overflow) return; /* grid-uniform (set before the last barrier): the batch is redone */
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
                u64 co = d.coff[c];
                x->loff = (u32)co;
                w = lim < 0 ? 0 : (u64)d.Lidx[co + d.cnu[c] + lim];
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
            Co la = d.L[x->loff + j];
            int u = d.Lidx[x->loff + j];
            i64 h = x_find<XF_GPU>(d.xtab, d.xmask, c, OK_ADD(x->key, d.ukey[u]), e + d.udeg[u], d.norder, d.order, d.xcap, d.direct, d.lvlcnt, d.lvlpool, d.lvlstride, d.xep, 2);
            if (h < 0) { *d.overflow = 4; continue; }
            u64 y[NL];
            for (int l = 0; l < NL; l++) y[l] = mmul(la.v[l], x->acc[l], l);
            acc_addv(d.xtab[h].acc, y);
        }
        if (gtid == 0) bsum[nb + 1] += W; /* work counter for the profile */
        grid.sync();
    }
}

/* emit: I[C + b] -= s_C f_C[b] for the finished entries [lo, hi) of the order with b != 0. Their slots need
 * no clearing (the next batch on this buffer has a new epoch, see x_find) */
__global__ void k_emit_range(Dev d, u32 lo, u32 hi, int *eflag) {
    u32 i = lo + blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= hi) return;
    XEnt *x = &d.xtab[d.order[i]];
    int nz = 0;
    for (int l = 0; l < NL; l++) nz |= x->acc[l] != 0;
    if (x->deg != 0 && nz) {
        u32 c = x->curve;
        i64 h = i_find(d.itab, d.imask, EMIT_KEY(d.ckey[c], x->key, x->deg), d.cdeg[c] + x->deg, d.icount, d.icap);
        if (h < 0) { *eflag = 5; } else
        { u64 y[NL]; for (int l = 0; l < NL; l++) y[l] = mneg(mmul(d.cs[c].v[l], x->acc[l], l), l); acc_addv(d.itab[h].acc, y); }
    }
}

/* ======================= class-pattern replay =======================
 * The curves of a layer with the same degree and the same nonzero coordinates have their L terms in one
 * set (the alpha supports of those coordinates), so their exps live on one point set with one recursion
 * pattern. The pattern is built once per class (patmode: one hashed run that records every scatter pair),
 * renumbered in level order and turned into per-target lists of (source, L term); then 32 curves at a
 * time (one per lane) replay it by plain indexing: no probing, no atomics, one reduction per target. */
/* class r of a build batch: its representative (the batch's curve c0 + r) and other curves (ost..), its points
 * [pbase, pbase + np) of the batch numbering, its L terms ljx[ljo..ljo + nlu) */
struct CbCls { u32 rep, ost, pbase, np, ljo; int nlu; };
/* a replay item: curves ci0..ci0 + nv - 1 of class r on S lanes (S | 32), values at F + Fo, L values at LV + Lo */
struct CbItem { u32 r, ci0, nv, S; u64 Fo, Lo; };
__device__ __forceinline__ u32 cb_curve(const CbCls &c, u32 ci) { return ci ? c.ost + ci - 1 : c.rep; }
/* last i < n with off[i] <= t */
__device__ __forceinline__ u32 cb_find(const u64 *off, u32 n, u64 t) {
    u32 lo = 0, hi = n;
    while (hi - lo > 1) { u32 mid = (lo + hi) / 2; if (off[mid] <= t) lo = mid; else hi = mid; }
    return lo;
}
/* points of the build batch per (class, level) */
__global__ void k_pat_hist(Dev d, u32 c0, int T, u32 nent, u32 *cnt) {
    u32 i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nent) return;
    const XEnt *x = &d.xtab[d.order[i]];
    atomicAdd(&cnt[(x->curve - c0) * (T + 1) + x->deg], 1u);
}
/* number the points class by class, level by level (any order within a level); index into the slot's loff */
__global__ void k_pat_assign(Dev d, u32 c0, int T, u32 nent, u32 *cur, u128 *pkey, int *plev) {
    u32 i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nent) return;
    XEnt *x = &d.xtab[d.order[i]];
    u32 q = atomicAdd(&cur[(x->curve - c0) * (T + 1) + x->deg], 1u);
    x->loff = q;   /* (loff is free once the levels are done) */ pkey[q] = x->key; plev[q] = x->deg;
}
/* pairs to (target number, source number within its class << 32 | L index) */
__global__ void k_pat_conv(Dev d, u32 c0, const CbCls *cl, u32 n) {
    u32 x = blockIdx.x * blockDim.x + threadIdx.x;
    if (x >= n) return;
    d.pk[x] = d.xtab[d.pk[x]].loff;
    u64 v = d.pv[x];
    const XEnt *s = &d.xtab[v >> 32];
    d.pv[x] = ((u64)(s->loff - cl[s->curve - c0].pbase) << 32) | (u32)v;
}
/* list offsets of the target-sorted pairs (points without pairs, e.g. each class's point 0, get empty lists) */
__global__ void k_pat_off(const u32 *key, u32 n, u32 npt, u32 *toff) {
    u32 x = blockIdx.x * blockDim.x + threadIdx.x;
    if (x > n) return;
    u32 a = x ? key[x - 1] + 1 : 0, b = x < n ? key[x] : npt;
    for (u32 t = a; t <= b; t++) toff[t] = x;
}
/* L values of the items' curves (thread per (item, L term, lane)); the extra row j = nlu sets the level-0 value 1 */
__global__ void k_cb_lv(Dev d, const CbItem *it, const CbCls *cl, const u64 *off, u32 nit, u64 W, const int *ljx, u64 *F, u64 *LV) {
    u64 t = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= W) return;
    u32 i = cb_find(off, nit, t);
    const CbItem &m = it[i]; const CbCls &c = cl[m.r];
    u64 tt = t - off[i];
    u32 j = (u32)(tt / m.S), k = (u32)(tt % m.S);
    int ok = k < m.nv;
    if (j == (u32)c.nlu) { for (int l = 0; l < NL; l++) F[m.Fo + (u64)l * m.S + k] = ok ? c_one[l] : 0; return; }
    u64 v[NL];
    for (int l = 0; l < NL; l++) v[l] = 0;
    if (ok) {
        int C[64];
        unpack_d(d.ckey[cb_curve(c, m.ci0 + k)], C);
        const Co *al = d.ualpha + (size_t)ljx[c.ljo + j] * c_h11;
        for (int q = 0; q < c_h11; q++) if (C[q]) {
            Co a = al[q];
            for (int l = 0; l < NL; l++) if (a.v[l]) v[l] = madd(v[l], mmul(mint(C[q], l), a.v[l], l), l);
        }
    }
    for (int l = 0; l < NL; l++) LV[m.Lo + ((u64)j * NL + l) * m.S + k] = v[l];
}
/* y < p 2^64 (y_hi < p): y 2^-64 mod p */
__device__ __forceinline__ u64 mredc2(u64 lo, u64 hi, int l) {
    u64 m = lo * c_pinv[l], mhi = __umul64hi(m, c_p[l]);
    u64 u = hi + mhi + (lo != 0);
    return u >= c_p[l] ? u - c_p[l] : u;
}
/* level g: thread per (item, target of level g, lane). The products are summed unreduced (128 bits + a count
 * of 2^128 wraps) and folded once: X = c 2^128 + hi 2^64 + lo == c K2 + hi K1 + lo (mod p) with K1 = 2^64 mod p,
 * K2 = 2^128 mod p both < 2^32 (checked on the host), so the folded value is < 2^97 and one REDC gives
 * sum mmul(L, f) mod p, the residue the scatter accumulates */
__global__ void k_cb_level(Dev d, int g, int T, const CbItem *it, const CbCls *cl, const u32 *clv, const u64 *off, u32 nit, u64 W,
                           const u64 *pat, const u32 *toff, u64 *F, const u64 *LV) {
    u64 th = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (th >= W) return;
    u32 i = cb_find(off, nit, th);
    const CbItem &m = it[i];
    u64 tt = th - off[i];
    u32 S = m.S, k = (u32)(tt % S);
    if (k >= m.nv) return;
    const CbCls &c = cl[m.r];
    u32 t = clv[m.r * (T + 2) + g] + (u32)(tt / S);
    const u64 *Fg = F + m.Fo + k, *Lg = LV + m.Lo + k;
    u64 lo[NL], hi[NL], cy[NL];
    for (int l = 0; l < NL; l++) lo[l] = hi[l] = cy[l] = 0;
    u32 x1 = toff[c.pbase + t + 1];
#pragma unroll 4
    for (u32 x = toff[c.pbase + t]; x < x1; x++) {
        u64 v = pat[x];
        const u64 *fs = Fg + (v >> 32) * NL * S, *ls = Lg + (u64)(u32)v * NL * S;
        for (int l = 0; l < NL; l++) {
            u64 a = fs[l * S], b = ls[l * S];
            u64 pl = a * b, ph = __umul64hi(a, b);
            lo[l] += pl; ph += lo[l] < pl;
            u64 h0 = hi[l]; hi[l] += ph; cy[l] += hi[l] < h0;
        }
    }
    u64 *Ft = F + m.Fo + (u64)t * NL * S + k;
    for (int l = 0; l < NL; l++) {
        u64 k1 = c_r64[l], k2 = c_r2[l];
        u64 ylo = hi[l] * k1, yhi = __umul64hi(hi[l], k1);
        ylo += lo[l]; yhi += ylo < lo[l];
        u64 c2 = cy[l] * k2;
        ylo += c2; yhi += ylo < c2;
        Ft[l * S] = mmul(mredc2(ylo, yhi, l), d.INV[g].v[l], l);
    }
}
/* a point of the I table, or -1 if absent. Plain loads, no acquire: only while no thread inserts (the I
 * table's keys and states do not change; the concurrent atomics touch acc only) */
__device__ __forceinline__ i64 i_look(const IEnt *tab, u64 mask, GKey k) {
    u64 h = gk_hash(k) & mask;
    for (u64 probes = 0; probes <= mask; probes++) {
        const IEnt *e = &tab[h];
        GKey ek = e->key; u32 st = e->state;
        if ((st == 0) | GK_EQ(ek, k)) return st ? (i64)h : -1;   /* one branch: both loads issue together */
        h = (h + 1) & mask;
    }
    return -1;
}
/* I[C + t] -= s_C f_C[t] for points t >= 1, emit thread th of (item, point, lane). With miss set (first
 * pass), the point is looked up without inserting, and a thread whose point is new sets its bit (th - lo)
 * instead; the second pass (miss NULL) inserts those */
__device__ __forceinline__ void cb_emit_t(Dev &d, const CbItem *it, const CbCls *cl, const u64 *off, u32 nit, u64 th, u64 lo,
                                          const u64 *F, const u128 *pkey, const int *plev, int *eflag, u32 *miss) {
    u32 i = cb_find(off, nit, th);
    const CbItem &m = it[i];
    u64 tt = th - off[i];
    u32 S = m.S, k = (u32)(tt % S), t = 1 + (u32)(tt / S);
    if (k >= m.nv) return;
    const u64 *f = F + m.Fo + (u64)t * NL * S + k;
    u64 fv[NL]; int nz = 0;
    for (int l = 0; l < NL; l++) { fv[l] = f[l * S]; nz |= fv[l] != 0; }
    if (!nz) return;
    const CbCls &c = cl[m.r];
    u32 cu = cb_curve(c, m.ci0 + k);
    GKey key = EMIT_KEY(d.ckey[cu], pkey[c.pbase + t], plev[c.pbase + t]);
    i64 h;
    if (miss) {
        h = i_look(d.itab, d.imask, key);
        if (h < 0) { atomicOr(&miss[(th - lo) >> 5], 1u << ((th - lo) & 31)); return; }
    } else {
        h = i_find(d.itab, d.imask, key, d.cdeg[cu] + plev[c.pbase + t], d.icount, d.icap);
        if (h < 0) { *eflag = 5; return; }
    }
    u64 y[NL];
    for (int l = 0; l < NL; l++) y[l] = mneg(mmul(d.cs[cu].v[l], fv[l], l), l);
    acc_addv(d.itab[h].acc, y);
}
/* emit threads [lo, hi); miss: see cb_emit_t (NULL: one pass that inserts) */
__global__ void k_cb_emit(Dev d, const CbItem *it, const CbCls *cl, const u64 *off, u32 nit, u64 lo, u64 hi,
                          const u64 *F, const u128 *pkey, const int *plev, int *eflag, u32 *miss) {
    u64 th = lo + (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (th < hi) cb_emit_t(d, it, cl, off, nit, th, lo, F, pkey, plev, eflag, miss);
}
/* second pass: thread per bitmap word, inserts the marked threads' points and clears the word */
__global__ void k_cb_emit_miss(Dev d, const CbItem *it, const CbCls *cl, const u64 *off, u32 nit, u64 lo, u32 nw,
                               const u64 *F, const u128 *pkey, const int *plev, int *eflag, u32 *miss) {
    u32 w = blockIdx.x * blockDim.x + threadIdx.x;
    if (w >= nw) return;
    u32 v = miss[w];
    if (!v) return;
    miss[w] = 0;
    for (; v; v &= v - 1) cb_emit_t(d, it, cl, off, nit, lo + (u64)w * 32 + __ffs(v) - 1, lo, F, pkey, plev, eflag, NULL);
}
/* reorder the layer's curves (ckey, cdeg, cs, cnu) by perm (new position -> old index) */
__global__ void k_permute(Dev d, const u32 *perm, u32 n, GKey *k2, int *g2, Co *s2, int *u2) {
    u32 i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    u32 o = perm[i];
    k2[i] = d.ckey[o]; g2[i] = d.cdeg[o]; s2[i] = d.cs[o]; u2[i] = d.cnu[o];
}

/* move every entry of an old I table into the (empty) current one */
__global__ void k_rehash_I(Dev d, const IEnt *old, u64 oldcap) {
    u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= oldcap) return;
    const IEnt *e = &old[i];
    if (e->state != 2) return;
    u64 h = gk_hash(e->key) & d.imask;
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
__global__ void k_collect(Dev d, GKey *okey, int *odeg, Co *oA, u32 *n) {
    u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    if (i > d.imask) return;
    IEnt *e = &d.itab[i];
    if (e->state != 2) return;
    u32 o = atomicAdd(n, 1u);
    okey[o] = e->key; odeg[o] = e->deg;
    for (int l = 0; l < NL; l++) oA[o].v[l] = e->acc[l];
}

static double wall(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }
enum { T_SETUP, T_EXTRACT, T_BUILDL, T_BUCKET, T_LEVEL, T_EMIT, T_WASTE, T_N };
static const char *TNAME[T_N] = {"setup", "extract", "buildL+seed", "bucket", "level", "emit", "retried(wasted)"};
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
template <class T> static T *talloc(size_t n) {   /* NULL if the allocation fails */
    T *p = NULL;
    if (cudaMalloc(&p, (n ? n : 1) * sizeof(T)) != cudaSuccess) { (void)cudaGetLastError(); return NULL; }
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
    u64 h = (hash128(k) ^ ((u64)c * 0xA24BAED4963EE407ull) XHDEG) & mask;
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
        if (*(volatile u32 *)&e->curve == c && *(volatile u128 *)&e->key == k XDEG(e)) return (i64)h;
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
        h = x_findB(xb, d.xmask, beta, OK_ADD(x->key, d.ukey[u]), e + d.udeg[u], d.norder, d.order, d.xcap, d.direct, d.lvlcnt, d.lvlpool, d.lvlstride);
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
        i64 h = i_find(d.itab, d.imask, EMIT_KEY(d.ckey[c], x->key, x->deg), d.cdeg[c] + x->deg, d.icount, d.icap);
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
    {
        int kbv = in->kbw != NULL;
        CK(cudaMemcpyToSymbol(c_kbv, &kbv, sizeof(int)));
        if (kbv) CK(cudaMemcpyToSymbol(c_kbw, in->kbw, in->h11 * sizeof(int)));
    }
    CK(cudaMemcpyToSymbol(c_maxdeg, &in->maxdeg, sizeof(int)));
#ifdef CGV_WIDE
    if (!in->dK || in->dn > 64) { fprintf(stderr, "cgv gpu: wide keys need the level-code decoder\n"); return 1; }
    {
        u128 plo = (u128)in->dplo[1] << 64 | in->dplo[0]; i64 wlo = in->dwlo;
        CK(cudaMemcpyToSymbol(c_d128, &in->d128, sizeof(int))); CK(cudaMemcpyToSymbol(c_dn, &in->dn, sizeof(int)));
        CK(cudaMemcpyToSymbol(c_dt0, &in->dt0, sizeof(int))); CK(cudaMemcpyToSymbol(c_dw0, &in->dw0, sizeof(int)));
        CK(cudaMemcpyToSymbol(c_dwlo, &wlo, sizeof(i64))); CK(cudaMemcpyToSymbol(c_dplo, &plo, sizeof(u128)));
        CK(cudaMemcpyToSymbol(c_doff, in->doff, sizeof(int) * 64)); CK(cudaMemcpyToSymbol(c_db, in->db, sizeof(int) * 64));
        CK(cudaMemcpyToSymbol(c_dw, in->dw, sizeof(int) * 64));
        if (in->dn) CK(cudaMemcpyToSymbol(c_dK, in->dK, sizeof(GKey) * in->dn));
        CK(cudaMemcpyToSymbol(c_dbase, in->dbase, sizeof(GKey))); CK(cudaMemcpyToSymbol(c_dKt0, in->dKt0, sizeof(GKey)));
    }
#endif
    int H = in->h11, D = in->maxdeg;

    Dev d; memset(&d, 0, sizeof(d));
    if (getenv("CGV_XEP0")) d.xep = (u32)atoi(getenv("CGV_XEP0"));   /* test knob: start near the epoch wrap */
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
        GKey *k = dalloc<GKey>(in->nI); int *g = dalloc<int>(in->nI); Co *v = dalloc<Co>(in->nI);
        h2d(k, (const GKey *)in->ikey, in->nI); h2d(g, in->ideg, in->nI); h2d(v, (const Co *)in->ival, in->nI);
        k_init_I<<<nblk(in->nI, 256), 256>>>(d, in->nI, k, g, v);
        SYNC();
        cudaFree(k); cudaFree(g); cudaFree(v);
    }

    /* curve buffers */
    d.ccap = getenv("CGV_CCAP_LOG2") ? 1u << atoi(getenv("CGV_CCAP_LOG2")) : 1u << 22;  /* test knob: start small to exercise growth */
    auto alloc_curves = [&](u32 n) {
        d.ccap = n;
        d.ckey = dalloc<GKey>(n); d.cdeg = dalloc<int>(n); d.cs = dalloc<Co>(n);
        d.cnu = dalloc<int>(n); d.coff = dalloc<u64>(n); d.cnl = dalloc<int>(n);
    };
    auto free_curves = [&]() { cudaFree(d.ckey); cudaFree(d.cdeg); cudaFree(d.cs); cudaFree(d.cnu); cudaFree(d.coff); cudaFree(d.cnl); };
    d.ckey = dalloc<GKey>(d.ccap); d.cdeg = dalloc<int>(d.ccap); d.cs = dalloc<Co>(d.ccap);
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
    /* double-buffered exp table (c028): batch k's emit runs on its own stream over one buffer
     * while batch k+1's levels run over the other; ev_emit[b] marks the end of the last emit on buffer b */
    XEnt *xtb[2] = {NULL, NULL}; u32 *ordb[2] = {NULL, NULL};
    int dbl = 0, dbl_off = getenv("CGV_NO_DBL") != NULL, cur = 0, ev_used[2] = {0, 0};
    const u64 fill_free8 = getenv("CGV_FILL_FREE8") ? atoi(getenv("CGV_FILL_FREE8")) : 2; /* batch sub-tables may fill to m - m*k/8 slots (k = 2: 3/4) before a batch counts as overflowing */
    u64 dbl_slots = 0; /* capacity of buffer 1 (slots; its order array has as many entries) */
    u64 xtb0_slots = 0; /* capacity of buffer 0 */
    auto alloc_exp = [&](u64 n) {
        d.xmask = n - 1; d.xcap = (u32)std::min<u64>(n / 2, 0xffffffffu);
        d.xtab = dalloc<XEnt>(n); CK(cudaMemset(d.xtab, 0, n * sizeof(XEnt))); xtb0_slots = n;
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
        CK(cub::DeviceScan::ExclusiveSum(NULL, d.scan_tmp_bytes, d.work, d.woff, (int)std::min<u64>((u64)d.xcap + 1, 0x7fffffff)));
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
    d.Lcap = std::min<u64>(std::max<u64>(Lbytes / (sizeof(Co) + sizeof(int)), 1 << 16), 0xffffffffu); /* XEnt.loff is 32-bit */
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
        n = std::min<u64>(n, 0xffffffffu); /* XEnt.loff is 32-bit */
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
    /* level scheduling (environment): CGV_DLEV=0|1 device-driven levels (default on with CUDA, off with HIP);
     * CGV_DLEV_SG resident blocks per SM for them (tuning; default from the occupancy); CGV_COOP=1 one cooperative
     * kernel for all levels; CGV_SMALL_LEVELS=0 no single-block kernel for small levels (both test/tuning) */
    /* class-pattern replay. Switches (environment):
     *   CGV_CB=0         class replay off (every curve takes the hashed path); CGV_MEM=low implies it (both
     *                    read in gv.c, passed as in->cb)
     *   CGV_CB_MIN       smallest class replayed (curves, default 4)
     *   CGV_CB_PAIRS     pattern pair buffer (pairs, 24 B each; default 2^25, at most 4% of the device's memory)
     *   CGV_CB_FMEM_GB   replay value buffer (default 0.5 GB, at most 2% of the device's memory)
     *   CGV_CB_XT        exp entries per pattern-build batch (default 1e6)
     * (the 5090: +1.3 GB, peak +23% on g32). The buffers are allocated at the first class batch and capped by
     * the memory then free (minus room for an I growth); an I or exp table growth that needs the room releases
     * them (they are reallocated at the next class batch), and classes that find too little memory fall back. */
    int cb_on = in->cb;
    int cb_min = getenv("CGV_CB_MIN") ? atoi(getenv("CGV_CB_MIN")) : 4;
    double cb_pairs_req = getenv("CGV_CB_PAIRS") ? atof(getenv("CGV_CB_PAIRS")) : (double)std::min<u64>((u64)1 << 25, totmem / 25 / 24);
    size_t cb_fmem_req = getenv("CGV_CB_FMEM_GB") ? (size_t)(atof(getenv("CGV_CB_FMEM_GB")) * 1e9) : std::min<size_t>((size_t)5e8, totmem / 50);
    for (int l = 0; l < NL; l++) if (r64[l] >> 32 || in->r2[l] >> 32) cb_on = 0;   /* the unreduced fold needs small 2^64, 2^128 mod p */
    struct HCls { u32 rep, ost, n; };
    std::vector<HCls> cls; u32 nrest = 0; std::vector<int> cdg;
    u32 cb_pcap = 0; size_t cb_fmem = 0;
    u32 *cb_pk2 = NULL; u64 *cb_pv2 = NULL; void *cb_stmp = NULL; size_t cb_stb = 0;
    u32 *cb_miss = NULL; int cb_miss_tried = 0;   /* class emits: bitmap of the first pass's new points (K bits) */
    u32 *cb_toff = NULL; u128 *cb_pkey = NULL; int *cb_plev = NULL, *cb_ljx = NULL; u64 cb_npcap = 0, cb_ljcap = 0;
    u64 *cb_F = NULL, *cb_LV = NULL; u64 cb_Fcap = 0;
    u32 *cb_cnt = NULL, *cb_clv = NULL; CbCls *cb_cl = NULL; u64 cb_clcap = 0;
    CbItem *cb_it = NULL; u64 *cb_off = NULL; u64 cb_itcap = 0, cb_offcap = 0;
    double cb_xt = getenv("CGV_CB_XT") ? atof(getenv("CGV_CB_XT")) : 1e6, est_p = 1;
    double t_cbb = 0, t_cbr = 0, t_cbe = 0; u64 cb_bW = 0; i64 cb_ncl = 0, cb_ncur = 0, cb_fail = 0;
    int cb_busy = 0;   /* a class batch is being replayed: its buffers stay */
    /* release the big class buffers (pairs, sort temp, values, points, L indices); 1 if that freed anything */
    auto cb_release = [&]() {
        if (cb_busy || !(d.pk || cb_F || cb_toff || cb_ljx)) return 0;
        SYNC();
        void *p[] = {d.pk, d.pv, cb_pk2, cb_pv2, cb_stmp, cb_toff, cb_pkey, cb_plev, cb_ljx, cb_F};
        for (void *q : p) if (q) cudaFree(q);
        d.pk = cb_pk2 = cb_toff = NULL; d.pv = cb_pv2 = cb_F = NULL; cb_stmp = NULL; cb_pkey = NULL; cb_plev = cb_ljx = NULL;
        cb_npcap = cb_ljcap = cb_Fcap = 0;
        if (in->verbose) fprintf(stderr, "  gpu: class buffers released\n");
        return 1;
    };
    /* (re)allocate the pair buffers and the sort temp; 0 if memory is too short (the classes fall back) */
    auto cb_alloc = [&]() {
        if (d.pk) return 1;
        size_t f = 0, t = 0; CK(cudaMemGetInfo(&f, &t)); f = budget_free(f);
        size_t keep = icap * sizeof(IEnt) + ((size_t)256 << 20), av = f > keep ? f - keep : 0;   /* room for an I growth, as ensure_L */
        cb_pcap = (u32)std::min<double>(std::min<double>(cb_pairs_req, 0x7fffffff), (double)(av / 4 / 24));
        cb_fmem = std::min<size_t>(cb_fmem_req, av / 8);
        if (cb_pcap < (1u << 16)) return 0;
        if (!d.npair) d.npair = dalloc<u32>(1);
        d.pcap = cb_pcap;
        d.pk = talloc<u32>(cb_pcap); d.pv = talloc<u64>(cb_pcap); cb_pk2 = talloc<u32>(cb_pcap); cb_pv2 = talloc<u64>(cb_pcap);
        cb_stb = 0;
        if (d.pk && d.pv && cb_pk2 && cb_pv2) {
            cub::DoubleBuffer<u32> kb(d.pk, cb_pk2); cub::DoubleBuffer<u64> vb(d.pv, cb_pv2);   /* no alternate buffers in the temp storage */
            CK(cub::DeviceRadixSort::SortPairs(NULL, cb_stb, kb, vb, (int)cb_pcap));
            cb_stmp = talloc<char>(cb_stb);
        }
        if (cb_stmp) return 1;
        void *p[] = {d.pk, d.pv, cb_pk2, cb_pv2};
        for (void *q : p) if (q) cudaFree(q);
        d.pk = cb_pk2 = NULL; d.pv = cb_pv2 = NULL;
        return 0;
    };
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
        auto grow_dev = [&](const char *how) {   /* room for both tables: rehash on the device */
            d.itab = dalloc<IEnt>(ncap); CK(cudaMemset(d.itab, 0, ncap * sizeof(IEnt)));
            icap = ncap; d.imask = icap - 1; d.icap = (u32)std::min<u64>((u64)(icap * i_capf), 0xffffffffu);
            k_rehash_I<<<nblk(ocap, 256), 256>>>(d, old, ocap);
            SYNC();
            cudaFree(old);
            if (in->verbose) fprintf(stderr, "  gpu: I table grown to %llu entries%s\n", (unsigned long long)icap, how);
        };
        if (nbytes + slack <= f) { grow_dev(""); return; }
        /* memory is short: release the class buffers (reallocated at the next class batch) */
        if (cb_release()) {
            CK(cudaMemGetInfo(&f, &t));
            if (nbytes + slack <= f) { grow_dev(" (class buffers released)"); return; }
        }
        /* staged: park the old table in host memory, so only the new one is on the device.
         * If that still does not fit, release the L pool, then (between layers only, when the
         * exp table is empty) halve the exp table until it does. */
        const char *how = "staged via host";
        if (nbytes + slack > f + obytes && dbl && cur == 0) {
            /* go single-buffered (buffer 1 is idle: the host synced the emit stream) */
            SYNC(); cudaFree(xtb[1]); cudaFree(ordb[1]); xtb[1] = NULL; ordb[1] = NULL; dbl = 0; dbl_off = 1;
            CK(cudaMemGetInfo(&f, &t));
            if (nbytes + slack <= f) { grow_dev(" (exp table single-buffered)"); return; }
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
#ifdef CGV_HIP
    int dlev = getenv("CGV_DLEV") && atoi(getenv("CGV_DLEV"));    /* device-driven levels: off by default on AMD (slower there) */
#else
    int dlev = !getenv("CGV_DLEV") || atoi(getenv("CGV_DLEV"));   /* device-driven levels (CGV_DLEV=0: off) */
#endif
    u64 *d_dw = dalloc<u64>(DL_G + 3), *d_dr = dalloc<u64>(5); CK(cudaMemset(d_dw, 0, (DL_G + 3) * sizeof(u64)));
    int dl_per = 0; CK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&dl_per, k_dscatter, 256, 0));
    unsigned dl_sg = (unsigned)prop.multiProcessorCount * (getenv("CGV_DLEV_SG") ? atoi(getenv("CGV_DLEV_SG")) : std::max(1, dl_per));
    cudaStream_t s_emit; CK(cudaStreamCreateWithFlags(&s_emit, cudaStreamNonBlocking));
    cudaEvent_t ev_emit[2], ev_lev;
    CK(cudaEventCreateWithFlags(&ev_emit[0], cudaEventDisableTiming)); CK(cudaEventCreateWithFlags(&ev_emit[1], cudaEventDisableTiming));
    CK(cudaEventCreateWithFlags(&ev_lev, cudaEventDisableTiming));
    int *d_eovf = dalloc<int>(1); CK(cudaMemset(d_eovf, 0, sizeof(int)));
    u32 *h_icnt = NULL; if (HOST_ALLOC((void **)&h_icnt, sizeof(u32)) != cudaSuccess) { h_icnt = NULL; (void)cudaGetLastError(); }
    int icnt_ok = 0; /* *h_icnt holds the I count after the last queued emit (valid once the emit stream is synced) */
    i64 iub = -1; /* upper bound on the I count once the queued emits are done (each entry adds at most one point); -1: unknown */
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
    /* progress lines (verbose runs): at most one per CGV_PROGRESS_S seconds (default 30), checked once per batch */
    double hb_every = getenv("CGV_PROGRESS_S") ? atof(getenv("CGV_PROGRESS_S")) : 30, hb_t0 = wall(), hb_last = hb_t0;
#define GPU_HB(DD, NC, B) do { double t_ = wall(); if (in->verbose && t_ - hb_last >= hb_every) { hb_last = t_; \
        fprintf(stderr, "    ... gpu extraction %.0fs: layer d=%d of %d (%u curves), batch %d of this layer; %lld curves through this layer\n", \
                t_ - hb_t0, (DD), D, (unsigned)(NC), (B), (long long)napplied); } } while (0)
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
            if (f + Lb < need && cb_release()) CK(cudaMemGetInfo(&f, &t));   /* the class buffers first */
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
        std::vector<int> hdeg(nc);
        cdg.resize(nc);
        {
            if (pat_cap < nc) { cudaFree(d_pat); pat_cap = nc; d_pat = dalloc<u64>(pat_cap); }
            k_pattern<<<nblk(nc, 256), 256>>>(d, nc, d_pat);
            std::vector<u64> hp(nc);
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
            /* curve classes (same degree, same nonzero coordinates) of >= cb_min curves go to the end of the
             * layer's curve arrays, one contiguous range per class; the rest stay in front for the hashed path */
            cls.clear(); nrest = nc;
            if (cb_on && !(bundle_min && nc >= (u32)bundle_min && super_delta == 1)) {
                std::vector<u32> idx(nc);
                for (u32 i = 0; i < nc; i++) idx[i] = i;
                std::sort(idx.begin(), idx.end(), [&](u32 a, u32 b) { return hdeg[a] != hdeg[b] ? hdeg[a] < hdeg[b] : hp[a] != hp[b] ? hp[a] < hp[b] : a < b; });
                std::vector<u32> perm, tail; std::vector<std::pair<u32, u32>> cr;
                for (u32 i = 0; i < nc;) {
                    u32 j = i;
                    while (j < nc && hdeg[idx[j]] == hdeg[idx[i]] && hp[idx[j]] == hp[idx[i]]) j++;
                    if (j - i >= (u32)cb_min) { cr.push_back({(u32)tail.size(), j - i}); for (u32 q = i; q < j; q++) tail.push_back(idx[q]); }
                    else for (u32 q = i; q < j; q++) perm.push_back(idx[q]);
                    i = j;
                }
                if (!cr.empty()) {
                    /* order: hashed curves (in their old order), the classes' first curves (the representatives,
                     * contiguous so that pattern builds batch them), the classes' other curves */
                    std::sort(perm.begin(), perm.end());
                    nrest = (u32)perm.size();
                    u32 R0 = nrest, nk = (u32)cr.size(), o = R0 + nk;
                    for (auto &c : cr) perm.push_back(tail[c.first]);
                    for (u32 k = 0; k < nk; k++) {
                        cls.push_back({R0 + k, o, cr[k].second});
                        for (u32 q = 1; q < cr[k].second; q++) perm.push_back(tail[cr[k].first + q]);
                        o += cr[k].second - 1;
                    }
                    u32 *dp = dalloc<u32>(nc); h2d(dp, perm.data(), nc);
                    GKey *k2 = dalloc<GKey>(nc); int *g2 = dalloc<int>(nc); Co *s2 = dalloc<Co>(nc); int *u2 = dalloc<int>(nc);
                    k_permute<<<nblk(nc, 256), 256>>>(d, dp, nc, k2, g2, s2, u2);
                    CK(cudaMemcpy(d.ckey, k2, nc * sizeof(GKey), cudaMemcpyDeviceToDevice)); CK(cudaMemcpy(d.cdeg, g2, nc * sizeof(int), cudaMemcpyDeviceToDevice));
                    CK(cudaMemcpy(d.cs, s2, nc * sizeof(Co), cudaMemcpyDeviceToDevice)); CK(cudaMemcpy(d.cnu, u2, nc * sizeof(int), cudaMemcpyDeviceToDevice));
                    SYNC();
                    cudaFree(dp); cudaFree(k2); cudaFree(g2); cudaFree(s2); cudaFree(u2);
                    std::vector<int> nu2(nc); std::vector<double> nl2(nc);
                    for (u32 i = 0; i < nc; i++) { nu2[i] = hnu[perm[i]]; nl2[i] = hnl[perm[i]]; cdg[i] = hdeg[perm[i]]; }
                    hnu.swap(nu2); hnl.swap(nl2); hdeg = cdg;
                }
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
                GPU_HB(dd, nc, layB);
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
        std::vector<std::pair<u32, u32>> ranges;   /* curve ranges for the hashed path */
        double lcb0[3] = {t_cbb, t_cbr, t_cbe};
        if (nrest) ranges.push_back({0, nrest});
        /* ---- classes: pattern builds for batches of representatives, then replay of the batch's classes ---- */
        u32 ncls = (u32)cls.size(), k0 = 0, kmax = ncls; u64 force_m = 0;
        auto cb_fallback = [&](u32 k) { ranges.push_back({cls[k].rep, cls[k].rep + 1}); if (cls[k].n > 1) ranges.push_back({cls[k].ost, cls[k].ost + cls[k].n - 1}); cb_fail++; };
        /* the hashed path for the not yet replayed curves ci0..ci0 + nv - 1 of a class (curve ci: rep for ci = 0, else ost + ci - 1) */
        auto cb_fallback_item = [&](const CbItem &m) {
            const HCls &c = cls[k0 + m.r];
            u32 a = m.ci0, b = m.ci0 + m.nv;
            if (a == 0) { ranges.push_back({c.rep, c.rep + 1}); a = 1; }
            if (a < b) ranges.push_back({c.ost + a - 1, c.ost + b - 1});
            cb_fail++;
        };
        i64 cib = -1;   /* upper bound on the I count while the class emits need no room check; -1: unknown */
        while (k0 < ncls) {
            double tc0 = wall();
            if (!cb_alloc()) { for (; k0 < ncls; k0++) cb_fallback(k0); break; }   /* memory is short: hashed path */
            int T = D - cdg[cls[k0].rep];    /* one degree per class batch */
            u32 nb = std::min(kmax, ncls - k0), c0 = cls[k0].rep;
            { u32 q = 1; while (q < nb && cdg[c0 + q] == cdg[c0]) q++; nb = q; }
            if (!force_m) { double a = 0; u32 q = 0; for (; q < nb; q++) { a += est_p * hnl[c0 + q] + 1; if (a > cb_xt && q > 0) break; } nb = std::max<u32>(1, q); }
            ensure_L();
            u64 tot = 0; u32 q = 0;
            for (; q < nb; q++) { u64 sz = hnu[c0 + q] + T + 1; if (tot + sz > d.Lcap) break; hoff[c0 + q] = tot; tot += sz; }   /* L terms + breakpoints */
            nb = q;
            if (!nb) { cb_fallback(k0); k0++; continue; }
            double sumnl = 0; for (q = 0; q < nb; q++) sumnl += hnl[c0 + q];
            u64 m = 1 << 16, need = (u64)((est_p * sumnl + nb) * 3) + 1024;
            while (m < need && m <= full_mask) m <<= 1;
            if (force_m) m = force_m;
            if (m > full_mask + 1) m = full_mask + 1;
            d.xmask = m - 1; d.xcap = (u32)std::min<u64>(m <= full_mask ? m - m * fill_free8 / 8 : m / 2, full_cap);
            d.lvlstride = d.xcap; d.direct = 1;   /* (T + 1) xcap <= (D + 2) full_cap = lvl_cap */
            d.xtab = xtb[0]; d.order = ordb[0]; d.patmode = 1;
            if (++d.xep >= 1u << 24) {   /* new epoch, as for a hashed batch (the wrap clears both buffers) */
                CK(cudaStreamSynchronize(s_emit));
                for (int q = 0; q < 2; q++) if (xtb[q]) CK(cudaMemset(xtb[q], 0, (q ? dbl_slots : xtb0_slots) * sizeof(XEnt)));
                d.xep = 1;
            }
            CK(cudaMemset(d.lvlcnt, 0, (D + 2) * sizeof(u32))); CK(cudaMemset(d.norder, 0, sizeof(u32))); CK(cudaMemset(d.npair, 0, sizeof(u32)));
            CK(cudaMemcpy(d.coff + c0, hoff.data() + c0, nb * sizeof(u64), cudaMemcpyHostToDevice));
            k_buildL<<<nb, 256>>>(d, c0);
            k_seed<<<nblk(nb, 256), 256>>>(d, c0, nb);
            for (int e = 0; e <= T; e++) {
                k_small_levels<<<1, 1024>>>(d, e, T, sl_nmax, sl_wmax, d_sl);
                u64 so[4];
                CK(cudaGetLastError()); CK(cudaMemcpy(so, d_sl, sizeof(so), cudaMemcpyDeviceToHost));
                if ((int)so[0] > T) break;
                e = (int)so[0];
                u32 n = (u32)so[3];
                if (n > d.lvlstride) { int ov = 6; CK(cudaMemcpy(d.overflow, &ov, sizeof(int), cudaMemcpyHostToDevice)); break; }
                const u32 *ent = d.lvlpool + (u64)e * d.lvlstride;
                u64 W = so[1];
                if (!so[2]) {
                    k_final<<<nblk((u64)n + 1, 256), 256>>>(d, e, ent, n);
                    CK(cub::DeviceScan::ExclusiveSum(d.scan_tmp, d.scan_tmp_bytes, d.work, d.woff, (int)(n + 1)));
                    W = d2h1(d.woff + n);
                }
                if (W) k_scatter<<<nblk(W, 256), 256>>>(d, e, ent, n, W);
            }
            SYNC();
            d.patmode = 0;
            int ov = d2h1(d.overflow);
            u32 npair = d2h1(d.npair), nent = std::min(d2h1(d.norder), d.xcap);
            if (ov) {
                CK(cudaMemset(d.overflow, 0, sizeof(int)));
                d.xmask = full_mask; d.xcap = full_cap;
                if (ov != 7) est_p *= 2;
                if (ov != 7 && m <= full_mask) { force_m = std::min<u64>(m * 4, full_mask + 1); kmax = nb; }
                else if (nb == 1) { cb_fallback(k0); k0++; force_m = 0; kmax = ncls; }
                else { kmax = std::max<u32>(1, nb / 2); force_m = ov == 7 ? 0 : full_mask + 1; }
                t_cbb += wall() - tc0;
                continue;
            }
            est_p = std::max((double)nent / (sumnl + nb), 0.7 * est_p);
            /* number the points per (class, level), convert and sort the pairs into per-target lists */
            if (cb_clcap < nb) {
                cudaFree(cb_cnt); cudaFree(cb_clv); cudaFree(cb_cl);
                cb_clcap = nb + nb / 4 + 16; cb_cnt = dalloc<u32>(cb_clcap * (D + 2)); cb_clv = dalloc<u32>(cb_clcap * (D + 3)); cb_cl = dalloc<CbCls>(cb_clcap);
            }
            CK(cudaMemset(cb_cnt, 0, (size_t)nb * (T + 1) * sizeof(u32)));
            k_pat_hist<<<nblk(nent, 256), 256>>>(d, c0, T, nent, cb_cnt);
            std::vector<u32> hc((size_t)nb * (T + 1)), hcl((size_t)nb * (T + 2)), hcur((size_t)nb * (T + 1));
            std::vector<int> hnlu(nb);
            CK(cudaMemcpy(hc.data(), cb_cnt, hc.size() * sizeof(u32), cudaMemcpyDeviceToHost));
            CK(cudaMemcpy(hnlu.data(), d.cnl + c0, nb * sizeof(int), cudaMemcpyDeviceToHost));
            std::vector<CbCls> hcls(nb);
            u32 npt = 0;
            for (u32 r = 0; r < nb; r++) {
                CbCls &c = hcls[r];
                c.rep = cls[k0 + r].rep; c.ost = cls[k0 + r].ost; c.pbase = npt; c.ljo = (u32)hoff[c0 + r]; c.nlu = hnlu[r];
                u32 a = 0;
                for (int e = 0; e <= T; e++) { hcl[(size_t)r * (T + 2) + e] = a; hcur[(size_t)r * (T + 1) + e] = npt + a; a += hc[(size_t)r * (T + 1) + e]; }
                hcl[(size_t)r * (T + 2) + T + 1] = a;
                c.np = a; npt += a;
            }
            u32 *skey = cb_pk2; u64 *spat = cb_pv2;   /* the sorted pairs (in whichever buffer cub left them) */
            if (cb_npcap < (u64)npt + 1) {
                cudaFree(cb_toff); cudaFree(cb_pkey); cudaFree(cb_plev); cb_npcap = (u64)npt + npt / 4 + 1;
                cb_toff = talloc<u32>(cb_npcap); cb_pkey = talloc<u128>(cb_npcap); cb_plev = talloc<int>(cb_npcap);
                if (!cb_toff || !cb_pkey || !cb_plev) cb_npcap = 0;
            }
            if (cb_ljcap < tot) { cudaFree(cb_ljx); cb_ljcap = tot + tot / 4 + 1; if (!(cb_ljx = talloc<int>(cb_ljcap))) cb_ljcap = 0; }
            if ((!npair && npt > nb) || !cb_npcap || !cb_ljcap) {
                /* no memory for the per-point lists (or, never expected, a build without pairs): the batch's classes fall back */
                if (!npair && npt > nb) fprintf(stderr, "cgv gpu: warning: pattern build recorded no pairs (layer %d); hashed path instead\n", dd);
                d.xmask = full_mask; d.xcap = full_cap;
                for (u32 r = 0; r < nb; r++) cb_fallback(k0 + r);
                k0 += nb; force_m = 0; kmax = ncls;
                t_cbb += wall() - tc0;
                continue;
            }
            h2d(cb_cl, hcls.data(), nb); h2d(cb_clv, hcl.data(), hcl.size()); h2d(cb_cnt, hcur.data(), hcur.size());
            CK(cudaMemcpy(cb_ljx, d.Lidx, tot * sizeof(int), cudaMemcpyDeviceToDevice));   /* the L pool may be released by an I growth */
            k_pat_assign<<<nblk(nent, 256), 256>>>(d, c0, T, nent, cb_cnt, cb_pkey, cb_plev);
            if (npair) {
                k_pat_conv<<<nblk(npair, 256), 256>>>(d, c0, cb_cl, npair);
                int bits = 1; while (bits < 32 && ((u64)1 << bits) < npt) bits++;
                cub::DoubleBuffer<u32> kb(d.pk, cb_pk2); cub::DoubleBuffer<u64> vb(d.pv, cb_pv2);
                CK(cub::DeviceRadixSort::SortPairs(cb_stmp, cb_stb, kb, vb, (int)npair, 0, bits));
                skey = kb.Current(); spat = vb.Current();
            }
            k_pat_off<<<nblk((u64)npair + 1, 256), 256>>>(skey, npair, npt, cb_toff);
            d.xmask = full_mask; d.xcap = full_cap;
            CK(cudaGetLastError());
            double tc1 = wall(); t_cbb += tc1 - tc0;
            cb_bW += npair;
            /* ---- replay: items (class, up to S curves), packed into chunks whose values fit cb_fmem ---- */
            std::vector<CbItem> its;
            for (u32 r = 0; r < nb; r++) {
                u32 n = cls[k0 + r].n, S = 32;
                if (n < 32) { S = 1; while (S < n) S *= 2; }
                while (S > 1 && ((u64)hcls[r].np + hcls[r].nlu) * NL * S * sizeof(u64) > cb_fmem) S /= 2;   /* one item within the value budget */
                for (u32 ci = 0; ci < n; ci += S) its.push_back({r, ci, std::min(S, n - ci), S, 0, 0});
            }
            for (u32 r = 0; r < nb; r++) cb_ncur += cls[k0 + r].n;
            cb_ncl += nb;
            cb_busy = 1;
            for (size_t i0 = 0; i0 < its.size();) {
                double tr0 = wall();
                u64 fo = 0, lo = 0; size_t i1 = i0;
                for (; i1 < its.size(); i1++) {
                    const CbCls &c = hcls[its[i1].r];
                    u64 fs = (u64)c.np * NL * its[i1].S, ls = (u64)c.nlu * NL * its[i1].S;
                    if (i1 > i0 && (fo + fs + lo + ls) * sizeof(u64) > cb_fmem) break;
                    its[i1].Fo = fo; its[i1].Lo = lo; fo += fs; lo += ls;
                }
                u32 nit = (u32)(i1 - i0);
                /* values and L values share one buffer: cb_fmem, or what a single item needs */
                if (cb_Fcap < fo + lo + 1) { cudaFree(cb_F); cb_Fcap = std::max<u64>(fo + lo + 1, cb_fmem / sizeof(u64)); if (!(cb_F = talloc<u64>(cb_Fcap))) cb_Fcap = 0; }
                /* thread offsets: L values, each level, emits */
                std::vector<u64> off((size_t)(T + 2) * (nit + 1));
                for (int g = 0; g <= T + 1; g++) {
                    u64 a = 0; u64 *o = off.data() + (size_t)g * (nit + 1);
                    for (u32 i = 0; i < nit; i++) {
                        const CbItem &mt = its[i0 + i]; const CbCls &c = hcls[mt.r];
                        o[i] = a;
                        if (g == 0) a += (u64)(c.nlu + 1) * mt.S;
                        else if (g <= T) a += (u64)(hcl[(size_t)mt.r * (T + 2) + g + 1] - hcl[(size_t)mt.r * (T + 2) + g]) * mt.S;
                        else a += (u64)(c.np - 1) * mt.S;
                    }
                    o[nit] = a;
                }
                if (cb_itcap < nit) { cudaFree(cb_it); cb_itcap = nit + nit / 4; if (!(cb_it = talloc<CbItem>(cb_itcap))) cb_itcap = 0; }
                if (cb_offcap < off.size()) { cudaFree(cb_off); cb_offcap = off.size() + off.size() / 4; if (!(cb_off = talloc<u64>(cb_offcap))) cb_offcap = 0; }
                if (!cb_Fcap || !cb_itcap || !cb_offcap) {   /* no memory: the items not yet replayed fall back */
                    for (size_t i = i0; i < its.size(); i++) cb_fallback_item(its[i]);
                    break;
                }
                cb_LV = cb_F + fo;
                h2d(cb_it, its.data() + i0, nit); h2d(cb_off, off.data(), off.size());
                u64 W0 = off[nit];
                k_cb_lv<<<nblk(W0, 256), 256>>>(d, cb_it, cb_cl, cb_off, nit, W0, cb_ljx, cb_F, cb_LV);
                for (int g = 1; g <= T; g++) {
                    const u64 *o = off.data() + (size_t)g * (nit + 1);
                    if (o[nit]) k_cb_level<<<nblk(o[nit], 256), 256>>>(d, g, T, cb_it, cb_cl, cb_clv, cb_off + (size_t)g * (nit + 1), nit, o[nit], spat, cb_toff, cb_F, cb_LV);
                }
                CK(cudaGetLastError());
                double tr1 = wall(); t_cbr += tr1 - tr0;
                /* emits in chunks of K threads (item, point, lane); each thread of a lane < nv adds at most one point. The
                 * I count is read back (a sync) for a room check only when the bound cib does not show that the chunk fits */
                const u64 K = 1u << 22, *oe = off.data() + (size_t)(T + 1) * (nit + 1), nem = oe[nit];
                u32 ie = 0;
                if (!cb_miss && !cb_miss_tried) { cb_miss_tried = 1; if ((cb_miss = talloc<u32>(K / 32))) CK(cudaMemset(cb_miss, 0, K / 32 * sizeof(u32))); }
                for (u64 a = 0; a < nem; a += K) {
                    u64 b = std::min(nem, a + K), nv = 0;
                    while (oe[ie + 1] <= a) ie++;
                    for (u32 i = ie; i < nit && oe[i] < b; i++) {   /* threads of [a, b) with a lane < nv */
                        u64 S = its[i0 + i].S, v = its[i0 + i].nv, r0 = std::max(a, oe[i]) - oe[i], r1 = std::min(b, oe[i + 1]) - oe[i];
                        nv += r1 / S * v + std::min(r1 % S, v) - (r0 / S * v + std::min(r0 % S, v));
                    }
                    if (nv && (cib < 0 || (u64)cib + nv > d.icap)) { cib = (i64)d2h1(d.icount); grow_I(nv, cib); }
                    cib += nv;
                    k_cb_emit<<<nblk(b - a, 256), 256>>>(d, cb_it, cb_cl, cb_off + (size_t)(T + 1) * (nit + 1), nit, a, b, cb_F, cb_pkey, cb_plev, d_eovf, cb_miss);
                    if (cb_miss) k_cb_emit_miss<<<nblk((b - a + 31) / 32, 256), 256>>>(d, cb_it, cb_cl, cb_off + (size_t)(T + 1) * (nit + 1), nit, a, (u32)((b - a + 31) / 32),
                                                                                    cb_F, cb_pkey, cb_plev, d_eovf, cb_miss);
                }
                t_cbe += wall() - tr1;
                i0 = i1;
            }
            cb_busy = 0;
            k0 += nb; force_m = 0; kmax = ncls;
        }
        if (!cls.empty()) SYNC();
        /* process the layer in batches of curves whose L fits and whose exp entries fit */
        ensure_L();
        for (auto &rg : ranges) {
        u32 c0 = rg.first, cend = rg.second;
        u32 bmax = cend - c0;
        while (c0 < cend) {
            ensure_L(); /* an I growth during the last emit may have released it */
            u32 nb = std::min(bmax, cend - c0);
            /* size the batch so its exp entries (estimated) stay near xtarget,
             * i.e. a table that mostly lives in L2 */
            if (!force_mask) {
                /* entries scale with the curves' L sizes: take curves until the estimate reaches xtarget */
                double acc_e = 0; u32 k = 0;
                for (; k < nb; k++) { acc_e += est_nu * hnl[c0 + k] + 1; if (acc_e > xtarget && k > 0) break; }
                nb = std::max<u32>(1, k);
            }
            /* L offsets within the pool */
            u64 tot = 0; u32 k = 0; int numax = 0;
            for (; k < nb; k++) { u64 sz = hnu[c0 + k] + D - hdeg[c0 + k] + 1; if (tot + sz > d.Lcap) break; hoff[c0 + k] = tot; tot += sz; numax = std::max(numax, hnu[c0 + k]); }
            nb = k;
            if (!nb) {
                /* one curve needs more L terms than the pool holds: grow the pool and retry */
                u64 sz = hnu[c0] + D - hdeg[c0] + 1;
                if (sz <= d.Lcap || L_want >= sz) { fprintf(stderr, "cgv gpu: L pool too small\n"); return 1; }
                L_want = (u64)(sz * 1.25) + 1024;
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
            /* the last emit over this buffer must have read its entries before the seed inserts */
            if (ev_used[cur]) CK(cudaStreamWaitEvent(0, ev_emit[cur], 0));
            /* new epoch: every slot of the earlier batches counts as empty. B lives in the state's top 24 bits
             * (see x_find) and skips 0 (bundled-layer states); at the wrap both buffers are cleared once */
            if (++d.xep >= 1u << 24) {
                CK(cudaStreamSynchronize(s_emit));
                for (int q = 0; q < 2; q++) if (xtb[q]) CK(cudaMemset(xtb[q], 0, (q ? dbl_slots : xtb0_slots) * sizeof(XEnt)));
                d.xep = 1;
            }
            CK(cudaMemset(d.norder, 0, sizeof(u32)));
            {
                u32 nch = (u32)std::max(1, (numax + BL_CH - 1) / BL_CH);
                if ((u64)nb * nch <= 2 * ((u64)full_cap + 1) && nch <= 65535) {
                    u32 *cnt = (u32 *)d.work; /* free until the levels */
                    k_buildL_cnt<<<dim3(nb, nch), 256>>>(d, c0, cnt);
                    k_buildL_put<<<dim3(nb, nch), 256>>>(d, c0, cnt);
                } else {
                    k_buildL<<<nb, 256>>>(d, c0);
                    k_seed<<<nblk(nb, 256), 256>>>(d, c0, nb);
                }
            }
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
                    if (dlev) {
                        /* the rest of the batch on the device: level e (unless k_small_levels finalized it), then e+1..T */
                        int e1 = e;
                        if (so[2]) { totW += W; layW += W; if (W) k_scatter<<<nblk(W, 256), 256>>>(d, e, ent, n, W); e1 = e + 1; }
                        for (int g = e1; g <= T; g++) {
                            k_dfin<<<DL_G, 256>>>(d, g, d_dw);
                            k_dscan<<<1, DL_G>>>(d, g, d_dw);
                            k_doff<<<DL_G, 256>>>(d, g, d_dw);
                            k_dscatter<<<dl_sg, 256>>>(d, g, d_dw);
                        }
                        k_dreport<<<1, 1>>>(d, T, d_dw, d_dr);
                        u64 so5[5];
                        CK(cudaGetLastError());
                        CK(cudaMemcpy(so5, d_dr, sizeof(so5), cudaMemcpyDeviceToHost));
                        totW += so5[4]; layW += so5[4];
                        tlev += wall() - tb;
                        fin = 1; for (int q = 0; q < 4; q++) fin_so[q] = so5[q];
                        if (so5[1]) ovf = 1;
                        break;
                    }
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
                    if (e > T) { if (d2h1(d.overflow)) ovf = 1; break; }
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
                    if (e == T) { if (d2h1(d.overflow)) ovf = 1; }
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
            if (!ovf && !fin) { if (d2h1(d.overflow)) ovf = 1; } /* d2h1 syncs stream 0 only: the last emit may still run */
            u32 nent = std::min(fin ? (u32)fin_so[2] : d2h1(d.norder), d.xcap);
            if (ovf) {
                /* table full: retry with a bigger table or fewer curves */
                CK(cudaMemset(d.overflow, 0, sizeof(int)));
                SYNC(); /* the table needs no clearing: the retry has a new epoch */
                if (getenv("CGV_PROF") && in->verbose) fprintf(stderr, "    retry: layer %d batch at %u, %u curves, mask 2^%d, entries %u\n", dd, c0, nb, __builtin_ctzll(d.xmask + 1), nent);
                est_nu *= 2; /* remember that curves here are bigger than estimated */
                if (d.xmask < full_mask) {
                    /* estimate was short: same curves, 4x the table */
                    force_mask = std::min<u64>((d.xmask + 1) * 4, full_mask + 1);
                    bmax = nb;
                } else {
                    if (nb == 1) {
                        size_t f = 0, t = 0; CK(cudaMemGetInfo(&f, &t));
                        u64 n = (full_mask + 1) * 2, nby = (n - (full_mask + 1)) * (per + (u64)(D + 2) * 2);
                        if (n <= xcap_max && nby > f && cb_release()) CK(cudaMemGetInfo(&f, &t));
                        if (n > xcap_max || nby > f) { fprintf(stderr, "cgv gpu: a single curve overflows the exp table (%llu slots, no memory to grow)\n", (unsigned long long)(full_mask + 1)); return 1; }
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
                    /* while the bound shows that no growth can be needed, the earlier emits are not waited for */
                    if (iub < 0 || (u64)iub + (hi - lo) > d.icap) {
                        CK(cudaStreamSynchronize(s_emit));
                        iub = lo == 0 && icnt_ok ? (i64)*h_icnt : (i64)d2h1(d.icount);
                        grow_I(hi - lo, iub);
                    }
                    iub += hi - lo;
                    if (lo == 0) { CK(cudaEventRecord(ev_lev, 0)); CK(cudaStreamWaitEvent(s_emit, ev_lev, 0)); }
                    k_emit_range<<<nblk(hi - lo, 256), 256, 0, s_emit>>>(d, lo, hi, d_eovf);
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
            GPU_HB(dd, nc, layB);
            force_mask = 0; bmax = cend - c0;
        }
        }
        /* layer end: wait for the last emit, back to buffer 0, check the emit failure flag */
        SYNC();
        cur = 0; d.xtab = xtb[0]; d.order = ordb[0]; ev_used[0] = ev_used[1] = 0; icnt_ok = 0; iub = -1;
        { int ov = d2h1(d_eovf); if (ov) { fprintf(stderr, "cgv gpu: overflow flag %d after emit (layer %d)\n", ov, dd); return 1; } }
        if (getenv("CGV_PROF") && in->verbose)
            fprintf(stderr, "    gpu layer %d: curves %u batches %d entries %.3g W %.3g  level %.3fs layer %.3fs  (%.2g W/s)  classes %zu (%u curves hashed) build %.3fs replay %.3fs emit %.3fs\n", dd, nc, layB,
                    (double)layEnt, (double)layW, tt[T_LEVEL] - layLev0, wall() - layT0, layW / (tt[T_LEVEL] - layLev0 + 1e-9),
                    cls.size(), nrest, t_cbb - lcb0[0], t_cbr - lcb0[1], t_cbe - lcb0[2]);
        dd = dhi;
    }
    SYNC();
    cudaFree(d_next); cudaFree(d_sl); cudaFree(d_dw); cudaFree(d_dr); cudaFree(d_bsum); cudaFree(d_eovf);
    if (xtb[1]) cudaFree(xtb[1]);
    if (ordb[1]) cudaFree(ordb[1]);
    if (h_icnt) HOST_FREE(h_icnt);
    cudaStreamDestroy(s_emit); cudaEventDestroy(ev_emit[0]); cudaEventDestroy(ev_emit[1]); cudaEventDestroy(ev_lev);
    if (test_hold) cudaFree(test_hold);
    cudaFree(d_pat); cudaFree(d_bcur); cudaFree(d_bofs); cudaFree(d_bcnl);
    { void *p[] = {d.npair, d.pk, d.pv, cb_pk2, cb_pv2, cb_stmp, cb_toff, cb_pkey, cb_plev, cb_ljx, cb_F, cb_cnt, cb_clv, cb_cl, cb_it, cb_off, cb_miss};
      for (void *q : p) if (q) cudaFree(q); }
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
    GKey *okey = dalloc<GKey>(npts); int *odeg = dalloc<int>(npts); Co *oA = dalloc<Co>(npts); u32 *on = dalloc<u32>(1);
    CK(cudaMemset(on, 0, sizeof(u32)));
    k_collect<<<nblk(icap, 256), 256>>>(d, okey, odeg, oA, on);
    SYNC();
    out->n = (int)d2h1(on);
    out->key = malloc((size_t)out->n * sizeof(GKey));
    out->deg = (int *)malloc((size_t)out->n * sizeof(int));
    out->A = malloc((size_t)out->n * sizeof(Co));
    CK(cudaMemcpy(out->key, okey, (size_t)out->n * sizeof(GKey), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(out->deg, odeg, (size_t)out->n * sizeof(int), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(out->A, oA, (size_t)out->n * sizeof(Co), cudaMemcpyDeviceToHost));
    out->napplied = napplied;
    if (in->verbose) {
        fprintf(stderr, "  gpu: %lld curves applied in %d batches (%d retries), %d points\n", (long long)napplied, nbatches, nretry, out->n);
        fprintf(stderr, "  gpu: %.3g scatter threads (%.3g/s in level)\n", (double)totW, totW / (tt[T_LEVEL] + 1e-9));
        if (cb_on) fprintf(stderr, "  gpu class buffers: pairs %.2f GB, sort tmp %.2f GB, values %.2f GB, points %.2f GB, L idx %.2f GB\n",
                           cb_pcap * 24.0 / 1e9, cb_stb / 1e9, cb_Fcap * 8.0 / 1e9, cb_npcap * 24.0 / 1e9, cb_ljcap * 4.0 / 1e9);
        if (cb_on) fprintf(stderr, "  gpu classes: %lld classes, %lld curves (%lld fell back); build %.2fs (%.3g pairs), replay %.2fs, emit %.2fs\n",
                           (long long)cb_ncl, (long long)cb_ncur, (long long)cb_fail, t_cbb, (double)cb_bW, t_cbr, t_cbe);
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
