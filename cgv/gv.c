/*
 * cgv: Gopakumar-Vafa invariants of CY threefold hypersurfaces via HKTY.
 *
 * Same mathematics as cygv (github.com/ariostas/cygv), reorganized for speed:
 *
 *  - Arithmetic is done modulo ~62-bit primes (Montgomery form) instead of
 *    exact rationals. Every denominator in HKTY is a product of small integers
 *    (factorials, harmonic numbers, n^3, degrees), so the pipeline is exact mod
 *    p, and the integer GVs are recovered by CRT. NL primes are carried as
 *    lanes of one coefficient, so all of them share every lattice lookup; more
 *    passes are run until the CRT lift is stable.
 *  - Lattice points are identified by 128-bit keys (256 or 512 when needed), packed linearly so that
 *    key(a+b) = key(a)+key(b). All polynomials are sparse (key, degree, value)
 *    lists; there is no global enumeration of the Mori cone.
 *  - The fundamental period and its derivatives only live on curves with at
 *    most two negative GLSM intersections (three, on a vex stratum). That set is
 *    the union over T of the lattice points of the cones
 *    K_T = Mori cone cap {Q_r.C >= 0 for r not in T}, each the nonnegative
 *    integer combinations of its Hilbert basis (computed by normaliz in the
 *    Python driver, from any description of the Mori cone).
 *  - The h11 instanton polynomials are contracted with the grading vector w,
 *    I = sum_t w_t inst_t, and cygv's beta products collapse to one:
 *        I = c0^{-1} S2 - sum_a alpha_a V_a,
 *    S2 = sum_{a<=b} K_ab c2_ab, V_a = 1/2 sum_b K_ab alpha_b, K_ab = sum_t w_t kappa_tab.
 *    After subtracting all lower-degree curves, I[C] = deg(C) A_C with
 *        A_C = sum_{n | C} GV_{C/n} / n^3.
 *  - cygv's Li2(q^N) power series are replaced by the multicover sum above,
 *    so each curve C costs one truncated exp(C.alpha).
 *
 * Input (text, whitespace separated), written by tools/cgv_run.py (write_input):
 *   h11 ndiv maxdeg             (older files start with the tag "cgv 2", skipped)
 *   w[h11]                      grading vector (positive on the Mori cone)
 *   Q[ndiv][h11]                GLSM charges (row r = divisor r)
 *   nint   (i j k val)[nint]    triple intersection numbers
 *   ncones { n k T_1..T_k  hb[n][h11] }   Hilbert basis of each cone K_T, |T| = k in 1..3. Only the union of
 *                               their lattice points is used (cones may repeat or contain one another), except
 *                               that a cone labeled by a pair T whose quadratic form vanishes is skipped (see
 *                               build_candidates). "Mori cone" may also be a face of it (then: the curves on that face).
 *   optional sections, each a keyword then its data:
 *   lightcone  m  p[m][h11]  nf  H[nf][h11]
 *                               keep only the backward lightcones of the points p: the curves C with
 *                               p - C in the Mori cone {x : H x >= 0} for some p (H: any inequalities
 *                               describing the Mori cone). Every curve feeding into such a C is again in
 *                               the set, so its GVs are exactly those of the full computation.
 *   vex  nvex { k i_1..i_k v[h11] }   vex strata (bgv; v entries may be "p/q")
 * Output: one line per nonzero GV: "c_0 ... c_{h11-1} gv".
 *
 * Environment switches (all optional):
 *   CGV_MEM=low         low-memory mode: no curve-class patterns (CPU and GPU) and, with glibc, big allocations
 *                       mapped separately (returned to the system when freed); CGV_LOW_MEM=1 is an alias,
 *                       CGV_LOW_MEM_MB the mapping threshold (default 1)
 *   CGV_CB=0            no curve-class patterns (CPU and GPU), nothing else changed
 *   CGV_LT_LAYOUT=packed|split   level-table layout (default: from the L2 size, see lt_layout)
 *   CGV_PRIME0, CGV_RESIDUES=1, CGV_MIN_PRIMES   distributed runs / prime count (see the pass loop)
 *   CGV_PROF=1          timing and memory details;  CGV_PROGRESS_S  seconds between progress lines (30)
 *   CGV_GPU_MIN, CGV_GPU_MAXDEGS, CGV_GPU_MIN_GB, CGV_ICAP, CGV_XCAP   GPU use and table sizes (-g builds)
 *   CGV_GPU_MEM_GB      GPU memory budget;  CGV_ALLOW_DISPLAY_GPU=1  allow a GPU that drives a display
 *                       (the other GPU switches: the block in gpu.cu's gpu_extract)
 *   CGV_LSYNC, CGV_CS_NU, CGV_NO_SUPER   extraction scheduling (tuning)
 *   CGV_MAJ_CAND, CGV_MAJ_ALL           majorant build only
 *   CGV_FORCE_WIDE=256|512              test hook: run the wide-key build
 */
#define _GNU_SOURCE
#include <stdint.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <sched.h>
#include <unistd.h>
#ifdef _WIN32
#include <windows.h>   /* GetActiveProcessorCount */
#endif
#ifdef __GLIBC__
#include <malloc.h>
#endif
#ifdef __APPLE__
#include <sys/sysctl.h>
#endif
#ifdef USE_GPU
#include "gpu.h"
#endif
#ifdef CGV_MAJ
#include <fenv.h>
#include <math.h> /* before the names below, which glibc also declares */
#define fadd cgv_fadd
#define fsub cgv_fsub
#define fmul cgv_fmul
#endif

typedef uint64_t u64;
typedef int64_t i64;
typedef unsigned __int128 u128;
typedef __int128 i128;
/* lattice keys: 128 bits, or (CGV_WIDE = 256 or 512 builds, used when the coordinates need more) wider */
#ifdef CGV_WIDE
/* KW little-endian 64-bit words; only the operations the key paths use */
#define KEYBITS CGV_WIDE
#define KW (CGV_WIDE / 64)
typedef struct { u64 w[KW]; } KEY;
static inline KEY k_fromi(i64 v) { KEY r; r.w[0] = (u64)v; for (int i = 1; i < KW; i++) r.w[i] = v < 0 ? ~(u64)0 : 0; return r; }
static inline KEY k_fromu(u64 v) { KEY r; r.w[0] = v; for (int i = 1; i < KW; i++) r.w[i] = 0; return r; }
static inline KEY k_add(KEY a, KEY b) { KEY r; u64 c = 0; for (int i = 0; i < KW; i++) { u128 t = (u128)a.w[i] + b.w[i] + c; r.w[i] = (u64)t; c = (u64)(t >> 64); } return r; }
static inline KEY k_sub(KEY a, KEY b) { KEY r; u64 c = 0; for (int i = 0; i < KW; i++) { u128 t = (u128)a.w[i] - b.w[i] - c; r.w[i] = (u64)t; c = (u64)(t >> 64) & 1; } return r; }
static inline KEY k_muli(KEY a, u64 j) { KEY r; u64 c = 0; for (int i = 0; i < KW; i++) { u128 t = (u128)a.w[i] * j + c; r.w[i] = (u64)t; c = (u64)(t >> 64); } return r; }
static inline int k_eq(KEY a, KEY b) { u64 x = 0; for (int i = 0; i < KW; i++) x |= a.w[i] ^ b.w[i]; return !x; }
static inline int k_le(KEY a, KEY b) { for (int i = KW - 1; i >= 0; i--) if (a.w[i] != b.w[i]) return a.w[i] < b.w[i]; return 1; }
static inline KEY k_shl(KEY a, int s) {
    KEY r; int q = s >> 6, b = s & 63;
    for (int i = KW - 1; i >= 0; i--) {
        u64 v = i - q >= 0 ? a.w[i - q] << b : 0;
        if (b && i - q - 1 >= 0) v |= a.w[i - q - 1] >> (64 - b);
        r.w[i] = v;
    }
    return r;
}
static inline KEY k_shr(KEY a, int s) {   /* logical */
    KEY r; int q = s >> 6, b = s & 63;
    for (int i = 0; i < KW; i++) {
        u64 v = i + q < KW ? a.w[i + q] >> b : 0;
        if (b && i + q + 1 < KW) v |= a.w[i + q + 1] << (64 - b);
        r.w[i] = v;
    }
    return r;
}
static inline u128 k_lo128(KEY a) { return (u128)a.w[1] << 64 | a.w[0]; }
#define K_ZERO k_fromu(0)
#else
typedef u128 KEY; typedef i128 SKEY;
#define KEYBITS 128
/* the 128-bit build uses the native operators (these expand to exactly the former expressions) */
#define k_add(a, b) ((a) + (b))
#define k_sub(a, b) ((a) - (b))
#define k_muli(a, j) ((a) * (KEY)(j))
#define k_eq(a, b) ((a) == (b))
#define k_le(a, b) ((a) <= (b))
#define k_lo128(a) (a)
#define K_ZERO ((KEY)0)
#endif

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}
static int verbose = 1;
static double rss_gb(void) { /* current resident set, GB */
    long pages = 0, rss = 0;
    FILE *f = fopen("/proc/self/statm", "r");
    if (f) { if (fscanf(f, "%ld %ld", &pages, &rss) != 2) rss = 0; fclose(f); }
    return rss * 4096.0 / 1e9;
}
#define MEMLOG(tag) do { if (getenv("CGV_PROF")) LOG("    [mem] %-28s %.2f GB\n", tag, rss_gb()); } while (0)
#define LOG(...) do { if (verbose) { fprintf(stderr, __VA_ARGS__); } } while (0)
static void die(const char *msg) { fprintf(stderr, "cgv: %s\n", msg); exit(1); }
static void *xmalloc(size_t n) { void *p = malloc(n ? n : 1); if (!p) die("out of memory"); return p; }
static void *xcalloc(size_t n, size_t s) { void *p = calloc(n ? n : 1, s ? s : 1); if (!p) die("out of memory"); return p; }
static void *xrealloc(void *p, size_t n) { p = realloc(p, n ? n : 1); if (!p) die("out of memory"); return p; }
/* low-memory mode: CGV_MEM=low, or its alias CGV_LOW_MEM=1 (see low_mem_setup) */
static int mem_low(void) {
    const char *m = getenv("CGV_MEM"), *v = getenv("CGV_LOW_MEM");
    return (m && !strcmp(m, "low")) || (v && *v && strcmp(v, "0"));
}
/* curve-class patterns allowed (CPU and GPU): not CGV_CB=0 (empty = unset), not low-memory mode */
static inline __attribute__((unused)) int cb_allowed(void) { const char *v = getenv("CGV_CB"); return (!v || !*v || atoi(v)) && !mem_low(); }
typedef struct { int *v; int n, cap; } IVec;
static void iv_push(IVec *a, int x) {
    if (a->n == a->cap) { a->cap = a->cap ? 2 * a->cap : 16; a->v = xrealloc(a->v, a->cap * sizeof(int)); }
    a->v[a->n++] = x;
}

/* ------------------------------------------------------------------------- */
/* Montgomery arithmetic modulo a prime p < 2^62                             */
/* ------------------------------------------------------------------------- */
typedef struct { u64 p, pinv, r2; } Field; /* pinv = -p^{-1} mod 2^64, r2 = 2^128 mod p */

static inline u64 mredc(const Field *F, u128 t) {
    u64 m = (u64)t * F->pinv;
    u64 u = (u64)((t + (u128)m * F->p) >> 64);
    return u >= F->p ? u - F->p : u;
}
static inline u64 fmul(const Field *F, u64 a, u64 b) { return mredc(F, (u128)a * b); }
static inline u64 fadd(const Field *F, u64 a, u64 b) { u64 c = a + b; return c >= F->p ? c - F->p : c; }
static inline u64 fsub(const Field *F, u64 a, u64 b) { return a >= b ? a - b : a + F->p - b; }
static inline u64 fneg(const Field *F, u64 a) { return a ? F->p - a : 0; }
static inline u64 fto(const Field *F, u64 x) { return fmul(F, x % F->p, F->r2); }
static inline u64 ffrom(const Field *F, u64 a) { return mredc(F, a); }
static __attribute__((unused)) u64 fint(const Field *F, i64 v) {
    u64 x = v >= 0 ? (u64)v % F->p : (F->p - (u64)(-(v + 1)) % F->p - 1);
    return fto(F, x);
}
static u64 fpow(const Field *F, u64 a, u64 e) {
    u64 r = fto(F, 1);
    while (e) { if (e & 1) r = fmul(F, r, a); a = fmul(F, a, a); e >>= 1; }
    return r;
}
static __attribute__((unused)) u64 finv(const Field *F, u64 a) { return fpow(F, a, F->p - 2); }

static u64 mulmod_plain(u64 a, u64 b, u64 m) { return (u64)((u128)a * b % m); }
static u64 powmod_plain(u64 a, u64 e, u64 m) {
    u64 r = 1; a %= m;
    while (e) { if (e & 1) r = mulmod_plain(r, a, m); a = mulmod_plain(a, a, m); e >>= 1; }
    return r;
}
static int is_prime(u64 n) {
    if (n < 2) return 0;
    static const u64 small[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int i = 0; i < 12; i++) if (n % small[i] == 0) return n == small[i];
    u64 d = n - 1; int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    for (int i = 0; i < 12; i++) {
        u64 x = powmod_plain(small[i], d, n);
        if (x == 1 || x == n - 1) continue;
        int ok = 0;
        for (int r = 1; r < s; r++) { x = mulmod_plain(x, x, n); if (x == n - 1) { ok = 1; break; } }
        if (!ok) return 0;
    }
    return 1;
}
static void field_init(Field *F, u64 p) {
    F->p = p;
    u64 inv = 1; /* Newton iteration for p^{-1} mod 2^64 */
    for (int i = 0; i < 7; i++) inv *= 2 - p * inv;
    F->pinv = (u64)0 - inv;
    u128 r = ((u128)1 << 64) % p;
    F->r2 = (u64)((r * r) % p);
}
/* The k-th prime below 2^62, counting down. */
static u64 nth_prime(int k) {
    u64 n = ((u64)1 << 62) - 1;
    for (;; n -= 2) if (is_prime(n) && k-- == 0) return n;
}

/* ------------------------------------------------------------------------- */
/* Multi-prime coefficients                                                  */
/* ------------------------------------------------------------------------- */
#ifndef NL
#define NL 2
#endif
#define MAX_PASSES 8
#ifdef CGV_MAJ
/* Majorant build: the same pipeline over nonnegative reals, rounded upward.
 * Subtraction and negation become addition and identity, so every value is an
 * upper bound on the absolute value of the true (integer/rational) one; the
 * run yields a proven bound on |GV| per curve, which fixes how many primes the
 * exact run needs. co_inv is only applied to exact values (n, n^3, c0's
 * constant term 1). */
typedef struct { long double v[NL]; } Co;
static Field FL[NL];
static inline Co co_zero(void) { Co r = {{0}}; return r; }
static inline int co_isz(Co a) { return a.v[0] == 0; }
static inline Co co_mul(Co a, Co b) { Co r; r.v[0] = a.v[0] * b.v[0]; return r; }
static inline Co co_add(Co a, Co b) { Co r; r.v[0] = a.v[0] + b.v[0]; return r; }
static inline Co co_sub(Co a, Co b) { return co_add(a, b); }
static inline Co co_neg(Co a) { return a; }
static inline Co co_int(i64 x) { Co r; r.v[0] = (long double)(x < 0 ? -x : x); return r; }
static inline Co co_inv(Co a) { Co r; r.v[0] = 1.0L / a.v[0]; return r; }
static inline void co_fma(Co *a, Co b, Co c) { a->v[0] += b.v[0] * c.v[0]; }
#else
typedef struct { u64 v[NL]; } Co;
static Field FL[NL];

static inline Co co_zero(void) { Co r; for (int l = 0; l < NL; l++) r.v[l] = 0; return r; }
static inline int co_isz(Co a) { u64 x = 0; for (int l = 0; l < NL; l++) x |= a.v[l]; return x == 0; }
static inline Co co_mul(Co a, Co b) { Co r; for (int l = 0; l < NL; l++) r.v[l] = fmul(&FL[l], a.v[l], b.v[l]); return r; }
static inline Co co_add(Co a, Co b) { Co r; for (int l = 0; l < NL; l++) r.v[l] = fadd(&FL[l], a.v[l], b.v[l]); return r; }
static inline Co co_sub(Co a, Co b) { Co r; for (int l = 0; l < NL; l++) r.v[l] = fsub(&FL[l], a.v[l], b.v[l]); return r; }
static inline Co co_neg(Co a) { Co r; for (int l = 0; l < NL; l++) r.v[l] = fneg(&FL[l], a.v[l]); return r; }
static inline Co co_int(i64 x) { Co r; for (int l = 0; l < NL; l++) r.v[l] = fint(&FL[l], x); return r; }
static inline Co co_inv(Co a) { Co r; for (int l = 0; l < NL; l++) r.v[l] = finv(&FL[l], a.v[l]); return r; }
static inline void co_fma(Co *a, Co b, Co c) { *a = co_add(*a, co_mul(b, c)); }
#endif

/* ------------------------------------------------------------------------- */
/* Problem data and lattice keys                                             */
/* ------------------------------------------------------------------------- */
static int h11, ndiv, maxdeg;
/* Curves are kept in HC coordinates: those of a basis LB (HC rows of length h11) of the lattice spanned by the cone
 * generators. HC = h11 and LB = NULL (identity) unless the generators span a lower-rank sublattice (e.g. the curves
 * on a face of the Mori cone): then keys pack HC coordinates, not h11. W, Q, Q0 are the curve-side (HC-column)
 * grading and charges; Wd, Qd, Q0d the divisor-side ones (h11 columns), equal to W, Q, Q0 when LB = NULL. */
static int HC; static i64 *LB;
static int *W, *Q, *Q0, *Wd, *Qd, *Q0d;
static int nlight, *light_pts, nlight_h, *light_h;   /* optional: the backward lightcones of light_pts (Mori cone: light_h x >= 0) */
static long long *light_hp, *light_hc;   /* H p for each point p (computed on first use); scratch for H C */
static int nint, (*intnums)[4];
static int ncones, *cone_n, **cone_g, (*cone_T)[3];
/* bgv: vex strata (cones S of the fan whose rays share no facet). |S| = 2: [V_S]_X = sum_a c_a J_a;
 * |S| = 3: d_t = J_t . V_S. Stored as exact rationals, converted per pass (per set of primes). */
static int nvex, *vex_k, (*vex_S)[3];
static i64 *vex_num, *vex_den;   /* [nvex][h11] */
static Co *vex_val;   /* [nvex][h11], filled by tables_init */
static int keybits;
/* KBV: per-coordinate field widths KBW[t] at offsets KBO[t] (used only when the uniform KEYBITS/HC bits do not fit) */
static int KBV, KBW[64], KBO[64];
#ifndef CGV_WIDE
static KEY fieldmask;
#endif
static KEY LT_XK;   /* level-table key marker: pack(-e_t) with W_t > 0 has negative degree, so it is never a level key */
static u128 LT_HXK;  /* lt_h(LT_XK): tables store lt_h(k) ^ LT_HXK, zero only for the marker */

static inline int dot(const int *a, const int *b) { int s = 0; for (int t = 0; t < HC; t++) s += a[t] * b[t]; return s; }
static inline long long dotll(const int *a, const int *b) { long long s = 0; for (int t = 0; t < HC; t++) s += (long long)a[t] * b[t]; return s; }
#ifdef CGV_WIDE
static inline KEY pack(const int *v) {
    KEY k = K_ZERO;
    for (int t = HC - 1; t >= 0; t--) k = k_add(k_shl(k, KBV ? KBW[t] : keybits), k_fromi(v[t]));
    return k;
}
static inline void unpack(KEY k, int *v) {
    for (int t = 0; t < HC; t++) {
        int w = KBV ? KBW[t] : keybits;   /* <= 40 */
        i64 f = (i64)(k.w[0] & (((u64)1 << w) - 1));
        if (f >= ((i64)1 << (w - 1))) f -= (i64)1 << w;
        v[t] = (int)f;
        k = k_shr(k_sub(k, k_fromi(f)), w);
    }
}
static inline u64 hash128(KEY k) {
    static const u64 C[8] = {1, 0x9E3779B97F4A7C15ull, 0xC2B2AE3D27D4EB4Full, 0x165667B19E3779F9ull,
                             0x27D4EB2F165667C5ull, 0x94D049BB133111EBull, 0xBF58476D1CE4E5B9ull, 0xFF51AFD7ED558CCDull};
    u64 x = 0;
    for (int i = 0; i < KW; i++) x ^= k.w[i] * C[i];
    x ^= x >> 31; x *= 0xD6E8FEB86659FD93ull; x ^= x >> 32;
    return x;
}
#else
static inline KEY pack(const int *v) {
    KEY k = 0;
    if (KBV) { for (int t = HC - 1; t >= 0; t--) k = (k << KBW[t]) + (KEY)(SKEY)v[t]; return k; }
    for (int t = HC - 1; t >= 0; t--) k = (k << keybits) + (KEY)(SKEY)v[t];
    return k;
}
static inline void unpack(KEY k, int *v) {
    if (KBV) {
        for (int t = 0; t < HC; t++) {
            int w = KBW[t];
            i64 f = (i64)(u64)(k & (((KEY)1 << w) - 1));
            if (f >= ((i64)1 << (w - 1))) f -= (i64)1 << w;
            v[t] = (int)f;
            k = (k - (KEY)(SKEY)f) >> w;
        }
        return;
    }
    for (int t = 0; t < HC; t++) {
        i64 f = (i64)(u64)(k & fieldmask);
        if (f >= ((i64)1 << (keybits - 1))) f -= (i64)1 << keybits;
        v[t] = (int)f;
        k = (k - (KEY)(SKEY)f) >> keybits;
    }
}
static inline u64 hash128(KEY k) {
    u64 x = (u64)k ^ ((u64)(k >> 64) * 0x9E3779B97F4A7C15ull);
    x ^= x >> 31; x *= 0xD6E8FEB86659FD93ull; x ^= x >> 32;
    return x;
}
#endif
/* level tables hash multiplicatively mod 2^128: H(k) = k M with M odd, so H(a + b) = H(a) + H(b) and H is invertible.
 * The tables store H(k) in place of k: the scatter loop forms H(kb + lkey[j]) = H(kb) + H(lkey[j]) with one add,
 * the slot is the top bits, and k = H(k) M^{-1} is recovered once per finalized slot */
#define LT_HM (((u128)0x9E3779B97F4A7C15ull << 64) | 0xD6E8FEB86659FD93ull)
static inline u128 lt_h(u128 k) { return k * LT_HM; }
static u128 LT_HMI;   /* M^{-1} mod 2^128 */
#define LT_SLOT(t, H) ((u64)((H) >> 64) >> (t)->sh)
/* 64-bit level-table keys: a level table holds offsets of one degree, all inside a box lo_t <= x_t <= hi_t (bounded
 * from the L terms), so phi(x) = sum_{t != t0} x_t 2^off_t mod 2^64 (fields of b_t bits, x_t0 fixed by W.x = degree) is
 * exact when the fields fit in 63 bits. phi is linear, so H64 = phi M mod 2^64 keeps the one-add scatter; it is kept in
 * the top half of the u128 hash (low half 0), so the slot and the scatter are as in the 128-bit fallback (LT_KS = 1). */
static int LT_PK;   /* level tables as packed records (lt_layout) */
static int LT_KS, L64_nf, L64_nlo, L64_t0, L64_p0, L64_w0;
/* field of coordinate t: bits M = m << off of phi, packed at bit ps >= off. Fields 0..nlo-1 move by d = ps - off < 64
 * (a 64x64 -> 128 multiply by 2^d), the others only within the high word */
#ifdef CGV_WIDE
static struct { int t, off, d, w; u64 M, mul; int ps; } L64_f[64];   /* ps: key position of the field */
#else
static struct { int t, off, d, w; u64 M, mul; } L64_f[64];
#endif
static u64 L64_plo, LT_HMI64; static KEY L64_klo; static i64 L64_wlo;
#define LT_HM64 0xD6E8FEB86659FD93ull
static u128 *ULH;     /* level-table hash of UKEY[u] */
static inline u64 l64_phi(const int *x) { u64 p = 0; for (int f = 0; f < L64_nf; f++) p += (u64)(i64)x[L64_f[f].t] << L64_f[f].off; return p; }
#ifdef CGV_WIDE
/* wide keys, level box of 64..127 bits: phi128(x) = sum_{t != t0} x_t 2^off_t mod 2^128, hashed H = phi128 M (LT_KS = 1) */
static struct { int t, off, b, w, ps; } L128_f[64];
static int L128_nf, L128_t0, L128_p0, L128_w0; static u128 L128_plo; static KEY L128_klo; static i64 L128_wlo;
static inline u128 l128_phi(const int *x) { u128 p = 0; for (int f = 0; f < L128_nf; f++) p += (u128)(i128)x[L128_f[f].t] << L128_f[f].off; return p; }
#endif
static inline u128 lt_hk(KEY k) {
#ifdef CGV_WIDE
    if (LT_KS) { int x[64]; unpack(k, x); return l128_phi(x) * LT_HM; }
#else
    if (LT_KS) return lt_h(k);
#endif
    int x[64]; unpack(k, x);
    return (u128)(l64_phi(x) * LT_HM64) << 64;
}
/* offset from its hash, in a level table of degree e: pack is linear, so the fields add straight into the key */
/* Level keys in the coordinates c of an LLL basis LBL_k of the lattice of the L terms (lt_setup uses them when the box in
 * key coordinates needs more bits: the L terms often span a lower-rank lattice): phi(c) packs the c_j like the key
 * coordinates above, and a key is decoded as sum_j c_j key(basis_j) */
static int LT_LB, LT_DEC, LBL_nf, LBL_t0, LBL_w0; static i64 LBL_wlo; static u128 LBL_plo;
static struct { int j, off, b, w; } LBL_f[64];
static KEY LBL_k[64], LBL_k0;
static inline KEY k_addmul(KEY k, KEY b, i64 c) { return c >= 0 ? k_add(k, k_muli(b, (u64)c)) : k_sub(k, k_muli(b, -(u64)c)); }
static __attribute__((noinline)) KEY lt_key_lb(u128 H, int e) {   /* (off the hot path: keeps lt_key small) */
    u128 y = LT_KS ? H * LT_HMI - LBL_plo : (u128)((u64)(H >> 64) * LT_HMI64 - (u64)LBL_plo);
    KEY k = LBL_k0; u64 ur = (u64)(e - LBL_wlo);   /* mod 2^64: the final value fits */
    for (int f = 0; f < LBL_nf; f++) {
        u64 v = (u64)(y >> LBL_f[f].off) & (((u64)1 << LBL_f[f].b) - 1);
        k = k_add(k, k_muli(LBL_k[LBL_f[f].j], v)); ur -= (u64)(i64)LBL_f[f].w * v;
    }
    if (LBL_t0 >= 0) { i64 r = (i64)ur; k = k_addmul(k, LBL_k[LBL_t0], LBL_w0 == 1 ? r : LBL_w0 == -1 ? -r : r / LBL_w0); }
    return k;
}
static inline KEY lt_key(u128 H, int e) {
#ifdef CGV_WIDE
    if (LT_LB) return lt_key_lb(H, e);
    if (LT_KS) {
        u128 y = H * LT_HMI - L128_plo; KEY k = L128_klo; u64 ur = (u64)(e - L128_wlo);   /* mod 2^64: the final value fits */
        for (int f = 0; f < L128_nf; f++) {
            u64 v = (u64)(y >> L128_f[f].off) & (((u64)1 << L128_f[f].b) - 1);
            k = k_add(k, k_shl(k_fromu(v), L128_f[f].ps)); ur -= (u64)(i64)L128_f[f].w * v;
        }
        i64 r = (i64)ur;
        if (L128_t0 >= 0) k = k_add(k, k_shl(k_fromi(L128_w0 == 1 ? r : L128_w0 == -1 ? -r : r / L128_w0), L128_p0));
        return k;
    }
    {   /* fields shifted to their key positions (anywhere in the wide key) */
        u64 y = (u64)(H >> 64) * LT_HMI64 - L64_plo; KEY k = L64_klo;
        for (int f = 0; f < L64_nf; f++) k = k_add(k, k_shl(k_fromu((y & L64_f[f].M) >> L64_f[f].off), L64_f[f].ps));
        if (L64_t0 >= 0) {
            u64 ur = (u64)(e - L64_wlo);
            for (int f = 0; f < L64_nf; f++) ur -= (u64)(i64)L64_f[f].w * ((y & L64_f[f].M) >> L64_f[f].off);
            i64 r = (i64)ur;
            k = k_add(k, k_shl(k_fromi(L64_w0 == 1 ? r : L64_w0 == -1 ? -r : r / L64_w0), L64_p0));
        }
        return k;
    }
#else
    if (LT_DEC) return LT_DEC == 1 ? H * LT_HMI : lt_key_lb(H, e);   /* LT_DEC = LT_KS | 2 LT_LB: one test on the common path */
    u64 y = (u64)(H >> 64) * LT_HMI64 - L64_plo, hi = 0; KEY k = L64_klo;
    for (int f = 0; f < L64_nlo; f++) k += (KEY)(u128)(y & L64_f[f].M) * L64_f[f].mul;
    for (int f = L64_nlo; f < L64_nf; f++) hi += (y & L64_f[f].M) << L64_f[f].d;
    k += (KEY)hi << 64;
    if (L64_t0 >= 0) {
        u64 ur = (u64)(e - L64_wlo);   /* mod 2^64: the final value fits */
        for (int f = 0; f < L64_nf; f++) ur -= (u64)(i64)L64_f[f].w * ((y & L64_f[f].M) >> L64_f[f].off);
        i64 r = (i64)ur;
        k += (KEY)(SKEY)(L64_w0 == 1 ? r : L64_w0 == -1 ? -r : r / L64_w0) << L64_p0;
    }
    return k;
#endif
}

static void read_input(FILE *f) {
#define RD(x) do { if (fscanf(f, "%d", &(x)) != 1) die("bad input"); } while (0)
    char tok[64], *e;
    int ver;
    if (fscanf(f, "%63s", tok) != 1) die("empty input");
    if (!strcmp(tok, "cgv")) {   /* files from before 2026-10 start with the tag "cgv 2" */
        if (fscanf(f, "%d", &ver) != 1 || ver != 2) die("unknown input format: write it with tools/cgv_run.py");
        RD(h11);
    } else {
        long v = strtol(tok, &e, 10);
        if (*e || e == tok || v < 1 || v > INT_MAX) die("bad input: write it with tools/cgv_run.py");
        h11 = (int)v;
    }
    RD(ndiv); RD(maxdeg);
#ifdef CGV_MAJ
    if (h11 < 1 || h11 > 64 || ndiv > 64) die("h11 or ndiv > 64 not supported");
#else
    if (h11 < 1 || ndiv < 1) die("bad input");
#endif
    W = xmalloc(h11 * sizeof(int));
    for (int a = 0; a < h11; a++) RD(W[a]);
    Q = xmalloc((size_t)ndiv * h11 * sizeof(int));
    for (int i = 0; i < ndiv * h11; i++) RD(Q[i]);
    RD(nint);
    intnums = xmalloc((size_t)nint * sizeof(*intnums));
    for (int i = 0; i < nint; i++) for (int j = 0; j < 4; j++) RD(intnums[i][j]);
    RD(ncones);
    if (ncones < 1) die("no K_T cones in the input");
    cone_n = xmalloc(ncones * sizeof(int));
    cone_g = xmalloc(ncones * sizeof(int *));
    cone_T = xmalloc(ncones * sizeof(*cone_T));
    for (int c = 0; c < ncones; c++) {
        int k;
        RD(cone_n[c]); RD(k);
        if (k < 1 || k > 3) die("a cone K_T needs 1 to 3 divisors in T");
        cone_T[c][0] = cone_T[c][1] = cone_T[c][2] = -1;
        for (int j = 0; j < k; j++) RD(cone_T[c][j]);
        cone_g[c] = xmalloc((size_t)(cone_n[c] ? cone_n[c] : 1) * h11 * sizeof(int));
        for (int i = 0; i < cone_n[c] * h11; i++) RD(cone_g[c][i]);
    }
    /* optional sections, each a keyword then its data */
    nvex = 0; nlight = 0;
    while (fscanf(f, "%63s", tok) == 1) {
        if (!strcmp(tok, "lightcone")) {
            RD(nlight);
            light_pts = xmalloc((size_t)(nlight ? nlight : 1) * h11 * sizeof(int));
            for (int i = 0; i < nlight * h11; i++) RD(light_pts[i]);
            RD(nlight_h);
            light_h = xmalloc((size_t)(nlight_h ? nlight_h : 1) * h11 * sizeof(int));
            for (int i = 0; i < nlight_h * h11; i++) RD(light_h[i]);
            if (!nlight) die("lightcone section with no points");
            continue;
        }
        if (strcmp(tok, "vex")) die("unknown input section");
        /* bgv: vex strata:  nvex  { k  i_1 .. i_k  v_1 .. v_h11 }  (k = 3: d_t; "p/q" allowed) */
        RD(nvex);
        if (nvex > 0) {
            vex_k = xmalloc(nvex * sizeof(int)); vex_S = xmalloc(nvex * sizeof(*vex_S));
            vex_num = xmalloc((size_t)nvex * h11 * sizeof(i64)); vex_den = xmalloc((size_t)nvex * h11 * sizeof(i64));
            for (int v = 0; v < nvex; v++) {
                RD(vex_k[v]);
                if (vex_k[v] == 2) die("vex 2-cone in the input: impossible for a fine fan on a reflexive polytope (MacFadden-Sheridan, arXiv:2512.14817, Prop. 5); check the fan");
                if (vex_k[v] != 3) die("vex stratum must have 3 divisors");
                vex_S[v][2] = -1;
                for (int j = 0; j < vex_k[v]; j++) RD(vex_S[v][j]);
                for (int a = 0; a < h11; a++) {
                    char vt[64]; long long nu, de = 1;
                    if (fscanf(f, "%63s", vt) != 1) die("bad vex stratum");
                    if (sscanf(vt, "%lld/%lld", &nu, &de) < 1 || de <= 0) die("bad vex rational");
                    vex_num[(size_t)v * h11 + a] = nu; vex_den[(size_t)v * h11 + a] = de;
                }
            }
        }
    }
#undef RD
    Q0 = xcalloc(h11, sizeof(int));
    for (int r = 0; r < ndiv; r++) for (int a = 0; a < h11; a++) Q0[a] += Q[r * h11 + a];
    Wd = W; Qd = Q; Q0d = Q0; HC = h11; LB = NULL;
}

/* Integer lattices spanned by vectors of n coordinates: a row-echelon basis by gcd row operations (exact, i64 with overflow
 * checks), coordinates in it, and LLL reduction with the matching coordinate change. Used for the curve lattice
 * (setup_lattice) and for the lattice of the L terms (lt_setup). */
/* The helpers return 0, or -1 with the reason in LAT_ERR: setup_lattice dies on it, lt_setup_lb falls back. */
static const char *LAT_ERR;
#define OVF(x) do { if (x) goto ovf; } while (0)
typedef struct { int n, rk; i64 *bas, *v; int *piv; } Lat;
static void lat_init(Lat *L, int n) {
    L->n = n; L->rk = 0; L->bas = xcalloc((size_t)n * n, sizeof(i64)); L->v = xmalloc(n * sizeof(i64)); L->piv = xmalloc(n * sizeof(int));
    for (int t = 0; t < n; t++) L->piv[t] = -1;
}
static void lat_free(Lat *L) { free(L->bas); free(L->v); free(L->piv); }
static int lat_add(Lat *L, const int *g) {
    int n = L->n; i64 *v = L->v;
    for (int t = 0; t < n; t++) v[t] = g[t];
    for (int col = 0; col < n; col++) {
        if (!v[col]) continue;
        if (L->piv[col] < 0) { L->piv[col] = L->rk; memcpy(L->bas + (size_t)L->rk * n, v, n * sizeof(i64)); L->rk++; break; }
        i64 *b = L->bas + (size_t)L->piv[col] * n, a0 = b[col], a1 = v[col];
        /* extended gcd: x a0 + y a1 = gg */
        i64 r0 = a0, r1 = a1, x0 = 1, x1 = 0, y0 = 0, y1 = 1;
        while (r1) { i64 q = r0 / r1, tt; tt = r0 - q * r1; r0 = r1; r1 = tt; tt = x0 - q * x1; x0 = x1; x1 = tt; tt = y0 - q * y1; y0 = y1; y1 = tt; }
        if (r0 < 0) { r0 = -r0; x0 = -x0; y0 = -y0; }
        i64 u0 = a0 / r0, u1 = a1 / r0;
        for (int t = col; t < n; t++) {
            i64 p, q, nb, nv;
            OVF(__builtin_mul_overflow(x0, b[t], &p)); OVF(__builtin_mul_overflow(y0, v[t], &q)); OVF(__builtin_add_overflow(p, q, &nb));
            OVF(__builtin_mul_overflow(u0, v[t], &p)); OVF(__builtin_mul_overflow(u1, b[t], &q)); OVF(__builtin_sub_overflow(p, q, &nv));
            b[t] = nb; v[t] = nv;
        }
    }
    return 0;
ovf:
    LAT_ERR = "lattice basis: integer overflow";
    return -1;
}
/* the basis rows in pivot order (rk x n), entries above each pivot reduced (smaller coordinates); pivot columns in pc */
static int lat_finish(Lat *L, i64 **Bp, int **pcp) {
    int n = L->n, rk = L->rk;
    i64 *B = xmalloc((size_t)(rk ? rk : 1) * n * sizeof(i64));
    int *pc = xmalloc((rk ? rk : 1) * sizeof(int));
    for (int col = 0, j = 0; col < n; col++) if (L->piv[col] >= 0) { memcpy(B + (size_t)j * n, L->bas + (size_t)L->piv[col] * n, n * sizeof(i64)); pc[j++] = col; }
    for (int j = 0; j < rk; j++) {
        i64 *b = B + (size_t)j * n;
        if (b[pc[j]] < 0) for (int t = 0; t < n; t++) b[t] = -b[t];
        for (int i = 0; i < j; i++) {
            i64 *a = B + (size_t)i * n, q = a[pc[j]] / b[pc[j]], p;
            if (a[pc[j]] - q * b[pc[j]] < 0) q--;
            if (2 * (a[pc[j]] - q * b[pc[j]]) > b[pc[j]]) q++;   /* symmetric remainder */
            if (q) for (int t = 0; t < n; t++) { OVF(__builtin_mul_overflow(q, b[t], &p)); OVF(__builtin_sub_overflow(a[t], p, &a[t])); }
        }
    }
    *Bp = B; *pcp = pc;
    return 0;
ovf:
    free(B); free(pc);
    LAT_ERR = "lattice basis: integer overflow";
    return -1;
}
/* coordinates c of a lattice vector g = sum_j c_j B_j in an echelon basis (one pivot at a time); v: n scratch */
static int lat_coords(const i64 *B, const int *pc, int k, int n, const int *g, int *c, i64 *v) {
    for (int t = 0; t < n; t++) v[t] = g[t];
    for (int j = 0; j < k; j++) {
        const i64 *b = B + (size_t)j * n; i64 q = v[pc[j]] / b[pc[j]], p;
        if (q * b[pc[j]] != v[pc[j]]) { LAT_ERR = "lattice: vector outside the basis"; return -1; }
        if (q) for (int t = pc[j]; t < n; t++) { OVF(__builtin_mul_overflow(q, b[t], &p)); OVF(__builtin_sub_overflow(v[t], p, &v[t])); }
        if (q > 0x3fffffff || q < -0x3fffffff) { LAT_ERR = "lattice: coordinates too large"; return -1; }
        c[j] = (int)q;
    }
    for (int t = 0; t < n; t++) if (v[t]) { LAT_ERR = "lattice: vector outside the basis"; return -1; }
    return 0;
ovf:
    LAT_ERR = "lattice basis: integer overflow";
    return -1;
}
/* LLL-reduce the k rows of B (short, near-orthogonal rows give small coordinates, i.e. fewer key bits); T (k x k) maps
 * the old coordinates to the new ones: c' = c T (b_i -= q b_j: column j of T += q column i; swap: swap columns) */
static int lat_lll(i64 *B, int k, int n, i64 **Tp) {
    i64 *T = xcalloc((size_t)k * k, sizeof(i64));
    for (int j = 0; j < k; j++) T[(size_t)j * k + j] = 1;
    /* double, not long double: the same reduced basis on every platform (long double is 64-bit on arm64) */
    double *bs = xmalloc((size_t)k * n * sizeof(double)), *nb = xmalloc(k * sizeof(double)), *mu = xmalloc((size_t)k * k * sizeof(double));
#define GS() do { for (int i = 0; i < k; i++) { \
        for (int t = 0; t < n; t++) bs[(size_t)i * n + t] = B[(size_t)i * n + t]; \
        for (int j = 0; j < i; j++) { double d = 0; \
            for (int t = 0; t < n; t++) d += (double)B[(size_t)i * n + t] * bs[(size_t)j * n + t]; \
            mu[(size_t)i * k + j] = d / nb[j]; \
            for (int t = 0; t < n; t++) bs[(size_t)i * n + t] -= mu[(size_t)i * k + j] * bs[(size_t)j * n + t]; } \
        nb[i] = 0; for (int t = 0; t < n; t++) nb[i] += bs[(size_t)i * n + t] * bs[(size_t)i * n + t]; } } while (0)
    GS();
    int i = 1, it = 0;
    while (i < k && it++ < 100000) {
        for (int j = i - 1; j >= 0; j--) {
            double m = mu[(size_t)i * k + j];
            if (m > 0.5 || m < -0.5) {
                OVF(!(m < 4e18 && m > -4e18));   /* (also NaN) */
                i64 q = (i64)(m >= 0 ? m + 0.5 : m - 0.5), p;
                for (int t = 0; t < n; t++) { OVF(__builtin_mul_overflow(q, B[(size_t)j * n + t], &p)); OVF(__builtin_sub_overflow(B[(size_t)i * n + t], p, &B[(size_t)i * n + t])); }
                for (int r = 0; r < k; r++) { OVF(__builtin_mul_overflow(q, T[(size_t)r * k + i], &p)); OVF(__builtin_add_overflow(T[(size_t)r * k + j], p, &T[(size_t)r * k + j])); }
                GS();
            }
        }
        double m = mu[(size_t)i * k + i - 1];
        if (nb[i] < (0.99 - m * m) * nb[i - 1]) {
            for (int t = 0; t < n; t++) { i64 x = B[(size_t)i * n + t]; B[(size_t)i * n + t] = B[(size_t)(i - 1) * n + t]; B[(size_t)(i - 1) * n + t] = x; }
            for (int r = 0; r < k; r++) { i64 x = T[(size_t)r * k + i]; T[(size_t)r * k + i] = T[(size_t)r * k + i - 1]; T[(size_t)r * k + i - 1] = x; }
            GS();
            if (i > 1) i--;
        } else i++;
    }
#undef GS
    free(bs); free(nb); free(mu);
    *Tp = T;
    return 0;
ovf:
    free(bs); free(nb); free(mu); free(T);
    LAT_ERR = "lattice basis: integer overflow";
    return -1;
}
#undef OVF
/* c <- c T (k coordinates) */
static int lat_apply(int *c, const i64 *T, int k) {
    int cn[64];
    for (int j = 0; j < k; j++) {
        i128 sum = 0;
        for (int r = 0; r < k; r++) sum += (i128)c[r] * T[(size_t)r * k + j];
        if (sum > 0x3fffffff || sum < -0x3fffffff) { LAT_ERR = "lattice: coordinates too large"; return -1; }
        cn[j] = (int)sum;
    }
    for (int j = 0; j < k; j++) c[j] = cn[j];
    return 0;
}
#define LAT_CK(x) do { if (x) die(LAT_ERR); } while (0)

/* The lattice spanned by the cone generators. If its rank HC < h11, every curve used is an integer combination of its
 * basis (sums of generators), so curves are re-expressed in its HC coordinates: the cones, the grading and the charges
 * (their curve side), and the lightcone inequalities. Exact, not a projection: the multicover gcd and C/n are computed
 * in the same lattice. */
static void setup_lattice(void) {
    HC = h11; LB = NULL;
#ifndef CGV_MAJ   /* the majorant build keeps the input coordinates (and h11 <= 64) */
    Lat L; lat_init(&L, h11);
    for (int c = 0; c < ncones; c++) for (int i = 0; i < cone_n[c]; i++) LAT_CK(lat_add(&L, cone_g[c] + (size_t)i * h11));
    int rk = L.rk;
    if (rk == 0 || rk == h11) {
        lat_free(&L);
        if (h11 > 64) die("h11 > 64 needs curves in a sublattice of rank <= 64 (e.g. a face of the Mori cone)");
        return;
    }
    if (rk > 64) die("the curves span a lattice of rank > 64: not supported");
    int *pc; LAT_CK(lat_finish(&L, &LB, &pc));
    lat_free(&L);
    HC = rk;
    i64 *v = xmalloc(h11 * sizeof(i64));
    for (int c = 0; c < ncones; c++) {
        int *gc = xmalloc((size_t)(cone_n[c] ? cone_n[c] : 1) * HC * sizeof(int));
        for (int i = 0; i < cone_n[c]; i++) LAT_CK(lat_coords(LB, pc, HC, h11, cone_g[c] + (size_t)i * h11, gc + (size_t)i * HC, v));
        free(cone_g[c]); cone_g[c] = gc;
    }
    {
        i64 *T; LAT_CK(lat_lll(LB, HC, h11, &T));
        for (int c = 0; c < ncones; c++) for (int g = 0; g < cone_n[c]; g++) LAT_CK(lat_apply(cone_g[c] + (size_t)g * HC, T, HC));
        free(T);
    }
    /* curve side: x . W = c . (LB W) etc. */
    W = xmalloc(HC * sizeof(int)); Q = xmalloc((size_t)ndiv * HC * sizeof(int)); Q0 = xmalloc(HC * sizeof(int));
    for (int j = 0; j < HC; j++) {
        const i64 *b = LB + (size_t)j * h11;
        i128 w = 0, q0 = 0;
        for (int t = 0; t < h11; t++) { w += (i128)b[t] * Wd[t]; q0 += (i128)b[t] * Q0d[t]; }
        for (int r = 0; r < ndiv; r++) {
            i128 q = 0;
            for (int t = 0; t < h11; t++) q += (i128)b[t] * Qd[(size_t)r * h11 + t];
            if (q > 0x3fffffff || q < -0x3fffffff) die("curve lattice: charges too large");
            Q[(size_t)r * HC + j] = (int)q;
        }
        if (w > 0x3fffffff || w < -0x3fffffff || q0 > 0x3fffffff || q0 < -0x3fffffff) die("curve lattice: charges too large");
        W[j] = (int)w; Q0[j] = (int)q0;
    }
    if (nlight) {
        /* H p in the input coordinates now; H itself on the curve side */
        light_hp = xmalloc((size_t)nlight * (nlight_h ? nlight_h : 1) * sizeof(long long));
        light_hc = xmalloc((nlight_h ? nlight_h : 1) * sizeof(long long));
        for (int i = 0; i < nlight; i++) for (int f = 0; f < nlight_h; f++) {
            long long s = 0;
            for (int t = 0; t < h11; t++) s += (long long)light_pts[(size_t)i * h11 + t] * light_h[(size_t)f * h11 + t];
            light_hp[(size_t)i * nlight_h + f] = s;
        }
        int *lh = xmalloc((size_t)(nlight_h ? nlight_h : 1) * HC * sizeof(int));
        for (int f = 0; f < nlight_h; f++) for (int j = 0; j < HC; j++) {
            i128 s = 0;
            for (int t = 0; t < h11; t++) s += (i128)LB[(size_t)j * h11 + t] * light_h[(size_t)f * h11 + t];
            if (s > 0x3fffffff || s < -0x3fffffff) die("curve lattice: lightcone inequalities too large");
            lh[(size_t)f * HC + j] = (int)s;
        }
        free(light_h); light_h = lh;
    }
    free(v); free(pc);
#endif
}
#undef OVF
/* curve x (HC coordinates) in the input coordinates */
static void to_input_coords(const int *x, i64 *y) {
    if (!LB) { for (int t = 0; t < h11; t++) y[t] = x[t]; return; }
    for (int t = 0; t < h11; t++) y[t] = 0;
    for (int j = 0; j < HC; j++) if (x[j]) for (int t = 0; t < h11; t++) {
        i64 p;
        if (__builtin_mul_overflow((i64)x[j], LB[(size_t)j * h11 + t], &p) || __builtin_add_overflow(y[t], p, &y[t])) die("curve lattice: output coordinates overflow");
    }
}

/* Key packing needs |coordinate| < 2^(keybits-1) for every point of degree
 * <= maxdeg. Every curve used is a nonnegative combination x = sum lambda_i g_i of the
 * K_T Hilbert basis elements, so
 * |x_t| <= maxdeg * max_i |g_it| / deg(g_i). This also checks that the grading is positive on them. */
static int setup_keys(void) {
    if (getenv("CGV_FORCE_WIDE") && KEYBITS < atoi(getenv("CGV_FORCE_WIDE"))) return atoi(getenv("CGV_FORCE_WIDE")) > 256 ? 2 : 1;   /* test hook: run the 256/512-bit build */
    keybits = KEYBITS / HC;
    if (keybits > 40) keybits = 40;
#ifndef CGV_WIDE
    fieldmask = ((KEY)1 << keybits) - 1;
#endif
    double ratio = 0, rt[64] = {0};
    for (int c = 0; c < ncones; c++) {
        for (int i = 0; i < cone_n[c]; i++) {
            const int *g = cone_g[c] + (size_t)i * HC;
            int d = dot(g, W), nz = 0;
            for (int t = 0; t < HC; t++) nz |= g[t] != 0;
            if (!nz) continue;
            if (d <= 0) die("grading vector is not positive on the Mori cone");
            for (int t = 0; t < HC; t++) { double r = (double)(g[t] < 0 ? -g[t] : g[t]) / d; if (r > ratio) ratio = r; if (r > rt[t]) rt[t] = r; }
        }
    }
    double bound = ratio * maxdeg;
    if (bound >= (double)(1 << 30)) die("coordinates too large for key packing (|coordinate| >= 2^30)");   /* unpack gives ints */
    KBV = 0;
    if (bound >= (double)((i64)1 << (keybits - 1)) - 1) {
        /* per-coordinate widths: |x_t| <= maxdeg rt_t */
        int tot = 0;
        for (int t = 0; t < HC; t++) {
            int w = 2;
            while (rt[t] * maxdeg >= (double)((i64)1 << (w - 1)) - 1) w++;   /* <= 32: |x_t| < 2^30 */
            KBO[t] = tot; KBW[t] = w; tot += w;
        }
        if (tot > KEYBITS) {
            LOG("keys: coordinates need %d bits (HC = %d, max degree %d)\n", tot, HC, maxdeg);
            if (KEYBITS < 512 && tot <= 512) { LOG("keys: switching to %d-bit keys\n", tot <= 256 ? 256 : 512); return tot <= 256 ? 1 : 2; }   /* the caller reruns with that build */
            die("coordinates too large for key packing");
        }
        KBV = 1;
        LOG("keys: per-coordinate field widths, %d of %d bits\n", tot, KEYBITS);
    }
    int v[64] = {0}, t0 = -1;
    for (int t = 0; t < HC && t0 < 0; t++) if (W[t] > 0) t0 = t;
    for (int t = 0; t < HC && t0 < 0; t++) if (W[t] < 0) t0 = t;   /* (sublattice coordinates) */
    if (t0 < 0) die("grading vector is zero on the curves");
    v[t0] = W[t0] > 0 ? -1 : 1; LT_XK = pack(v);
    LT_HMI = LT_HM; for (int i = 0; i < 8; i++) LT_HMI *= 2 - LT_HM * LT_HMI; /* Newton: doubles the correct bits */
    return 0;
}

/* ------------------------------------------------------------------------- */
/* Slot maps: key -> slot, with per-slot key, degree and value          */
/* ------------------------------------------------------------------------- */
/* hash slots are one u64 each: low 32 bits = slot index + 1 (0 = empty), high 32 bits = a tag from the hash,
 * so a probe compares the tag and reads the stored key only on a tag match (8 B per slot instead of 20) */
typedef struct {
    u64 *hs; u64 mask;
    KEY *key; int *deg; Co *val; int n, cap;
} SM;
static void sm_init(SM *m, int cap) {
    u64 c = 16;
    while (c < 2 * (u64)cap) c <<= 1;
    m->hs = xcalloc(c, sizeof(u64)); m->mask = c - 1;
    m->cap = cap > 16 ? cap : 16; m->n = 0;
    m->key = xmalloc(m->cap * sizeof(KEY)); m->deg = xmalloc(m->cap * sizeof(int)); m->val = xmalloc(m->cap * sizeof(Co));
}
static void sm_free(SM *m) { free(m->hs); free(m->key); free(m->deg); free(m->val); }
static inline int sm_get(const SM *m, KEY k) {
    u64 x = hash128(k), h = x & m->mask, tg = x >> 32;
    for (;;) {
        u64 e = m->hs[h];
        if (!e) return -1;
        if (e >> 32 == tg) { int s = (int)(uint32_t)e - 1; if (k_eq(m->key[s], k)) return s; }
        h = (h + 1) & m->mask;
    }
}
static void sm_rehash(SM *m) {
    u64 c = (m->mask + 1) * 2;
    free(m->hs);
    m->hs = xcalloc(c, sizeof(u64)); m->mask = c - 1;
    for (int s = 0; s < m->n; s++) {
        u64 x = hash128(m->key[s]), h = x & m->mask;
        while (m->hs[h]) h = (h + 1) & m->mask;
        m->hs[h] = (x >> 32 << 32) | (u64)(s + 1);
    }
}
/* get or create; *created set to 1 if new (value zeroed) */
static inline int sm_add_x(SM *m, KEY k, u64 x, int deg, int *created) {
    u64 h = x & m->mask, tg = x >> 32;
    for (;;) {
        u64 e = m->hs[h];
        if (!e) break;
        if (e >> 32 == tg) { int s = (int)(uint32_t)e - 1; if (k_eq(m->key[s], k)) { *created = 0; return s; } }
        h = (h + 1) & m->mask;
    }
    if (m->n == m->cap) {
        m->cap *= 2;
        m->key = xrealloc(m->key, m->cap * sizeof(KEY)); m->deg = xrealloc(m->deg, m->cap * sizeof(int)); m->val = xrealloc(m->val, m->cap * sizeof(Co));
    }
    int s = m->n++;
    m->key[s] = k; m->deg[s] = deg; m->val[s] = co_zero();
    m->hs[h] = (tg << 32) | (u64)(s + 1);
    *created = 1;
    if (2 * (u64)m->n > m->mask) sm_rehash(m);
    return s;
}
static inline int sm_add(SM *m, KEY k, int deg, int *created) { return sm_add_x(m, k, hash128(k), deg, created); }
/* remove everything; cost proportional to the number of slots */
static void sm_clear(SM *m) {
    if ((u64)m->n * 8 > m->mask) { memset(m->hs, 0, (m->mask + 1) * sizeof(u64)); m->n = 0; return; }
    /* locate all positions first, then blank, so probe chains stay intact */
    static __thread u64 *pos = NULL; static __thread int pcap = 0;
    if (pcap < m->n) { pcap = m->n * 2; pos = xrealloc(pos, pcap * sizeof(u64)); }
    for (int s = 0; s < m->n; s++) {
        u64 h = hash128(m->key[s]) & m->mask;
        while ((int)(uint32_t)m->hs[h] != s + 1) h = (h + 1) & m->mask;
        pos[s] = h;
    }
    for (int s = 0; s < m->n; s++) m->hs[pos[s]] = 0;
    m->n = 0;
}

/* Sparse polynomial, sorted by (degree, key). */
typedef struct { int n; KEY *key; int *deg; Co *val; } SPoly;
static void sp_free(SPoly *p) { free(p->key); free(p->deg); free(p->val); memset(p, 0, sizeof(*p)); }
static void sp_alloc(SPoly *p, int n) { p->n = 0; p->key = xmalloc((n ? n : 1) * sizeof(KEY)); p->deg = xmalloc((n ? n : 1) * sizeof(int)); p->val = xmalloc((n ? n : 1) * sizeof(Co)); }

/* sort slot indices by (degree, key) */
static void sort_slots(const SM *m, int *ord, int n) {
    /* insertion-merge sort on (deg, key); qsort needs a global for the map */
    if (n < 2) return;
    int *tmp = xmalloc(n * sizeof(int));
    for (int w = 1; w < n; w *= 2) {
        for (int lo = 0; lo < n; lo += 2 * w) {
            int mid = lo + w < n ? lo + w : n, hi = lo + 2 * w < n ? lo + 2 * w : n;
            int i = lo, j = mid, k = lo;
            while (i < mid && j < hi) {
                int a = ord[i], b = ord[j];
                int le = m->deg[a] < m->deg[b] || (m->deg[a] == m->deg[b] && k_le(m->key[a], m->key[b]));
                tmp[k++] = le ? ord[i++] : ord[j++];
            }
            while (i < mid) tmp[k++] = ord[i++];
            while (j < hi) tmp[k++] = ord[j++];
        }
        memcpy(ord, tmp, n * sizeof(int));
    }
    free(tmp);
}
/* Nonzero slots as a sorted sparse polynomial; clears the map. */
static SPoly sm_take(SM *m) {
    int *ord = xmalloc((m->n ? m->n : 1) * sizeof(int));
    int k = 0;
    for (int s = 0; s < m->n; s++) if (!co_isz(m->val[s])) ord[k++] = s;
    sort_slots(m, ord, k);
    SPoly p; sp_alloc(&p, k);
    for (int z = 0; z < k; z++) { int s = ord[z]; p.key[z] = m->key[s]; p.deg[z] = m->deg[s]; p.val[z] = m->val[s]; }
    p.n = k;
    free(ord);
    sm_clear(m);
    return p;
}

/* Binary min-heap of slots keyed by degree. */
typedef struct { int *s; int n, cap; } Heap;
static inline void hp_push(Heap *h, const int *deg, int s) {
    if (h->n == h->cap) { h->cap = h->cap ? 2 * h->cap : 256; h->s = xrealloc(h->s, h->cap * sizeof(int)); }
    int i = h->n++;
    while (i > 0) { int p = (i - 1) / 2; if (deg[h->s[p]] <= deg[s]) break; h->s[i] = h->s[p]; i = p; }
    h->s[i] = s;
}
static inline int hp_pop(Heap *h, const int *deg) {
    int top = h->s[0], last = h->s[--h->n];
    int i = 0;
    for (;;) {
        int c = 2 * i + 1;
        if (c >= h->n) break;
        if (c + 1 < h->n && deg[h->s[c + 1]] < deg[h->s[c]]) c++;
        if (deg[h->s[c]] >= deg[last]) break;
        h->s[i] = h->s[c]; i = c;
    }
    if (h->n) h->s[i] = last;
    return top;
}

/* ------------------------------------------------------------------------- */
/* Threads                                                                   */
/* ------------------------------------------------------------------------- */
static int nthreads = 1;
static int use_gpu = 0, gpu_device = 0;
/* Persistent thread pool: parallel_for posts a job, the nthreads-1 workers and
 * the caller take indices from a shared counter. Creating threads per call was
 * the dominant cost for gradings with many tiny layers. */
typedef struct {
    pthread_mutex_t lock; pthread_cond_t wake, done;
    void (*fn)(void *, int, int); void *ctx;
    int n, next, active, gen, quit, nw, started;
    pthread_t *th;
} Pool;
static Pool POOL = {PTHREAD_MUTEX_INITIALIZER, PTHREAD_COND_INITIALIZER, PTHREAD_COND_INITIALIZER, NULL, NULL, 0, 0, 0, 0, 0, 0, 0, NULL};
typedef struct { int tid; } PoolArg;
static void pool_run(int tid) {
    for (;;) {
        pthread_mutex_lock(&POOL.lock);
        int i = POOL.next < POOL.n ? POOL.next++ : -1;
        pthread_mutex_unlock(&POOL.lock);
        if (i < 0) break;
        POOL.fn(POOL.ctx, i, tid);
    }
}
static void *pool_worker(void *arg) {
    int tid = ((PoolArg *)arg)->tid;
#ifdef CGV_MAJ
    fesetround(FE_UPWARD);
#endif
    int seen = 0;
    for (;;) {
        pthread_mutex_lock(&POOL.lock);
        while (POOL.gen == seen && !POOL.quit) pthread_cond_wait(&POOL.wake, &POOL.lock);
        if (POOL.quit) { pthread_mutex_unlock(&POOL.lock); break; }
        seen = POOL.gen;
        POOL.active++;
        pthread_mutex_unlock(&POOL.lock);
        pool_run(tid);
        pthread_mutex_lock(&POOL.lock);
        if (--POOL.active == 0) pthread_cond_signal(&POOL.done);
        pthread_mutex_unlock(&POOL.lock);
    }
    return NULL;
}
static void pool_start(int nth) {
    if (POOL.started) return;
    POOL.nw = nth - 1; POOL.started = 1;
    POOL.th = xmalloc((POOL.nw ? POOL.nw : 1) * sizeof(pthread_t));
    PoolArg *args = xmalloc((POOL.nw ? POOL.nw : 1) * sizeof(PoolArg));
    for (int t = 0; t < POOL.nw; t++) { args[t].tid = t + 1; pthread_create(&POOL.th[t], NULL, pool_worker, &args[t]); }
}
static void parallel_for(int n, int nth, void (*fn)(void *, int, int), void *ctx) {
    if (nth <= 1 || n <= 1) { for (int i = 0; i < n; i++) fn(ctx, i, 0); return; }
    pool_start(nthreads);
    pthread_mutex_lock(&POOL.lock);
    POOL.fn = fn; POOL.ctx = ctx; POOL.n = n; POOL.next = 0; POOL.gen++;
    pthread_cond_broadcast(&POOL.wake);
    pthread_mutex_unlock(&POOL.lock);
    pool_run(0);
    /* wait until every worker that joined this job has finished */
    pthread_mutex_lock(&POOL.lock);
    while (POOL.active > 0) pthread_cond_wait(&POOL.done, &POOL.lock);
    POOL.fn = NULL;
    pthread_mutex_unlock(&POOL.lock);
}

/* ------------------------------------------------------------------------- */
/* Candidate curves: lattice points of the cones K_T (<= 2 negative GLSM      */
/* intersections and degree <= maxdeg                                         */
/* ------------------------------------------------------------------------- */
static SM CAND;

/* backward lightcones: C is kept if H (p - C) >= 0 for some point p, i.e. H C <= H p row by row */
static int in_lightcone(const int *C) {
    if (!light_hp) {
        light_hc = xmalloc((nlight_h ? nlight_h : 1) * sizeof(long long));
        light_hp = xmalloc((size_t)nlight * (nlight_h ? nlight_h : 1) * sizeof(long long));
        for (int i = 0; i < nlight; i++)
            for (int f = 0; f < nlight_h; f++) light_hp[(size_t)i * nlight_h + f] = dotll(light_pts + (size_t)i * h11, light_h + (size_t)f * h11);
    }
    long long *hc = light_hc;   /* single-threaded callers only */
    for (int f = 0; f < nlight_h; f++) hc[f] = dotll(C, light_h + (size_t)f * HC);
    for (int i = 0; i < nlight; i++) {
        const long long *hp = light_hp + (size_t)i * nlight_h;
        int f = 0;
        while (f < nlight_h && hc[f] <= hp[f]) f++;
        if (f == nlight_h) return 1;
    }
    return 0;
}

/* Enumerate the nonnegative integer sums of gens up to maxdeg (BFS in the
 * scratch map vis), inserting them into out. */
static void enumerate_sums(const int *gens, int ng, SM *vis, SM *out) {
    KEY *gk = xmalloc((ng ? ng : 1) * sizeof(KEY));
    int *gd = xmalloc((ng ? ng : 1) * sizeof(int));
    int m = 0;
    for (int i = 0; i < ng; i++) {
        const int *g = gens + (size_t)i * HC;
        int d = dot(g, W);
        if (d <= 0 || d > maxdeg) continue;
        gk[m] = pack(g); gd[m] = d; m++;
    }
    int cr;
    sm_add(vis, K_ZERO, 0, &cr);
    for (int s = 0; s < vis->n; s++) {
        KEY k = vis->key[s]; int d = vis->deg[s];
        for (int j = 0; j < m; j++) {
            int dd = d + gd[j];
            if (dd > maxdeg) continue;
            sm_add(vis, k_add(k, gk[j]), dd, &cr);
        }
    }
    for (int s = 0; s < vis->n; s++) sm_add(out, vis->key[s], vis->deg[s], &cr);
    sm_clear(vis);
    free(gk); free(gd);
}

typedef struct { SM *vis, *uni; const int *use; } ConeCtx;
static void cone_fn(void *vctx, int i, int tid) {
    ConeCtx *c = vctx;
    if (!c->use[i]) return;
    enumerate_sums(cone_g[i], cone_n[i], &c->vis[tid], &c->uni[tid]);
}

typedef struct { SM *uni; int nth; SM *sh; int nsh; } CandMerge;
static inline int shard_n(KEY k, int nsh) { return (int)((hash128(k) >> 40) % (u64)nsh); }
static void cand_merge_fn(void *vctx, int k, int tid) {
    CandMerge *c = vctx; int cr;
    for (int t = 0; t < c->nth; t++) {
        const SM *u = &c->uni[t];
        for (int s = 0; s < u->n; s++) if (shard_n(u->key[s], c->nsh) == k) sm_add(&c->sh[k], u->key[s], u->deg[s], &cr);
    }
}

static void build_candidates(void) {
    double t0 = now();
    sm_init(&CAND, 1 << 16);
    {
        /* Skip cones contained in another. The negatives of any point of K_T lie
         * in the union of the negatives of its Hilbert basis, negmask[T], which
         * is itself realized inside K_T. So K_c is inside K_c2 whenever
         * negmask[c] is a subset of negmask[c2]; equal masks keep the first. */
        int *use = xmalloc(ncones * sizeof(int));
        int nw = (ndiv + 63) / 64;   /* negmask[c]: nw words of divisor bits */
        u64 *negmask = xcalloc((size_t)ncones * nw, sizeof(u64));
        for (int c = 0; c < ncones; c++) {
            u64 *m = negmask + (size_t)c * nw;
            for (int i = 0; i < cone_n[c]; i++) {
                const int *g = cone_g[c] + (size_t)i * HC;
                for (int r = 0; r < ndiv; r++) if (dot(g, Q + (size_t)r * HC) < 0) m[r >> 6] |= (u64)1 << (r & 63);
            }
        }
        /* A curve with exactly two negative intersections s1, s2 only enters through
         * S2 = (factorial) * form(s1, s2), form(s1, s2) = sum_{a<=b} K_ab (Q_s1a Q_s2b + Q_s1b Q_s2a)
         * (half on the diagonal): if that integer vanishes the whole cone K_{s1,s2}
         * adds nothing beyond its one-negative points, which the |T| = 1 cones cover. */
        int nzero = 0;
        int *skip = xcalloc(ncones, sizeof(int));
        for (int c = 0; c < ncones; c++) {
            if (cone_T[c][0] < 0 || cone_T[c][1] < 0 || cone_T[c][2] >= 0) continue;   /* |T| = 2 cones only */
            int s1 = cone_T[c][0], s2 = cone_T[c][1];
            i64 form2 = 0; /* 2 * form, over the integers */
            /* recompute K_ab exactly as tables_init does (distinct permutations) */
            i64 *Kab2 = xcalloc((size_t)h11 * h11, sizeof(i64));
            for (int z = 0; z < nint; z++) {
                int idx[3] = {intnums[z][0], intnums[z][1], intnums[z][2]};
                static const int perms[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
                int seen[6][3]; int ns = 0;
                for (int p = 0; p < 6; p++) {
                    int t = idx[perms[p][0]], a = idx[perms[p][1]], b = idx[perms[p][2]];
                    int dup = 0;
                    for (int q = 0; q < ns; q++) if (seen[q][0] == t && seen[q][1] == a && seen[q][2] == b) dup = 1;
                    if (dup) continue;
                    seen[ns][0] = t; seen[ns][1] = a; seen[ns][2] = b; ns++;
                    Kab2[a * h11 + b] += (i64)Wd[t] * intnums[z][3];
                }
            }
            for (int a = 0; a < h11; a++) for (int b = a; b < h11; b++) {
                i64 x = (i64)Qd[s1 * h11 + a] * Qd[s2 * h11 + b] + (i64)Qd[s1 * h11 + b] * Qd[s2 * h11 + a];
                form2 += (a == b ? 1 : 2) * Kab2[a * h11 + b] * x;
            }
            free(Kab2);
            if (form2 == 0) { skip[c] = 1; nzero++; }
        }
        int nuse = 0;
        for (int c = 0; c < ncones; c++) {
            use[c] = 1;
            if (skip[c]) { use[c] = 0; continue; }
            const u64 *mc = negmask + (size_t)c * nw;
            for (int c2 = 0; c2 < ncones && use[c]; c2++) {
                if (c2 == c || skip[c2]) continue;
                const u64 *m2 = negmask + (size_t)c2 * nw;
                int sub = 1, eq = 1;
                for (int k = 0; k < nw; k++) { sub &= (mc[k] & ~m2[k]) == 0; eq &= mc[k] == m2[k]; }
                if (sub && (!eq || c2 < c)) use[c] = 0;
            }
            nuse += use[c];
        }
        int nth = nthreads;
        SM *vis = xmalloc(nth * sizeof(SM)), *uni = xmalloc(nth * sizeof(SM));
        for (int t = 0; t < nth; t++) { sm_init(&vis[t], 1024); sm_init(&uni[t], 1024); }
        ConeCtx cc = {vis, uni, use};
        parallel_for(ncones, nth, cone_fn, &cc);
        /* merge the per-thread unions in parallel, one hash shard per task */
        double tm = now();
        int nsh = nth * 4;
        SM *sh = xmalloc(nsh * sizeof(SM));
        for (int k = 0; k < nsh; k++) sm_init(&sh[k], 1024);
        CandMerge cm = {uni, nth, sh, nsh};
        parallel_for(nsh, nth, cand_merge_fn, &cm);
        int tot = 0;
        for (int k = 0; k < nsh; k++) tot += sh[k].n;
        CAND.n = 0;
        CAND.key = xrealloc(CAND.key, (tot ? tot : 1) * sizeof(KEY)); CAND.deg = xrealloc(CAND.deg, (tot ? tot : 1) * sizeof(int));
        CAND.val = xrealloc(CAND.val, sizeof(Co)); CAND.cap = tot;
        for (int k = 0; k < nsh; k++) {
            memcpy(CAND.key + CAND.n, sh[k].key, sh[k].n * sizeof(KEY));
            memcpy(CAND.deg + CAND.n, sh[k].deg, sh[k].n * sizeof(int));
            CAND.n += sh[k].n;
            sm_free(&sh[k]);
        }
        free(sh);
        for (int t = 0; t < nth; t++) { sm_free(&vis[t]); sm_free(&uni[t]); }
        if (getenv("CGV_PROF")) LOG("    candidates: cones %.2fs, merge %.2fs\n", tm - t0, now() - tm);
        free(vis); free(uni); free(use); free(negmask); free(skip);
        LOG("candidates: %d points with <=2 negatives from %d of %d cones (%d with vanishing form), %.2fs\n", CAND.n, nuse, ncones, nzero, now() - t0);
    }
    if (nlight) {
        /* every curve feeding into a curve of the lightcones is in them, so the other candidates are dropped */
        double tr = now();
        int m = 0, before = CAND.n, C[64];
        for (int s = 0; s < CAND.n; s++) {
            unpack(CAND.key[s], C);
            if (in_lightcone(C)) { CAND.key[m] = CAND.key[s]; CAND.deg[m] = CAND.deg[s]; m++; }
        }
        CAND.n = m;
        LOG("backward lightcones of %d point%s: %d of %d candidates kept, %.2fs\n", nlight, nlight == 1 ? "" : "s", m, before, now() - tr);
    }
}

/* ------------------------------------------------------------------------- */
/* Tables and fundamental period                                             */
/* ------------------------------------------------------------------------- */
static int tabn;
static Co *INV, *H1, *H2, *FACT, *IFACT, *DEGC;
static Co *Kab, *Qm, *Q0m;
/* Row mode (RMODE): c1 and alpha in the coordinates of the charge rows rho = 0 (Q0) and the divisors whose charge is
 * nonzero on the curves (NR of them) instead of the h11 divisor coordinates; used when NR < h11 (curves on a face
 * meet few divisors). c1 = sum_rho c1_rho q_rho with q_rho the full charge row, so alpha = sum_rho beta_rho q_rho,
 * C.alpha = sum_rho beta_rho (C.q_rho), and every K-contraction becomes one of M = q_rho^T K q_sigma (integers). */
static int RMODE, NR;
#ifndef CGV_MAJ
static int *RROW, *ROWOF;   /* RROW[rho]: divisor (-1: Q0); ROWOF[r]: row of divisor r or -1 */
#endif
static i64 *MINT; static Co *Mm;
static void setup_rows(void) {
    RMODE = 0;
#ifndef CGV_MAJ
    RROW = xmalloc((ndiv + 1) * sizeof(int)); ROWOF = xmalloc(ndiv * sizeof(int));
    NR = 0; RROW[NR++] = -1;
    for (int r = 0; r < ndiv; r++) {
        int nz = 0;
        for (int j = 0; j < HC; j++) nz |= Q[(size_t)r * HC + j] != 0;
        ROWOF[r] = nz ? NR : -1;
        if (nz) RROW[NR++] = r;
    }
    if (NR >= h11) return;
    RMODE = 1;
    /* K_ab = sum_t W_t kappa_tab as a sparse list (distinct permutations, as tables_init) */
    int nk = 0; int (*ke)[2] = xmalloc((size_t)(nint ? nint : 1) * 6 * sizeof(*ke)); i64 *kv = xmalloc((size_t)(nint ? nint : 1) * 6 * sizeof(i64));
    for (int z = 0; z < nint; z++) {
        int idx[3] = {intnums[z][0], intnums[z][1], intnums[z][2]};
        static const int perms[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
        int seen[6][3]; int ns = 0;
        for (int p = 0; p < 6; p++) {
            int t = idx[perms[p][0]], a = idx[perms[p][1]], b = idx[perms[p][2]];
            int dup = 0;
            for (int q = 0; q < ns; q++) if (seen[q][0] == t && seen[q][1] == a && seen[q][2] == b) dup = 1;
            if (dup) continue;
            seen[ns][0] = t; seen[ns][1] = a; seen[ns][2] = b; ns++;
            if (Wd[t] && intnums[z][3]) { ke[nk][0] = a; ke[nk][1] = b; kv[nk] = (i64)Wd[t] * intnums[z][3]; nk++; }
        }
    }
    MINT = xmalloc((size_t)NR * NR * sizeof(i64));
    i128 *kq = xmalloc(h11 * sizeof(i128));
    for (int sg = 0; sg < NR; sg++) {
        const int *qs = RROW[sg] < 0 ? Q0d : Qd + (size_t)RROW[sg] * h11;
        for (int a = 0; a < h11; a++) kq[a] = 0;
        for (int e = 0; e < nk; e++) kq[ke[e][0]] += (i128)kv[e] * qs[ke[e][1]];
        for (int rh = 0; rh < NR; rh++) {
            const int *qr = RROW[rh] < 0 ? Q0d : Qd + (size_t)RROW[rh] * h11;
            i128 m = 0;
            for (int a = 0; a < h11; a++) if (qr[a]) m += kq[a] * qr[a];
            if (m > ((i128)1 << 62) || m < -((i128)1 << 62)) die("row mode: q K q too large");
            MINT[(size_t)rh * NR + sg] = (i64)m;
        }
    }
    free(kq); free(ke); free(kv);
#endif
}

static void tables_init(void) {
    int maxv = 1;
    int C[64];
    for (int s = 0; s < CAND.n; s++) {
        unpack(CAND.key[s], C);
        int s0 = dot(C, Q0);
        if (s0 > maxv) maxv = s0;
        if (-s0 > maxv) maxv = -s0;   /* bgv: |m| for negative anticanonical degree */
        for (int r = 0; r < ndiv; r++) { int v = dot(C, Q + (size_t)r * HC); if (v < 0) v = -v; if (v > maxv) maxv = v; }
    }
    tabn = (maxv > maxdeg ? maxv : maxdeg) + 2;
    INV = xmalloc(tabn * sizeof(Co)); H1 = xmalloc(tabn * sizeof(Co)); H2 = xmalloc(tabn * sizeof(Co));
    FACT = xmalloc(tabn * sizeof(Co)); IFACT = xmalloc(tabn * sizeof(Co)); DEGC = xmalloc(tabn * sizeof(Co));
    FACT[0] = co_int(1); INV[0] = H1[0] = H2[0] = co_zero(); DEGC[0] = co_zero();
    for (int i = 1; i < tabn; i++) { DEGC[i] = co_int(i); FACT[i] = co_mul(FACT[i - 1], DEGC[i]); }
#ifdef CGV_MAJ
    IFACT[0] = co_int(1); /* 1/FACT would round the wrong way: build it from 1/i (rounded up) */
    for (int i = 1; i < tabn; i++) IFACT[i] = co_mul(IFACT[i - 1], co_inv(DEGC[i]));
    if (!isfinite(FACT[tabn - 1].v[0])) die("majorant: factorial table overflows");
#else
    IFACT[tabn - 1] = co_inv(FACT[tabn - 1]);
    for (int i = tabn - 1; i > 0; i--) IFACT[i - 1] = co_mul(IFACT[i], DEGC[i]);
#endif
    for (int i = 1; i < tabn; i++) {
#ifdef CGV_MAJ
        INV[i] = co_inv(DEGC[i]);
#else
        INV[i] = co_mul(IFACT[i], FACT[i - 1]);
#endif
        H1[i] = co_add(H1[i - 1], INV[i]);
        H2[i] = co_add(H2[i - 1], co_mul(INV[i], INV[i]));
    }
    int H = h11;
    Kab = xmalloc((size_t)H * H * sizeof(Co));
    for (int i = 0; i < H * H; i++) Kab[i] = co_zero();
    for (int z = 0; z < nint; z++) {
        int idx[3] = {intnums[z][0], intnums[z][1], intnums[z][2]};
        int v = intnums[z][3];
        if (!v) continue;
        static const int perms[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
        int seen[6][3]; int ns = 0;
        for (int p = 0; p < 6; p++) {
            int t = idx[perms[p][0]], a = idx[perms[p][1]], b = idx[perms[p][2]];
            int dup = 0;
            for (int q = 0; q < ns; q++) if (seen[q][0] == t && seen[q][1] == a && seen[q][2] == b) dup = 1;
            if (dup) continue;
            seen[ns][0] = t; seen[ns][1] = a; seen[ns][2] = b; ns++;
            Kab[a * H + b] = co_add(Kab[a * H + b], co_int((i64)Wd[t] * v));
        }
    }
    Qm = xmalloc((size_t)ndiv * H * sizeof(Co)); Q0m = xmalloc(H * sizeof(Co));
    if (nvex) {
        vex_val = xmalloc((size_t)nvex * H * sizeof(Co));
        for (int v = 0; v < nvex; v++) for (int a = 0; a < H; a++) {
            Co d = co_int(vex_den[(size_t)v * h11 + a]);
            if (co_isz(d)) die("vex stratum denominator divisible by a prime");
            vex_val[(size_t)v * H + a] = co_mul(co_int(vex_num[(size_t)v * h11 + a]), co_inv(d));
        }
    }
    for (int i = 0; i < ndiv * H; i++) Qm[i] = co_int(Qd[i]);
    for (int a = 0; a < H; a++) Q0m[a] = co_int(Q0d[a]);
    if (RMODE) { Mm = xmalloc((size_t)NR * NR * sizeof(Co)); for (int i = 0; i < NR * NR; i++) Mm[i] = co_int(MINT[i]); }
}
static void tables_free(void) {
    free(INV); free(H1); free(H2); free(FACT); free(IFACT); free(DEGC); free(Kab); free(Qm); free(Q0m);
    if (RMODE) free(Mm);
    if (nvex) free(vex_val);
}

/* coordinates of c1 and alpha: h11 of them, or (majorant build) one per row of the charge
 * matrix [Q0; Q_1..Q_ndiv]. Every c1 is a combination of those integer rows, so there
 * alpha_u = sum_r beta_ur Q_r and |C.alpha_u| <= sum_r |beta_ur| |C.Q_r| with C.Q_r exact:
 * tighter than sum_t |C_t| |alpha_ut|, and zero exactly where C meets none of the rows. */
static int NA;
static inline const int *qrow(int r) { return r == 0 ? Q0 : Q + (size_t)(r - 1) * HC; }
static inline int curve_coords(const int *C, Co *cm, int *nzt) {
    int n = 0;
#ifdef CGV_MAJ
    for (int r = 0; r < NA; r++) { int x = dot(C, qrow(r)); if (x) { cm[n] = co_int(x); nzt[n++] = r; } }
#else
    for (int t = 0; t < HC; t++) if (C[t]) { cm[n] = co_int(C[t]); nzt[n++] = t; }
#endif
    return n;
}

/* Fundamental period data at curve C: c0, c1[a], S2 = sum_{a<=b} K_ab c2_ab
 * (1/2 on the diagonal). Formulas follow cygv's compute_c_{0,1,2}neg. */
static void fp_point(const int *C, Co *c0, Co *c1, Co *s2) {
    int H = h11;
    int dq[ndiv];
    int dq0 = dot(C, Q0);
    int nneg = 0, s1 = -1, sb = -1;
    for (int r = 0; r < ndiv; r++) {
        dq[r] = dot(C, Q + (size_t)r * HC);
        if (dq[r] < 0) { nneg++; if (s1 < 0) s1 = r; else sb = r; }
    }
    *c0 = co_zero(); *s2 = co_zero();
    for (int a = 0; a < H; a++) c1[a] = co_zero();
    if (dq0 < 0) {
        /* bgv: m = q0.n < 0 (vex phase). S = negative set is a vex cone; the 1/L pole cancels against the
         * stratum class [V_S]_X (Liam's note / Wang's quasimap I-function):
         *   c = [V_S]_X * U0 * exp(A.rho + ...),  U0 = (-1)^(sum_S (k_i-1) + |m|-1) prod_S (k_i-1)! / (prod_notS b_i! (|m|-1)!)
         *   A_a = q0_a H_{|m|-1} - sum_notS q_ia H_{b_i} - sum_S q_ia H_{k_i-1}
         * |S| = 3: S2 = U0 (w . d); |S| = 2 cannot occur for fine fans (MacFadden-Sheridan Prop. 5); else 0. */
        if (nneg < 2) die("negative anticanonical degree with fewer than two negative intersections");
        /* |S| = 2: S cannot be a cone of a fine fan (MacFadden-Sheridan Prop. 5; the driver rejects vex 2-cones),
         * so S is not a cone and the term vanishes (Stanley-Reisner). |S| >= 4: invisible to GVs. */
        if (nneg != 3) return;
#ifdef CGV_MAJ
        die("majorant build: vex strata not supported");
#endif
        int neg[3], nn = 0;
        for (int r = 0; r < ndiv; r++) if (dq[r] < 0) neg[nn++] = r;
        int vs = -1;
        for (int v = 0; v < nvex && vs < 0; v++) {
            if (vex_k[v] != nneg) continue;
            int ok = 1;
            for (int j = 0; j < nneg; j++) if (vex_S[v][j] != neg[j]) ok = 0;
            if (ok) vs = v;
        }
        if (vs < 0) return;   /* S is not a cone of the fan: the term vanishes (Stanley-Reisner) */
        int mm = -dq0;
        Co u = IFACT[mm - 1];
        int e = mm - 1;
        for (int r = 0; r < ndiv; r++) {
            if (dq[r] < 0) { u = co_mul(u, FACT[-dq[r] - 1]); e += -dq[r] - 1; }
            else u = co_mul(u, IFACT[dq[r]]);
        }
        if (e % 2) u = co_neg(u);
        const Co *cv = vex_val + (size_t)vs * H;
        {
            Co wd = co_zero();
            for (int t = 0; t < H; t++) if (Wd[t]) wd = co_add(wd, co_mul(co_int(Wd[t]), cv[t]));
            *s2 = co_mul(u, wd);
        }
        return;
    }
    if (nneg > 2) return;
    Co half = INV[2];
    Co fact = FACT[dq0];
    for (int r = 0; r < ndiv; r++) fact = co_mul(fact, dq[r] >= 0 ? IFACT[dq[r]] : FACT[-dq[r] - 1]);
    Co Av[H];
    Co s = co_zero();
    if (nneg == 0) {
        *c0 = fact;
        for (int a = 0; a < H; a++) {
            Co v = co_mul(H1[dq0], Q0m[a]);
            for (int r = 0; r < ndiv; r++) if (Qd[r * H + a]) v = co_sub(v, co_mul(H1[dq[r]], Qm[r * H + a]));
            Av[a] = v;
            c1[a] = co_mul(fact, v);
        }
        for (int a = 0; a < H; a++) for (int b = a; b < H; b++) {
            Co K = Kab[a * H + b];
            if (co_isz(K)) continue;
            Co v = co_mul(Av[a], Av[b]);
            v = co_sub(v, co_mul(H2[dq0], co_mul(Q0m[a], Q0m[b])));
            for (int r = 0; r < ndiv; r++)
                if (Qd[r * H + a] && Qd[r * H + b]) v = co_add(v, co_mul(H2[dq[r]], co_mul(Qm[r * H + a], Qm[r * H + b])));
            if (a == b) K = co_mul(K, half);
            s = co_add(s, co_mul(K, v));
        }
        *s2 = co_mul(fact, s);
    } else if (nneg == 1) {
        int v = dq[s1];
        Co fs = (v % 2 == 0) ? co_neg(fact) : fact; /* sn = -1 if v even */
        for (int a = 0; a < H; a++) {
            Co x = co_mul(H1[dq0], Q0m[a]);
            for (int r = 0; r < ndiv; r++) {
                if (!Qd[r * H + a]) continue;
                int n = dq[r] < 0 ? -dq[r] - 1 : dq[r];
                x = co_sub(x, co_mul(H1[n], Qm[r * H + a]));
            }
            Av[a] = x;
            c1[a] = co_mul(fs, Qm[s1 * H + a]);
        }
        for (int a = 0; a < H; a++) for (int b = a; b < H; b++) {
            Co K = Kab[a * H + b];
            if (co_isz(K)) continue;
            Co x = co_add(co_mul(Av[a], Qm[s1 * H + b]), co_mul(Av[b], Qm[s1 * H + a]));
            if (a == b) K = co_mul(K, half);
            s = co_add(s, co_mul(K, x));
        }
        *s2 = co_mul(fs, s);
    } else {
        int v1 = dq[s1], v2 = dq[sb];
        Co fs = ((v1 + v2) % 2 == 0) ? fact : co_neg(fact);
        for (int a = 0; a < H; a++) for (int b = a; b < H; b++) {
            Co K = Kab[a * H + b];
            if (co_isz(K)) continue;
            Co x = co_add(co_mul(Qm[s1 * H + a], Qm[sb * H + b]), co_mul(Qm[s1 * H + b], Qm[sb * H + a]));
            if (a == b) K = co_mul(K, half);
            s = co_add(s, co_mul(K, x));
        }
        *s2 = co_mul(fs, s);
    }
#ifdef CGV_MAJ
    for (int r = 0; r < NA; r++) c1[r] = co_zero();
    if (nneg == 0) { c1[0] = co_mul(fact, H1[dq0]); for (int r = 0; r < ndiv; r++) c1[1 + r] = co_mul(fact, H1[dq[r]]); }
    else if (nneg == 1) c1[1 + s1] = fact;
#endif
}

#ifndef CGV_MAJ
/* fp_point in row mode: c1row[rho] (NR entries); same S2 through M = q^T K q (sum_{a<=b} K'_ab x_ab = 1/2 x^T K x
 * for the symmetric products, and = u^T K v for x_ab = u_a v_b + u_b v_a) */
static void fp_point_r(const int *C, Co *c0, Co *c1, Co *s2) {
    int dq[ndiv];
    int dq0 = dot(C, Q0);
    int nneg = 0, s1 = -1, sb = -1;
    for (int r = 0; r < ndiv; r++) {
        dq[r] = ROWOF[r] < 0 ? 0 : dot(C, Q + (size_t)r * HC);
        if (dq[r] < 0) { nneg++; if (s1 < 0) s1 = r; else sb = r; }
    }
    *c0 = co_zero(); *s2 = co_zero();
    for (int a = 0; a < NR; a++) c1[a] = co_zero();
    if (dq0 < 0) {   /* vex phase: c1 = 0; S2 as fp_point */
        Co c1d[h11];
        fp_point(C, c0, c1d, s2);
        return;
    }
    if (nneg > 2) return;
    Co half = INV[2];
    Co fact = FACT[dq0];
    for (int r = 0; r < ndiv; r++) fact = co_mul(fact, dq[r] >= 0 ? IFACT[dq[r]] : FACT[-dq[r] - 1]);
    if (nneg == 2) {
        int v1 = dq[s1], v2 = dq[sb];
        Co fs = ((v1 + v2) % 2 == 0) ? fact : co_neg(fact);
        *s2 = co_mul(fs, Mm[(size_t)ROWOF[s1] * NR + ROWOF[sb]]);
        return;
    }
    /* gamma: Av = sum_rho gamma_rho q_rho over the rows S where it is nonzero */
    int S[NR], ns = 0; Co g[NR];
    S[ns] = 0; g[ns++] = H1[dq0];
    for (int r = 0; r < ndiv; r++) {
        if (!dq[r]) continue;
        int n = dq[r] < 0 ? -dq[r] - 1 : dq[r];
        if (!n) continue;
        S[ns] = ROWOF[r]; g[ns++] = co_neg(H1[n]);
    }
    Co s = co_zero();
    if (nneg == 0) {
        *c0 = fact;
        for (int i = 0; i < ns; i++) c1[S[i]] = co_mul(fact, g[i]);
        for (int i = 0; i < ns; i++) {
            const Co *Mr = Mm + (size_t)S[i] * NR;
            Co x = co_zero();
            for (int j = 0; j < ns; j++) co_fma(&x, g[j], Mr[S[j]]);
            co_fma(&s, g[i], x);
        }
        s = co_sub(s, co_mul(H2[dq0], Mm[0]));
        for (int r = 0; r < ndiv; r++) if (dq[r]) { int q = ROWOF[r]; s = co_add(s, co_mul(H2[dq[r]], Mm[(size_t)q * NR + q])); }
        *s2 = co_mul(fact, co_mul(half, s));
    } else {
        int v = dq[s1], q1 = ROWOF[s1];
        Co fs = (v % 2 == 0) ? co_neg(fact) : fact;
        c1[q1] = fs;
        for (int i = 0; i < ns; i++) co_fma(&s, g[i], Mm[(size_t)S[i] * NR + q1]);
        *s2 = co_mul(fs, s);
    }
}
#endif

typedef struct { int nchunk; SM *c0, *s2, **c1; } FPCtx;
static void fp_chunk(void *vctx, int ch, int tid) {
    FPCtx *ctx = vctx;
    int lo = (int)((i64)CAND.n * ch / ctx->nchunk), hi = (int)((i64)CAND.n * (ch + 1) / ctx->nchunk);
    int C[64]; Co c0, c1[NA], s2; int cr;
    for (int s = lo; s < hi; s++) {
        unpack(CAND.key[s], C);
#ifndef CGV_MAJ
        if (RMODE) fp_point_r(C, &c0, c1, &s2); else
#endif
        fp_point(C, &c0, c1, &s2);
        KEY k = CAND.key[s]; int d = CAND.deg[s];
        if (!co_isz(c0)) { int z = sm_add(&ctx->c0[ch], k, d, &cr); ctx->c0[ch].val[z] = c0; }
        if (!co_isz(s2)) { int z = sm_add(&ctx->s2[ch], k, d, &cr); ctx->s2[ch].val[z] = s2; }
        for (int a = 0; a < NA; a++) if (!co_isz(c1[a])) { int z = sm_add(&ctx->c1[a][ch], k, d, &cr); ctx->c1[a][ch].val[z] = c1[a]; }
    }
}

/* ------------------------------------------------------------------------- */
/* Sparse algebra                                                            */
/* ------------------------------------------------------------------------- */
/* a*b truncated at maxdeg, accumulated in the scratch map A. */
static SPoly sp_mul(const SPoly *a, const SPoly *b, SM *A) {
    int cr;
    for (int x = 0; x < a->n; x++) {
        int lim = maxdeg - a->deg[x];
        KEY ki = a->key[x]; Co ai = a->val[x]; int di = a->deg[x];
        for (int y = 0; y < b->n; y++) {
            if (b->deg[y] > lim) break;
            int s = sm_add(A, k_add(ki, b->key[y]), di + b->deg[y], &cr);
            co_fma(&A->val[s], ai, b->val[y]);
        }
    }
    return sm_take(A);
}

/* Reciprocal of a series whose first entry is the (nonzero) constant term. */
static SPoly sp_recip(const SPoly *c0, SM *A) {
    if (c0->n == 0 || !k_eq(c0->key[0], K_ZERO)) die("series has no constant term");
    Co inv0 = co_inv(c0->val[0]);
    Heap H = {0};
    int cr;
    int s0 = sm_add(A, K_ZERO, 0, &cr); hp_push(&H, A->deg, s0);
    while (H.n) {
        int sb = hp_pop(&H, A->deg);
        int db = A->deg[sb];
        Co rb = db == 0 ? inv0 : co_neg(co_mul(A->val[sb], inv0));
        A->val[sb] = rb;
        if (co_isz(rb)) continue;
        KEY kb = A->key[sb];
        for (int y = 1; y < c0->n; y++) {
            int d = db + c0->deg[y];
            if (d > maxdeg) break;
            int s = sm_add(A, k_add(kb, c0->key[y]), d, &cr);
            if (cr) hp_push(&H, A->deg, s);
            co_fma(&A->val[s], c0->val[y], rb);
        }
    }
    free(H.s);
    return sm_take(A);
}

/* ------------------------------------------------------------------------- */
/* Extraction                                                                */
/* ------------------------------------------------------------------------- */
static int NU;          /* union of alpha supports, sorted by (deg, key) */
static KEY *UKEY;
static int *UDEG;
static Co *UALPHA;      /* [NU][h11] */

static int NSH = 1; /* shards of the instanton polynomial, folded in parallel */
static inline int shard_of(KEY k) { return (int)((hash128(k) >> 40) % (u64)NSH); }
typedef struct { KEY *key; int *deg; Co *val; int n, cap; } OutBuf;
static inline void ob_push(OutBuf *o, KEY k, int d, Co v) {
    if (o->n == o->cap) { o->cap = o->cap ? 2 * o->cap : 1024; o->key = xrealloc(o->key, o->cap * sizeof(KEY)); o->deg = xrealloc(o->deg, o->cap * sizeof(int)); o->val = xrealloc(o->val, o->cap * sizeof(Co)); }
    o->key[o->n] = k; o->deg[o->n] = d; o->val[o->n] = v; o->n++;
}
static void ob_free(OutBuf *o) { free(o->key); free(o->deg); free(o->val); }

/* Extraction outputs are subtracted from the sharded instanton polynomial.
 * Buffers are flushed under a per-shard lock once they reach EMIT_FLUSH
 * entries, so memory stays bounded; only degrees above the layer being
 * extracted change, so folding early is safe. */
#ifndef EMIT_FLUSH
#define EMIT_FLUSH 1024
#endif
static SM *EX_IS; static Heap *EX_IH; static pthread_mutex_t *EX_LK;
#ifndef FOLD_B
#define FOLD_B 256
#endif
/* I -= o, in chunks: hash the chunk and prefetch its hash slots, then (slots now cached) prefetch the key behind
 * each first-probe tag match, then apply in order. One entry at a time, hash slot -> key is a dependent miss chain;
 * this keeps a whole chunk's misses in flight per pass. Prefetches are hints only (the apply pass probes as before). */
static void sm_sub_buf(SM *I, Heap *IH, OutBuf *o) {
    u64 x[FOLD_B]; int cr;
    for (int z0 = 0; z0 < o->n; z0 += FOLD_B) {
        int nb = o->n - z0 < FOLD_B ? o->n - z0 : FOLD_B;
        const KEY *k = o->key + z0;
        for (int z = 0; z < nb; z++) { x[z] = hash128(k[z]); __builtin_prefetch(&I->hs[x[z] & I->mask], 0, 3); }
        for (int z = 0; z < nb; z++) {
            u64 e = I->hs[x[z] & I->mask];
            if (e >> 32 == x[z] >> 32 && e) { uint32_t s = (uint32_t)e - 1; __builtin_prefetch(&I->key[s], 0, 3); }
        }
        for (int z = 0; z < nb; z++) {
            int s = sm_add_x(I, k[z], x[z], o->deg[z0 + z], &cr);
            if (cr) hp_push(IH, I->deg, s);
            I->val[s] = co_sub(I->val[s], o->val[z0 + z]);
        }
    }
    o->n = 0;
}
static void ex_fold(OutBuf *o, int sh) { sm_sub_buf(&EX_IS[sh], &EX_IH[sh], o); }
static inline void emit_out(OutBuf *out, KEY k, int d, Co v) {
    int sh = shard_of(k);
    OutBuf *o = &out[sh];
    ob_push(o, k, d, v);
    if (o->n >= EMIT_FLUSH && EX_LK) {
        pthread_mutex_lock(&EX_LK[sh]);
        ex_fold(o, sh);
        pthread_mutex_unlock(&EX_LK[sh]);
    }
}

/* L2 of cpu0 (a performance core; on Apple silicon the per-cluster value) in bytes, 0 if unknown */
static long l2_size(void) {
#if defined(__APPLE__)
    long v = 0; size_t n = sizeof(v);
    if (!sysctlbyname("hw.perflevel0.l2cachesize", &v, &n, NULL, 0) && v > 0) return v;
    n = sizeof(v); if (!sysctlbyname("hw.l2cachesize", &v, &n, NULL, 0) && v > 0) return v;
#elif defined(__linux__)
    for (int i = 0; i < 8; i++) {
        char p[96], b[32]; int lev = 0; long sz = 0; char u = 0; FILE *f;
        snprintf(p, sizeof p, "/sys/devices/system/cpu/cpu0/cache/index%d/level", i);
        if (!(f = fopen(p, "r"))) break;
        if (fscanf(f, "%d", &lev) != 1) lev = 0;
        fclose(f);
        snprintf(p, sizeof p, "/sys/devices/system/cpu/cpu0/cache/index%d/type", i);
        if (lev != 2 || !(f = fopen(p, "r"))) continue;
        if (fscanf(f, "%31s", b) != 1 || !strcmp(b, "Instruction")) { fclose(f); continue; }
        fclose(f);
        snprintf(p, sizeof p, "/sys/devices/system/cpu/cpu0/cache/index%d/size", i);
        if (!(f = fopen(p, "r"))) continue;
        if (fscanf(f, "%ld%c", &sz, &u) >= 1) { fclose(f); return u == 'K' ? sz << 10 : u == 'M' ? sz << 20 : sz; }
        fclose(f);
    }
#endif
    return 0;
}
/* packed records when L2 < 1 MB: there even the split key array of a common 2^15-slot table (256 KB) overflows L2, so
 * the separate value line is a second L3 access (i5-10600K, 256 KB: -16% at 8 threads); with a larger L2 the key
 * arrays stay resident and records (4x the bytes per slot) lose (M1, Core Ultra: +4 to +6%). Unknown: split.
 * CGV_LT_LAYOUT=packed|split overrides */
static int said_ks, said_pk;   /* last logged key width / layout: logged once per run, again only if it changes */
static void lt_layout(void) {
    const char *e = getenv("CGV_LT_LAYOUT");
    if (e && strcmp(e, "packed") && strcmp(e, "split")) {
        if (*e && said_pk < 0) LOG("cgv: note: CGV_LT_LAYOUT=%s is not packed or split: ignored\n", e);
        e = NULL;
    }
    long l2 = e ? 0 : l2_size();
    LT_PK = e ? !strcmp(e, "packed") : l2 > 0 && l2 < (1L << 20);
    if (said_pk != LT_PK) {
        if (e) LOG("  level-table layout: %s (CGV_LT_LAYOUT)\n", LT_PK ? "packed" : "split");
        else LOG("  level-table layout: %s (L2 %ld KB)\n", LT_PK ? "packed" : "split", l2 >> 10);
    }
    said_pk = LT_PK;
}

/* level-table key encoding for this pass (needs UKEY): 64-bit keys when they fit, else 128-bit */
/* level keys in a basis of the L terms' lattice (see lt_key_lb), if their box fits maxbits (63 or 127) bits: sets the
 * LBL_* fields, LT_KS, LT_HXK and ULH; 0 if it does not fit */
/* the lattice part depends only on the L terms, the same in every prime pass: computed once per run and maxbits
 * (LBC[maxbits > 63]; a signature of UKEY/UDEG guards the reuse). Any overflow: not used (0) */
typedef struct { int done, ok, nu, k, t0, b[64]; u64 sig; i64 wl[64], lw[64], wlo; int *bj, *cu; } LbCache;
static LbCache LBC[2];
static void lbc_reset(void) { for (int i = 0; i < 2; i++) { free(LBC[i].bj); free(LBC[i].cu); memset(&LBC[i], 0, sizeof(LBC[i])); } }
static u64 lbc_sig(void) {
    u64 h = (u64)NU * 0x9E3779B97F4A7C15ull ^ (u64)maxdeg;
    for (int u = 0; u < NU; u++) {
        u64 w[sizeof(KEY) / 8]; memcpy(w, &UKEY[u], sizeof(KEY));
        for (size_t i = 0; i < sizeof(KEY) / 8; i++) h = (h ^ w[i]) * 0xD6E8FEB86659FD93ull;
        h = (h ^ (u64)(unsigned)UDEG[u]) * 0xD6E8FEB86659FD93ull;
    }
    return h;
}
static void lbc_compute(int maxbits, u64 sig) {
    LbCache *C = &LBC[maxbits > 63];
    free(C->bj); free(C->cu); memset(C, 0, sizeof(*C));
    C->done = 1; C->nu = NU; C->sig = sig;
    Lat L; lat_init(&L, HC);
    int x[64], ok = 1;
    for (int u = 0; u < NU && ok; u++) { unpack(UKEY[u], x); ok = !lat_add(&L, x); }
    int k = L.rk, *pc = NULL;
    i64 *B = NULL, *T = NULL, *v = xmalloc(HC * sizeof(i64));
    ok = ok && k > 0 && !lat_finish(&L, &B, &pc);
    lat_free(&L);
    int *cu = xmalloc((size_t)(NU ? NU : 1) * (k ? k : 1) * sizeof(int));
    for (int u = 0; u < NU && ok; u++) { unpack(UKEY[u], x); ok = !lat_coords(B, pc, k, HC, x, cu + (size_t)u * k, v); }
    ok = ok && !lat_lll(B, k, HC, &T);
    i64 hi[64], *lw = C->lw, *wl = C->wl; int *b = C->b, tot = 0, t0 = -1;
    for (int j = 0; j < k && ok; j++) {   /* basis entries and degrees must fit int (pack, LBL_f[].w) */
        hi[j] = lw[j] = 0; wl[j] = 0;
        for (int t = 0; t < HC; t++) ok &= B[(size_t)j * HC + t] <= 0x3fffffff && B[(size_t)j * HC + t] >= -0x3fffffff;
        for (int t = 0; t < HC && ok; t++) { i64 p; ok = !__builtin_mul_overflow(B[(size_t)j * HC + t], (i64)W[t], &p) && !__builtin_add_overflow(wl[j], p, &wl[j]); }
        ok &= wl[j] <= INT_MAX && wl[j] >= -INT_MAX;
    }
    for (int u = 0; u < NU && ok; u++) {   /* a sum of L terms of total degree <= maxdeg has c_j in [lw, hi] */
        int *c = cu + (size_t)u * k, d = UDEG[u];
        if (lat_apply(c, T, k)) { ok = 0; break; }
        for (int j = 0; j < k; j++) {
            i64 m = (i64)maxdeg * c[j], f = m >= 0 ? m / d : -((-m + d - 1) / d), cc = m >= 0 ? (m + d - 1) / d : -(-m / d);
            if (f > hi[j]) hi[j] = f;
            if (cc < lw[j]) lw[j] = cc;
        }
    }
    for (int j = 0; j < k && ok; j++) {
        b[j] = 0; while (b[j] < 63 && ((i64)1 << b[j]) < hi[j] - lw[j] + 1) b[j]++;
        tot += b[j];
        if (wl[j] && (t0 < 0 || b[j] > b[t0] || (b[j] == b[t0] && (wl[j] == 1 || wl[j] == -1)))) t0 = j;
    }
    if (tot <= maxbits) t0 = -1;
    if (ok && (t0 < 0 ? tot : tot - b[t0]) > maxbits) {
        if (maxbits > 63) LOG("  level box: %d bits in a basis of the L-term lattice (rank %d; %d without the degree-fixed coordinate)\n", tot, k, t0 < 0 ? tot : tot - b[t0]);
        ok = 0;
    }
    for (int j = 0; j < k && ok; j++) {   /* sum_j wl lw, |.| < 2^62 (lt_key_lb's e - LBL_wlo) */
        i64 p;
        if (j == t0) continue;
        ok = !__builtin_mul_overflow(wl[j], lw[j], &p) && !__builtin_add_overflow(C->wlo, p, &C->wlo) && C->wlo < ((i64)1 << 62) && C->wlo > -((i64)1 << 62);
    }
    if (ok) {
        C->ok = 1; C->k = k; C->t0 = t0; C->cu = cu; cu = NULL;
        C->bj = xmalloc((size_t)k * HC * sizeof(int));
        for (size_t i = 0; i < (size_t)k * HC; i++) C->bj[i] = (int)B[i];
    }
    free(B); free(pc); free(v); free(cu); free(T);
}
/* level keys in a basis of the L terms' lattice (see lt_key_lb), if their box fits maxbits (63 or 127) bits: sets the
 * LBL_* fields, LT_KS, LT_HXK and ULH; 0 if it does not fit */
static int lt_setup_lb(int maxbits) {
    LbCache *C = &LBC[maxbits > 63];
    u64 sig = lbc_sig();
    if (!C->done || C->nu != NU || C->sig != sig) lbc_compute(maxbits, sig);
    if (!C->ok) return 0;
    int k = C->k, t0 = C->t0;
    LT_KS = maxbits > 63; LBL_t0 = t0; LBL_w0 = t0 >= 0 ? (int)C->wl[t0] : 1; LBL_nf = 0; LBL_wlo = C->wlo;
    for (int j = 0; j < k; j++) LBL_k[j] = pack(C->bj + (size_t)j * HC);
    LBL_k0 = K_ZERO; LBL_plo = 0;
    for (int j = 0, o = 0; j < k; j++) {
        if (j == t0) continue;
        LBL_f[LBL_nf].j = j; LBL_f[LBL_nf].off = o; LBL_f[LBL_nf].b = C->b[j]; LBL_f[LBL_nf].w = (int)C->wl[j]; LBL_nf++;
        LBL_k0 = k_addmul(LBL_k0, LBL_k[j], C->lw[j]); LBL_plo += (u128)(i128)C->lw[j] << o;
        o += C->b[j];
    }
    ULH = xmalloc((NU ? NU : 1) * sizeof(u128));
    LT_HMI64 = LT_HM64; for (int i = 0; i < 6; i++) LT_HMI64 *= 2 - LT_HM64 * LT_HMI64;
    for (int u = 0; u < NU; u++) {
        const int *c = C->cu + (size_t)u * k; u128 p = 0;
        for (int f = 0; f < LBL_nf; f++) p += (u128)(i128)c[LBL_f[f].j] << LBL_f[f].off;
        ULH[u] = LT_KS ? p * LT_HM : (u128)((u64)p * LT_HM64) << 64;
    }
    /* phi = plo + 2^(maxbits) is outside the box: never a key */
    LT_HXK = LT_KS ? (LBL_plo + ((u128)1 << 127)) * LT_HM : (u128)(((u64)LBL_plo + ((u64)1 << 63)) * LT_HM64) << 64;
    return 1;
}
static void lt_setup(void) {
    int x[64], b[64], lo[64], ok = 1, t0 = -1, B = 0;
    i64 hi[64], lw[64];
    for (int t = 0; t < HC; t++) hi[t] = lw[t] = 0;
    for (int u = 0; u < NU && ok; u++) {
        int d = UDEG[u];
        unpack(UKEY[u], x);
        if (d < 1 || dot(x, W) != d) { ok = 0; break; }
        for (int t = 0; t < HC; t++) {   /* a sum of L terms of total degree <= maxdeg has x_t in [lw, hi] */
            i64 m = (i64)maxdeg * x[t], f = m >= 0 ? m / d : -((-m + d - 1) / d), c = m >= 0 ? (m + d - 1) / d : -(-m / d);
            if (f > hi[t]) hi[t] = f;
            if (c < lw[t]) lw[t] = c;
        }
    }
    for (int t = 0; t < HC && ok; t++) {
        b[t] = 0; while (b[t] < 63 && ((i64)1 << b[t]) < hi[t] - lw[t] + 1) b[t]++;
        lo[t] = (int)lw[t]; B += b[t]; ok &= b[t] <= (KBV ? KBW[t] : keybits);
        if (W[t] && (t0 < 0 || b[t] > b[t0] || (b[t] == b[t0] && abs(W[t]) == 1))) t0 = t;
    }
    if (B <= 63) t0 = -1;   /* x_t0 from the degree only when needed: it costs a division per decode */
    LT_LB = 0;
    if (!(ok && (t0 < 0 ? B : B - b[t0]) <= 63) && lt_setup_lb(63)) LT_LB = 1;
    else if (ok && (t0 < 0 ? B : B - b[t0]) <= 63) {
        LT_KS = 0; L64_t0 = t0; L64_nf = L64_nlo = 0; L64_wlo = 0;
        for (int t = 0, o = 0; t < HC; t++) {   /* b_t <= the key field width, so its position ps >= off; ps < 64 come first */
            int ps = KBV ? KBO[t] : t * keybits;
            if (t == t0) { L64_p0 = ps; L64_w0 = W[t]; x[t] = 0; continue; }
            int d = ps - o;
            L64_f[L64_nf].t = t; L64_f[L64_nf].off = o; L64_f[L64_nf].w = W[t]; L64_f[L64_nf].M = (((u64)1 << b[t]) - 1) << o;
#ifdef CGV_WIDE
            L64_f[L64_nf].ps = ps;
#endif
            L64_f[L64_nf].d = d < 64 ? d : d - 64; L64_f[L64_nf].mul = d < 64 ? (u64)1 << d : 0;
            if (d < 64) L64_nlo = L64_nf + 1;
            L64_nf++;
            o += b[t]; x[t] = lo[t]; L64_wlo += (i64)W[t] * lo[t];
        }
        L64_klo = pack(x);
        L64_plo = l64_phi(lo);
        LT_HMI64 = LT_HM64; for (int i = 0; i < 6; i++) LT_HMI64 *= 2 - LT_HM64 * LT_HMI64;
        LT_HXK = (u128)((L64_plo + ((u64)1 << 63)) * LT_HM64) << 64;  /* phi = plo + 2^63 is outside the box: never a key */
    } else {
        LT_KS = 1; LT_HXK = lt_h(k_lo128(LT_XK));
#ifdef CGV_WIDE
        /* the wide key itself does not fit the 128-bit hash: pack the level box instead (as the 64-bit keys, in 127 bits) */
        if (B <= 127) t0 = -1;
        if (!ok || (t0 < 0 ? B : B - b[t0]) > 127) {
            if (lt_setup_lb(127)) { LT_LB = 1; goto ulh; }
            LOG("  level box: %d bits in key coordinates (%d without the degree-fixed coordinate)\n", B, t0 < 0 ? B : B - b[t0]);
            die("wide keys: a level table's offsets need more than 127 bits");
        }
        L128_t0 = t0; L128_nf = 0; L128_wlo = 0;
        for (int t = 0, o = 0; t < HC; t++) {
            int ps = KBV ? KBO[t] : t * keybits;
            if (t == t0) { L128_p0 = ps; L128_w0 = W[t]; x[t] = 0; continue; }
            L128_f[L128_nf].t = t; L128_f[L128_nf].off = o; L128_f[L128_nf].b = b[t]; L128_f[L128_nf].w = W[t]; L128_f[L128_nf].ps = ps; L128_nf++;
            o += b[t]; x[t] = lo[t]; L128_wlo += (i64)W[t] * lo[t];
        }
        L128_klo = pack(x);
        L128_plo = l128_phi(lo);
        LT_HXK = (L128_plo + ((u128)1 << 127)) * LT_HM;   /* outside the box: never a key */
#endif
    }
#ifdef CGV_WIDE
ulh:
#endif
    if (!LT_LB) {   /* (lt_setup_lb sets ULH) */
        ULH = xmalloc((NU ? NU : 1) * sizeof(u128));
        for (int u = 0; u < NU; u++) ULH[u] = lt_hk(UKEY[u]);
    }
    if (said_ks != LT_KS + 2 * LT_LB) LOG("  level-table keys: %d bits%s\n", LT_KS ? 128 : 64, LT_LB ? " (in a basis of the L-term lattice)" : "");
    said_ks = LT_DEC = LT_KS + 2 * LT_LB;
    lt_layout();
}

/* Per-thread workspace for exp(C.alpha), organized by degree: one small open-
 * addressing table per degree level. All sources of degree e are final at the
 * same time, so they are scattered together, one target degree g at a time;
 * each (e, g) sweep then works inside the single level-g table, which (unlike
 * the whole curve's support) mostly fits in L2. */
/* keys and values in separate arrays: a probe reads only keys (8 or 16 B each, no padding), the value once on a hit.
 * Keys store lt_h(key) ^ LT_HXK, so an empty (zeroed) slot is xk == 0 and needs no occupied flag; 2^LT_KS words per slot,
 * the high word first (packed layout: xk is the record array, aligned in its allocation, and vl is unused) */
typedef struct { u64 *xk; Co *vl; u64 mask; int n; int *used; int ucap; int sh; } LT;
typedef struct {
    LT **lt;                          /* per degree 0..maxdeg, allocated on first use */
    int *act; int nact, actcap;       /* levels used by the current curve */
    int *lvh; int nlvh, lvhcap;       /* min-heap of levels holding pending entries */
    u128 *skey; Co *sval; int scap;   /* nonzero sources of the current level, as lt_h(key) */
    u128 *lkey; int *ldeg; Co *lval;  /* nonzero terms of deg(u) (C.alpha)[u], degree-sorted; lkey holds lt_h(u) */
    int *gdeg, *gst; int ng;          /* degree groups of L: degree gdeg[q] is lkey[gst[q] .. gst[q+1]) */
    OutBuf *out;                      /* per shard */
    i64 n_inner, n_slots;             /* profiling counters */
    int nopush;                       /* CS mode: levels are registered by the owner, not in lt_find */
    unsigned char *inh;               /* CS mode: level already in the heap */
} WS;
static void lt_release(LT *t) { free(LT_PK && t->xk ? ((void **)t->xk)[-1] : (void *)t->xk); free(t->vl); free(t->used); }
static void lt_alloc(LT *t, u64 cap) {
    t->xk = xcalloc(cap << LT_KS, sizeof(u64)); t->vl = xmalloc(cap * sizeof(Co)); t->mask = cap - 1; t->sh = 64 - __builtin_ctzll(cap); /* xk = 0: empty */
}
static void ws_init_sized(WS *w, int nlterms) {
    memset(w, 0, sizeof(*w));
    w->lt = xcalloc(maxdeg + 1, sizeof(LT *));
    w->actcap = 64; w->act = xmalloc(w->actcap * sizeof(int));
    w->lvhcap = 64; w->lvh = xmalloc(w->lvhcap * sizeof(int));
    w->scap = 64; w->skey = xmalloc(w->scap * sizeof(u128)); w->sval = xmalloc(w->scap * sizeof(Co));
    int n = nlterms ? nlterms : 1;
    w->lkey = xmalloc(n * sizeof(u128)); w->ldeg = xmalloc(n * sizeof(int)); w->lval = xmalloc(n * sizeof(Co));
    w->gdeg = xmalloc((n + 1) * sizeof(int)); w->gst = xmalloc((n + 2) * sizeof(int));
    w->out = xcalloc(NSH, sizeof(OutBuf));
}
static void ws_init(WS *w) { ws_init_sized(w, NU); }
static inline void lvh_push(WS *w, int v) {
    if (w->nlvh == w->lvhcap) { w->lvhcap *= 2; w->lvh = xrealloc(w->lvh, w->lvhcap * sizeof(int)); }
    int i = w->nlvh++;
    while (i > 0) { int p = (i - 1) >> 1; if (w->lvh[p] <= v) break; w->lvh[i] = w->lvh[p]; i = p; }
    w->lvh[i] = v;
}
static inline int lvh_pop(WS *w) {
    int top = w->lvh[0], last = w->lvh[--w->nlvh], i = 0, n = w->nlvh;
    for (;;) {
        int c = 2 * i + 1;
        if (c >= n) break;
        if (c + 1 < n && w->lvh[c + 1] < w->lvh[c]) c++;
        if (w->lvh[c] >= last) break;
        w->lvh[i] = w->lvh[c]; i = c;
    }
    if (n) w->lvh[i] = last;
    return top;
}
static void ws_free(WS *w) {
    for (int d = 0; d <= maxdeg; d++) if (w->lt[d]) { lt_release(w->lt[d]); free(w->lt[d]); }
    free(w->lt); free(w->act); free(w->lvh); free(w->skey); free(w->sval); free(w->lkey); free(w->ldeg); free(w->lval);
    free(w->gdeg); free(w->gst); free(w->inh);
    for (int s = 0; s < NSH; s++) ob_free(&w->out[s]);
    free(w->out);
}
static inline u128 lt_xk(const u64 *xk, u64 h) { return LT_KS ? (u128)xk[2 * h] << 64 | xk[2 * h + 1] : (u128)xk[h] << 64; }
static inline void lt_setxk(u64 *xk, u64 h, u128 x) { if (LT_KS) { xk[2 * h] = (u64)(x >> 64); xk[2 * h + 1] = (u64)x; } else xk[h] = (u64)(x >> 64); }
static void lt_grow(LT *t) {
    u64 *ok = t->xk; Co *ov = t->vl; u64 ocap = t->mask + 1;
    lt_alloc(t, 2 * ocap);
    for (int z = 0; z < t->n; z++) {
        int o = t->used[z];
        u128 x = lt_xk(ok, o);
        u64 h = LT_SLOT(t, x ^ LT_HXK);
        while (lt_xk(t->xk, h)) h = (h + 1) & t->mask;
        lt_setxk(t->xk, h, x); t->vl[h] = ov[o];
        t->used[z] = (int)h;
    }
    free(ok); free(ov);
}
static inline Co *lt_find_h(WS *w, int d, u128 H) {
    LT *t = w->lt[d];
    if (!t) { t = w->lt[d] = xcalloc(1, sizeof(LT)); }
    if (!t->xk) { lt_alloc(t, 256); t->ucap = 128; t->used = xmalloc(t->ucap * sizeof(int)); }
    u64 h = LT_SLOT(t, H); u128 xk = H ^ LT_HXK; u64 xh = (u64)(xk >> 64), xl = (u64)xk;
    if (LT_KS) for (;;) {
        u64 e0 = t->xk[2 * h], e1 = t->xk[2 * h + 1];
        if (e0 == xh && e1 == xl) return &t->vl[h];
        if (!(e0 | e1)) break;
        h = (h + 1) & t->mask;
    } else for (;;) {
        u64 e = t->xk[h];
        if (e == xh) return &t->vl[h];
        if (!e) break;
        h = (h + 1) & t->mask;
    }
#ifndef LT_LFN
#define LT_LFN 3   /* level tables grow above LT_LFN/LT_LFD full: 3/4 (was 1/2; fuller tables stay in cache) */
#define LT_LFD 4
#endif
    if (LT_LFD * (u64)(t->n + 1) > LT_LFN * (u64)t->mask) { lt_grow(t); return lt_find_h(w, d, H); }
    lt_setxk(t->xk, h, xk); t->vl[h] = co_zero();
    Co *e = &t->vl[h];
    if (t->n == 0 && !w->nopush) {
        if (w->nact == w->actcap) { w->actcap *= 2; w->act = xrealloc(w->act, w->actcap * sizeof(int)); }
        w->act[w->nact++] = d;
        lvh_push(w, d);
    }
    if (t->n == t->ucap) { t->ucap *= 2; t->used = xrealloc(t->used, t->ucap * sizeof(int)); }
    t->used[t->n++] = (int)h;
    return e;
}
#ifndef LT_PF
#define LT_PF 4   /* prefetch distance in the scatter loop (0 = off) */
#endif
/* the ns sources of w times the terms j0..j1 into level g. The table's fields live in locals: the stores of co_fma
 * (u64) may alias t->mask, so through t they would be reloaded for every term. Only a new key (lt_find_h, which may
 * grow the table) refreshes them. */
static inline __attribute__((always_inline)) void scatter_ks(WS *w, int g, int j0, int j1, int ns, const int ks) {
    const u128 *lkey = w->lkey, *skey = w->skey; const Co *lval = w->lval, *sval = w->sval;
    if (!w->lt[g] || !w->lt[g]->xk) lt_find_h(w, g, skey[0] + lkey[j0]);   /* allocate (the key is added anyway) */
    LT *t = w->lt[g];
    u64 *xk = t->xk; Co *vl = t->vl; u64 mask = t->mask; int sh = t->sh;
    u64 xkh = (u64)(LT_HXK >> 64), xkl = (u64)LT_HXK;
    for (int z = 0; z < ns; z++) {
        u128 kb = skey[z]; Co fb = sval[z]; u64 kh = (u64)(kb >> 64);
        for (int j = j0; j < j1; j++) {
            /* 64-bit keys: the hash is the high word only (low words are 0) */
#if LT_PF > 0
            if (j + LT_PF < j1) __builtin_prefetch(&xk[(ks ? (u64)((kb + lkey[j + LT_PF]) >> 64) : kh + (u64)(lkey[j + LT_PF] >> 64)) >> sh << ks], 1, 1);
#endif
            u128 H; u64 hh, xl = 0;
            if (ks) { H = kb + lkey[j]; hh = (u64)(H >> 64); xl = (u64)H ^ xkl; }
            else { hh = kh + (u64)(lkey[j] >> 64); H = (u128)hh << 64; }
            u64 xh = hh ^ xkh, h = hh >> sh;
            Co *c;
            for (;;) {
                if (ks) {
                    u64 e0 = xk[2 * h], e1 = xk[2 * h + 1];
                    if (e0 == xh && e1 == xl) { c = &vl[h]; break; }
                    if (!(e0 | e1)) { c = lt_find_h(w, g, H); xk = t->xk; vl = t->vl; mask = t->mask; sh = t->sh; break; }
                } else {
                    u64 e = xk[h];
                    if (e == xh) { c = &vl[h]; break; }
                    if (!e) { c = lt_find_h(w, g, H); xk = t->xk; vl = t->vl; mask = t->mask; sh = t->sh; break; }
                }
                h = (h + 1) & mask;
            }
            co_fma(c, lval[j], fb);
        }
    }
}
static void scatter_group(WS *w, int g, int j0, int j1, int ns) { if (LT_KS) scatter_ks(w, g, j0, j1, ns, 1); else scatter_ks(w, g, j0, j1, ns, 0); }
#ifndef LT_KEEP
#define LT_KEEP 32768 /* level tables larger than this are released after each curve */
#endif
static void ws_clear(WS *w) {
    for (int a = 0; a < w->nact; a++) {
        LT *t = w->lt[w->act[a]];
        if (!t) continue;
        if (t->mask + 1 > LT_KEEP) { lt_release(t); memset(t, 0, sizeof(*t)); continue; }
        for (int z = 0; z < t->n; z++) lt_setxk(t->xk, t->used[z], 0);
        t->n = 0;
    }
    w->nact = 0; w->nlvh = 0;
}
/* degree groups of the L terms in w (ldeg sorted) */
static void ws_groups(WS *w, int nl) {
    w->ng = 0;
    for (int j = 0; j < nl; j++)
        if (j == 0 || w->ldeg[j] != w->ldeg[j - 1]) { w->gdeg[w->ng] = w->ldeg[j]; w->gst[w->ng] = j; w->ng++; }
    w->gst[w->ng] = nl;
}

/* Subtract s * z^C exp(C.alpha) from I (recorded in w's output), truncated
 * at maxdeg. exp via the Euler recurrence deg(m) f[m] = sum_{a+b=m} deg(a) L[a] f[b]. */
static void apply_curve(WS *w, KEY kc, int dC, Co s) {
    int T = maxdeg - dC;
    int C[64];
    unpack(kc, C);
    Co cm[64]; int nzt[64], nnz = curve_coords(C, cm, nzt);
    int nl = 0;
    for (int u = 0; u < NU && UDEG[u] <= T; u++) {
        Co v = co_zero();
        const Co *al = UALPHA + (size_t)u * NA;
        for (int q = 0; q < nnz; q++) co_fma(&v, cm[q], al[nzt[q]]);
        if (co_isz(v)) continue;
        w->lkey[nl] = ULH[u]; w->ldeg[nl] = UDEG[u]; w->lval[nl] = co_mul(v, DEGC[UDEG[u]]); nl++;
    }
    ws_groups(w, nl);
    *lt_find_h(w, 0, 0) = co_int(1);
    i64 inner = 0, slots = 0;
    while (w->nlvh) {
        int e = lvh_pop(w);
        LT *le = w->lt[e];
        slots += le->n;
        /* finalize this level; collect its nonzero sources */
        if (w->scap < le->n) { w->scap = 2 * le->n; w->skey = xrealloc(w->skey, w->scap * sizeof(u128)); w->sval = xrealloc(w->sval, w->scap * sizeof(Co)); }
        int ns = 0;
        for (int z = 0; z < le->n; z++) {
            int xi = le->used[z];
            Co f = e == 0 ? le->vl[xi] : co_mul(le->vl[xi], INV[e]);
            if (co_isz(f)) continue;
            u128 hk0 = lt_xk(le->xk, xi) ^ LT_HXK; KEY ko = k_add(kc, lt_key(hk0, e));
            emit_out(w->out, ko, dC + e, co_mul(s, f));
            w->skey[ns] = hk0; w->sval[ns] = f; ns++;
        }
        /* scatter into each higher level g = e + k, over the degrees k present in L */
        for (int q = 0; q < w->ng && w->gdeg[q] <= T - e; q++) {
            int j0 = w->gst[q], j1 = w->gst[q + 1];
            if (ns) scatter_group(w, e + w->gdeg[q], j0, j1, ns);
            inner += (i64)ns * (j1 - j0);
        }
    }
    w->n_inner += inner;
    w->n_slots += slots;
    ws_clear(w);
}

/* ---- small layers: all curves advanced level-synchronously, each with its own
 * per-degree tables, so (curve, target level) pairs are independent tasks ---- */
typedef struct {
    KEY kc; int dC, T; Co s;
    WS w;                             /* reuses the table machinery; w.out = emissions */
    int ns;
} CS;
static Co *lt_find_hp(WS *w, int d, u128 H);
static void cs_setup(CS *c) {
    WS *w = &c->w;
    int nuT;
    { int lo = 0, hi = NU; while (lo < hi) { int mid = (lo + hi) / 2; if (UDEG[mid] <= c->T) lo = mid + 1; else hi = mid; } nuT = lo; }
    ws_init_sized(w, nuT);
    w->nopush = 1;
    w->inh = xcalloc(maxdeg + 1, 1);
    int C[64];
    unpack(c->kc, C);
    Co cm[64]; int nzt[64], nnz = curve_coords(C, cm, nzt);
    int nl = 0;
    for (int u = 0; u < nuT; u++) {
        Co v = co_zero();
        const Co *al = UALPHA + (size_t)u * NA;
        for (int q = 0; q < nnz; q++) co_fma(&v, cm[q], al[nzt[q]]);
        if (co_isz(v)) continue;
        w->lkey[nl] = ULH[u]; w->ldeg[nl] = UDEG[u]; w->lval[nl] = co_mul(v, DEGC[UDEG[u]]); nl++;
    }
    ws_groups(w, nl);
    if (LT_PK) *lt_find_hp(w, 0, 0) = co_int(1);
    else *lt_find_h(w, 0, 0) = co_int(1);
    w->inh[0] = 1; lvh_push(w, 0);
}
static inline int cs_next(const CS *c) { return c->w.nlvh ? c->w.lvh[0] : -1; }
static void cs_finalize(CS *c, int e) {
    WS *w = &c->w;
    c->ns = 0;
    if (cs_next(c) != e) return;
    lvh_pop(w);
    LT *le = w->lt[e];
    if (!le || !le->n) return;
    if (w->scap < le->n) { w->scap = 2 * le->n; w->skey = xrealloc(w->skey, w->scap * sizeof(u128)); w->sval = xrealloc(w->sval, w->scap * sizeof(Co)); }
    int ns = 0;
    for (int z = 0; z < le->n; z++) {
        int xi = le->used[z];
        Co f = e == 0 ? le->vl[xi] : co_mul(le->vl[xi], INV[e]);
        if (co_isz(f)) continue;
        u128 hk0 = lt_xk(le->xk, xi) ^ LT_HXK; KEY ko = k_add(c->kc, lt_key(hk0, e));
        emit_out(w->out, ko, c->dC + e, co_mul(c->s, f));
        w->skey[ns] = hk0; w->sval[ns] = f; ns++;
    }
    c->ns = ns;
    /* this level is done: free its table (it is never touched again) */
    w->n_slots += le->n;
    lt_release(le); free(le); w->lt[e] = NULL;
}
/* sources of level e times the L terms of degree group q, into level e + gdeg[q] */
static void cs_scatter(CS *c, int e, int q) {
    WS *w = &c->w;
    if (c->ns) scatter_group(w, w->gdeg[q] + e, w->gst[q], w->gst[q + 1], c->ns);
}
/* ---- packed layout (LT_PK, chosen in lt_layout): one record {key words, value} per slot, so a hit finds its value
 * on the line it probed (32 B with 64-bit keys and NL = 3, two per 64-B aligned line). Copies of the split functions
 * above, kept separate so the split code compiles as before (a version folded into one implementation with a constant
 * layout flag compiled the split hot loop differently); LT.xk is the aligned record array, with the allocation's
 * address in the word before it. ---- */
#define LT_AL(x) (((x) + _Alignof(Co) - 1) / _Alignof(Co) * _Alignof(Co))
#define LT_VOK(ks) LT_AL(sizeof(u64) << (ks))         /* value offset in a record */
#define LT_RSK(ks) LT_AL(LT_VOK(ks) + sizeof(Co))     /* record size */
static inline u64 *ltp_k(const LT *t, u64 h, const int ks) { return (u64 *)((unsigned char *)t->xk + h * LT_RSK(ks)); }
static inline Co *ltp_v(const LT *t, u64 h, const int ks) { return (Co *)((unsigned char *)t->xk + h * LT_RSK(ks) + LT_VOK(ks)); }
static inline u128 ltp_xk(const LT *t, u64 h) { const u64 *k = ltp_k(t, h, LT_KS); return LT_KS ? (u128)k[0] << 64 | k[1] : (u128)k[0] << 64; }
static inline void ltp_setxk(LT *t, u64 h, u128 x) { u64 *k = ltp_k(t, h, LT_KS); k[0] = (u64)(x >> 64); if (LT_KS) k[1] = (u64)x; }
static void lt_alloc_p(LT *t, u64 cap) {
    unsigned char *m = xcalloc(cap * LT_RSK(LT_KS) + 64 + sizeof(void *), 1), *a = m + sizeof(void *);
    t->vl = NULL; t->xk = (u64 *)(a + (-(uintptr_t)a & 63)); ((void **)t->xk)[-1] = m; t->mask = cap - 1; t->sh = 64 - __builtin_ctzll(cap);
}
static void lt_grow_p(LT *t) {
    LT o = *t;
    lt_alloc_p(t, 2 * (o.mask + 1));
    for (int z = 0; z < t->n; z++) {
        int s = t->used[z];
        u128 x = ltp_xk(&o, s);
        u64 h = LT_SLOT(t, x ^ LT_HXK);
        while (ltp_xk(t, h)) h = (h + 1) & t->mask;
        memcpy(ltp_k(t, h, LT_KS), ltp_k(&o, s, LT_KS), LT_RSK(LT_KS));
        t->used[z] = (int)h;
    }
    free(((void **)o.xk)[-1]);
}
static Co *lt_find_hp(WS *w, int d, u128 H) {
    LT *t = w->lt[d];
    if (!t) { t = w->lt[d] = xcalloc(1, sizeof(LT)); }
    if (!t->xk) { lt_alloc_p(t, 256); t->ucap = 128; t->used = xmalloc(t->ucap * sizeof(int)); }
    u64 h = LT_SLOT(t, H); u128 xk = H ^ LT_HXK; u64 xh = (u64)(xk >> 64), xl = (u64)xk;
    if (LT_KS) for (;;) {
        u64 e0 = ltp_k(t, h, 1)[0], e1 = ltp_k(t, h, 1)[1];
        if (e0 == xh && e1 == xl) return ltp_v(t, h, 1);
        if (!(e0 | e1)) break;
        h = (h + 1) & t->mask;
    } else for (;;) {
        u64 e = ltp_k(t, h, 0)[0];
        if (e == xh) return ltp_v(t, h, 0);
        if (!e) break;
        h = (h + 1) & t->mask;
    }
    if (LT_LFD * (u64)(t->n + 1) > LT_LFN * (u64)t->mask) { lt_grow_p(t); return lt_find_hp(w, d, H); }
    ltp_setxk(t, h, xk);
    Co *e = ltp_v(t, h, LT_KS); *e = co_zero();
    if (t->n == 0 && !w->nopush) {
        if (w->nact == w->actcap) { w->actcap *= 2; w->act = xrealloc(w->act, w->actcap * sizeof(int)); }
        w->act[w->nact++] = d;
        lvh_push(w, d);
    }
    if (t->n == t->ucap) { t->ucap *= 2; t->used = xrealloc(t->used, t->ucap * sizeof(int)); }
    t->used[t->n++] = (int)h;
    return e;
}
static inline __attribute__((always_inline)) void scatter_ks_p(WS *w, int g, int j0, int j1, int ns, const int ks) {
    const u128 *lkey = w->lkey, *skey = w->skey; const Co *lval = w->lval, *sval = w->sval;
    if (!w->lt[g] || !w->lt[g]->xk) lt_find_hp(w, g, skey[0] + lkey[j0]);   /* allocate (the key is added anyway) */
    LT *t = w->lt[g];
    const size_t RS = LT_RSK(ks), VO = LT_VOK(ks);
    unsigned char *r = (unsigned char *)t->xk; u64 mask = t->mask; int sh = t->sh;
    u64 xkh = (u64)(LT_HXK >> 64), xkl = (u64)LT_HXK;
    for (int z = 0; z < ns; z++) {
        u128 kb = skey[z]; Co fb = sval[z]; u64 kh = (u64)(kb >> 64);
        for (int j = j0; j < j1; j++) {
#if LT_PF > 0
            if (j + LT_PF < j1) __builtin_prefetch(r + ((ks ? (u64)((kb + lkey[j + LT_PF]) >> 64) : kh + (u64)(lkey[j + LT_PF] >> 64)) >> sh) * RS, 1, 1);
#endif
            u128 H; u64 hh, xl = 0;
            if (ks) { H = kb + lkey[j]; hh = (u64)(H >> 64); xl = (u64)H ^ xkl; }
            else { hh = kh + (u64)(lkey[j] >> 64); H = (u128)hh << 64; }
            u64 xh = hh ^ xkh, h = hh >> sh;
            Co *c;
            for (;;) {
                unsigned char *p = r + h * RS;
                if (ks) {
                    u64 e0 = ((u64 *)p)[0], e1 = ((u64 *)p)[1];
                    if (e0 == xh && e1 == xl) { c = (Co *)(p + VO); break; }
                    if (!(e0 | e1)) { c = lt_find_hp(w, g, H); r = (unsigned char *)t->xk; mask = t->mask; sh = t->sh; break; }
                } else {
                    u64 e = *(u64 *)p;
                    if (e == xh) { c = (Co *)(p + VO); break; }
                    if (!e) { c = lt_find_hp(w, g, H); r = (unsigned char *)t->xk; mask = t->mask; sh = t->sh; break; }
                }
                h = (h + 1) & mask;
            }
            co_fma(c, lval[j], fb);
        }
    }
}
static void scatter_group_p(WS *w, int g, int j0, int j1, int ns) { if (LT_KS) scatter_ks_p(w, g, j0, j1, ns, 1); else scatter_ks_p(w, g, j0, j1, ns, 0); }
static void ws_clear_p(WS *w) {
    for (int a = 0; a < w->nact; a++) {
        LT *t = w->lt[w->act[a]];
        if (!t) continue;
        if (t->mask + 1 > LT_KEEP) { lt_release(t); memset(t, 0, sizeof(*t)); continue; }
        for (int z = 0; z < t->n; z++) ltp_setxk(t, t->used[z], 0);
        t->n = 0;
    }
    w->nact = 0; w->nlvh = 0;
}
static void apply_curve_p(WS *w, KEY kc, int dC, Co s) {
    int T = maxdeg - dC;
    int C[64];
    unpack(kc, C);
    Co cm[64]; int nzt[64], nnz = curve_coords(C, cm, nzt);
    int nl = 0;
    for (int u = 0; u < NU && UDEG[u] <= T; u++) {
        Co v = co_zero();
        const Co *al = UALPHA + (size_t)u * NA;
        for (int q = 0; q < nnz; q++) co_fma(&v, cm[q], al[nzt[q]]);
        if (co_isz(v)) continue;
        w->lkey[nl] = ULH[u]; w->ldeg[nl] = UDEG[u]; w->lval[nl] = co_mul(v, DEGC[UDEG[u]]); nl++;
    }
    ws_groups(w, nl);
    *lt_find_hp(w, 0, 0) = co_int(1);
    i64 inner = 0, slots = 0;
    while (w->nlvh) {
        int e = lvh_pop(w);
        LT *le = w->lt[e];
        slots += le->n;
        if (w->scap < le->n) { w->scap = 2 * le->n; w->skey = xrealloc(w->skey, w->scap * sizeof(u128)); w->sval = xrealloc(w->sval, w->scap * sizeof(Co)); }
        int ns = 0;
        for (int z = 0; z < le->n; z++) {
            int xi = le->used[z];
            Co f = e == 0 ? *ltp_v(le, xi, LT_KS) : co_mul(*ltp_v(le, xi, LT_KS), INV[e]);
            if (co_isz(f)) continue;
            u128 hk0 = ltp_xk(le, xi) ^ LT_HXK; KEY ko = k_add(kc, lt_key(hk0, e));
            emit_out(w->out, ko, dC + e, co_mul(s, f));
            w->skey[ns] = hk0; w->sval[ns] = f; ns++;
        }
        for (int q = 0; q < w->ng && w->gdeg[q] <= T - e; q++) {
            int j0 = w->gst[q], j1 = w->gst[q + 1];
            if (ns) scatter_group_p(w, e + w->gdeg[q], j0, j1, ns);
            inner += (i64)ns * (j1 - j0);
        }
    }
    w->n_inner += inner;
    w->n_slots += slots;
    ws_clear_p(w);
}
static void cs_finalize_p(CS *c, int e) {
    WS *w = &c->w;
    c->ns = 0;
    if (cs_next(c) != e) return;
    lvh_pop(w);
    LT *le = w->lt[e];
    if (!le || !le->n) return;
    if (w->scap < le->n) { w->scap = 2 * le->n; w->skey = xrealloc(w->skey, w->scap * sizeof(u128)); w->sval = xrealloc(w->sval, w->scap * sizeof(Co)); }
    int ns = 0;
    for (int z = 0; z < le->n; z++) {
        int xi = le->used[z];
        Co f = e == 0 ? *ltp_v(le, xi, LT_KS) : co_mul(*ltp_v(le, xi, LT_KS), INV[e]);
        if (co_isz(f)) continue;
        u128 hk0 = ltp_xk(le, xi) ^ LT_HXK; KEY ko = k_add(c->kc, lt_key(hk0, e));
        emit_out(w->out, ko, c->dC + e, co_mul(c->s, f));
        w->skey[ns] = hk0; w->sval[ns] = f; ns++;
    }
    c->ns = ns;
    w->n_slots += le->n;
    lt_release(le); free(le); w->lt[e] = NULL;
}
static void cs_scatter_p(CS *c, int e, int q) {
    WS *w = &c->w;
    if (c->ns) scatter_group_p(w, w->gdeg[q] + e, w->gst[q], w->gst[q + 1], c->ns);
}
/* after a scatter step: register the curve's newly used levels (single-threaded per curve) */
static void cs_register(CS *c, int e) {
    WS *w = &c->w;
    for (int q = 0; q < w->ng && w->gdeg[q] <= c->T - e; q++) {
        int g = e + w->gdeg[q];
        if (!w->inh[g] && w->lt[g] && w->lt[g]->n) { w->inh[g] = 1; lvh_push(w, g); }
    }
}
typedef struct { CS *cs; int e; int *ti, *tk; } CSCtx;
static void cs_setup_fn(void *v, int i, int tid) { cs_setup(&((CSCtx *)v)->cs[i]); }
static void cs_final_fn(void *v, int i, int tid) { CSCtx *c = v; if (LT_PK) cs_finalize_p(&c->cs[i], c->e); else cs_finalize(&c->cs[i], c->e); }
static void cs_scat_fn(void *v, int i, int tid) { CSCtx *c = v; if (LT_PK) cs_scatter_p(&c->cs[c->ti[i]], c->e, c->tk[i]); else cs_scatter(&c->cs[c->ti[i]], c->e, c->tk[i]); }
static void cs_reg_fn(void *v, int i, int tid) { CSCtx *c = v; cs_register(&c->cs[i], c->e); }
static i64 *cs_sort_w;
static int cmp_task(const void *a, const void *b) { i64 x = cs_sort_w[*(const int *)a], y = cs_sort_w[*(const int *)b]; return (x < y) - (x > y); }

typedef struct { const KEY *key; const int *deg; const Co *rc; WS *ws; } LayerCtx;
/* progress lines during extraction (verbose runs): at most one per HB_EVERY seconds, from thread 0 only, so the cost
 * is one clock read per curve on one thread */
static double HB_EVERY = 30, HB_T0, HB_LAST; static int HB_D, HB_DMAX, HB_NL; static i64 HB_DONE;
static void hb_line(const char *what, int i, int n) {
    double t = now();
    if (!verbose || t - HB_LAST < HB_EVERY) return;
    HB_LAST = t;
    LOG("    ... extraction %.0fs: layer d=%d of %d (%d curves; %lld applied in earlier layers), %s %d of %d\n",
        t - HB_T0, HB_D, HB_DMAX, HB_NL, (long long)HB_DONE, what, i, n);
}
static void layer_fn(void *vctx, int i, int tid) {
    LayerCtx *ctx = vctx;
    if (tid == 0) hb_line("curve", i, HB_NL);
    if (LT_PK) apply_curve_p(&ctx->ws[tid], ctx->key[i], ctx->deg[i], ctx->rc[i]);
    else apply_curve(&ctx->ws[tid], ctx->key[i], ctx->deg[i], ctx->rc[i]);
}

#if !defined(CGV_MAJ) && !defined(CGV_WIDE)   /* the pattern hash is lt_h of 128-bit keys: with wide keys the classes are off */
/* ---- curve classes (off with CGV_MEM=low or CGV_CB=0). The curves of a layer with the same degree and the same nonzero
 * coordinates have their L terms in one union of alpha supports, so exp(C.alpha) of all of them lives on one point set
 * with one pattern of (source, L term) -> target pairs. The pattern is built once per class by hashing, then replayed per
 * curve with direct indexing (no keys, no probing), each target gathering its pairs unreduced and reducing once.
 * Residues are exact, so the different summation order changes no bit. The patterns cost memory (8 B per pair), so
 * a class is batched only while its pattern fits: all live patterns and replay buffers share the budget CB_TOTAL, set per
 * layer as a fraction of the global and level tables; builds reserve their allocated bytes as they grow. */
#ifndef CB_KMIN
#define CB_KMIN 4        /* smallest class worth a pattern */
#endif
#ifndef CB_FRAC
#define CB_FRAC 0.25    /* pattern budget, as a fraction of the global + level tables */
#endif
#ifndef CB_MININ
#define CB_MININ 256  /* smallest per-curve work (pairs) for which a pattern pays */
#endif
static size_t CB_TOTAL, CB_USED;   /* pattern budget (bytes) and bytes reserved now, shared by the threads */
static int cb_resv_add(size_t *have, size_t want) {   /* move a reservation to want bytes; 0 (unchanged) if over budget */
    if (want <= *have) { __atomic_sub_fetch(&CB_USED, *have - want, __ATOMIC_RELAXED); *have = want; return 1; }
    size_t d = want - *have;
    if (__atomic_add_fetch(&CB_USED, d, __ATOMIC_RELAXED) > CB_TOTAL) { __atomic_sub_fetch(&CB_USED, d, __ATOMIC_RELAXED); return 0; }
    *have = want; return 1;
}
static int cb_on;
static u128 CB_HXK;      /* lt_h(LT_XK): the hash stores lt_h(key) ^ CB_HXK, zero only for the marker */
static u64 CB_K1[NL], CB_K2[NL];
typedef struct {
    int nlu, *lu, *lnj;                 /* union L terms (indices into UKEY), degree-sorted; lnj[e] = # of degree <= T - e */
    u128 *lh;                           /* lt_h(UKEY[lu[j]]) */
    int np, pcap; u128 *ph; int *plev, *pnx, *pnew;   /* points: hashed key, level, next in its level list, final index */
    int *lhd, *ltl;                     /* level lists */
    u128 *hk; int *hv; u64 hmask; int hsh;
    int *toff, *ps, *pj; size_t pcap2;  /* pairs (source, L term) of each target, CSR in level order */
    size_t tcap;                        /* pairs room during the build */
    u128 *pk; int *pl;                  /* real keys and levels, final order */
    Co *F, *LV; size_t fcap, lvcap;
    unsigned char *um;
    size_t hint;                        /* expected pairs (from the class's first curve) */
    size_t resv;                        /* bytes reserved in CB_USED */
    i64 n_built, n_run, n_cls, n_crv;
} CB;
static int cb_setup(void) {
    for (int l = 0; l < NL; l++) {
        CB_K1[l] = (u64)(((u128)1 << 64) % FL[l].p); CB_K2[l] = FL[l].r2;
        if (CB_K1[l] >> 32 || CB_K2[l] >> 32) return 0;   /* the one-step fold below needs both small */
    }
    CB_HXK = lt_h(LT_XK);
    return 1;
}
static void cb_init(CB *c) {
    memset(c, 0, sizeof(*c));
    c->um = xcalloc(NU ? NU : 1, 1);
    c->lu = xmalloc((NU ? NU : 1) * sizeof(int)); c->lh = xmalloc((NU ? NU : 1) * sizeof(u128));
    c->lnj = xmalloc((maxdeg + 2) * sizeof(int)); c->lhd = xmalloc((maxdeg + 2) * sizeof(int)); c->ltl = xmalloc((maxdeg + 2) * sizeof(int));
}
static void cb_trim(CB *c, int all) {
    cb_resv_add(&c->resv, 0);
    free(c->ps); free(c->pj); c->ps = c->pj = NULL; c->pcap2 = 0;
    free(c->F); c->F = NULL; c->fcap = 0;
    free(c->hk); free(c->hv); c->hk = NULL; c->hv = NULL; c->hmask = 0;
    free(c->ph); free(c->plev); free(c->pnx); free(c->pnew); free(c->pk); free(c->pl); free(c->toff);
    c->ph = c->pk = NULL; c->plev = c->pnx = c->pnew = c->pl = c->toff = NULL; c->pcap = 0;
    free(c->LV); c->LV = NULL; c->lvcap = 0;
    if (all) { free(c->um); free(c->lu); free(c->lh); free(c->lnj); free(c->lhd); free(c->ltl); }
}
static void cb_hash_alloc(CB *c, u64 cap) {
    free(c->hk); free(c->hv);
    c->hk = xcalloc(cap, sizeof(u128)); c->hv = xmalloc(cap * sizeof(int)); c->hmask = cap - 1; c->hsh = 64 - __builtin_ctzll(cap);
}
/* bytes a build holds: targets (+ the 8 B per pair still to come), point arrays, hash */
static inline size_t cb_bbytes(size_t tcap, size_t pcap, u64 hcap) { return 12 * tcap + 24 * pcap + 20 * hcap; }
/* index of the point with hash H, added if new; -1 if growing would exceed the budget */
static inline int cb_point(CB *c, u128 H, int lev) {
    u128 xk = H ^ CB_HXK;
    u64 h = (u64)(H >> 64) >> c->hsh;
    for (;;) {
        u128 e = c->hk[h];
        if (!e) break;
        if (e == xk) return c->hv[h];
        h = (h + 1) & c->hmask;
    }
    if (4 * (u64)(c->np + 1) > 3 * c->hmask) {
        if (!cb_resv_add(&c->resv, cb_bbytes(c->tcap, c->pcap, 2 * (c->hmask + 1)))) return -1;
        cb_hash_alloc(c, 2 * (c->hmask + 1));
        for (int p = 0; p < c->np; p++) {
            u64 g = (u64)(c->ph[p] >> 64) >> c->hsh;
            while (c->hk[g]) g = (g + 1) & c->hmask;
            c->hk[g] = c->ph[p] ^ CB_HXK; c->hv[g] = p;
        }
        return cb_point(c, H, lev);
    }
    if (c->np == c->pcap) {
        int nc = c->pcap ? 2 * c->pcap : 4096;
        if (!cb_resv_add(&c->resv, cb_bbytes(c->tcap, nc, c->hmask + 1))) return -1;
        c->pcap = nc;
        c->ph = xrealloc(c->ph, c->pcap * sizeof(u128)); c->plev = xrealloc(c->plev, c->pcap * sizeof(int)); c->pnx = xrealloc(c->pnx, c->pcap * sizeof(int));
    }
    int p = c->np++;
    c->hk[h] = xk; c->hv[h] = p;
    c->ph[p] = H; c->plev[p] = lev; c->pnx[p] = -1;
    if (c->lhd[lev] < 0) c->lhd[lev] = p; else c->pnx[c->ltl[lev]] = p;
    c->ltl[lev] = p;
    return p;
}
/* mark the nonzero L terms (degree <= T) of curve kc */
static void cb_mark(CB *c, u128 kc, int T) {
    int C[64]; unpack(kc, C);
    Co cm[64]; int nzt[64], nnz = curve_coords(C, cm, nzt);
    for (int u = 0; u < NU && UDEG[u] <= T; u++) {
        if (c->um[u]) continue;
        Co v = co_zero();
        const Co *al = UALPHA + (size_t)u * NA;
        for (int q = 0; q < nnz; q++) co_fma(&v, cm[q], al[nzt[q]]);
        if (!co_isz(v)) c->um[u] = 1;
    }
}
/* pattern over the marked L terms: one lookup pass records each pair's target in source order (pj), then a
 * counting sort moves the pairs into per-target lists. 0 if it would exceed the budget (12 B per pair at peak) or
 * int indexing (more than INT_MAX pairs) */
static int cb_build(CB *c, int T) {
    c->nlu = 0;
    for (int u = 0; u < NU && UDEG[u] <= T; u++) if (c->um[u]) { c->um[u] = 0; c->lu[c->nlu] = u; c->lh[c->nlu] = lt_h(UKEY[u]); c->nlu++; }
    for (int e = 0, j = c->nlu; e <= T; e++) { while (j > 0 && UDEG[c->lu[j - 1]] > T - e) j--; c->lnj[e] = j; }
    for (int e = 0; e <= T; e++) c->lhd[e] = -1;
    c->np = 0; c->tcap = 0;
    if (!cb_resv_add(&c->resv, cb_bbytes(0, c->pcap, 4096))) return 0;
    cb_hash_alloc(c, 4096);
    if (cb_point(c, 0, 0) < 0) return 0;
    size_t x = 0, cap = 0;
    int *tg = NULL;
    for (int e = 0; e <= T; e++) {
        int nj = c->lnj[e];
        for (int p = c->lhd[e]; p >= 0; p = c->pnx[p]) {
            if (x + nj > cap) {
                if (x + nj > INT_MAX) { free(tg); return 0; }
                size_t nc = cap ? cap + cap / 2 : (c->hint > 4096 ? c->hint : 4096);
                while (nc < x + nj) nc += nc / 2;
                if (!cb_resv_add(&c->resv, cb_bbytes(nc, c->pcap, c->hmask + 1))) { free(tg); return 0; }
                cap = c->tcap = nc; tg = xrealloc(tg, cap * sizeof(int));
            }
            u128 H = c->ph[p];
            for (int j = 0; j < nj; j++) { int q = cb_point(c, H + c->lh[j], e + UDEG[c->lu[j]]); if (q < 0) { free(tg); return 0; } tg[x + j] = q; }
            x += nj;
        }
    }
    free(c->hk); free(c->hv); c->hk = NULL; c->hv = NULL; c->hmask = 0;   /* the replay needs no hash */
    int np = c->np;
    if (!cb_resv_add(&c->resv, 4 * cap + 8 * x + 52 * (size_t)c->pcap)) { free(tg); return 0; }   /* peak of the sort below */
    c->pnew = xrealloc(c->pnew, c->pcap * sizeof(int));
    c->pk = xrealloc(c->pk, c->pcap * sizeof(u128)); c->pl = xrealloc(c->pl, c->pcap * sizeof(int));
    c->toff = xrealloc(c->toff, (c->pcap + 1) * sizeof(int));
    int k = 0;
    for (int e = 0; e <= T; e++) for (int p = c->lhd[e]; p >= 0; p = c->pnx[p]) { c->pnew[p] = k; c->pk[k] = c->ph[p] * LT_HMI; c->pl[k] = e; k++; }
    memset(c->toff, 0, (np + 1) * sizeof(int));
    for (size_t z = 0; z < x; z++) { tg[z] = c->pnew[tg[z]]; c->toff[tg[z] + 1]++; }
    for (int t = 0; t < np; t++) c->toff[t + 1] += c->toff[t];
    c->pcap2 = x; c->ps = xmalloc((x ? x : 1) * sizeof(int)); c->pj = xmalloc((x ? x : 1) * sizeof(int));
    int *cur = c->pnew; memcpy(cur, c->toff, np * sizeof(int));   /* pnew is done: reuse as fill cursors */
    size_t z = 0;
    for (int p = 0; p < np; p++) { int nj = c->lnj[c->pl[p]]; for (int j = 0; j < nj; j++, z++) { int q = cur[tg[z]]++; c->ps[q] = p; c->pj[q] = j; } }
    free(tg);
    free(c->ph); free(c->plev); free(c->pnx); c->ph = NULL; c->plev = c->pnx = NULL;   /* the replay needs pk, pl only */
    cb_resv_add(&c->resv, 8 * x + 28 * (size_t)c->pcap);
    c->n_built += x;
    return 1;
}
/* replay for curve kc: target t of level e gets X = sum over its pairs of L[j] F[s], kept unreduced per lane (128 bits
 * + a carry count). X = c 2^128 + hi 2^64 + lo == c K2 + hi K1 + lo (mod p) with K1 = 2^64 mod p and K2 = 2^128 mod p
 * both < 2^32, so the folded value is < p 2^64 and one Montgomery step gives F[t] in Montgomery form */
static int cb_run(const CB *c, CB *r, WS *w, u128 kc, int dC, Co s) {   /* 0: buffers over budget, nothing done */
    int np = c->np, nlu = c->nlu;
    if (r->lvcap < (size_t)nlu || r->fcap < (size_t)np) {
        size_t lc = r->lvcap > (size_t)nlu ? r->lvcap : (size_t)nlu, fc = r->fcap > (size_t)np ? r->fcap : (size_t)np;
        if (!cb_resv_add(&r->resv, (lc + fc) * sizeof(Co))) return 0;
        if (r->lvcap < lc) { r->lvcap = lc; free(r->LV); r->LV = xmalloc(r->lvcap * sizeof(Co)); }
        if (r->fcap < fc) { r->fcap = fc; free(r->F); r->F = xmalloc(r->fcap * sizeof(Co)); }
    }
    Co *LV = r->LV, *F = r->F;
    int C[64], nzt[64]; Co cm[64];
    unpack(kc, C); int nnz = curve_coords(C, cm, nzt);
    for (int j = 0; j < nlu; j++) {
        Co v = co_zero();
        const Co *al = UALPHA + (size_t)c->lu[j] * NA;
        for (int q = 0; q < nnz; q++) co_fma(&v, cm[q], al[nzt[q]]);
        LV[j] = co_mul(v, DEGC[UDEG[c->lu[j]]]);
    }
    F[0] = co_int(1);
    emit_out(w->out, kc, dC, s);
    const int *toff = c->toff, *ps = c->ps, *pj = c->pj;
    for (int t = 1; t < np; t++) {
        int e = c->pl[t], x0 = toff[t], x1 = toff[t + 1];
        u128 acc[NL]; u64 cy[NL];
        for (int l = 0; l < NL; l++) { acc[l] = 0; cy[l] = 0; }
        for (int x = x0; x < x1;) {   /* blocks of CB_BLK pairs: F and LV are reduced (< p < 2^62, nth_prime), so products are
                                        * < 2^124 and a block sum < 16 * 2^124 = 2^128 fits; only the block needs a carry check
                                        * (on x86 the sums then stay in registers) */
            enum { CB_BLK = 16 };
            _Static_assert(CB_BLK <= 16, "block sums of products < 2^124 must fit 128 bits");
            int xe = x1 - x > CB_BLK ? x + CB_BLK : x1;
            u128 b[NL];
            for (int l = 0; l < NL; l++) b[l] = 0;
            for (; x < xe; x++) {
                const Co *fs = &F[ps[x]], *lv = &LV[pj[x]];
                for (int l = 0; l < NL; l++) b[l] += (u128)lv->v[l] * fs->v[l];
            }
            for (int l = 0; l < NL; l++) { acc[l] += b[l]; cy[l] += acc[l] < b[l]; }
        }
        Co f;
        for (int l = 0; l < NL; l++) f.v[l] = mredc(&FL[l], (u128)cy[l] * CB_K2[l] + (u128)(u64)(acc[l] >> 64) * CB_K1[l] + (u64)acc[l]);
        f = co_mul(f, INV[e]);
        F[t] = f;
        if (!co_isz(f)) emit_out(w->out, kc + c->pk[t], dC + e, co_mul(s, f));
    }
    r->n_run += toff[np] - toff[1];
    w->n_slots += np;
    return 1;
}
/* A layer is one task list, classes largest first: B(c) runs class c's first curve the hash way (its work decides whether
 * a pattern pays and fits) and builds the pattern into slot c mod CB_NS; the replays of c's other curves (one task each,
 * so a big class spreads over all threads) are listed after B(c + nthreads), so builds run ahead of the replays that
 * need them and threads rarely wait. A replay whose build is still running spins. Slot s serves classes s, s + CB_NS, ...
 * in turn: turn[s] is the class allowed in, passed on to c + CB_NS by c's failed build or by c's last replay. Tasks are
 * taken in list order and c's replays precede B(c + CB_NS), so a build waiting for its turn only waits for tasks
 * already taken. The
 * layer's single curves go first (they can be big); the fine-grained replays make the tail. */
typedef struct { u64 m; int d, i; } CBKey;
static int cmp_cbkey(const void *a, const void *b) {
    const CBKey *x = a, *y = b;
    if (x->d != y->d) return x->d - y->d;
    return (x->m > y->m) - (x->m < y->m);
}
/* cost model, learned at run time per thread (seconds per pair): h = hash way per inner op, b = pattern build, r = replay;
 * q = pattern pairs / first-curve inner ops. A class of k curves gets a pattern when (k-1) in1 h > P b + (k-1) P r,
 * P = q in1, with a 10% margin */
typedef struct { double ht, hn, bt, bn, rt, rn, qp, qi; char pad[64]; } CBCost;
static CBCost *CBC; static double CB_h = 1, CB_b = 1.5, CB_r = 0.4, CB_q = 1.3;
static void cb_cost_update(void) {
    CBCost a = {0};
    for (int t = 0; t < nthreads; t++) { a.ht += CBC[t].ht; a.hn += CBC[t].hn; a.bt += CBC[t].bt; a.bn += CBC[t].bn; a.rt += CBC[t].rt; a.rn += CBC[t].rn; a.qp += CBC[t].qp; a.qi += CBC[t].qi; }
    if (a.hn > 1e5 && a.bn > 1e5 && a.rn > 1e5) { CB_h = a.ht / a.hn; CB_b = a.bt / a.bn; CB_r = a.rt / a.rn; }
    if (a.qi > 1e5) CB_q = a.qp / a.qi;
}
typedef struct { const u128 *key; const int *deg; const Co *rc; WS *ws; CB *cb, *rb; CBKey *ck;
                 int *cst, *cn;             /* classes (largest first): start in ck, size */
                 int *st, *left, *turn;     /* per class: 0 building, 1 pattern, 2 hash way; replays left; per slot: class in turn */
                 int *tc, *tz, nt; } CBCtx; /* tasks: class (or -1 - curve for singles), member index (0 = build) */
#define CB_NS (2 * nthreads)
static void cb_build_task(CBCtx *x, int c, int tid) {
    WS *w = &x->ws[tid];
    const CBKey *k = x->ck + x->cst[c]; int n = x->cn[c], sl = c % CB_NS;
    CB *cb = &x->cb[sl];
    while (__atomic_load_n(&x->turn[sl], __ATOMIC_ACQUIRE) != c) sched_yield();
    i64 in0 = w->n_inner; double t0 = now();
    if (LT_PK) apply_curve_p(w, x->key[k[0].i], x->deg[k[0].i], x->rc[k[0].i]);
    else apply_curve(w, x->key[k[0].i], x->deg[k[0].i], x->rc[k[0].i]);
    i64 in1 = w->n_inner - in0;
    double t1 = now(); CBCost *cc = &CBC[tid];
    cc->ht += t1 - t0; cc->hn += in1;
    double P = CB_q * in1, k1 = n - 1;
    int okb = 0;
    if (in1 >= CB_MININ && __atomic_load_n(&CB_USED, __ATOMIC_RELAXED) + 16 * (size_t)P <= CB_TOTAL && 0.9 * k1 * in1 * CB_h > P * CB_b + k1 * P * CB_r) {
        int T = maxdeg - k[0].d;
        for (int z = 1; z < n; z++) cb_mark(cb, x->key[k[z].i], T);
        cb->hint = (size_t)P + (size_t)P / 4;
        okb = cb_build(cb, T);
        cc->bt += now() - t1;
        if (okb) { cc->bn += cb->toff[cb->np]; cc->qp += cb->toff[cb->np]; cc->qi += in1; cb->n_cls++; cb->n_crv += n - 1; }
        else { for (int u = 0; u < NU; u++) cb->um[u] = 0; cb_trim(cb, 0); }
    }
    if (!okb) __atomic_store_n(&x->turn[sl], c + CB_NS, __ATOMIC_RELEASE);
    __atomic_store_n(&x->st[c], okb ? 1 : 2, __ATOMIC_RELEASE);
}
static void cb_task_fn(void *vctx, int t, int tid) {
    CBCtx *x = vctx;
    int c = x->tc[t], z = x->tz[t];
    if (tid == 0) hb_line("task", t, x->nt);
    if (c < 0) { int i = -1 - c; if (LT_PK) apply_curve_p(&x->ws[tid], x->key[i], x->deg[i], x->rc[i]); else apply_curve(&x->ws[tid], x->key[i], x->deg[i], x->rc[i]); return; }
    if (z == 0) { cb_build_task(x, c, tid); return; }
    int i = x->ck[x->cst[c] + z].i, st;
    while (!(st = __atomic_load_n(&x->st[c], __ATOMIC_ACQUIRE))) sched_yield();
    if (st == 2) { if (LT_PK) apply_curve_p(&x->ws[tid], x->key[i], x->deg[i], x->rc[i]); else apply_curve(&x->ws[tid], x->key[i], x->deg[i], x->rc[i]); return; }
    CB *cb = &x->cb[c % CB_NS];
    double t0 = now();
    if (cb_run(cb, &x->rb[tid], &x->ws[tid], x->key[i], x->deg[i], x->rc[i])) { CBC[tid].rt += now() - t0; CBC[tid].rn += cb->toff[cb->np]; }
    else if (LT_PK) apply_curve_p(&x->ws[tid], x->key[i], x->deg[i], x->rc[i]); else apply_curve(&x->ws[tid], x->key[i], x->deg[i], x->rc[i]);
    if (__atomic_sub_fetch(&x->left[c], 1, __ATOMIC_ACQ_REL) == 0) { cb_trim(cb, 0); __atomic_store_n(&x->turn[c % CB_NS], c + CB_NS, __ATOMIC_RELEASE); }
}
static i64 *cb_sortw;
static int cmp_cbw(const void *a, const void *b) { i64 x = cb_sortw[*(const int *)a], y = cb_sortw[*(const int *)b]; return (x < y) - (x > y); }
static void cb_layer(const u128 *key, const int *deg, const Co *rc, int nl, WS *ws, CB *cb, CB *rb) {
    CBKey *ck = xmalloc(nl * sizeof(CBKey));
    for (int i = 0; i < nl; i++) {
        int C[64]; unpack(key[i], C);
        u64 m = 0;
        for (int t = 0; t < HC; t++) if (C[t]) m |= 1ull << t;   /* curve coordinates (HC <= 64) */
        ck[i].m = m; ck[i].d = deg[i]; ck[i].i = i;
    }
    cb_cost_update();
    qsort(ck, nl, sizeof(CBKey), cmp_cbkey);
    int *cs0 = xmalloc(nl * sizeof(int)), *cn0 = xmalloc(nl * sizeof(int)), *sg = xmalloc(nl * sizeof(int)), ncl = 0, ns = 0;
    for (int i = 0; i < nl;) {
        int j = i;
        while (j < nl && ck[j].d == ck[i].d && ck[j].m == ck[i].m) j++;
        if (j - i >= CB_KMIN) { cs0[ncl] = i; cn0[ncl++] = j - i; }
        else for (int z = i; z < j; z++) sg[ns++] = ck[z].i;
        i = j;
    }
    /* largest classes first (lower degree = deeper exponentials, then more curves) */
    int *ord = xmalloc((ncl ? ncl : 1) * sizeof(int)); i64 *wt = xmalloc((ncl ? ncl : 1) * sizeof(i64));
    for (int c = 0; c < ncl; c++) { ord[c] = c; wt[c] = (i64)cn0[c] * 64 - ck[cs0[c]].d; }
    cb_sortw = wt; qsort(ord, ncl, sizeof(int), cmp_cbw);
    int *cst = xmalloc((ncl ? ncl : 1) * sizeof(int)), *cn = xmalloc((ncl ? ncl : 1) * sizeof(int));
    int *st = xcalloc(ncl ? ncl : 1, sizeof(int)), *left = xmalloc((ncl ? ncl : 1) * sizeof(int)), *turn = xmalloc(CB_NS * sizeof(int));
    for (int c = 0; c < ncl; c++) { cst[c] = cs0[ord[c]]; cn[c] = cn0[ord[c]]; left[c] = cn[c] - 1; }
    for (int z = 0; z < CB_NS; z++) turn[z] = z;
    int *tc = xmalloc(nl * sizeof(int)), *tz = xmalloc(nl * sizeof(int)), nt = 0, R = nthreads;
    for (int z = 0; z < ns; z++) { tc[nt] = -1 - sg[z]; tz[nt++] = 0; }
    for (int c = 0; c < ncl + R; c++) {
        if (c < ncl) { tc[nt] = c; tz[nt++] = 0; }
        int r = c - R;
        if (r >= 0 && r < ncl) for (int z = 1; z < cn[r]; z++) { tc[nt] = r; tz[nt++] = z; }
    }
    CBCtx x = {key, deg, rc, ws, cb, rb, ck, cst, cn, st, left, turn, tc, tz, nt};
    parallel_for(nt, nthreads, cb_task_fn, &x);
    for (int t = 0; t < nthreads; t++) { free(rb[t].F); free(rb[t].LV); rb[t].F = rb[t].LV = NULL; rb[t].fcap = rb[t].lvcap = 0; cb_resv_add(&rb[t].resv, 0); }
    free(ck); free(cs0); free(cn0); free(sg); free(ord); free(wt); free(cst); free(cn); free(st); free(left); free(turn); free(tc); free(tz);
}
#endif

typedef struct { SM *IS; Heap *IH; OutBuf **outs; int nouts; } FoldCtx;
static void fold_fn(void *vctx, int sh, int tid) {
    FoldCtx *c = vctx;
    SM *I = &c->IS[sh]; Heap *IH = &c->IH[sh];
    for (int t = 0; t < c->nouts; t++) sm_sub_buf(I, IH, &c->outs[t][sh]);
}

typedef struct { const SPoly *a, *b; SPoly *out; SM *acc; int recip; } MulJob;
static void mul_fn(void *vctx, int i, int tid) {
    MulJob *j = &((MulJob *)vctx)[i];
    *j->out = j->recip ? sp_recip(j->a, &j->acc[tid]) : sp_mul(j->a, j->b, &j->acc[tid]);
}

/* merge per-chunk slot maps into one sorted sparse polynomial */
static SPoly merge_chunks(SM *ch, int nch) {
    SM all; sm_init(&all, 1024); int cr;
    for (int c = 0; c < nch; c++) {
        for (int s = 0; s < ch[c].n; s++) { int z = sm_add(&all, ch[c].key[s], ch[c].deg[s], &cr); all.val[z] = co_add(all.val[z], ch[c].val[s]); }
        sm_free(&ch[c]);
    }
    SPoly p = sm_take(&all);
    sm_free(&all);
    return p;
}

/* ---- chunked sparse products whose results go straight into I shards ---- */
typedef struct { const SPoly *a, *b; int xlo, xhi, neg; OutBuf *out; } PJob;
typedef struct { PJob *jobs; SM *acc; OutBuf **out; } PCtx;
/* number of terms of b with degree <= lim (b sorted by degree) */
static int count_le(const SPoly *b, int lim) {
    int lo = 0, hi = b->n;
    while (lo < hi) { int mid = (lo + hi) / 2; if (b->deg[mid] <= lim) lo = mid + 1; else hi = mid; }
    return lo;
}
static void make_prod_jobs(const SPoly **pa, const SPoly **pb, const int *pneg, int np, int target_jobs, PJob **jobs, int *nj, int *jcap) {
    double tot = 0;
    for (int p = 0; p < np; p++) {
        const SPoly *a = pa[p]->n >= pb[p]->n ? pa[p] : pb[p], *b = a == pa[p] ? pb[p] : pa[p];
        for (int x = 0; x < a->n; x++) tot += count_le(b, maxdeg - a->deg[x]);
    }
    double chunk = tot / target_jobs + 1;
    for (int p = 0; p < np; p++) {
        const SPoly *a = pa[p]->n >= pb[p]->n ? pa[p] : pb[p], *b = a == pa[p] ? pb[p] : pa[p];
        int x0 = 0; double acc = 0;
        for (int x = 0; x < a->n; x++) {
            acc += count_le(b, maxdeg - a->deg[x]);
            if (acc >= chunk || x == a->n - 1) {
                if (*nj == *jcap) { *jcap = *jcap ? 2 * *jcap : 256; *jobs = xrealloc(*jobs, *jcap * sizeof(PJob)); }
                (*jobs)[(*nj)++] = (PJob){a, b, x0, x + 1, pneg[p], NULL};
                x0 = x + 1; acc = 0;
            }
        }
    }
}
static void pjob_fn(void *vctx, int i, int tid) {
    PCtx *c = vctx; PJob *j = &c->jobs[i]; SM *A = &c->acc[tid];
    int cr;
    for (int x = j->xlo; x < j->xhi; x++) {
        int lim = maxdeg - j->a->deg[x], di = j->a->deg[x];
        KEY ki = j->a->key[x]; Co ai = j->a->val[x];
        for (int y = 0; y < j->b->n; y++) {
            if (j->b->deg[y] > lim) break;
            int s = sm_add(A, k_add(ki, j->b->key[y]), di + j->b->deg[y], &cr);
            co_fma(&A->val[s], ai, j->b->val[y]);
        }
    }
    /* subtract the chunk's result (negated when it adds) from I via the locked flush */
    for (int s = 0; s < A->n; s++) {
        Co v = A->val[s];
        if (co_isz(v)) continue;
        emit_out(c->out[tid], A->key[s], A->deg[s], j->neg ? v : co_neg(v));
    }
    sm_clear(A);
}
typedef struct { SM *IS; Heap *IH; OutBuf **outs; int nouts; } IFoldCtx;
static void ifold_fn(void *vctx, int sh, int tid) {
    IFoldCtx *c = vctx;
    for (int t = 0; t < c->nouts; t++) ex_fold(&c->outs[t][sh], sh);
}

/* Result of one pass: points with nonzero GV residues (NL lanes). */
typedef struct { int n; KEY *key; int *deg; Co *gv; } PassResult;

/* ---- parallel multicover ---- */
/* index map (host-memory step 2): key -> point index, 4 bytes per slot; the keys and values stay in the
 * MC arrays (key[i], gv[i]), so the map holds no second copy of them */
typedef struct { int *hv; u64 mask; int n; } IM;
static void im_init(IM *t, int cap) {
    u64 c = 16;
    while (c < 2 * (u64)cap) c <<= 1;
    t->hv = xmalloc(c * sizeof(int)); memset(t->hv, 0xff, c * sizeof(int)); t->mask = c - 1; t->n = 0;
}
static inline int im_get(const IM *t, const KEY *key, KEY k) {
    u64 h = hash128(k) & t->mask;
    for (;;) {
        int v = t->hv[h];
        if (v < 0) return -1;
        if (k_eq(key[v], k)) return v;
        h = (h + 1) & t->mask;
    }
}
/* insert index i (key[i] must already hold its key); the key must not be in the map yet */
static void im_put(IM *t, const KEY *key, int i) {
    if (2 * (u64)(t->n + 1) > t->mask) {
        u64 c = (t->mask + 1) * 2; int *old = t->hv; u64 om = t->mask;
        t->hv = xmalloc(c * sizeof(int)); memset(t->hv, 0xff, c * sizeof(int)); t->mask = c - 1;
        for (u64 h = 0; h <= om; h++) if (old[h] >= 0) {
            u64 g = hash128(key[old[h]]) & t->mask;
            while (t->hv[g] >= 0) g = (g + 1) & t->mask;
            t->hv[g] = old[h];
        }
        free(old);
    }
    u64 h = hash128(key[i]) & t->mask;
    while (t->hv[h] >= 0) h = (h + 1) & t->mask;
    t->hv[h] = i; t->n++;
}
typedef struct {
    KEY *key; int *deg; Co *gv; int *gcd; int *sh; int n, cap;
    IM *GS; int ng;              /* sharded key -> point index (values: gv[i] = A, then GV once its gcd level is done) */
    Co *inv3;
    const int *lst; int nl;      /* current gcd level */
    int g;
    KEY **nk; int **nd; int *nn, *ncap; /* per-thread new multiples */
} MC;
static void mc_prep(void *vctx, int ch, int tid) {
    MC *m = vctx;
    int lo = (int)((i64)m->n * ch / 256), hi = (int)((i64)m->n * (ch + 1) / 256);
    int C[64];
    for (int i = lo; i < hi; i++) {
        unpack(m->key[i], C);
        int gg = 0;
        for (int t = 0; t < HC; t++) { int x = C[t] < 0 ? -C[t] : C[t]; int a = gg, b = x; while (b) { int r = a % b; a = b; b = r; } gg = a; }
        m->gcd[i] = gg ? gg : 1;
        m->sh[i] = shard_n(m->key[i], m->ng);
    }
}
typedef struct { MC *m; int **idx; int *cnt; } MCIns;
static void mc_insert(void *vctx, int k, int tid) {
    MCIns *c = vctx; MC *m = c->m;
    for (int z = 0; z < c->cnt[k]; z++) im_put(&m->GS[k], m->key, c->idx[k][z]);   /* points are distinct */
}
static void mc_level(void *vctx, int ch, int tid) {
    MC *m = vctx;
    int lo = (int)((i64)m->nl * ch / 256), hi = (int)((i64)m->nl * (ch + 1) / 256);
    int C[64], Cn[64];
    for (int z = lo; z < hi; z++) {
        int i = m->lst[z];
        Co g = m->gv[i];
        if (m->g > 1) {
            unpack(m->key[i], C);
            for (int n = 2; n <= m->g; n++) {
                if (m->g % n) continue;
                for (int t = 0; t < HC; t++) Cn[t] = C[t] / n;
                KEY kn = pack(Cn);
                int k = im_get(&m->GS[shard_n(kn, m->ng)], m->key, kn);
                if (k < 0) continue;
                Co v = m->gv[k];   /* a divisor: lower gcd level, already final */
                if (!co_isz(v)) g = co_sub(g, co_mul(v, m->inv3[n]));
            }
        }
        m->gv[i] = g;
        if (co_isz(g)) continue;
        for (int j = 2; (i64)j * m->deg[i] <= maxdeg; j++) {
            KEY kj = k_muli(m->key[i], j);
            if (im_get(&m->GS[shard_n(kj, m->ng)], m->key, kj) >= 0) continue;
            if (m->nn[tid] == m->ncap[tid]) { m->ncap[tid] = m->ncap[tid] ? 2 * m->ncap[tid] : 256; m->nk[tid] = xrealloc(m->nk[tid], m->ncap[tid] * sizeof(KEY)); m->nd[tid] = xrealloc(m->nd[tid], m->ncap[tid] * 2 * sizeof(int)); }
            m->nk[tid][m->nn[tid]] = kj;
            m->nd[tid][2 * m->nn[tid]] = j * m->deg[i];
            m->nd[tid][2 * m->nn[tid] + 1] = j * m->g;
            m->nn[tid]++;
        }
    }
}
static PassResult multicover(KEY *rkey, int *rdeg, Co *rA, int nord) {
    MC m; memset(&m, 0, sizeof(m));
    m.key = rkey; m.deg = rdeg; m.gv = rA; m.n = nord; m.cap = nord;
    m.gcd = xmalloc((nord ? nord : 1) * sizeof(int)); m.sh = xmalloc((nord ? nord : 1) * sizeof(int));
    m.ng = 4 * nthreads;
    parallel_for(256, nthreads, mc_prep, &m);
    /* bucket by shard, insert in parallel */
    int *cnt = xcalloc(m.ng, sizeof(int)), **idx = xmalloc(m.ng * sizeof(int *));
    for (int i = 0; i < nord; i++) cnt[m.sh[i]]++;
    for (int k = 0; k < m.ng; k++) { idx[k] = xmalloc((cnt[k] ? cnt[k] : 1) * sizeof(int)); cnt[k] = 0; }
    for (int i = 0; i < nord; i++) idx[m.sh[i]][cnt[m.sh[i]]++] = i;
    m.GS = xmalloc(m.ng * sizeof(IM));
    for (int k = 0; k < m.ng; k++) im_init(&m.GS[k], cnt[k] + 16);
    MCIns ins = {&m, idx, cnt};
    parallel_for(m.ng, nthreads, mc_insert, &ins);
    for (int k = 0; k < m.ng; k++) free(idx[k]);
    free(idx); free(cnt);
    /* gcd levels */
    int maxg = maxdeg + 1;
    m.inv3 = xmalloc((maxg + 1) * sizeof(Co));
    for (int n = 1; n <= maxg; n++) m.inv3[n] = co_inv(co_int((i64)n * n * n));
    IVec *lv = xcalloc(maxg + 1, sizeof(IVec));
    for (int i = 0; i < nord; i++) if (m.deg[i] > 0) iv_push(&lv[m.gcd[i] <= maxg ? m.gcd[i] : maxg], i);
    m.nk = xcalloc(nthreads, sizeof(KEY *)); m.nd = xcalloc(nthreads, sizeof(int *));
    m.nn = xcalloc(nthreads, sizeof(int)); m.ncap = xcalloc(nthreads, sizeof(int));
    for (int g = 1; g <= maxg; g++) {
        if (!lv[g].n) continue;
        m.g = g; m.lst = lv[g].v; m.nl = lv[g].n;
        parallel_for(256, nthreads, mc_level, &m);
        /* untouched multiples found at this level: add with A = 0 */
        for (int t = 0; t < nthreads; t++) {
            for (int z = 0; z < m.nn[t]; z++) {
                KEY kj = m.nk[t][z]; int dj = m.nd[t][2 * z], gj = m.nd[t][2 * z + 1];
                int k = shard_n(kj, m.ng);
                if (im_get(&m.GS[k], m.key, kj) >= 0) continue;   /* found by another point or thread */
                if (m.n == m.cap) {
                    m.cap += m.cap / 2 + 1024;   /* 1.5x (host-memory step 3; was 2x) */
                    m.key = xrealloc(m.key, m.cap * sizeof(KEY)); m.deg = xrealloc(m.deg, m.cap * sizeof(int)); m.gv = xrealloc(m.gv, m.cap * sizeof(Co));
                    m.gcd = xrealloc(m.gcd, m.cap * sizeof(int)); m.sh = xrealloc(m.sh, m.cap * sizeof(int));
                }
                int i = m.n++;
                m.key[i] = kj; m.deg[i] = dj; m.gv[i] = co_zero(); m.gcd[i] = gj; m.sh[i] = k;
                im_put(&m.GS[k], m.key, i);
                iv_push(&lv[gj <= maxg ? gj : maxg], i);
            }
            m.nn[t] = 0;
        }
    }
    for (int k = 0; k < m.ng; k++) free(m.GS[k].hv);
    free(m.gcd); free(m.sh); m.gcd = NULL; m.sh = NULL;
    /* nonzero GVs, compacted in place (host-memory step 3; was a copy into new arrays) */
    PassResult R; R.n = 0;
    for (int i = 0; i < m.n; i++) if (m.deg[i] > 0 && !co_isz(m.gv[i])) { m.key[R.n] = m.key[i]; m.deg[R.n] = m.deg[i]; m.gv[R.n] = m.gv[i]; R.n++; }
    int nz = R.n ? R.n : 1;
    R.key = xrealloc(m.key, nz * sizeof(KEY)); R.deg = xrealloc(m.deg, nz * sizeof(int)); R.gv = xrealloc(m.gv, nz * sizeof(Co));
    m.key = NULL; m.deg = NULL; m.gv = NULL;
    for (int g = 0; g <= maxg; g++) free(lv[g].v);
    for (int t = 0; t < nthreads; t++) { free(m.nk[t]); free(m.nd[t]); }
    free(lv); free(m.GS); free(m.inv3); free(m.nk); free(m.nd); free(m.nn); free(m.ncap);
    free(m.key); free(m.deg); free(m.gv); free(m.gcd); free(m.sh);
    return R;
}

#ifdef CGV_MAJ
/* candidate GVs (from a modular run), as |GV| rounded up */
static SM CANDGV; static int have_cand;
static void read_cand(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) die("cannot open candidates");
    sm_init(&CANDGV, 1 << 16);
    int C[64], cr; char num[1100];
    for (;;) {
        int t;
        for (t = 0; t < h11; t++) if (fscanf(f, "%d", &C[t]) != 1) break;
        if (t < h11) break;
        if (fscanf(f, "%1099s", num) != 1) die("bad candidates");
        long double g = strtold(num[0] == '-' ? num + 1 : num, NULL); /* rounded up (FE_UPWARD) */
        int s = sm_add(&CANDGV, pack(C), dot(C, W), &cr);
        CANDGV.val[s].v[0] = g;
    }
    fclose(f);
    have_cand = 1;
}
/* |I[C]| = deg(C) |A_C| <= deg(C) sum_{n | C} |GV_{C/n}| / n^3 */
static Co cand_I(KEY k, int d) {
    int C[64], Cn[64], g = 0;
    unpack(k, C);
    for (int t = 0; t < h11; t++) { int a = C[t] < 0 ? -C[t] : C[t]; while (a) { int r = g % a; g = a; a = r; } }
    long double A = 0;
    for (int n = 1; n <= g; n++) {
        if (g % n) continue;
        for (int t = 0; t < h11; t++) Cn[t] = C[t] / n;
        int z = sm_get(&CANDGV, pack(Cn));
        if (z >= 0) A += CANDGV.val[z].v[0] / ((long double)n * n * n);
    }
    Co r; r.v[0] = A * d; return r;
}
#endif
static PassResult run_pass(void) {
    int H = h11;
    double t0 = now();
    tables_init();
#ifdef CGV_MAJ
    NA = ndiv + 1;
    if (NA > 64) die("majorant: too many divisors");
#else
    NA = RMODE ? NR : h11;
#endif
    /* the K-contraction for V and I: K_ab on the divisor coordinates, or M on the charge rows */
    const Co *KF = RMODE ? Mm : Kab; int NF = RMODE ? NR : H;

    /* fundamental period on the candidates */
    int nch = 64;
    FPCtx fc; fc.nchunk = nch;
    fc.c0 = xmalloc(nch * sizeof(SM)); fc.s2 = xmalloc(nch * sizeof(SM)); fc.c1 = xmalloc(NA * sizeof(SM *));
    for (int c = 0; c < nch; c++) { sm_init(&fc.c0[c], 64); sm_init(&fc.s2[c], 256); }
    for (int a = 0; a < NA; a++) { fc.c1[a] = xmalloc(nch * sizeof(SM)); for (int c = 0; c < nch; c++) sm_init(&fc.c1[a][c], 64); }
    parallel_for(nch, nthreads, fp_chunk, &fc);
    /* per-coordinate arrays on the stack up to 64 entries (as before the h11 > 64 support), else on the heap */
    SPoly c0 = merge_chunks(fc.c0, nch), S2 = merge_chunks(fc.s2, nch), c1s[64], *c1 = NA <= 64 ? c1s : xmalloc(NA * sizeof(SPoly));
    for (int a = 0; a < NA; a++) { c1[a] = merge_chunks(fc.c1[a], nch); free(fc.c1[a]); }
    free(fc.c0); free(fc.s2); free(fc.c1);
    double t1 = now();
    LOG("  fundamental period: nnz c0 %d, S2 %d  %.2fs\n", c0.n, S2.n, t1 - t0);
    MEMLOG("after fundamental period");

    SM *acc = xmalloc(nthreads * sizeof(SM));
    for (int t = 0; t < nthreads; t++) sm_init(&acc[t], 1024);
    double tp0 = now();
    SPoly c0inv;
    { MulJob j = {&c0, NULL, &c0inv, acc, 1}; mul_fn(&j, 0, 0); }
    double tp1 = now();

    /* alpha_a = c0inv * c1_a and P = c0inv * S2 */
    SPoly alphas[64], *alpha = (NA > H ? NA : H) <= 64 ? alphas : xmalloc((NA > H ? NA : H) * sizeof(SPoly)), P;
    {
        MulJob jobss[64], *jobs = NA <= 64 ? jobss : xmalloc(NA * sizeof(MulJob));
        for (int a = 0; a < NA; a++) jobs[a] = (MulJob){&c0inv, &c1[a], &alpha[a], acc, 0};
        parallel_for(NA, nthreads, mul_fn, jobs);
        if (jobs != jobss) free(jobs);
    }
    for (int a = 0; a < NA; a++) sp_free(&c1[a]);
    if (c1 != c1s) free(c1);
    sp_free(&c0);
    memset(&P, 0, sizeof(P));
    double tp2 = now();

    /* union of alpha supports */
    {
        SM U; sm_init(&U, 1024); int cr;
        for (int a = 0; a < NA; a++) for (int x = 0; x < alpha[a].n; x++) sm_add(&U, alpha[a].key[x], alpha[a].deg[x], &cr);
        int nu = U.n;
        int *ord = xmalloc((nu ? nu : 1) * sizeof(int));
        for (int s = 0; s < nu; s++) ord[s] = s;
        sort_slots(&U, ord, nu);
        int *pos = xmalloc((nu ? nu : 1) * sizeof(int));
        NU = nu;
        UKEY = xmalloc((nu ? nu : 1) * sizeof(KEY)); UDEG = xmalloc((nu ? nu : 1) * sizeof(int));
        UALPHA = xmalloc((size_t)(nu ? nu : 1) * NA * sizeof(Co));
        for (size_t z = 0; z < (size_t)(nu ? nu : 1) * NA; z++) UALPHA[z] = co_zero();
        for (int u = 0; u < nu; u++) { UKEY[u] = U.key[ord[u]]; UDEG[u] = U.deg[ord[u]]; pos[ord[u]] = u; }
        for (int a = 0; a < NA; a++) for (int x = 0; x < alpha[a].n; x++) UALPHA[(size_t)pos[sm_get(&U, alpha[a].key[x])] * NA + a] = alpha[a].val[x];
        free(ord); free(pos); sm_free(&U);
    }
    Co *UAH = UALPHA; /* alpha in the h11 coordinates, for V and I */
#ifdef CGV_MAJ
    /* |alpha_ua| <= sum_r |beta_ur| |Q_ra| */
    for (int a = 0; a < NA; a++) sp_free(&alpha[a]);
    UAH = xmalloc((size_t)(NU ? NU : 1) * H * sizeof(Co));
    for (int a = 0; a < H; a++) sp_alloc(&alpha[a], NU);
    for (int u = 0; u < NU; u++) for (int a = 0; a < H; a++) {
        Co v = co_zero();
        for (int r = 0; r < NA; r++) if (qrow(r)[a]) co_fma(&v, UALPHA[(size_t)u * NA + r], co_int(qrow(r)[a]));
        UAH[(size_t)u * H + a] = v;
        if (!co_isz(v)) { alpha[a].key[alpha[a].n] = UKEY[u]; alpha[a].deg[alpha[a].n] = UDEG[u]; alpha[a].val[alpha[a].n++] = v; }
    }
#endif

    /* V_a = 1/2 sum_b K_ab alpha_b ; I = P - sum_a alpha_a V_a */
    SPoly Vs[64], *V = NF <= 64 ? Vs : xmalloc(NF * sizeof(SPoly));
    Co half = INV[2];
    for (int a = 0; a < NF; a++) {
        sp_alloc(&V[a], NU);
        for (int u = 0; u < NU; u++) {
            Co v = co_zero();
            for (int b = 0; b < NF; b++) co_fma(&v, KF[(size_t)a * NF + b], UAH[(size_t)u * NF + b]);
            v = co_mul(v, half);
            if (!co_isz(v)) { V[a].key[V[a].n] = UKEY[u]; V[a].deg[V[a].n] = UDEG[u]; V[a].val[V[a].n++] = v; }
        }
    }
    double tp3 = now();
    /* I = P - sum_a alpha_a V_a with P = c0inv * S2: all products split into
     * chunks of similar work, each chunk's result pre-sharded by key, then the
     * shards folded in parallel */
    NSH = 4 * nthreads;
    SM *IS = xmalloc(NSH * sizeof(SM));
    Heap *IH = xcalloc(NSH, sizeof(Heap));
    for (int sh = 0; sh < NSH; sh++) sm_init(&IS[sh], 1024);
    int nI = 0;
    {
        const SPoly *pas[65], *pbs[65], **pa = NF < 65 ? pas : xmalloc((NF + 1) * sizeof(SPoly *)), **pb = NF < 65 ? pbs : xmalloc((NF + 1) * sizeof(SPoly *)); int pnegs[65], *pneg = NF < 65 ? pnegs : xmalloc((NF + 1) * sizeof(int)), np_ = 0;
        pa[np_] = &c0inv; pb[np_] = &S2; pneg[np_++] = 0;
        for (int a = 0; a < NF; a++) { pa[np_] = &alpha[a]; pb[np_] = &V[a]; pneg[np_++] = 1; }
        PJob *jobs = NULL; int nj = 0, jcap = 0;
        make_prod_jobs(pa, pb, pneg, np_, 16 * nthreads, &jobs, &nj, &jcap);
        EX_IS = IS; EX_IH = IH; EX_LK = xmalloc(NSH * sizeof(pthread_mutex_t));
        for (int sh = 0; sh < NSH; sh++) pthread_mutex_init(&EX_LK[sh], NULL);
        OutBuf **pout = xmalloc(nthreads * sizeof(OutBuf *));
        for (int t = 0; t < nthreads; t++) pout[t] = xcalloc(NSH, sizeof(OutBuf));
        PCtx pc = {jobs, acc, pout};
        parallel_for(nj, nthreads, pjob_fn, &pc);
        double tp35 = now();
        /* remaining buffered entries */
        IFoldCtx fc = {IS, IH, pout, nthreads};
        parallel_for(NSH, nthreads, ifold_fn, &fc);
        for (int t = 0; t < nthreads; t++) { for (int sh = 0; sh < NSH; sh++) ob_free(&pout[t][sh]); free(pout[t]); }
        free(pout); free(jobs); if (pa != pas) { free(pa); free(pb); free(pneg); }
        for (int sh = 0; sh < NSH; sh++) pthread_mutex_destroy(&EX_LK[sh]);
        free(EX_LK); EX_LK = NULL;
        for (int sh = 0; sh < NSH; sh++) nI += IS[sh].n;
        if (getenv("CGV_PROF")) LOG("    instanton: %d product chunks %.2fs, fold %.2fs\n", nj, tp35 - tp3, now() - tp35);
    }
    for (int a = 0; a < NF; a++) { sp_free(&V[a]); sp_free(&alpha[a]); }
    if (V != Vs) free(V);
    if (alpha != alphas) free(alpha);
    if (UAH != UALPHA) free(UAH);
#ifndef CGV_MAJ
    /* extraction pairs curves (HC coordinates) with alpha: C.alpha = sum_j c_j alpha'_j, alpha'_j = sum_t LB_jt alpha_t
     * (divisor coordinates) or sum_rho beta_rho (q_rho . LB_j) (row mode) */
    if (LB || RMODE) {
        Co *X = xmalloc((size_t)(NU ? NU : 1) * HC * sizeof(Co)), *cf = xmalloc((size_t)NA * HC * sizeof(Co));
        for (int a = 0; a < NA; a++) for (int j = 0; j < HC; j++)
            cf[(size_t)a * HC + j] = co_int(RMODE ? (RROW[a] < 0 ? Q0[j] : Q[(size_t)RROW[a] * HC + j]) : LB[(size_t)j * h11 + a]);
        for (int u = 0; u < NU; u++) for (int j = 0; j < HC; j++) {
            Co v = co_zero();
            for (int a = 0; a < NA; a++) co_fma(&v, UALPHA[(size_t)u * NA + a], cf[(size_t)a * HC + j]);
            X[(size_t)u * HC + j] = v;
        }
        free(UALPHA); free(cf); UALPHA = X; NA = HC;
    }
#endif
    sp_free(&P); sp_free(&c0inv); sp_free(&S2);
    for (int t = 0; t < nthreads; t++) sm_free(&acc[t]);
    free(acc);
    if (getenv("CGV_PROF")) LOG("    instanton: c0inv %.2fs (nnz %d), alpha %.2fs, union+V %.2fs\n", tp1 - tp0, c0inv.n, tp2 - tp1, tp3 - tp2);
    double t2 = now();
    if (getenv("CGV_PROF")) LOG("    instanton: build I %.2fs\n", t2 - tp3);
    MEMLOG("after instanton");
    LOG("  instanton polynomial: nnz I %d, |union supp alpha| %d  %.2fs\n", nI, NU, t2 - t1);

    KEY *rkey = NULL; int *rdeg = NULL; Co *rA = NULL; int nord = 0; /* A per extracted point */
    double t3;
#ifdef USE_GPU
    /* the GPU has ~0.5 s of fixed costs; small problems finish sooner on the CPU */
    int gpu_min = getenv("CGV_GPU_MIN") ? atoi(getenv("CGV_GPU_MIN")) : 200000;  /* measured: the GPU wins from ~200k instanton points (5090 and AMD, right-sized tables) */
    /* gradings with many distinct degrees (e.g. p-vectors) mean many small levels per
     * curve; the GPU pays a synchronization per level and loses to the CPU there */
    int ndeg = 0;
    if (use_gpu) {
        unsigned char *seen = xcalloc(maxdeg + 1, 1);
        for (int sh = 0; sh < NSH; sh++) for (int s = 0; s < IS[sh].n; s++) if (!seen[IS[sh].deg[s]]) { seen[IS[sh].deg[s]] = 1; ndeg++; }
        free(seen);
    }
    int gpu_maxdeg = getenv("CGV_GPU_MAXDEGS") ? atoi(getenv("CGV_GPU_MAXDEGS")) : 128;
    int gpu_ok = nI >= gpu_min && ndeg <= gpu_maxdeg && !nlight;   /* lightcones: the CPU skips the points outside */
    if (use_gpu && nlight) LOG("  (backward lightcones: extraction on the CPU)\n");
    else if (use_gpu && !gpu_ok) LOG("  (%d instanton points, %d distinct degrees: extraction on the CPU)\n", nI, ndeg);
    if (use_gpu && gpu_ok) {
        /* a busy or missing GPU is not an error: fall back to the CPU */
        double fgb = gpu_free_gb(gpu_device), need = getenv("CGV_GPU_MIN_GB") ? atof(getenv("CGV_GPU_MIN_GB")) : 4.0;
        if (fgb < need) { LOG("  (GPU %d has %.1f GB free < %.1f GB: extraction on the CPU)\n", gpu_device, fgb, need); gpu_ok = 0; }
    }
    if (use_gpu && gpu_ok) {
        GpuIn gin; memset(&gin, 0, sizeof(gin));
        gin.device = gpu_device; gin.verbose = verbose;
        for (int l = 0; l < NL; l++) { gin.p[l] = FL[l].p; gin.pinv[l] = FL[l].pinv; gin.r2[l] = FL[l].r2; }
        gin.h11 = HC; gin.keybits = keybits; gin.maxdeg = maxdeg; gin.kbw = KBV ? KBW : NULL;   /* curve coordinates; UALPHA is [NU][HC] */
        gin.NU = NU; gin.ukey = UKEY; gin.udeg = UDEG; gin.ualpha = UALPHA;
        KEY *ik = xmalloc((nI ? nI : 1) * sizeof(KEY)); int *ig = xmalloc((nI ? nI : 1) * sizeof(int)); Co *iv = xmalloc((nI ? nI : 1) * sizeof(Co));
        int k = 0;
        for (int sh = 0; sh < NSH; sh++) {
            for (int s = 0; s < IS[sh].n; s++) { ik[k] = IS[sh].key[s]; ig[k] = IS[sh].deg[s]; iv[k] = IS[sh].val[s]; k++; }
            sm_free(&IS[sh]); free(IH[sh].s);
        }
        free(IS); free(IH);
        gin.nI = k; gin.ikey = ik; gin.ideg = ig; gin.ival = iv;
        gin.tabn = tabn; gin.INV = INV; gin.DEGC = DEGC;
        gin.icap_log2 = getenv("CGV_ICAP") ? atoi(getenv("CGV_ICAP")) : 0;
        gin.xcap_log2 = getenv("CGV_XCAP") ? atoi(getenv("CGV_XCAP")) : 0;
        gin.cb = cb_allowed();
#ifdef CGV_WIDE
        /* wide keys: the GPU's exp tables key by lt_setup's level codes (phi of the offset, unhashed), decoded at emit */
        lt_setup();
        u128 *uc = xmalloc((NU ? NU : 1) * sizeof(u128));
        for (int u = 0; u < NU; u++) uc[u] = LT_KS ? ULH[u] * LT_HMI : (u128)((u64)(ULH[u] >> 64) * LT_HMI64);
        free(ULH); ULH = NULL;
        KEY dK[64], dbase, dKt0 = K_ZERO; u128 plo;
        if (LT_LB) {
            gin.dn = LBL_nf; gin.d128 = LT_KS;
            for (int f = 0; f < LBL_nf; f++) { gin.doff[f] = LBL_f[f].off; gin.db[f] = LBL_f[f].b; gin.dw[f] = LBL_f[f].w; dK[f] = LBL_k[LBL_f[f].j]; }
            dbase = LBL_k0; gin.dt0 = LBL_t0 >= 0; gin.dw0 = LBL_w0; if (gin.dt0) dKt0 = LBL_k[LBL_t0]; gin.dwlo = LBL_wlo; plo = LBL_plo;
        } else if (!LT_KS) {
            gin.dn = L64_nf; gin.d128 = 0;
            for (int f = 0; f < L64_nf; f++) { gin.doff[f] = L64_f[f].off; gin.db[f] = __builtin_popcountll(L64_f[f].M); gin.dw[f] = L64_f[f].w; dK[f] = k_shl(k_fromu(1), L64_f[f].ps); }
            dbase = L64_klo; gin.dt0 = L64_t0 >= 0; gin.dw0 = L64_w0; if (gin.dt0) dKt0 = k_shl(k_fromu(1), L64_p0); gin.dwlo = L64_wlo; plo = L64_plo;
        } else {
            gin.dn = L128_nf; gin.d128 = 1;
            for (int f = 0; f < L128_nf; f++) { gin.doff[f] = L128_f[f].off; gin.db[f] = L128_f[f].b; gin.dw[f] = L128_f[f].w; dK[f] = k_shl(k_fromu(1), L128_f[f].ps); }
            dbase = L128_klo; gin.dt0 = L128_t0 >= 0; gin.dw0 = L128_w0; if (gin.dt0) dKt0 = k_shl(k_fromu(1), L128_p0); gin.dwlo = L128_wlo; plo = L128_plo;
        }
        gin.dplo[0] = (u64)plo; gin.dplo[1] = (u64)(plo >> 64);
        gin.dK = dK; gin.dbase = &dbase; gin.dKt0 = &dKt0; gin.ukey = uc;
#endif
        GpuOut gout;
        if (gpu_extract(&gin, &gout)) die("gpu extraction failed");
        free(ik); free(ig); free(iv);
#ifdef CGV_WIDE
        free(uc);
#endif
        rkey = gout.key; rdeg = gout.deg; rA = gout.A; nord = gout.n;
        free(UKEY); free(UDEG); free(UALPHA);
        t3 = now();
        LOG("  extraction (gpu): %lld curves applied, %d points  %.2fs\n", gout.napplied, nord, t3 - t2);
    } else
#endif
    {
    /* extraction by degree layers */
    EX_IS = IS; EX_IH = IH; EX_LK = xmalloc(NSH * sizeof(pthread_mutex_t));
    for (int sh = 0; sh < NSH; sh++) pthread_mutex_init(&EX_LK[sh], NULL);
    lt_setup();
    WS *ws = xmalloc(nthreads * sizeof(WS));
    for (int t = 0; t < nthreads; t++) ws_init(&ws[t]);
    int prof = getenv("CGV_PROF") != NULL;
#if !defined(CGV_MAJ) && !defined(CGV_WIDE)
    /* curve-class patterns (see cb_build); off in the low-memory mode */
    CB *cbs = NULL;
    cb_on = cb_allowed() && cb_setup();
    CB *rbs = NULL;   /* per-thread replay buffers */
    if (cb_on) { cbs = xmalloc(2 * nthreads * sizeof(CB)); rbs = xcalloc(nthreads, sizeof(CB)); CBC = xcalloc(nthreads, sizeof(CBCost)); for (int t = 0; t < 2 * nthreads; t++) cb_init(&cbs[t]); }
#endif
    KEY *lkey = NULL; int *ldeg = NULL; Co *rc = NULL; int lcap = 0;
    int ocap = 0;
    i64 napplied = 0;
    double t_apply = 0, t_fold = 0;
    int lsync_max = getenv("CGV_LSYNC") ? atoi(getenv("CGV_LSYNC")) : 4 * nthreads;
    double cs_min_nu = getenv("CGV_CS_NU") ? atof(getenv("CGV_CS_NU")) : 20000;
    int super_delta = NU ? UDEG[0] : maxdeg + 1;
    if (super_delta < 1) super_delta = 1;
    if (getenv("CGV_NO_SUPER")) super_delta = 1;
    LOG("  layers span %d degree%s (smallest degree in the alpha supports)\n", super_delta, super_delta == 1 ? "" : "s");
    if (getenv("CGV_PROGRESS_S")) HB_EVERY = atof(getenv("CGV_PROGRESS_S"));
    HB_T0 = HB_LAST = now(); HB_DMAX = maxdeg;
    for (;;) {
        int d = -1;
        for (int sh = 0; sh < NSH; sh++) if (IH[sh].n) { int dd = IS[sh].deg[IH[sh].s[0]]; if (d < 0 || dd < d) d = dd; }
        if (d < 0) break;
        /* A curve M only changes points of degree >= deg(M) + delta, delta = the
         * smallest degree in the alpha supports (every nonzero term of exp(M.alpha)
         * has at least that degree). So all points with degree in [d, d + delta)
         * are final together, and their curves can be applied as one layer. */
        int dhi = d + super_delta; /* exclusive */
        int nl = 0;
        for (int sh = 0; sh < NSH; sh++) {
            SM *I = &IS[sh];
            while (IH[sh].n && I->deg[IH[sh].s[0]] < dhi) {
                int s = hp_pop(&IH[sh], I->deg);
                int dd = I->deg[s];
                if (nlight) {   /* points outside the lightcones feed only points outside them: skip */
                    int C[64]; unpack(I->key[s], C);
                    if (!in_lightcone(C)) continue;
                }
                if (nord == ocap) { ocap = ocap ? 2 * ocap : 4096; rkey = xrealloc(rkey, ocap * sizeof(KEY)); rdeg = xrealloc(rdeg, ocap * sizeof(int)); rA = xrealloc(rA, ocap * sizeof(Co)); }
                rkey[nord] = I->key[s]; rdeg[nord] = dd; rA[nord] = dd ? co_mul(I->val[s], INV[dd]) : co_zero(); nord++;
#ifdef CGV_MAJ
                /* certification: sources are the candidates' exact |I| values (correct by
                 * induction on degree), not their bounds, so bounds do not compound */
                if (have_cand && dd > 0 && dd < maxdeg) I->val[s] = cand_I(I->key[s], dd);
#endif
                if (dd == 0 || dd >= maxdeg || co_isz(I->val[s])) continue;
                if (nl == lcap) { lcap = lcap ? 2 * lcap : 4096; lkey = xrealloc(lkey, lcap * sizeof(KEY)); ldeg = xrealloc(ldeg, lcap * sizeof(int)); rc = xrealloc(rc, lcap * sizeof(Co)); }
                lkey[nl] = I->key[s]; ldeg[nl] = dd; rc[nl] = I->val[s]; nl++;
            }
        }
        if (!nl) continue;
        HB_D = d; HB_NL = nl; HB_DONE = napplied;
        napplied += nl;
        LayerCtx lc = {lkey, ldeg, rc, ws};
        double tl = now();
        i64 in0 = 0, s0 = 0;
        if (prof) for (int t = 0; t < nthreads; t++) { in0 += ws[t].n_inner; s0 += ws[t].n_slots; }
        OutBuf **outs; int nouts;
        CS *cs = NULL;
        /* few but big curves: advance them level-synchronously (for small curves
         * the per-level synchronization costs more than it saves) */
        double avg_nu = 0;
        if (nl < lsync_max) {
            int T = maxdeg - d, lo = 0, hi = NU;
            while (lo < hi) { int mid = (lo + hi) / 2; if (UDEG[mid] <= T) lo = mid + 1; else hi = mid; }
            avg_nu = lo;
        }
        if (nl < lsync_max && nthreads > 1 && avg_nu >= cs_min_nu) {
            cs = xcalloc(nl, sizeof(CS));
            for (int i = 0; i < nl; i++) { cs[i].kc = lkey[i]; cs[i].dC = ldeg[i]; cs[i].T = maxdeg - ldeg[i]; cs[i].s = rc[i]; }
            parallel_for(nl, nthreads, cs_setup_fn, &(CSCtx){cs, 0, NULL, NULL});
            int maxg = 1;
            for (int i = 0; i < nl; i++) if (cs[i].w.ng > maxg) maxg = cs[i].w.ng;
            CSCtx cc = {cs, 0, xmalloc((size_t)nl * maxg * sizeof(int)), xmalloc((size_t)nl * maxg * sizeof(int))};
            i64 *wt = xmalloc((size_t)nl * maxg * sizeof(i64));
            int *ord = xmalloc((size_t)nl * maxg * sizeof(int));
            for (;;) {
                /* next level used by any curve */
                int e = -1;
                for (int i = 0; i < nl; i++) { int x = cs_next(&cs[i]); if (x >= 0 && (e < 0 || x < e)) e = x; }
                if (e < 0) break;
                cc.e = e;
                hb_line("level-synchronous level", e, HB_DMAX - HB_D);
                parallel_for(nl, nthreads, cs_final_fn, &cc);
                int nt = 0;
                for (int i = 0; i < nl; i++) {
                    if (!cs[i].ns) continue;
                    for (int q = 0; q < cs[i].w.ng && cs[i].w.gdeg[q] <= cs[i].T - e; q++) {
                        int gsz = cs[i].w.gst[q + 1] - cs[i].w.gst[q];
                        cc.ti[nt] = i; cc.tk[nt] = q; wt[nt] = (i64)cs[i].ns * gsz; nt++;
                        cs[i].w.n_inner += (i64)cs[i].ns * gsz;
                    }
                }
                if (!nt) continue;
                /* biggest tasks first */
                for (int z = 0; z < nt; z++) ord[z] = z;
                cs_sort_w = wt; qsort(ord, nt, sizeof(int), cmp_task);
                int *ti2 = xmalloc(nt * sizeof(int)), *tk2 = xmalloc(nt * sizeof(int));
                for (int z = 0; z < nt; z++) { ti2[z] = cc.ti[ord[z]]; tk2[z] = cc.tk[ord[z]]; }
                memcpy(cc.ti, ti2, nt * sizeof(int)); memcpy(cc.tk, tk2, nt * sizeof(int));
                free(ti2); free(tk2);
                parallel_for(nt, nthreads, cs_scat_fn, &cc);
                parallel_for(nl, nthreads, cs_reg_fn, &cc);
            }
            free(cc.ti); free(cc.tk); free(wt); free(ord);
            outs = xmalloc(nl * sizeof(OutBuf *)); nouts = nl;
            for (int i = 0; i < nl; i++) { outs[i] = cs[i].w.out; ws[0].n_inner += cs[i].w.n_inner; ws[0].n_slots += cs[i].w.n_slots; }
        } else {
#if !defined(CGV_MAJ) && !defined(CGV_WIDE)
            if (cb_on) {
                /* pattern budget: a quarter of the global and level tables' footprint, split over the threads */
                double fp = 0;
                for (int sh = 0; sh < NSH; sh++) fp += (IS[sh].mask + 1) * 8.0 + IS[sh].cap * (double)(sizeof(u128) + sizeof(int) + sizeof(Co)) + IH[sh].cap * 4.0;
                for (int t = 0; t < nthreads; t++) for (int g = 0; g <= maxdeg; g++) if (ws[t].lt[g] && ws[t].lt[g]->xk) fp += (ws[t].lt[g]->mask + 1) * (double)(LT_PK ? LT_RSK(LT_KS) : (sizeof(u64) << LT_KS) + sizeof(Co));
                CB_TOTAL = (size_t)(fp * CB_FRAC);
                cb_layer(lkey, ldeg, rc, nl, ws, cbs, rbs);
            } else
#endif
            parallel_for(nl, nthreads, layer_fn, &lc);
            outs = xmalloc(nthreads * sizeof(OutBuf *)); nouts = nthreads;
            for (int t = 0; t < nthreads; t++) outs[t] = ws[t].out;
        }
        double tm = now();
        i64 no = 0;
        if (prof) for (int t = 0; t < nouts; t++) for (int sh = 0; sh < NSH; sh++) no += outs[t][sh].n;
        FoldCtx fcx = {IS, IH, outs, nouts};
        parallel_for(NSH, nthreads, fold_fn, &fcx);
        if (cs) { for (int i = 0; i < nl; i++) ws_free(&cs[i].w); free(cs); }
        free(outs);
        double tf = now();
        t_apply += tm - tl; t_fold += tf - tm;
        if (!prof && tf - tl >= HB_EVERY) {
            i64 pend = 0;
            for (int sh = 0; sh < NSH; sh++) pend += IH[sh].n;
            LOG("    layer d=%d of %d done: %d curves in %.0fs (apply %.0fs, fold %.0fs); extraction %.0fs, %lld points pending\n",
                d, maxdeg, nl, tf - tl, tm - tl, tf - tm, tf - HB_T0, (long long)pend);
        }
        if (prof) {
            i64 in1 = 0, s1 = 0;
            for (int t = 0; t < nthreads; t++) { in1 += ws[t].n_inner; s1 += ws[t].n_slots; }
            MEMLOG("layer");
            LOG("    layer d=%d curves %d: inner %.3g slots %.3g outputs %.3g  apply %.3fs fold %.3fs\n", d, nl,
                (double)(in1 - in0), (double)(s1 - s0), (double)no, tm - tl, tf - tm);
        }
    }
    if (getenv("CGV_PROF")) {
        double gI = 0, gW = 0, gO = 0;
        for (int sh = 0; sh < NSH; sh++) gI += (IS[sh].mask + 1) * (double)(sizeof(KEY) + sizeof(int)) + IS[sh].cap * (double)(sizeof(KEY) + sizeof(int) + sizeof(Co)) + IH[sh].cap * 4.0;
        for (int t = 0; t < nthreads; t++) {
            for (int d = 0; d <= maxdeg; d++) if (ws[t].lt[d] && ws[t].lt[d]->xk) gW += (ws[t].lt[d]->mask + 1) * (double)(LT_PK ? LT_RSK(LT_KS) : (sizeof(u64) << LT_KS) + sizeof(Co)) + ws[t].lt[d]->ucap * 4.0;
            for (int sh = 0; sh < NSH; sh++) gO += ws[t].out[sh].cap * (double)(sizeof(KEY) + sizeof(int) + sizeof(Co));
        }
        LOG("    [mem] I shards %.2f GB, thread level tables %.2f GB, thread out buffers %.2f GB, CAND %.2f GB\n", gI / 1e9, gW / 1e9, gO / 1e9,
            CAND.cap * (double)(sizeof(KEY) + sizeof(int)) / 1e9);
    }
#if !defined(CGV_MAJ) && !defined(CGV_WIDE)
    if (cb_on) {
        i64 nb = 0, nr = 0, ncl = 0, ncr = 0;
        for (int t = 0; t < 2 * nthreads; t++) { nb += cbs[t].n_built; ncl += cbs[t].n_cls; ncr += cbs[t].n_crv; cb_trim(&cbs[t], 1); }
        for (int t = 0; t < nthreads; t++) { nr += rbs[t].n_run; free(rbs[t].F); free(rbs[t].LV); }
        free(cbs); free(rbs); free(CBC);
        if (getenv("CGV_PROF")) LOG("    curve classes: cost per pair hash %.3g, build %.3g, replay %.3g ns; pattern/curve %.2f\n", CB_h * 1e9, CB_b * 1e9, CB_r * 1e9, CB_q);
        if (getenv("CGV_PROF")) LOG("    curve classes: %lld curves replayed from %lld class patterns (%.3g pairs built, %.3g replayed)\n", (long long)ncr, (long long)ncl, (double)nb, (double)nr);
    }
#endif
    for (int t = 0; t < nthreads; t++) ws_free(&ws[t]);
    for (int sh = 0; sh < NSH; sh++) pthread_mutex_destroy(&EX_LK[sh]);
    free(EX_LK); EX_LK = NULL; EX_IS = NULL; EX_IH = NULL;
    free(ws); free(lkey); free(ldeg); free(rc);
    int npts = 0;
    for (int sh = 0; sh < NSH; sh++) { npts += IS[sh].n; sm_free(&IS[sh]); free(IH[sh].s); }
    free(IS); free(IH);
    free(UKEY); free(UDEG); free(UALPHA); free(ULH);
    t3 = now();
    LOG("  extraction: %lld curves applied, %d points  %.2fs (apply %.2fs, fold %.2fs)\n", (long long)napplied, npts, t3 - t2, t_apply, t_fold);
    }

    /* GV_C = A_C - sum_{n>=2, n | C} GV_{C/n}/n^3, including multiples of
     * curves with nonzero GV that were never touched (A = 0). Points are
     * processed by gcd of their coordinates: a point only depends on its
     * divisors, which have smaller gcd, so each gcd level runs in parallel. */
    PassResult R = multicover(rkey, rdeg, rA, nord);
    if (nlight) {
        /* products of lightcone curves reach curves outside the lightcones; their values miss the curves
         * that were dropped (and are not integers), so only the lightcones are kept */
        int m = 0, C[64];
        for (int z = 0; z < R.n; z++) {
            unpack(R.key[z], C);
            if (in_lightcone(C)) { R.key[m] = R.key[z]; R.deg[m] = R.deg[z]; R.gv[m] = R.gv[z]; m++; }
        }
        R.n = m;
    }
    tables_free();
    LOG("  multicover: %d nonzero  %.2fs (pass total %.2fs)\n", R.n, now() - t3, now() - t0);
    return R;
}

/* ------------------------------------------------------------------------- */
/* CRT reconstruction                                                        */
/* ------------------------------------------------------------------------- */
typedef struct { uint32_t *d; int n; } Big;
static void big_mul_small_add(Big *x, u64 m, u64 add) {
    u128 carry = add;
    for (int i = 0; i < x->n; i++) { u128 t = (u128)x->d[i] * m + carry; x->d[i] = (uint32_t)t; carry = t >> 32; }
    while (carry) { x->d[x->n++] = (uint32_t)carry; carry >>= 32; }
}
static int big_cmp(const Big *a, const Big *b) {
    int na = a->n, nb = b->n;
    while (na && !a->d[na - 1]) na--;
    while (nb && !b->d[nb - 1]) nb--;
    if (na != nb) return na < nb ? -1 : 1;
    for (int i = na - 1; i >= 0; i--) if (a->d[i] != b->d[i]) return a->d[i] < b->d[i] ? -1 : 1;
    return 0;
}
static void big_sub(Big *a, const Big *b) {
    i64 br = 0;
    for (int i = 0; i < a->n; i++) {
        i64 t = (i64)a->d[i] - (i < b->n ? b->d[i] : 0) - br;
        br = t < 0; a->d[i] = (uint32_t)(t + (br << 32));
    }
}
static void big_to_dec(Big x, char *out) {
    uint32_t d[256]; memcpy(d, x.d, x.n * sizeof(uint32_t));
    int n = x.n; char buf[4096]; int len = 0;
    while (n && !d[n - 1]) n--;
    if (!n) { strcpy(out, "0"); return; }
    while (n) {
        u64 r = 0;
        for (int i = n - 1; i >= 0; i--) { u64 t = (r << 32) | d[i]; d[i] = (uint32_t)(t / 1000000000u); r = t % 1000000000u; }
        while (n && !d[n - 1]) n--;
        for (int k = 0; k < 9 && (n || r); k++) { buf[len++] = '0' + r % 10; r /= 10; }
    }
    for (int i = 0; i < len; i++) out[i] = buf[len - 1 - i];
    out[len] = 0;
}
/* Garner: residues (plain) -> symmetric integer as decimal. */
static void crt_symmetric(int np, const u64 *primes, const u64 *res, char *out) {
    u64 v[64];
    for (int i = 0; i < np; i++) {
        u64 x = res[i] % primes[i];
        for (int j = 0; j < i; j++) {
            u64 inv = powmod_plain(primes[j] % primes[i], primes[i] - 2, primes[i]);
            x = mulmod_plain((x + primes[i] - v[j] % primes[i]) % primes[i], inv, primes[i]);
        }
        v[i] = x;
    }
    uint32_t bx[256] = {0}, bm[256] = {0};
    Big X = {bx, 1}, M = {bm, 1};
    bm[0] = 1;
    for (int i = np - 1; i >= 0; i--) big_mul_small_add(&X, primes[i], v[i]);
    for (int i = 0; i < np; i++) big_mul_small_add(&M, primes[i], 0);
    uint32_t b2[256] = {0}; Big X2 = {b2, X.n};
    memcpy(b2, bx, X.n * sizeof(uint32_t));
    big_mul_small_add(&X2, 2, 0);
    if (big_cmp(&X2, &M) > 0) {
        uint32_t bmx[256] = {0}; Big MX = {bmx, M.n};
        memcpy(bmx, bm, M.n * sizeof(uint32_t));
        big_sub(&MX, &X);
        out[0] = '-';
        big_to_dec(MX, out + 1);
    } else {
        big_to_dec(X, out);
    }
}

/* ---- fast CRT: Garner digits with precomputed inverses; 128-bit fast path ---- */
#define CRT_CHUNKS 256
static u64 GINV[MAX_PASSES * NL][MAX_PASSES * NL]; /* GINV[i][j] = p_j^{-1} mod p_i */
static void crt_setup(int np, const u64 *primes) {
    for (int i = 0; i < np; i++) for (int j = 0; j < i; j++) GINV[i][j] = powmod_plain(primes[j] % primes[i], primes[i] - 2, primes[i]);
}
static void garner(int k, const u64 *primes, const u64 *r, u64 *v) {
    for (int i = 0; i < k; i++) {
        u64 x = r[i] % primes[i];
        for (int j = 0; j < i; j++) x = mulmod_plain((x + primes[i] - v[j] % primes[i]) % primes[i], GINV[i][j], primes[i]);
        v[i] = x;
    }
}
/* symmetric lift from the first k <= 2 primes as a signed 128-bit integer */
static i128 lift128(int k, const u64 *primes, const u64 *r) {
    u64 v[2] = {0, 0};
    garner(k, primes, r, v);
    u128 X = v[0], M = primes[0];
    if (k == 2) { X += (u128)v[1] * primes[0]; M *= primes[1]; }
    return X > M / 2 ? (i128)X - (i128)M : (i128)X;
}
static inline u64 i128_mod(i128 x, u64 p) { i128 m = x % (i128)p; return (u64)(m < 0 ? m + p : m); }
static int fmt_i128(i128 x, char *out) {
    char b[48]; int n = 0, neg = x < 0;
    u128 u = neg ? (u128)(-x) : (u128)x;
    do { b[n++] = '0' + (int)(u % 10); u /= 10; } while (u);
    int k = 0;
    if (neg) out[k++] = '-';
    while (n) out[k++] = b[--n];
    out[k] = 0;
    return k;
}
/* residues are stored by column: resv[j][s] = residue of slot s modulo prime j (host-memory step 1) */
typedef struct { u64 **resv; int n, np; const u64 *primes; int *bad; } CrtCtx;
static void crt_check_chunk(void *vctx, int ch, int tid) {
    CrtCtx *c = vctx;
    int lo = (int)((i64)c->n * ch / CRT_CHUNKS), hi = (int)((i64)c->n * (ch + 1) / CRT_CHUNKS);
    int k = c->np - 1;
    char a[1024], b[1024];
    u64 r[MAX_PASSES * NL];
    for (int s = lo; s < hi && !c->bad[ch]; s++) {
        for (int j = 0; j < c->np; j++) r[j] = c->resv[j][s];
        if (k <= 2) {
            i128 S = lift128(k, c->primes, r);
            if (i128_mod(S, c->primes[k]) != r[k] % c->primes[k]) c->bad[ch] = 1;
        } else {
            crt_symmetric(c->np, c->primes, r, a);
            crt_symmetric(k, c->primes, r, b);
            if (strcmp(a, b)) c->bad[ch] = 1;
        }
    }
}
typedef struct { u64 **resv; const KEY *key; int n; int np; const u64 *primes; char **buf; size_t *len; int *cnt; } OutCtx;
static void out_chunk(void *vctx, int ch, int tid) {
    OutCtx *c = vctx;
    int n = c->n;
    int lo = (int)((i64)n * ch / CRT_CHUNKS), hi = (int)((i64)n * (ch + 1) / CRT_CHUNKS);
    size_t cap = 4096, len = 0;
    char *buf = xmalloc(cap);
    int C[64]; char num[1024]; i64 Cf[h11];
    int k = c->np - 1;
    u64 r[MAX_PASSES * NL];
    for (int s = lo; s < hi; s++) {
        for (int j = 0; j < c->np; j++) r[j] = c->resv[j][s];
        if (k <= 2) { i128 S = lift128(k, c->primes, r); if (S == 0) continue; fmt_i128(S, num); }
        else { crt_symmetric(c->np, c->primes, r, num); if (!strcmp(num, "0")) continue; }
        while (len + 21 * (size_t)h11 + 1100 > cap) { cap *= 2; buf = xrealloc(buf, cap); }
        unpack(c->key[s], C); to_input_coords(C, Cf);
        for (int t = 0; t < h11; t++) len += fmt_i128(Cf[t], buf + len), buf[len++] = ' ';
        size_t m = strlen(num); memcpy(buf + len, num, m); len += m; buf[len++] = '\n';
        c->cnt[ch]++;
    }
    c->buf[ch] = buf; c->len[ch] = len;
}
/* one group of output chunks: task k of the group is chunk c0 + k */
typedef struct { OutCtx *oc; int c0; } OutGrp;
static void out_chunk_grp(void *vctx, int k, int tid) { OutGrp *g = vctx; out_chunk(g->oc, g->c0 + k, tid); }

#include "cgv.h"
#ifndef CGV_ENTRY
#define CGV_ENTRY cgv_entry_single
#define CGV_STANDALONE
#endif
/* Low-memory mode (CGV_MEM=low, or CGV_LOW_MEM=1): no curve-class patterns, and with glibc every allocation of at
 * least LOW_MEM_MB gets its own memory mapping, returned to the system as soon as it is freed. By default glibc raises this threshold on its own
 * after big blocks are freed, so the arrays cgv grows and frees in the first layers end up in shared pools that
 * keep the freed space (measured: ~2.5-2.9 GB held at max_deg 28 on 24 threads, CPU path). Costs some time
 * (the kernel zero-fills fresh pages); no effect with other allocators (e.g. macOS). CGV_LOW_MEM_MB sets the
 * threshold. */
#define LOW_MEM_MB 1
static double low_mem_mb;   /* threshold set (MB); 0: not in low-memory mode; -1: could not be set */
static void low_mem_setup(void) {
    low_mem_mb = 0;
    if (!mem_low()) return;
    double mb = getenv("CGV_LOW_MEM_MB") ? atof(getenv("CGV_LOW_MEM_MB")) : LOW_MEM_MB;
    if (mb < 0.125) mb = 0.125;
    if (mb > 32) mb = 32;   /* glibc's upper limit on 64-bit */
    low_mem_mb = -1;
#ifdef __GLIBC__
    if (mallopt(M_MMAP_THRESHOLD, (int)(mb * (1 << 20))) == 1) low_mem_mb = mb;
#endif
}
/* the setup notes, once the key width is settled (a run that switches to wider keys logs them from the wide entry only) */
static void setup_log(void) {
    if (low_mem_mb > 0) LOG("low-memory mode: no curve-class patterns, allocations >= %g MB mapped separately\n", low_mem_mb);
#ifdef __GLIBC__
    else if (low_mem_mb < 0) LOG("low-memory mode: no curve-class patterns\ncgv: note: mallopt failed: allocations not mapped separately\n");
#else
    else if (low_mem_mb < 0) LOG("low-memory mode: no curve-class patterns (separate mappings need glibc)\n");
#endif
    if (LB) LOG("curves span a rank-%d sublattice of Z^%d: keys in %d coordinates\n", HC, h11, HC);
#ifndef CGV_MAJ
    if (RMODE) LOG("row mode: c1 and alpha in %d charge-row coordinates (h11 = %d)\n", NR, h11);
#endif
}

int CGV_ENTRY(int argc, char **argv, CgvProbe *probe) {
    const char *in = NULL;
    verbose = 1; nthreads = 0; use_gpu = 0; gpu_device = 0; said_ks = said_pk = -1; lbc_reset(); /* settings are per call */
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-t") && i + 1 < argc) nthreads = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-q")) verbose = 0;
        else if (!strcmp(argv[i], "-g")) { use_gpu = 1; if (i + 1 < argc && argv[i + 1][0] >= '0' && argv[i + 1][0] <= '9' && strlen(argv[i + 1]) == 1) gpu_device = atoi(argv[++i]); }
        else in = argv[i];
    }
#ifdef _WIN32
    if (nthreads < 1) { long nc = (long)GetActiveProcessorCount(ALL_PROCESSOR_GROUPS); nthreads = nc > 0 ? (int)nc : 1; }   /* default: all cores */
#else
    if (nthreads < 1) { long nc = sysconf(_SC_NPROCESSORS_ONLN); nthreads = nc > 0 ? (int)nc : 1; }   /* default: all cores */
#endif
    low_mem_setup();
    FILE *f = in ? fopen(in, "r") : stdin;
    if (!f) die("cannot open input");
    double t0 = now();
    read_input(f);
    if (in) fclose(f);
    if (probe && probe->maxdeg_override > 0 && probe->maxdeg_override < maxdeg) maxdeg = probe->maxdeg_override;
    if (probe) probe->maxdeg = maxdeg;
    setup_lattice();
    setup_rows();
    { int w = setup_keys(); if (w) return CGV_RC_WIDE + w - 1; }
    setup_log();
    build_candidates();
    MEMLOG("after candidates");

#ifdef CGV_MAJ
    {
        /* one pass over the reals: per-curve upper bounds on |GV| */
        int rm = fegetround();
        fesetround(FE_UPWARD);
        if (getenv("CGV_MAJ_CAND")) read_cand(getenv("CGV_MAJ_CAND"));
        PassResult R = run_pass();
        fesetround(rm);
        double *mb = xmalloc((maxdeg + 1) * sizeof(double)), top = -1;
        for (int d = 0; d <= maxdeg; d++) mb[d] = -1;
        int C[64];
        for (int z = 0; z < R.n; z++) {
            long double v = R.gv[z].v[0];
            if (!isfinite(v)) die("majorant overflowed");
            double b = (double)log2l(v);
            if (b > mb[R.deg[z]]) mb[R.deg[z]] = b;
            if (b > top) top = b;
            if (!probe && getenv("CGV_MAJ_ALL")) {
                unpack(R.key[z], C);
                for (int t = 0; t < h11; t++) printf("%d ", C[t]);
                printf("%.3f\n", b);
            }
        }
        LOG("majorant: %d curves, max log2 bound on |GV| %.2f, total %.2fs\n", R.n, top, now() - t0);
        if (probe) { probe->maxdeg = maxdeg; probe->maxbits = mb; }
        else { for (int d = 1; d <= maxdeg; d++) if (mb[d] >= 0 && !getenv("CGV_MAJ_ALL")) printf("%d %.3f\n", d, mb[d]); free(mb); }
        free(R.key); free(R.deg); free(R.gv);
        return 0;
    }
#endif
    /* Passes of NL primes until the CRT lift is stable: the symmetric lift S
     * from all but the last prime must be congruent to the last residue (then
     * S is also the lift from all primes), for every curve. */
    u64 primes[MAX_PASSES * NL];
    int np = 0, stable = 0;
    /* distributed runs: CGV_PRIME0 = index of the first prime; CGV_RESIDUES=1 = one pass,
     * write raw residues (combined by tools/crt_combine.py) */
    int prime0 = getenv("CGV_PRIME0") ? atoi(getenv("CGV_PRIME0")) : 0, resid = getenv("CGV_RESIDUES") != NULL && !probe;
    /* slots: key rk[s], degree rd[s], residues resv[j][s] (0 for primes whose pass had GV = 0 there).
     * The first pass's result becomes the slot list as is; the key -> slot map RES is only built if a
     * later pass is needed (it can add keys that were 0 modulo the earlier primes). */
    SM RES; int res_map = 0;
    KEY *rk = NULL; int *rd = NULL; int rn = 0, rcap = 0;
    u64 *resv[MAX_PASSES * NL]; memset(resv, 0, sizeof(resv));
    int cr;
    int minp = getenv("CGV_MIN_PRIMES") ? atoi(getenv("CGV_MIN_PRIMES")) : 0; /* e.g. from certify.py */
    for (int pass = 0; pass < MAX_PASSES && (!stable || np < minp); pass++) {
        for (int l = 0; l < NL; l++) { primes[np + l] = nth_prime(prime0 + np + l); field_init(&FL[l], primes[np + l]); }
        LOG("pass %d: primes %d..%d\n", pass, np, np + NL - 1);
        PassResult R = run_pass();
        if (np == 0) {
            rk = R.key; rd = R.deg; rn = R.n; rcap = R.n ? R.n : 1;
            for (int l = 0; l < NL; l++) {
                resv[l] = xmalloc((size_t)rcap * sizeof(u64));
                for (int z = 0; z < R.n; z++) resv[l][z] = ffrom(&FL[l], R.gv[z].v[l]);
            }
            free(R.gv);
        } else {
            if (!res_map) {
                sm_init(&RES, rn + 16);
                for (int s = 0; s < rn; s++) sm_add(&RES, rk[s], rd[s], &cr);
                res_map = 1;
            }
            for (int l = 0; l < NL; l++) resv[np + l] = xcalloc(rcap, sizeof(u64));
            for (int z = 0; z < R.n; z++) {
                int s = sm_add(&RES, R.key[z], R.deg[z], &cr);
                if (cr) {
                    if (s >= rcap) {
                        int nc = rcap + rcap / 2 + 16;
                        rk = xrealloc(rk, (size_t)nc * sizeof(KEY)); rd = xrealloc(rd, (size_t)nc * sizeof(int));
                        for (int j = 0; j < np + NL; j++) {
                            resv[j] = xrealloc(resv[j], (size_t)nc * sizeof(u64));
                            memset(resv[j] + rcap, 0, (size_t)(nc - rcap) * sizeof(u64));
                        }
                        rcap = nc;
                    }
                    rk[s] = R.key[z]; rd[s] = R.deg[z]; rn = s + 1;
                }
                for (int l = 0; l < NL; l++) resv[np + l][s] = ffrom(&FL[l], R.gv[z].v[l]);
            }
            free(R.key); free(R.deg); free(R.gv);
        }
        np += NL;
        if (resid) break;
        double tc = now();
        crt_setup(np, primes);
        CrtCtx cc = {resv, rn, np, primes, xcalloc(CRT_CHUNKS, sizeof(int))};
        parallel_for(CRT_CHUNKS, nthreads, crt_check_chunk, &cc);
        stable = 1;
        for (int c = 0; c < CRT_CHUNKS; c++) stable &= !cc.bad[c];
        free(cc.bad);
        LOG("primes used: %d, stable: %d  (crt %.2fs)\n", np, stable, now() - tc);
    }
    if (resid) {
        printf("# primes");
        for (int l = 0; l < np; l++) printf(" %llu", (unsigned long long)primes[l]);
        printf("\n");
        int C[64]; i64 Cf[h11];
        for (int s = 0; s < rn; s++) {
            unpack(rk[s], C); to_input_coords(C, Cf);
            for (int t = 0; t < h11; t++) printf("%lld ", (long long)Cf[t]);
            for (int l = 0; l < np; l++) printf(l ? " %llu" : "%llu", (unsigned long long)resv[l][s]);
            printf("\n");
        }
        fflush(stdout);
        return 0;
    }
    if (!stable)
        die("CRT did not stabilize: the extracted GVs are not integers. Usually the Mori cone given is too small:\n"
            "curves outside it with nonzero GVs feed into curves inside it (e.g. it is a subcone that is not a face of\n"
            "the Mori cone). For GVs of chosen curves use their backward lightcones (the lightcone input section).");
    if (probe) {
        /* largest log2|GV| per degree (lift from all but the last prime) */
        probe->maxbits = xmalloc((maxdeg + 1) * sizeof(double));
        for (int d = 0; d <= maxdeg; d++) probe->maxbits[d] = -1;
        char num[1024];
        u64 row[MAX_PASSES * NL];
        for (int s = 0; s < rn; s++) {
            double b;
            for (int j = 0; j < np; j++) row[j] = resv[j][s];
            if (np - 1 <= 2) { i128 S = lift128(np - 1, primes, row); if (!S) continue; u128 a = S < 0 ? (u128)(-S) : (u128)S; b = 0; while (a >>= 1) b++; }
            else { crt_symmetric(np, primes, row, num); if (!strcmp(num, "0")) continue; b = (strlen(num) - (num[0] == '-')) * 3.3219; }
            int d = rd[s];
            if (b > probe->maxbits[d]) probe->maxbits[d] = b;
        }
        return 0;
    }
    double tw = now();
    OutCtx oc = {resv, rk, rn, np, primes, xcalloc(CRT_CHUNKS, sizeof(char *)), xcalloc(CRT_CHUNKS, sizeof(size_t)), xcalloc(CRT_CHUNKS, sizeof(int))};
    /* in groups of chunks, written as each group finishes, so at most one group's text is in memory
     * (host-memory step 3; was all 256 chunks) */
    int nout = 0, grp = 4 * nthreads;
    for (int c0 = 0; c0 < CRT_CHUNKS; c0 += grp) {
        int c1 = c0 + grp < CRT_CHUNKS ? c0 + grp : CRT_CHUNKS;
        OutGrp gg = {&oc, c0};
        parallel_for(c1 - c0, nthreads, out_chunk_grp, &gg);
        for (int c = c0; c < c1; c++) { fwrite(oc.buf[c], 1, oc.len[c], stdout); free(oc.buf[c]); nout += oc.cnt[c]; }
    }
    fflush(stdout);
    LOG("%d nonzero GVs, output %.2fs, total %.2fs\n", nout, now() - tw, now() - t0);
    return 0;
}

#ifdef CGV_STANDALONE
int main(int argc, char **argv) { return CGV_ENTRY(argc, argv, NULL); }
#endif
