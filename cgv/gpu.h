/* Interface between gv.c (host pipeline) and gpu.cu (extraction on a CUDA GPU). */
#ifndef CGV_GPU_H
#define CGV_GPU_H
#include <stdint.h>
typedef struct {
    int device, verbose;
    uint64_t p[8], pinv[8], r2[8];        /* per lane (first NL used) */
    int h11, keybits, maxdeg;
    const int *kbw;                       /* per-coordinate key field widths [h11], or NULL: keybits each */
    /* CGV_WIDE builds: ukey holds the L terms' level codes (u128), ikey and the output keys are 256/512-bit; a code o of
     * level e decodes to dbase + sum_f v_f dK[f] (+ r dKt0), v_f = bits [doff, doff + db) of o - dplo (mod 2^64, or
     * 2^128 if d128), r = (e - dwlo - sum_f dw[f] v_f) / dw0 (if dt0) */
    int d128, dn, dt0, dw0; long long dwlo; uint64_t dplo[2]; int doff[64], db[64], dw[64]; const void *dK, *dbase, *dKt0;
    int NU; const void *ukey; const int *udeg; const void *ualpha;   /* union of alpha supports */
    int nI; const void *ikey; const int *ideg; const void *ival;     /* initial instanton polynomial */
    int tabn; const void *INV; const void *DEGC;                     /* 1/n and n, Montgomery */
    int icap_log2, xcap_log2;             /* optional table size overrides (0 = auto) */
    int cb;                               /* class replay allowed (0: CGV_CB=0 or low-memory mode) */
} GpuIn;
typedef struct { int n; void *key; int *deg; void *A; long long napplied; } GpuOut;
#ifdef __cplusplus
extern "C"
#endif
int gpu_extract(const GpuIn *in, GpuOut *out);
/* free memory on a CUDA device in GB (-1 if the device is unusable) */
#ifdef __cplusplus
extern "C"
#endif
double gpu_free_gb(int device);
#endif
