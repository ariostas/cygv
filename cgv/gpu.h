/* Interface between gv.c (host pipeline) and gpu.cu (extraction on a CUDA GPU). */
#ifndef CGV_GPU_H
#define CGV_GPU_H
#include <stdint.h>
typedef struct {
    int device, verbose;
    uint64_t p[8], pinv[8], r2[8];        /* per lane (first NL used) */
    int h11, keybits, maxdeg;
    int NU; const void *ukey; const int *udeg; const void *ualpha;   /* union of alpha supports */
    int nI; const void *ikey; const int *ideg; const void *ival;     /* initial instanton polynomial */
    int tabn; const void *INV; const void *DEGC;                     /* 1/n and n, Montgomery */
    int icap_log2, xcap_log2;             /* optional table size overrides (0 = auto) */
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
