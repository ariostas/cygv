/* One source for NVIDIA (nvcc) and AMD (hipcc): the CUDA names used in gpu.cu map to HIP here,
 * and the few device operations that differ get small portable wrappers.
 * AMD support assumes 32-lane waves (RDNA, e.g. gfx1151); checked at run time in gpu_extract. */
#ifndef GPU_COMPAT_H
#define GPU_COMPAT_H

#if defined(__HIPCC__) || defined(__HIP_PLATFORM_AMD__)
#define CGV_HIP 1
#include <hip/hip_runtime.h>
#include <hip/hip_cooperative_groups.h>
#include <hipcub/hipcub.hpp>
namespace cub = hipcub;
namespace cg = cooperative_groups;
/* runtime API */
#define cudaError_t hipError_t
#define cudaSuccess hipSuccess
#define cudaGetErrorString hipGetErrorString
#define cudaGetLastError hipGetLastError
#define cudaDeviceSynchronize hipDeviceSynchronize
#define cudaSetDevice hipSetDevice
#define cudaDeviceProp hipDeviceProp_t
#define cudaGetDeviceProperties hipGetDeviceProperties
#define cudaDeviceGetPCIBusId hipDeviceGetPCIBusId
#define cudaMalloc hipMalloc
#define cudaFree hipFree
#define cudaMemset hipMemset
#define cudaMemcpy hipMemcpy
#define cudaMemcpyAsync hipMemcpyAsync
#define cudaMemcpyHostToDevice hipMemcpyHostToDevice
#define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
#define cudaMemcpyDeviceToDevice hipMemcpyDeviceToDevice
#define cudaMemGetInfo hipMemGetInfo
#define cudaMemcpyToSymbol(sym, src, n) hipMemcpyToSymbol(HIP_SYMBOL(sym), src, n)
#define cudaLaunchCooperativeKernel hipLaunchCooperativeKernel
#define cudaOccupancyMaxActiveBlocksPerMultiprocessor hipOccupancyMaxActiveBlocksPerMultiprocessor
#define cudaFuncSetAttribute hipFuncSetAttribute
#define cudaFuncAttributeMaxDynamicSharedMemorySize hipFuncAttributeMaxDynamicSharedMemorySize

/* device: acquire/release on a 32-bit flag (agent = whole-GPU scope) */
__device__ __forceinline__ unsigned gpu_ld_acquire(const unsigned *s) { return __hip_atomic_load(s, __ATOMIC_ACQUIRE, __HIP_MEMORY_SCOPE_AGENT); }
__device__ __forceinline__ void gpu_st_release(unsigned *s, unsigned v) { __hip_atomic_store(s, v, __ATOMIC_RELEASE, __HIP_MEMORY_SCOPE_AGENT); }
/* warp operations on 32-lane waves; masks are the low 32 bits */
__device__ __forceinline__ unsigned warp_active(void) { return (unsigned)__ballot(1); }
template <class T> __device__ __forceinline__ T warp_shfl(unsigned mask, T v, int src) { (void)mask; return __shfl(v, src, 32); }
__device__ __forceinline__ unsigned warp_match_any(unsigned mask, int v) {
    /* lanes of `mask` holding the same v as this lane (emulation of __match_any_sync) */
    unsigned rem = mask, mine = 0;
    int lane = threadIdx.x & 31;
    while (rem) {
        int leader = __ffs(rem) - 1;
        int lv = __shfl(v, leader, 32);
        unsigned same = (unsigned)__ballot(v == lv) & rem;
        if (same & (1u << lane)) mine = same;
        rem &= ~same;
    }
    return mine;
}
__device__ __forceinline__ void warp_sync(void) { __builtin_amdgcn_wave_barrier(); }
__device__ __forceinline__ int warp_any(unsigned mask, int p) { (void)mask; return __any(p); }

#else /* CUDA */
#include <cuda_runtime.h>
#include <cub/device/device_scan.cuh>
#include <cub/device/device_radix_sort.cuh>
#include <cub/block/block_scan.cuh>
#include <cooperative_groups.h>
namespace cg = cooperative_groups;
__device__ __forceinline__ unsigned gpu_ld_acquire(const unsigned *s) {
    unsigned v;
    asm volatile("ld.acquire.gpu.global.u32 %0, [%1];" : "=r"(v) : "l"(s) : "memory");
    return v;
}
__device__ __forceinline__ void gpu_st_release(unsigned *s, unsigned v) {
    asm volatile("st.release.gpu.global.u32 [%0], %1;" :: "l"(s), "r"(v) : "memory");
}
__device__ __forceinline__ unsigned warp_active(void) { return __activemask(); }
template <class T> __device__ __forceinline__ T warp_shfl(unsigned mask, T v, int src) { return __shfl_sync(mask, v, src); }
__device__ __forceinline__ unsigned warp_match_any(unsigned mask, int v) { return __match_any_sync(mask, v); }
__device__ __forceinline__ void warp_sync(void) { __syncwarp(); }
__device__ __forceinline__ int warp_any(unsigned mask, int p) { return __any_sync(mask, p); }
#endif

#endif
