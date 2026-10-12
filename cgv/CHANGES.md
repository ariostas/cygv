# Changes since the 2026-09-23 snapshot

## Speed and capacity (2026-10): optimization pass
Every change below was A/B-tested against the previous state on three machines (M1 Pro, Core Ultra 7 + RTX 5090,
i5-10600K + RX 6700 XT) with identical output required (sha of the sorted GVs), then kept only if it was clearly faster in
most cases and used at most ~25% more memory.

Overall, production -> this version, same machines, identical output:
- CPU: -35 to -50% on a Core Ultra 7, -30 to -40% on an i5-10600K, -30 to -35% on an M1 Pro (measured side by side
  before the last two CPU changes, which were neutral on the M1).
- GPU extraction at max_deg 32: RTX 5090 39.7 -> 12.4 s (-69%), RX 6700 XT 176 -> 77.7 s (-56%).
- Memory: lower on the CPU; +13 to +23% GPU memory from class replay (`CGV_MEM=low` turns it off).
- High h11: faces of up to 40 dimensions at h11 = 491 to max degree 40, on CPU and GPU (e.g. a 32-dimensional face at
  degree 40: 1.95M GVs in 226 s on 8 CPU threads, ~46 s on the RTX 5090).

CPU:
- **Level tables**: structure-of-arrays slots, one 8-byte hash slot (32-bit tag + index) in the main point table
  instead of a key copy (1-thread memory -17 to -23%), level tables that store H(k) so no key is rebuilt in the
  scatter (-7.5 to -15%), and 64-bit level keys whenever they fit (decoded only for emitted slots; -7 to -18%).
- **Flush and scatter**: staged flush into the point table (hash, prefetch, then apply; -2.5 to -7%); scatter with
  the table fields held in locals (no alias reloads); packed key+value records on CPUs with an L2 below 1 MB, chosen
  at run time (i5-10600K -2 to -18%; split layout elsewhere; `CGV_LT_LAYOUT=packed|split` overrides).
- **Curve-class pattern replay** (default on): curves of the same degree and nonzero-coordinate pattern share one
  recorded scatter pattern, replayed with plain gathers; a run-time cost model skips it where it does not pay, and its
  memory comes from one shared budget. -3 to -23% (deep runs at 8 threads -16 to -23%), memory +1 to +23%.
  Off with `CGV_CB=0`, or in low-memory mode (`CGV_MEM=low`).
- **Replay inner loop** sums products in blocks of 16 with one carry fold per block (x86: -0.1 to -6%).

GPU (CUDA and HIP from the same source unless noted):
- **Level kernels** with shorter dependent-load chains (per-curve work breakpoints instead of a binary search, curve
  offsets in the entry header, 16-ary warp search): -1 to -2%. Parallel L build and fewer device syncs: -2.4%, GPU
  memory -3.6%.
- **Batch-epoch exp-table states**: the emit no longer clears the exp table (-2 to -4.5%); CUDA sm_100+ adds a
  relaxed-head fast path (-6.5%).
- **GPU curve-class replay** (default on): -52 to -60% extraction time at max_deg 32 (RTX 5090 38.7 -> 15.2 s,
  RX 6700 XT 170 -> 81 s) for +13 to +23% GPU memory (buffers scale with the card's memory). Two-pass class emit
  (probe, then insert the misses): a further -3 to -4%. Off with `CGV_CB=0` or `CGV_MEM=low`.
- **Device-driven levels** after the first big level (CUDA only; -3.5% on an RTX 5090, neutral-to-slower on AMD, so
  off there; `CGV_DLEV=0` turns it off).

Capacity (high h11):
- Curves are kept in coordinates of an LLL-reduced basis of the lattice they span, keys use per-coordinate field
  widths when needed, and level keys switch to a basis of the L-term lattice when that is smaller. Inputs on faces at
  h11 = 20-491 run (e.g. an h11 = 150 face to degree 34: 24,196 GVs in 5 s).
- **Wide keys, automatically**: when the coordinates need more than 128-bit keys, cgv reruns that entry with 256- or
  512-bit keys (CPU and GPU; input from a pipe is spooled to a temporary file so it can be reread). Nothing to set;
  the standard inputs keep 128/64-bit keys and run at the same speed. h11 = 491 faces of rank 32-40 to degree 40:
  CPU 159-226 s on 8 threads, RTX 5090 45 s.

Fixes from the review of this pass: a race on CPU pattern slots, an int overflow in the pattern build above 2^31
pairs, unchecked 64-bit overflows in coordinate conversion, gpu.cu compiling for NVIDIA GPUs older than sm_100, GPU
class-replay buffers counted against the GPU memory budget (with a fallback to the hashed path instead of an exit when
an allocation fails), wide reruns keeping the GPU (`-g`), and lattice reduction in double with `-ffp-contract=off` so
reduced bases match across architectures (x86 and Apple arm64).

Build: gv.c is now compiled for 3 lane counts x 3 key widths (and gpu.cu likewise), so a full build takes longer
(`make -j8 cgv cgv_gpu`: ~25 s) and the binaries are larger; `build.rs` (the Python wheel) compiles the same set in
parallel.

**Low-memory mode** is now `CGV_MEM=low` (`CGV_LOW_MEM=1` and `compute_gvs(..., low_memory=True)` still work as
aliases): it turns off curve-class replay on CPU and GPU, and on Linux/glibc also maps allocations of 1 MB or more
separately, as before.

## Inputs (2026-10): only what the computation needs
- **The Mori cone is given by any description**: `mori_rays` (any vectors spanning it) or `mori_hyperplanes`
  (any inequalities cutting it out), redundant either way. No lattice points of the Mori cone and no Mori cone
  Hilbert basis are needed. `tools/cgv_run.py` computes the Hilbert bases of the cones
  K_T = Mori cone ∩ {Q_r·C >= 0 for r not in T} with normaliz (cached in `~/.cache/cgv`), and cgv lists the
  candidate curves from them. The input dict key `generators` is now `mori_rays`, and `cgv_phase.compute_gv`'s
  `generators=` argument is now `mori_rays=`. `cygv_compat.compute_gv` keeps cygv's signature.
- **One candidate path**: the check that the given points cover every Mori cone lattice point up to the max degree
  (and the slow full enumeration it fell back to) is gone. It decided only agreement with cygv, never cgv's
  output, and cost up to 770 s / 79 GB on one region (now 91 s / 1 GB, identical output). `CGV_STRICT_SEMIGROUP`
  and the `cgv: note` about missing Hilbert basis elements are gone with it.
- **Input file format** (written by `cgv_run.write_input`): grading, Q, intersection numbers, the K_T
  Hilbert bases, then optional `lightcone` and `vex` sections. The original format is no longer read.
- **Lightcone GVs**: `lightcone=[p, ...]` computes only the backward lightcones of the curves p (the curves C with
  p - C in the Mori cone: exactly the curves whose GVs enter those of p), so the GVs are exact; `max_deg` then
  defaults to the largest degree of the p. It replaces `generators=pts, saturate=False` (curves that are
  nonnegative integer sums of pts), which dropped curves feeding in from outside and stopped with "CRT did not
  stabilize" unless the cone of pts was a face of the Mori cone. That error now explains its usual cause.
- **Removed**: `cgv_run.run_cygv` and the `--compare` option (`tests/make_refs.py` calls cygv itself).
- **Fewer K_T cones; faces through the plain interface** (`cgv_run.cone_data`): only divisors negative somewhere on
  the cone are put in T, and each distinct cone is written once (K_T = K_M, M the divisors negative on its Hilbert
  basis). Identical output; e.g. P(1,1,2,2,2)[8] 21 -> 1 cones, the h11 = 10 tests 91-105 -> 55-91. `mori_rays` may
  span a face of the Mori cone: normaliz then runs in a basis of the face's lattice (a 12-dimensional face at
  h11 = 491: 123 normaliz runs, 1.7 s; before: ~123k runs in Z^491, a sampled one unfinished after 10 min, and the face's equations were dropped).
  `grading_vector` defaults to CYTools' `find_grading_vector` and is checked to be positive on the cone; `intnums`
  may be a dict. The cache files are now `kt2_*.json`.
- **No format tag**: input files start with `h11 ndiv maxdeg`; a leading `cgv 2` (older files) is still accepted.

## Behaviour
- **Vex phases** (new, additive): cgv accepts an optional strata section (one exact curve class per vex 3-cone)
  and computes curves of negative anticanonical degree with the pole-free coefficient; before, it stopped with
  "negative anticanonical degree". A vex 2-cone in the input stops with a pointer to MacFadden–Sheridan Prop. 5.
  Inputs without the section behave exactly as before (identical output on all regression inputs, CPU and GPU).
  `tools/cgv_phase.py` builds the input from a fan, `q` and `kappa`; `tests/vex_regress.py` tests it.
- **Incomplete generator sets**: superseded by the input changes above (2026-10).
- **GPU hand-off threshold** 300k -> 200k instanton points (measured on NVIDIA and AMD): mid-size jobs
  now use the GPU (e.g. a D22 job 3.8 s on CPU -> 1.2 s on an RTX 5090).

- **`cgv` without `-t` uses all cores** (was 1 thread; `compute_gvs` already defaulted to all cores). Fewer
  threads are the simplest way to use less memory (README).

## GPU bug fixes (all present in the 2026-09-23 snapshot)
- **Intermittent segfault**: the batch-size estimator read past the end of an array after each batch
  (crash in ~1 of 5 runs on some inputs; otherwise a NaN estimate and poor batching). Fixing it also made
  GPU runs ~10% faster.
- **Silent kernel launch failures**: launch errors were not checked; a failed launch skipped work that
  only the CRT check caught. Now checked at every sync point.
- **Out of memory on deep runs**: I-table growth now stages through host memory, releases the L pool
  and shrinks the exp table when needed (e.g. a D34 run that failed now completes).
- **Layers with more than 4.2M curves** (error `overflow 2`): curve buffers now grow, sized by a
  counting pass before extraction.

## Performance
- **Right-sized GPU tables**: the exp table starts at 1 GB and grows only when one curve needs more
  (also fixes "a single curve overflows the exp table"); the L pool starts at <= 2 GB and grows; level
  lists have their own allocation. Integrated GPUs: small jobs 2.2 s -> 0.7 s; discrete GPUs equal or
  slightly faster (D32: 219 s -> 213 s on an RTX 5090).

- **Faster hash lookup in the GPU scatter** (v2; found with ShinkaEvolve): the common case, an
  occupied slot, is compared and skipped without re-reading the entry through volatile loads.
  Identical output; RTX 5090: +3-5% at max_deg 24-28, +7% at 32 (e.g. 213 s -> 199 s); AMD: +6% at
  D24, +9% at D28. Gains grow with degree.

- **GPU batch loop** (found by a Claude Code kernel evolution, then tuned by hand): fewer blocking
  host-device round trips per batch, one fused emit+clear pass, early abort of batches that overflow
  their table, batch tables may fill to 3/4 before a retry, and the emit of batch k overlapped with
  batch k+1's levels on a second stream (second table capped at 2^22 slots / 256 MB). Identical
  output; RTX 5090 7-14% faster at max_deg 28-32, 19% at 36 (18.7 -> 15.2 min).

- **CPU extraction ~12-13% faster** (the CPU analogue of the GPU table lessons): the per-degree level tables grow at
  3/4 full instead of 1/2 (smaller tables stay in cache), their entries shrink from 48 to 40 bytes (the key is
  stored XOR a marker key of negative degree, so an empty slot is all zeros and no occupied flag is needed), and the
  scatter loop prefetches 4 terms ahead. Identical output. 2c076428 (h11 = 10), i5-10600K, 10 threads:
  max_deg 24 35.5 -> 30.9 s, 26 136.7 -> 119.0 s, 30 1557 -> 1366 s.

## Memory
- **Deep runs no longer fail at the very end**: all GPU extraction buffers are freed before the
  results are collected (a max_deg 36 run of an h11 = 10 geometry that needed ~3 GB more now
  completes on a 32 GB GPU).
- **I table (the running answer) kept 70-85% full instead of 40-50%**: same speed, one fewer table
  doubling on deep runs (max_deg 36: peak GPU memory 29.6 -> 21.1 GB; max_deg 32: 13.7-16.8 ->
  10.6 GB). `CGV_I_TGT` / `CGV_I_CAPF` override.
- **Host memory at the end of a run about halved** (CPU and GPU paths; the peak is in the final stages:
  multicover, CRT, output). Residues are stored per prime in flat arrays (was a 256-byte row per GV
  plus a second map of all keys); the multicover maps hold 4-byte point indices instead of a second
  copy of every key and value; the results are compacted in place; the output is written in groups of
  chunks. Identical output. Peak host RSS (h11 = 10 geometry, GPU path): max_deg 30 3.7 -> 2.1 GB,
  32 7.0 -> 3.6 GB, 36 18.2 -> 9.1 GB; run time unchanged.
- **CPU peak memory ~20% lower** as a side effect of the fuller, smaller level tables: max_deg 26 2.1 -> 1.7 GB,
  max_deg 30 7.0 -> 5.6 GB (same runs as above).
- **Low-memory switch** `CGV_LOW_MEM=1` / `compute_gvs(..., low_memory=True)` (Linux/glibc): allocations
  of 1 MB or more get their own mapping and are returned to the system when freed, instead of staying in
  glibc's pools as fragmentation. CPU path, max_deg 26-28, 24 threads: peak host memory -22% to -30% for
  about +10% time. Off by default; no effect on macOS.

## Portability
- **AMD GPUs** via HIP from the same `gpu.cu` (`gpu_compat.h`, `make cgv_hip HIP_ARCH=...`). Tested on
  an RX 6700 XT (RDNA2, ROCm 5.7 with `HIP_ARCH=gfx1030 HSA_OVERRIDE_GFX_VERSION=10.3.0`) and a
  Radeon 8060S (RDNA3.5, ROCm 7.1). Integrated GPUs size their tables from a 24 GB budget
  (`CGV_GPU_MEM_GB` sets it on any GPU). Needs 32-lane waves.
- **macOS / Apple Silicon**: the CPU path builds and passes the stress tests (M1 Pro). There is no
  Metal GPU backend: a prototype was exact on the regression set but 3.8x slower than the CPU path on an
  M1 Pro (hash-table memory latency), so it is not included.

## New tools
- `tools/cygv_compat.py`: `compute_gv` / `compute_gw` with cygv's exact signatures; options cgv does
  not support (`nefpart`, `min_points`, `target_points`, Gromov–Witten) raise `NotImplementedError`.

## Validation of this snapshot
CPU regression 38/38; GPU regression in every mode (default, forced, tiny tables, memory pressure,
2 GB budget, tiny curve buffers); stress tests clean on NVIDIA and AMD; CPU stress clean on macOS.
Exactness certified (tools/certify.py) for several PFVs up to max_deg 28-34.
Update of 2026-09-28 (CPU extraction speedup, thread default): CPU regression 38/38 and CPU+HIP stress clean
(i5-10600K + RX 6700 XT, HIP_ARCH=gfx1030); GPU regression forced onto the GPU exact on every case that fit in the
memory left by another job (the 4 largest p-like cases were out of memory there); sorted output identical to the
previous version at max_deg 24, 26 and 30 (2c076428).
