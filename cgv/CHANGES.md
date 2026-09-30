# Changes since the 2026-09-23 snapshot

## Behaviour
- **Vex phases** (new, additive): cgv accepts an optional strata section (one exact curve class per vex 3-cone)
  and computes curves of negative anticanonical degree with the pole-free coefficient; before, it stopped with
  "negative anticanonical degree". A vex 2-cone in the input stops with a pointer to MacFadden–Sheridan Prop. 5.
  Inputs without the section behave exactly as before (identical output on all regression inputs, CPU and GPU).
  `tools/cgv_phase.py` builds the input from a fan, `q` and `kappa`; `tests/vex_regress.py` tests it.
- **Incomplete generator sets** (reported by the group): when the generators miss a Hilbert basis
  element of their cone, cgv now prints `cgv: note: ...` and stays on the fast path instead of
  falling back to slow full enumeration (178 s / 20 GB before). `compute_gvs` raises it as a Python
  warning. Semantics: cgv computes over all lattice points of the cone (the saturated semigroup);
  cygv uses the generated semigroup. Identical GVs in every case compared. `CGV_STRICT_SEMIGROUP=1`
  restores the old behaviour.
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
