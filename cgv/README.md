# cgv

Gopakumar–Vafa invariants of Calabi–Yau threefold hypersurfaces (HKTY), in C with an optional GPU
extraction stage (NVIDIA CUDA and AMD HIP from one source, `gpu_compat.h`). Results are exact and identical to cygv 0.2.3 (regression-tested), at roughly
1000x the speed and a fraction of the memory on deep degrees.

## Use

From Python (cytools):

```python
import sys; sys.path.insert(0, "path/to/cgv/tools")
from cgv_run import compute_gvs
gvs = compute_gvs(cy, max_deg=24)                       # {curve: GV}, like cy.compute_gvs(...).dok
gvs = compute_gvs(cy, max_deg=24, grading_vec=p)        # custom grading vector
gvs = compute_gvs(cy, max_deg=24, device="cpu")         # or "gpu", "gpu:1", "auto" (default)
```

From the shell (input written by `tools/cgv_run.py`'s `write_input`):

```
./cgv     input.txt                # CPU, all cores (-t N: N threads)
./cgv_gpu -t 24 -g 0 input.txt     # CUDA GPU 0 (small, sparse-degree, or GPU-busy cases stay on the CPU)
          -l N                     # force N prime lanes (default: chosen by a probe run)
          -q                       # quiet
```

Output: one line per nonzero GV, `c_1 ... c_h11 gv`.

Build: `make cgv cgv_gpu` (the GPU build needs nvcc; `CUDA_ARCH` defaults to sm_120);
`make cgv_hip HIP_ARCH=gfx1100` for AMD (hipcc, 32-lane waves). Tests: `python tests/regress.py`
(38 frozen cygv 0.2.3 outputs). Changes since the 2026-09-23 snapshot: `CHANGES.md`.

This directory is part of [cygv](https://github.com/ariostas/cygv) and is
licensed like cygv (GPL-3.0-or-later, see `../LICENSE`). It was developed with Claude (Anthropic)
through Claude Code, including two rounds of agent-driven kernel evolution on the GPU code.
AMD: `make cgv_hip HIP_ARCH=gfx1151` (hipcc; same source through `gpu_compat.h`; needs 32-lane waves, i.e. RDNA).
GPUs that ROCm does not officially support: build for the nearest supported target and set
`HSA_OVERRIDE_GFX_VERSION` when running (RX 6700 XT, gfx1031: `HIP_ARCH=gfx1030`, `HSA_OVERRIDE_GFX_VERSION=10.3.0`).
Integrated GPUs (shared memory) size their tables from a 24 GB budget; `CGV_GPU_MEM_GB` sets it on any GPU.

Low on memory? Two levers, which stack:
- **Fewer threads** (`-t N`, or `compute_gvs(..., threads=N)`; default: all cores). Each thread keeps its own
  working tables, so peak host memory drops with the thread count (on top of a fixed part), at the cost of time.
  CPU path, 24 -> 6 threads: max_deg 28 (h11 = 10) 6.1 -> 3.0 GB, 118 -> 192 s; max_deg 26 (h11 = 11) 5.5 -> 2.8 GB,
  62 -> 95 s. Works on every OS (MacBook Pro, 6 threads: the same cases in 3.2 GB).
- **`CGV_LOW_MEM=1`** (or `compute_gvs(..., low_memory=True)`), Linux/glibc only: every allocation of 1 MB or more
  gets its own mapping, returned to the system as soon as it is freed (`CGV_LOW_MEM_MB` sets the size).
  CPU path, max_deg 26-28, 24 threads: peak host memory -22% to -30% for about +10% time. No effect on macOS.

## Any phase: FRSTs and vex fans

`tools/cgv_phase.py` runs cgv in a given phase of the ambient toric variety, including *vex* phases: fine regular
fans on the polytope's boundary points (not interior to facets) that do not refine its face fan (MacFadden–Sheridan,
[arXiv:2512.14817](https://arxiv.org/abs/2512.14817)). There, curves of negative anticanonical degree exist and
HKTY's coefficient has a pole; it is removed exactly by the class of the stratum those curves force into the CY
(Liam McAllister's prescription, following J. Wang's quasimap I-function,
[arXiv:1910.14440](https://arxiv.org/abs/1910.14440)). For fine fans the stratum is always a toric curve
(MacFadden–Sheridan Prop. 5), so the extra input is one curve class per vex 3-cone, which the module derives itself.

```python
import sys; sys.path.insert(0, "path/to/cgv/tools")
from cgv_phase import compute_gv, phase_from_polytope
gvs = compute_gv(cones, q, kappa, max_deg=20)      # cones: maximal cones as column indices of q
cones, q, kappa = phase_from_polytope(p, fan)      # optional: from a CYTools polytope and its fan (point labels)
gvs = compute_gv(cones, q, kappa, 20, generators=subset, saturate=False)   # lightcone GVs: a subset of generators
```

Inputs are cygv's (`q`, `kappa`, optionally `generators` and `grading`) plus the fan. Everything else is derived with
numpy and normaliz: the rays (integer Gale dual of `q`), the fan's intersection ring (`kappa` is checked against it),
the Mori cone (wall curves) and its Hilbert basis, and the vex cones with their exact stratum classes. With no vex
cones the run is plain cgv. Validation: the FRST and vex phases of Liam's two examples agree with his exact-rational
reference and across every flop; and vex phases of one polytope reproduce FRSTs of *another* polytope with the same CY
(identified by a basis change of kappa and c2) class by class, e.g. about 3 million classes over 247 such pairs at
h11 = 3-5 to degree 50, with no disagreement.

## Exactness certificate

The default run checks its CRT lift with one extra prime, which is strong evidence but not a
proof. To prove an output exact:

```
make cgv_maj
./cgv -t 24 in.txt > gvs.txt              # stderr: "primes used: N"
python tools/certify.py in.txt gvs.txt N  # exit 0 = CERTIFIED
```

`cgv_maj` reruns the pipeline over nonnegative reals with upward rounding (subtraction becomes
addition), feeding in the candidates' |GV| for lower-degree curves. That gives a proven bound
B_C >= |GV_C| assuming everything of lower degree is right. If every B_C < M/2 (M = product of
the primes), induction on degree proves every GV exact. The cost is about one CPU run with
more memory. If it fails, it reports how many primes are needed; rerun with `CGV_MIN_PRIMES`
or a distributed run (below).

## Many jobs, many machines

```
python tools/dispatch.py slots.json --deploy           # copy sources to ssh hosts, build
python tools/dispatch.py slots.json OUT job*.txt       # queue: GPU slots take the largest jobs
```

One run split by primes across devices (e.g. primes 0-1 on the GPU, 2-3 on another host):

```
CGV_RESIDUES=1 CGV_PRIME0=0 ./cgv_gpu -l 2 -g 0 in.txt > a.res
CGV_RESIDUES=1 CGV_PRIME0=2 ssh host cgv -l 2 in.txt   > b.res
python tools/crt_combine.py a.res b.res > gvs.txt      # checks the lift is stable
python tools/certify.py in.txt gvs.txt 4
```

## How it works (same mathematics as cygv, reorganized)

1. Arithmetic mod ~62-bit primes (Montgomery), 2–4 primes carried as lanes; integer GVs by CRT,
   with one extra prime as a check (more passes run automatically if the check fails).
2. Only curves with <= 2 negative GLSM intersections carry fundamental-period data; they are the
   lattice points of a few cones (Hilbert bases by normaliz, cached in `~/.cache/cgv`), so the
   Mori cone is never enumerated (checked saturated; otherwise falls back to full enumeration).
3. One instanton polynomial I = sum_t w_t inst_t (w = grading vector) instead of h11 of them.
   After subtracting lower curves, I[C] = deg(C) A_C, A_C = sum_{n|C} GV_{C/n}/n^3.
4. Each curve's correction is s_C z^C exp(C.alpha), computed with the Euler recurrence per degree
   level, scattered one target degree at a time (CPU) or as flat (entry, term) work lists (GPU).

## Tests

```
python tests/regress.py [--bin ./cgv_gpu --extra "-g 0"] [--no-cones]   # vs frozen cygv outputs
python tests/vex_regress.py [--bin ./cgv_gpu --extra "-g 0"]            # phases incl. vex (tools/cgv_phase.py)
python tests/stress.py 300                                               # races: modes x thread counts
CGV_GPU_MIN=0 CGV_XCAP=19 CGV_ICAP=10 python tests/regress.py --bin ./cgv_gpu --extra "-g 0"   # GPU overflow/retry paths
CGV_ICAP=10 CGV_TEST_LEAVE_MB=300 ./cgv_gpu -g 0 in.txt   # GPU memory pressure: staged I growth, L release, exp shrink
```

References in `tests/refs_vex` (16, regenerated by `tests/make_vex_refs.py`, which needs CYTools): every fan of
Liam's h11 = 2 and 3 examples, and four vex/FRST pairs of different polytopes with the same CY (h11 = 4, 5), checked
against each other through their basis map.

References in `tests/refs` (35): quintic, ten h11=10 CYs (default grading; four also with p-like
gradings), and the h11=9 case pfv 2d7b127a at max_deg 10..20 (deg 20 from a 67-minute cygv run).

## Tuning knobs (environment)

`CGV_PROF=1` phase/layer profile; `CGV_GPU_MIN`, `CGV_GPU_MAXDEGS`, `CGV_GPU_MIN_GB` (when cgv_gpu uses the GPU);
`CGV_XTARGET` (GPU batch size); `CGV_LSYNC`, `CGV_CS_NU` (CPU small-layer mode).

Run one GPU job per device at a time (each sizes its tables to the free memory).
