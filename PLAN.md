# Plan: a Rust implementation of the fast GV algorithm from PR #94

Status: **investigation done, phases 0 and 0d done, nothing implemented yet.** Progress is
tracked by the check boxes in section 5. Where sections 1–5 mention `rug` for the new code, read
"big integers": decision 6.10 makes `rug` optional. Branch `gv-modular-rust` is PR #94
(`natemacfadden:cgv-cc-cpu`, 11 commits) rebased onto `origin/main` (6fe293a, v0.3.0), plus our
commits. The PR branch itself (`cgv-cc-cpu`) is untouched. The follow-up PRs #95, #97 and #100
are covered in section 1.1.

Section 6 lists the decisions that need an answer before implementation starts.

## 1. What PR #94 does

It adds `cgv/`, a standalone C program (`gv.c`, 2226 lines in #94, 2279 after the #100 sync, plus a 88-line driver), and wires it
in as `compute_gv(..., backend="cgv")`: `build.rs` compiles it with `cc`, the binary is embedded
in the crate with `include_bytes!`, written to a cache directory on first use, and run as a
subprocess through a text file. A Python helper (`_cgv_run.py`) shells out to the `normaliz`
binary for Hilbert bases. Scope: threefold hypersurfaces, `max_deg` only, GV (GW derived from GV
in Python). `cgv/gpu.cu` (CUDA/HIP) is in the tree but not built by the PR.

The mathematics is still HKTY, but reorganised. In pipeline order:

| # | Step | What is different from the current code |
|---|------|------------------------------------------|
| 1 | **Arithmetic mod ~62-bit primes** (Montgomery), 2–4 primes carried as "lanes" of one coefficient; integers recovered by CRT, one extra prime as a check, more passes if the check fails | replaces `rug::Rational`/`Float` everywhere in the hot path |
| 2 | **Curve classes are packed `u128` keys** with `key(a+b) = key(a)+key(b)`; polynomials are sparse `(key, degree, value)` lists and open-addressing tables | replaces `Semigroup` (all points enumerated up front) and `Polynomial` (`HashMap<index, T>` over that enumeration) |
| 3 | **Candidates**: only curves with ≤ 2 negative GLSM intersections carry period data. They are the lattice points of the cones $K_T = \text{Mori} \cap \{Q_r\cdot C \ge 0,\ r \notin T\}$, $\lvert T\rvert \le 2$, enumerated from each cone's Hilbert basis. The Mori cone is never enumerated. Fallback without cone data: enumerate the whole semigroup and filter | new; the fallback is the analogue of `Semigroup::with_max_degree` |
| 4 | **Fundamental period** at each candidate: `c0`, `c1[a]`, and $S_2 = \sum_{a\le b} K_{ab}\, c_{2,ab}$ with $K_{ab} = \sum_t w_t \kappa_{tab}$ | same formulas as `compute_c_{0,1,2}neg`, but factorials/harmonic numbers come from tables mod p, and `c2` is contracted immediately instead of kept per index pair |
| 5 | **One instanton series** $I = \sum_t w_t\, \mathrm{inst}_t = c_0^{-1} S_2 - \sum_a \alpha_a V_a$, $V_a = \tfrac12 \sum_b K_{ab}\alpha_b$ | replaces the $h^{1,1}$ `inst` polynomials and the per-pair beta/F products |
| 6 | **Extraction** by degree layers: after subtracting lower curves, $I[C] = \deg(C)\, A_C$; each curve's correction is $s_C\, z^C \exp(C\cdot\alpha)$, computed with the Euler recurrence into per-degree hash tables | replaces `invert_series` (qN windows, `exp_pos_neg`, `pow`, `li_2`) |
| 7 | **Multicover**: $GV_C = A_C - \sum_{n\ge2,\ n \mid C} GV_{C/n}/n^3$, by gcd level | replaces the $\mathrm{Li}_2(q^N)$ series. Note $A_C$ *is* the GW invariant |
| 8 | CRT lift, lane-count probe (a cheap run at 3/4 of the degree to predict the size of the GVs) | new |

Extras in `cgv/` that the PR does not expose through cygv: a certification build (`CGV_MAJ`: the
same pipeline over upward-rounded reals to prove the CRT lift exact), a distributed residue mode,
a glibc low-memory switch, the GPU stage.

### 1.1 The follow-up PRs: #95, #97, #100

All three are stacked on #94 (each contains the previous ones; only the last commit or two are
new). Read from their diffs against our copy of `cgv/`.

- **#95, GPU variant.** Builds `cgv/gpu.cu` (CUDA or HIP) in `build.rs` when `nvcc`/`hipcc` is
  found and adds `device=` to the Python API. Already a future option (6.8); it only confirms
  that the device belongs in the options struct (4.1.6) and that a backend must be able to
  decline (4.1.3): the C keeps lightcone runs on the CPU.
- **#97, any phase of the ambient toric variety, including vex fans** (MacFadden–Sheridan,
  arXiv:2512.14817). In a vex phase there are curves of negative anticanonical degree
  $m = q_0 \cdot C < 0$, where the HKTY coefficient has a pole; it cancels against the class of
  the stratum $V_S$ those curves force into the CY (McAllister's prescription via Wang,
  arXiv:1910.14440). For fine fans on reflexive polytopes $S$ is a 3-cone of the fan and $V_S$ a
  toric curve, so the extra input is one rational class $d_S$ per vex 3-cone. Two parts:
  - *Core* (`gv.c`, ~90 lines): curves with $m < 0$ and exactly three negative intersections
    $S$ that match a vex stratum get $S_2 = U_0\,(w \cdot d_S)$ and no `c0`/`c1`, with
    $U_0 = (-1)^{e}\, \prod_{S}(k_i-1)! \,/\, \big(\prod_{\notin S} b_i!\,(\lvert m\rvert-1)!\big)$;
    any other $m < 0$ curve contributes nothing. The $K_T$ list gains one cone per vex 3-cone,
    the factorial tables cover $\lvert m\rvert$, the strata values are rationals reduced mod p.
    Threefold hypersurfaces only; the C's certifier build rejects vex input.
  - *Preprocessing* (`cgv/tools/cgv_phase.py`, ~300 lines of Python with exact rational linear
    algebra and normaliz): from the fan's maximal cones, `q` and `kappa`, derive the rays
    (integer Gale dual of `q`), the toric intersection ring (to check `kappa`), the wall curves
    (default Mori rays), the vex cones (rays sharing no facet of the hull of the rays) and the
    exact stratum classes.
  - Validated by the author against an independent exact-rational implementation and across
    flops, with 16 frozen references in `cgv/tests/refs_vex`.
  - It also exposes a latent bug in the reference pipeline: `compute_c_1neg`/`_2neg` cast
    $q_0 \cdot C$ to `u32`, so a negative anticanonical degree silently wraps (section 7).
- **#100, minimal inputs, one candidate path, lightcone GVs.**
  - The Mori cone is given by any description: spanning vectors (`mori_rays`, redundant
    allowed) or inequalities (`mori_hyperplanes`). Lattice points and the Mori cone's own
    Hilbert basis are no longer inputs; the `generators` key is gone.
  - Candidates always come from the $K_T$ Hilbert bases. The saturation check and the full
    enumeration it fell back to are deleted ("it decided only agreement with cygv, never cgv's
    output", and cost up to 770 s / 79 GB on one region, now 91 s / 1 GB). The C therefore
    always needs normaliz.
  - *Lightcone*: given curves $p$, keep only their backward lightcones $\{C : p - C \in
    \text{Mori cone}\}$, tested as $H C \le H p$ row by row with $H$ the Mori facets. Applied
    at three points: when listing candidates, when choosing the curves to apply in the
    extraction, and on the output (products of kept curves can land outside, with wrong
    values). `max_deg` defaults to the largest degree of the $p$. This is exactly the
    target-point trimming discussed under 6.3, and it replaces the old "subset of generators"
    mode, which only worked for faces.
  - New input file format ("cgv 2"), clearer CRT error, progress lines.

What this changes in the plan: 6.1 (what cone data is), 6.2 (no separate Mori Hilbert basis),
6.3 (target points decided; phases added), 6.11 (the same features in the reference), phases
4b, 6–6c and 7a–7d, section 5.2, and section 7.

## 2. Measurements

This machine (8 threads), `benches/data/h11_9.yaml`, GV, results identical (12404 GVs compared at
degree 14):

| Case | reference (`cygv --file`) | C, no cone data | C, with cone data |
|------|---------------------------|-----------------|-------------------|
| `h11_9`, max_deg 10 | 5.9 s | – | – |
| `h11_9`, max_deg 14 | 222 s | 3.1 s (2.9 s enumerating the semigroup) | 0.15 s |
| `h11_9`, max_deg 18 | not run (bench table: 800 s on 16 cores) | 23.6 s (21.0 s enumerating) | 1.5 s |
| `h11_8`, max_deg 20 | not run (bench table: 1600 s on 16 cores) | not run | 7.9 s |

What this says about where the speed comes from:

- Steps 1, 2, 4–7 alone (no cone data) are worth roughly 70× at degree 14. Once they are in
  place the *pass itself* is tiny (0.2 s at degree 14) and enumerating the semigroup becomes
  ~90 % of the run.
- The cone trick (step 3) is worth another ~15–20× on top, and is the only part that needs
  Hilbert bases. `normaliz` took ~1.7 s per geometry for the 91 cones of `h11_9` (cached after).
- With cone data, extraction (step 6) dominates deep runs (7.6 of 7.9 s for `h11_8` at 20).
- Since the oracle sync (0d) the C has no path without cone data any more (#100 deleted it), so
  the middle column can no longer be reproduced with it. The performance target for tuning is
  the with-cone-data column (phase 5).

## 3. How much of the existing machinery can be shared

**The computational core: essentially nothing.** Module by module:

| Existing | Usable by the new algorithm? |
|----------|------------------------------|
| `semigroup` | No. The point of the new approach is not to enumerate it, and the fallback needs packed keys, not `DMatrix` columns in a `HashSet<DVector>`. (The reverse direction is possible: `reduce_generators` + packed-key BFS could later speed up `Semigroup` itself.) |
| `polynomial`, `PolynomialProperties` | No. Indexed by position in the enumerated semigroup; `mul` goes through `monomial_map`. |
| `PolynomialCoeff<T>` | No. It requires ordering, comparison with floats, rounding, `abs`, `zero_cutoff`; a residue has none of these. It is also built around in-place `*Assign` on heap bignums, while a residue lane is a `Copy` array of `u64`. |
| `factorial` | No. Per-call bignum products vs. precomputed tables mod p. |
| `fundamental_period` | Formulas only. A generic version over "ring element" would mean rewriting `compute_c_*neg` around a new trait and table lookups for ~150 lines of shared text, and the new code wants $S_2$, not `c2` per pair. |
| `instanton`, `series_inversion` | No. Different formulas and data flow. |
| `misc::process_int_nums` | Barely: the new code needs the symmetrised $K_{ab}$, a dozen lines. |

**The shell around the core: all of it.**

- `io::Input` (YAML parsing, validation, result sorting and formatting) and the CLI.
- `hkty::compute_gvgw_strings`, already the single runtime dispatch for the CLI and Python, with
  the right return type (`Vec<((DVector<i32>, usize), String)>`).
- `src/python.rs` and `python/cygv/hkty.py` (the subprocess/Ctrl-C handling, input normalisation).
- Dependencies: `rayon` (the C thread pool maps onto the private-pool pattern `run_hkty` now
  uses), `rug` (`Integer` for the CRT and `Rational` for GW replaces ~150 lines of hand-rolled
  bignum in C; `Float` with directed rounding could replace `long double` + `fesetround` if the
  certifier is ported), `nalgebra` at the API boundary only.
- Tests and benches: `benches/common` scenarios, `benches/data`, `tests/full_hkty.rs`.

**The reference implementation as an oracle.** Having both in one crate allows differential
tests at stage level, not only end to end: e.g. $\sum_t w_t\,\mathrm{inst}_t$ from
`compute_instanton_data` reduced mod p must equal the new $I$; `c0`/`c1` likewise.

Things Rust makes simpler than the C: `u128` is native on every target (the C needs gcc/clang,
no MSVC); the three compilations of `gv.c` with `-DNL=2,3,4` become one `const NL: usize`
generic, which is the pattern `run_hkty` already uses; no embedded executable, temp files or
text protocol.

## 4. Recommended structure

**One crate, one new self-contained module tree, shared front end. No intertwining with the
existing stage modules, and no separate crate for now.**

- Not intertwined: there is no code to share in the core, and forcing `Polynomial` or
  `PolynomialCoeff` to cover residues would damage the reference implementation's readability,
  which is the reason to keep it.
- Not a separate crate: the things that *are* shared (io, dispatch, Python, benches, the
  stage-level oracle tests) all live in this crate, and a workspace means touching maturin
  config, both CI cache keys, `deploy.yml` and the release process for two crates. The module
  boundary will be clean enough to split out later if there is a reason (e.g. a GMP-free build).

Layout (6.7):

```
src/lib.rs                  InvariantKind, CYKind, Method; re-exports
src/dispatch.rs             compute_gvgw_strings and the typed entry points: pick the method, call it
src/io.rs, src/python.rs    unchanged in place (shared front end)
src/cone.rs, src/cone/      shared by both methods (6.11): double description (7a), Hilbert bases (7b)
src/phase.rs                shared: fan -> strata, default Mori rays and grading, kappa check (7c)

src/reference.rs            run_hkty (today's src/hkty.rs minus the dispatch)
src/reference/semigroup.rs, polynomial.rs (+ polynomial/), factorial.rs,
              fundamental_period.rs, instanton.rs, series_inversion.rs, misc.rs
                            moved as they are, with their error modules and tests

src/modular.rs              entry point: run(...), options and stats structs, pass loop, lane probe
src/modular/field.rs        Montgomery arithmetic, prime search, Lanes<const NL: usize>
src/modular/key.rs          Key trait; packed u64/u128 with per-coordinate fields, inline arrays; bound check, hash
src/modular/table.rs        slot map, sparse polynomial, per-degree level tables, heap
src/modular/candidates.rs   generator reduction, semigroup enumeration, cone union, saturation check
src/modular/period.rs       tables (factorials, harmonic numbers), c0 / c1 / S2 per curve
src/modular/instanton.rs    sparse mul / reciprocal, alpha, V, I
src/modular/extraction.rs   the extraction interface (4.1) and its CPU implementation: layers, exp via Euler recurrence, sharded fold
src/modular/multicover.rs   A -> GV, GV -> GW
src/modular/crt.rs          Garner + big integers (pure-Rust crate, 6.10), stability check
```

`src/hkty.rs` currently mixes three things that end up in different places: `run_hkty` (goes to
`reference`), `compute_gvgw_strings` (the runtime dispatch used by the CLI and Python, goes to
`dispatch`), and the eight typed `compute_g{v,w}_{rat,float}_{threefold,nfold}` wrappers (also
`dispatch`; whether all eight survive unchanged is a detail to settle when the `Method` argument
is added). Per CLAUDE.md, adding an option there means touching every dispatch site.

### 4.1 Design constraints from future options (6.8)

These cost little now and are what keeps the later additions local. They come from how the C
connects the pieces (`cgv/gpu.h`, `CGV_RESIDUES` / `tools/crt_combine.py`, the `CGV_MAJ` build).

1. **Extraction sits behind a narrow interface.** In the C the GPU replaces exactly one stage,
   through one call: in go the field constants per lane, the key layout, `max_deg`, the union
   of the alpha supports (keys, degrees, values), the initial $I$ (keys, degrees, values) and
   the `1/n` and `n` tables; out comes a list of (key, degree, $A$). Everything before it
   (candidates, period, instanton) and after it (multicover, CRT) is untouched. So `extraction`
   exposes that signature as a trait or function type with the CPU code as one implementation,
   and nothing outside the module reaches into its tables.
2. **Plain data at that boundary.** Flat arrays of fixed-size `Copy` values with a defined
   layout (`#[repr(C)]`): residues as `[u64; NL]` in Montgomery form with the constants passed
   alongside, keys as plain integers or inline arrays. No heap-owning element types, no
   `HashMap` in the interface. A device kernel then does bit-identical arithmetic, which is
   what lets the CPU path serve as its test oracle.
3. **A backend may support a subset and decline.** The `Key` trait (6.4) and the
   multi-component $I$ for n-folds (6.3) multiply the layouts a device kernel would have to
   handle. The interface therefore lets an implementation refuse a configuration, and the
   caller falls back to the CPU, as the C already does for small or busy cases.
4. **One pass = residues; lifting is a separate step.** The pass loop takes the index of its
   first prime in a fixed, documented prime sequence and returns residues per curve. The CRT
   lift and the stability check take residue sets from any source. That is the distributed
   mode's whole requirement (plus a serialisation of the residue set), and it is also what
   `check_primes` and "at least n primes" from the certifier need.
5. **Stages are generic over the coefficient, with hooks for the certifier.** The certifier
   (6.5) is the pipeline over upward-rounded non-negative reals, but in the C it is more than a
   type swap (`CGV_MAJ` in `gv.c`): alpha gets one coordinate per charge row (`ndiv + 1`
   instead of $h^{1,1}$) to keep bounds tight, the factorial and inverse tables are built
   differently, subtraction and negation become addition and identity, and during extraction
   each point's source value is replaced by the bound implied by the candidate GVs. So:
   - the coefficient trait has zero, from-integer, add, sub, *neg*, mul, fused multiply-add,
     exact inverse, is-zero, and table construction belongs to the coefficient type;
   - the period stage produces alpha in a coordinate system chosen by the caller (the
     $h^{1,1}$ basis, or the charge rows), not hard-wired to $h^{1,1}$;
   - the extraction interface (constraint 1) accepts an optional per-point source override.
   Only the lane implementer needs to satisfy constraint 2. The alternative, a separate code
   path for the certifier, would drop this constraint; the hooks are cheap if designed in from
   phases 2–3 and expensive to retrofit.
6. **Options travel in a struct.** `check_primes`, the lane override, the thread count, the
   saturation opt-out and anything added later (device, first prime index, minimum primes,
   tuning thresholds for the level-synchronous mode and table load factors) are fields of one
   options struct with defaults, passed through `dispatch` as a unit. A new option is then a
   new field, not a new positional argument at every dispatch site. The C's environment
   variables map to fields; nothing reads the environment inside the library.
7. **Reporting is data.** Phase timings, candidate and point counts, primes used, missing
   Hilbert basis elements come back in a stats struct instead of being printed. The CLI and
   Python decide what to show; the benches and a future profiler read the same struct.
8. **Memory in few large allocations.** The low-memory switch in the C is an allocator
   setting (glibc returns large blocks to the system when freed), not an algorithm. It needs
   nothing structural beyond what the tables already do, and the library does not install a
   global allocator, so a binary or the Python module can choose one.

## 5. Phases

Tick the boxes as the work lands. The phases are listed in execution order; their numbers are
kept for cross-references, which is why phase 6 comes before phase 5. The order is 0–4b, then
6–6c, then 5–5b, then 7 onwards (7a can start any time). Tuning without cone data would
optimise the wrong profile: there the enumeration is ~90 % of the run, with cone data the
extraction is (section 2). Phases 1–4 need no decision except 6.5–6.7.

- [x] **0. Branch hygiene** (6.6): the PR's integration glue is removed and `cgv/*.c` stays as
  a development oracle and performance baseline.
- [ ] **0b. Move the existing pipeline under `src/reference/`** (6.7), as a commit of its own
  with no functional change. Done before any new code so that the diff stays a pure move.
  - [ ] `git mv` the modules and fix the `use` paths in `src`, `tests/`, `benches/` and the
    CLI.
  - [ ] Split `src/hkty.rs` into `reference.rs` and `dispatch.rs`.
  - [ ] Update CLAUDE.md.
- [ ] **0c. Drop floating-point computation** (6.9), as its own commit.
  - [ ] Check whether CYTools or other callers pass `prec` (before starting).
  - [ ] Remove `prec` and the `Float` variants from the Rust code, `io::Input` and the CLI.
  - [ ] Remove the `gv-float` / `gw-float` bench variants.
  - [ ] Python: `prec` deprecated for one release (`DeprecationWarning`, exact values
    returned), `mpmath` dropped.
- [x] **0d. Sync the oracle** (6.6): `cgv/` is the #100 version, without the integration glue.
- [ ] **1. Primitives.**
  - [ ] `field`, with unit tests of its operations against the big-number crate.
  - [ ] `key`, with pack/unpack round trips and overflow detection.
  - [ ] `table`.
- [ ] **2. Period and instanton, enumeration path only.**
  - [ ] Candidates by full enumeration.
  - [ ] `period`.
  - [ ] `instanton`.
  - [ ] $I$, `c0` and `c1` tested against the reference stages mod p.
- [ ] **3. Extraction, multicover, CRT**, single-threaded and simple (`apply_curve` only), for
  threefold hypersurfaces.
  - [ ] Provisional choice of the pure-Rust big-number crate, confined to `crt` and
    `multicover` (6.10, phase 5b).
  - [ ] `extraction`.
  - [ ] `multicover`.
  - [ ] `crt`, with the test knob that forces several CRT passes (5.2).
  - [ ] End-to-end differential tests against `run_hkty` on consistent geometries only (5.2):
    P(1,1,2,2,2)[8] and `benches/data`, not the toy threefold of `examples/`.
- [ ] **3b. Complete intersections and n-folds** (6.3). Done before tuning so that phase 5
  optimises the final data layout.
  - [ ] Multi-row `q0`.
  - [ ] The $M$-component $I$.
  - [ ] The exponent-2 multicover.
  - [ ] Check which geometry `examples/fourfold.yaml` and the `fourfold` bench model each are
    (5.2).
  - [ ] Differential tests against the reference on the fourfold models.
- [ ] **4. Front end.**
  - [ ] The `Method { Auto, Modular, Reference }` argument through `dispatch`, `io::Input`,
    the CLI and Python. `auto` falls back to `reference` for requests `modular` does not
    support; an explicit `modular` returns an error for them.
  - [ ] GW.
  - [ ] `n_threads` exposed in the Python API (it is accepted by `src/python.rs` but
    `python/cygv/hkty.py` always passes `None`).
  - [ ] `generators` optional in `io::Input` when `mori_hyperplanes` is given, with validation
    reworked accordingly.
  - [ ] Error messages for non-integer results that name the likely cause (an inconsistent
    input, or a cone that is not a face) instead of only "CRT did not stabilize".
  - [ ] `auto` made the default, once the differential tests of phases 3 and 3b pass.
- [ ] **4b. Vex phases in the reference pipeline** (6.3, 6.11).
  - [ ] The `strata` input through the front end.
  - [ ] The three-negative bucket and the extra term in the reference instanton step.
  - [ ] $m < 0$ kept out of the 1/2-negative routines (this also fixes the `as u32` wrap of
    section 7), and the error for $m < 0$ with fewer than two negatives.
  - [ ] Tests against `cgv/tests/refs_vex` and the C oracle.
- [ ] **6. Cone fast path through optional input.** Runs before phase 5, and fixes the
  interface that phase 7 plugs into.
  - [ ] Cone data (6.1) accepted by the Rust API and YAML. The Mori cone may be described by
    spanning vectors or by inequalities, plus equations for faces (section 7).
  - [ ] `candidates` enumerates each $K_T$ from its Hilbert basis, with the containment and
    vanishing-form skips. The candidate source is an interface, so phase 7b can plug in a
    native one (Hilbert bases or direct enumeration) without touching the rest.
  - [ ] The Python wrapper may produce the cone data with `normaliz`.
  - [ ] Tests: candidate sets with and without the skips, face inputs (5.2).
- [ ] **6b. Target points for `modular`** (6.3).
  - [ ] The facet test at the three filter points (listing candidates, choosing the curves to
    apply in the extraction, the output).
  - [ ] Optional `max_deg`, defaulting to the largest degree of the targets.
  - [ ] Differential tests against the full computation and the C's lightcone mode.
- [ ] **6c. Vex phases in `modular`** (6.3).
  - [ ] The $m < 0$ period branch.
  - [ ] The extra $K_T$ per stratum.
  - [ ] Tests against the reference (phase 4b), the C oracle and `cgv/tests/refs_vex`.
- [ ] **5. Parallelism and tuning.** Target parity with the C *with cone data* (section 2).
  - [ ] Rayon over candidates / product chunks / curves per layer.
  - [ ] Sharded $I$ with per-shard locks.
  - [ ] The level-synchronous mode for few large curves.
  - [ ] Lane probe.
  - [ ] Bench scenarios.
  - [ ] Stress test (5.2).
  - [ ] Compile time and binary size measured: the code is generic over the key type, the lane
    count and the coefficient type, and the combinations multiply.
- [ ] **5b. Make `rug` optional** (6.10). After phase 4, so that `auto` already routes most
  requests to `modular` when the reference pipeline's default build gets slower.
  - [ ] The pure-Rust crate chosen from benchmarks (if it is `malachite`, confirm that
    LGPL-3.0-only fits GPL-3.0-or-later).
  - [ ] `PolynomialCoeff` implemented for it.
  - [ ] Backend-neutral result types.
  - [ ] The `rug` feature.
  - [ ] CI workflows and CLAUDE.md simplified (including whether Windows still needs the
    MSYS2/MINGW64 route).
- [ ] **7a. Double description (committed, not a prototype).** Can start any time. Support
  hyperplanes from spanning vectors, in any dimension, with exact integers. Needed regardless
  of how 7b turns out: it gives the Mori cone's facets (target points, `mori_rays` input, the
  default grading) and the hull facets that phase 7c needs.
  - [ ] Decide: our own implementation or `howzat`.
  - [ ] Implementation and tests.
- [ ] **7b. Candidate lattice points without normaliz: time-boxed comparison, then integrate
  or drop** (5.1). Independent of phases 1–5 and can run alongside them.
  - [ ] Direct enumeration of polytope lattice points.
  - [ ] A correctly specified Pottier algorithm.
  - [ ] Comparison on the benchmark geometries (5.1).
  - [ ] The winner integrated, or the phase dropped.
- [ ] **7c. Phases from a fan** (6.3): derive `strata`, the default Mori rays and grading, and
  optionally check `kappa`, from the fan's maximal cones and `q`, without normaliz. Port of
  `cgv_phase.py`, with exact rational linear algebra on small matrices using the big-number
  crate of 6.10; estimated 400–600 lines with tests. A directly supplied `strata` input stays
  accepted. The vex path is only fully normaliz-free if 7b succeeds: each stratum adds a $K_T$
  whose lattice points are still needed.
  - [ ] Rays as the integer Gale dual of `q`.
  - [ ] Walls (3-faces shared by two maximal cones) and their classes by Cox–Little–Schenck
    6.4.4.
  - [ ] Vex cones: the fan's 2- and 3-faces whose rays share no facet of the hull of the rays
    (facets from 7a).
  - [ ] The toric intersection ring for the `kappa` check.
  - [ ] Rejection of inputs outside the paper's assumptions (a fine fan on the boundary points
    of a reflexive 4-polytope not interior to facets), e.g. vex 2-cones (MacFadden–Sheridan
    Prop. 5).
  - [ ] Tests piece by piece against `cgv_phase.py` (normaliz needed on the development
    machine only) and end to end against `cgv/tests/refs_vex`.
- [ ] **7d. Saturation for the reference pipeline** (6.11). Only if 7b is integrated;
  afterwards `auto` has one saturation semantics (6.2).
  - [ ] The lattice points of the Mori cone up to `max_deg`, which is exactly the saturated
    semigroup the reference enumerates, from the 7b machinery (directly if 7b uses polytope
    enumeration, from the Mori cone's Hilbert basis if it uses Pottier).
  - [ ] `mori_hyperplanes` accepted by both methods.
- [ ] **8. Certifier** (6.5).
  - [ ] The majorant run, fed with the candidate GVs, behind its own switch.
  - [ ] Tests against `cgv_maj` on the benchmark geometries.
  - [ ] Optional: the generator-reduction back-port to `Semigroup`.
- [ ] **9. Wrap up.**
  - [ ] Confirm with the PR author that the provenance of `refs_vex` allows shipping it as
    fixtures (5.2).
  - [ ] The selected oracle references converted into Rust fixtures (5.2), *before* `cgv/` is
    removed.
  - [ ] `cgv/` removed.
  - [ ] README, CLAUDE.md and docs updated.
  - [ ] Decide what happens to PR #94 and credit its author (section 7).

Outside the phases:

- [ ] Optional: a minimal fix (an error) for the `as u32` wrap in `compute_c_1neg` /
  `compute_c_2neg` on `main`, independently of this work (section 7).

Follow-ups after the main plan, not tracked here: GPU and the other C features (6.8).

### 5.1 Phase 7 in detail: cone computations without normaliz

Goal: produce the cone data of 6.1 inside the crate, so the fast path does not depend on an
external program. No suitable crate was found on crates.io (no Hilbert basis crate, no normaliz
or 4ti2 bindings; `howzat` does double description only and has not been evaluated).

What is actually needed is narrower than a Hilbert basis. `modular` needs the *lattice points*
of each $K_T$ up to `max_deg`, and the reference (7d) needs those of the Mori cone up to
`max_deg`. Each is the set of integer points of a bounded polytope
$P_T = K_T \cap \{w \cdot x \le \texttt{max\_deg}\}$, whose inequalities are all known once
7a has produced the Mori cone's facets: those facets, the rows $Q_r$ for $r \notin T$, and the
degree bound. Hilbert bases are one route to those points (the C's: enumerate non-negative
integer sums of the basis), not the only one. Two candidates, compared in a time-boxed
prototype:

**A. Direct enumeration of the polytope's lattice points** (suggested in review).
- Recursive coordinate-by-coordinate enumeration: at each level, bounds on the next coordinate
  given the ones already fixed. Exact bounds come from linear programming over the rationals, or
  cheaper floating-point LP with exact checks at the leaves; Fourier–Motzkin projection is the
  LP-free alternative but can blow up in dimension ~10.
- Risk: it is not automatically output-sensitive. Weak bounds make it explore empty branches,
  and an LP per node is costly at $10^5$–$10^6$ points. A change of lattice basis aligned with
  the cone (e.g. from a simplicial subcone) can tighten the bounds.
- Gives the reference's saturated semigroup (7d) with no Hilbert basis at all.

**B. Pottier's algorithm, correctly specified.** Cutting a cone with a known Hilbert basis by
one halfspace: split the basis by the sign of the new inequality, keep adding irreducible sums
of a positive and a negative element until nothing new appears, keep the non-negative ones.
- Start from a basis of the lattice ($\pm e_i$, or a lattice basis of the linear span for a
  face), then cut by the Mori facets *without* truncation: during this stage the cones are not
  inside the Mori cone, the grading is not positive on them, and truncating by degree would be
  unsound. The result is the Mori cone's Hilbert basis (which 7d can use).
- Then cut by the rows $Q_r$, $r \notin T$, *with* truncation at `max_deg`: from here every cone
  lies inside the Mori cone, the grading is positive, and sums only increase degree. The nested
  cones ($K_\emptyset \subset K_{\{r\}} \subset K_{\{r,s\}}$) share intermediate bases, and
  two-divisor cones with vanishing form are skipped before computing them.
- Risk: the untruncated first stage, i.e. the size of the Mori cone's Hilbert basis and of the
  intermediate bases, at $h^{1,1}$ = 10–11.

Both need the facets (7a), handle faces through the equations of their span, and must produce
saturated results whatever lattice points the input lists (decision 6.2).

Prototype protocol:

- Run on the four `benches/data` geometries (and `h11_9` under its sparse grading) at their
  Heavy and Huge degrees, plus a face input and the vex references.
- Correctness oracle: identical candidate sets to the C with normaliz-produced cone data, and
  identical GVs end to end; for B also normaliz's Hilbert bases truncated to `max_deg`.
- Unsaturated input (`h11_11.yaml`): the candidates must come out saturated.

Exit criteria:

- *Integrate* whichever approach takes at most a few seconds on every benchmark geometry
  (normaliz: ~1.7 s for the 91 cones of `h11_9`) and matches the oracle. It becomes
  `src/cone/` (shared, 6.11), feeding the phase 6 interface; the optional-input path stays for
  callers that already have cone data.
- *Drop* if neither does. The optional-input path of phase 6 remains the only fast path, and
  the findings are recorded here.

### 5.2 Testing strategy

Layers, from cheapest to most independent:

1. *Unit tests* per module (field arithmetic against the big-number crate, key round trips,
   table operations, Hilbert bases against brute force on small cones, fan derivations against
   `cgv_phase.py` output recorded as constants).
2. *Stage-level comparisons* with the reference pipeline mod p: `c0`, `c1`, the contracted
   $S_2$ and $I$ against the reference's period and instanton data reduced modulo a prime. This
   pins a mismatch to a stage.
3. *Differential tests in CI* (`tests/`): `modular` and `reference` on the same inputs must give
   identical GVs and GWs. Inputs must be *consistent geometries*: the toy threefold in
   `examples/threefold.yaml` (the same data as `python/tests/test_hkty.py`) is not one, and
   there the two methods legitimately disagree, because the reference reads each GV from the
   first nonzero coordinate (`src/series_inversion.rs:161`) while `modular` uses the
   grading-weighted sum. Replace it with P(1,1,2,2,2)[8] for $h^{1,1} = 2$, and add a cheap
   consistency check to the reference ($\mathrm{inst}_t[C]/C_t$ equal across $t$) so that
   inconsistent input is an error in both methods. Inputs then: P(1,1,2,2,2)[8], the realistic
   models at low degree including `h11_9_plike` (hundreds of degree levels, which stresses the
   level tables most), a face of the Mori cone, target points, complete intersections, n-folds,
   vex phases. Kept fast enough for every CI run.
   Paths that need forcing, because they never trigger at CI-sized degrees:
   - *Several CRT passes*: a test knob that forces a small lane count and a minimum number of
     primes, so the multi-pass merge (curves that appear only in later passes) and the
     big-integer lift for three or more primes run.
   - *"CRT did not stabilize"*: a non-face subcone input must fail with the explanatory error.
   - *Candidate skips*: the containment and vanishing-form skips (`gv.c`, `build_candidates`)
     silently drop curves if wrong; compare candidate sets with and without them.
   - *Keys*: the coordinate-bound abort, and the fallback from packed keys to array keys.
4. *Frozen fixtures*: expected values that come from neither implementation. A selection of
   `cgv/tests/refs` (outputs of cygv 0.2.3) and all or most of `refs_vex` (checked against an
   independent exact-rational implementation and across flops), converted to small compressed
   files under `tests/data/` with a size budget of about 1 MB, and read by an integration test.
   Converted in phase 9 at the latest, before `cgv/` is removed. They also catch a bug shared by
   both implementations, which layer 3 cannot. Excluded from the published crate (`exclude` in
   `Cargo.toml`). Gap: there is no independent reference for complete intersections or n-folds;
   look for published tables to use as fixtures. Also open: `examples/fourfold.yaml` describes
   itself as a hypersurface while the `fourfold` bench model has a nef partition; check which
   geometry each one is before relying on them. And before shipping `refs_vex` as fixtures,
   confirm with the PR author that their provenance allows it.
5. *Stress test*, in the spirit of `cgv/tests/stress.py`: repeated runs at several thread
   counts, outputs compared with each other, to catch races in phase 5's parallel code. Marked
   `#[ignore]` and run before releases and on parallelism changes, not on every CI run.
6. *The C oracle*, until phase 9: development-time comparisons through `cgv/tools`, never in
   CI (it needs a C toolchain and normaliz). Same prime sequence, so residues can be compared
   one for one when narrowing a mismatch down.

## 6. Decisions needed

1. **Hilbert bases (the main one).** The last ~15–20× needs the Hilbert basis of the Mori cone
   and of every $K_T$ (91 cones for `h11_9`), plus the Mori cone's support hyperplanes. The PR
   gets them from the `normaliz` binary in Python.

   "Cone data", following #100: (i) the Mori cone's support hyperplanes (facets), needed for
   target points and the coordinate bound; (ii) the Hilbert basis of each $K_T$: every $T$ with
   $\lvert T\rvert \le 2$ (13 + 78 = 91 cones for `h11_9`, of which 19 are enumerated after the
   skips) plus one per vex 3-cone (6.3); (iii) the divisor set $T$ of each cone, which enables
   the vanishing-form skip. Independent of the grading vector and `max_deg` when computed
   exactly. `modular` does *not* need the Mori cone's own Hilbert basis: the $K_T$ already contain
   every candidate, and the basis was only ever used for the saturation check, which #100
   deletes. (The reference pipeline needs it to saturate, 6.11.)
   The data is trusted: a wrong basis for some $K_T$ silently drops candidate curves. If phase
   7b settles on enumerating lattice points directly (5.1), the native path needs only the
   facets; the Hilbert bases remain the form in which callers can supply cone data.

   **Current direction**: optional input first (phase 6), then a time-boxed prototype of a
   native implementation (phase 7b, section 5.1) that either replaces the need for normaliz or
   is dropped. Linking libnormaliz (C++) is ruled out as contrary to the goal of the rewrite.
   The Mori cone itself can be described either way, as in #100: spanning vectors (the
   existing `generators` field, now read as "any vectors spanning the cone", lattice points or
   just rays) or inequalities (a new `mori_hyperplanes` field). Still open: whether the Python
   wrapper should call `normaliz` when it is installed.
2. **Saturation semantics. Direction: saturate by default.** The reference uses the semigroup
   the generators *generate* (non-negative integer combinations); the correct object is *all
   lattice points of the cone they span*, $\sigma \cap \mathbb{Z}^n$. The two differ whenever
   the input misses a Hilbert basis element of the cone, e.g. generators (1,0), (1,2) miss (1,1).
   This is a common problem for cygv users today: CYTools passes the rays plus a finite sample
   of lattice points, which can miss an element. `benches/data/h11_11.yaml` does (its 2032
   generators miss `[2,0,0,2,1,0,0,-1,-1,0,0]`, degree 6; not an artefact of the file's
   trimming at degree 16), and the reference then silently computes over a semigroup with holes.

   So the new method **works over the saturated semigroup by default**, rather than treating
   saturation as something to check. Since #100 this needs no separate step: the candidates come
   from the $K_T$ Hilbert bases, which are computed from the cone's description and so are
   saturated by construction.
   It is skipped only when (i) the dimension is too large for the computation to be practical,
   or (ii) the user disables it; in both cases the input generators are used as given, as the
   reference does. Consequences:
   - This raises the stakes of phase 7: computing the candidates natively (7b) is no longer
     only a speed-up, it is what makes the default correct without an external program.
   - Unlike the C after #100, we keep the enumeration of the *generated* semigroup, but only as
     the path for the opt-out cases above (and as the phase 2 stepping stone). Given only rays it
     still runs, but generates the semigroup of the rays, which is usually far from saturated;
     that combination gets a warning. With inequalities only it is unavailable.
   - `auto` must not mix semantics silently. `modular` saturates; the reference saturates only
     once phase 7d exists. Until then, a request that falls back to the reference reports (in
     the stats struct, and as a warning in the CLI and Python) that it used the generated
     semigroup. After 7d both methods saturate the same way.
   - Outputs can differ from the reference on unsaturated input, by design. Differential tests
     need saturated inputs, or the reference fed the completed generator set. Add `h11_11` as
     an explicit test of both.
   - Reporting which Hilbert basis elements the input missed (the C's old note) would need the
     Mori cone's basis again. Dropped, as in #100; it can come back as a diagnostic if wanted.

   Still open: what "too large" means (a dimension threshold, a time budget, or both; to be
   settled from the phase 7 measurements), whether skipping for size should warn or error, the
   name of the opt-out, and whether the completed basis should also be offered to the
   reference method as an input fix.
3. **Scope. Direction: the same geometries as the reference** — threefolds and n-folds,
   hypersurfaces and complete intersections. "Complete intersection" here means what cygv
   already supports through `nefpart`: a complete intersection in a toric variety given by a
   nef partition of the divisors. The configuration-matrix CICYs in products of projective
   spaces are a special case of that, not a separate input format. The C covers only threefold
   hypersurfaces, so the other three combinations are new work, not a port:
   - *Complete intersections*: `q0` gets one row per part of the nef partition, as in
     `compute_omega`. The period step becomes $\prod_j (C\cdot q_{0,j})! / \prod_r (C\cdot
     q_r)!$ with the harmonic sums running over the `q0` rows, and the "negative anticanonical
     degree" check is per row. Candidates, instanton, extraction are unchanged.
   - *N-folds*: there is no grading-vector contraction. Each reference surface $k$ has its own
     form $K^{(k)}_{ab} = \kappa_{kab}$ and its own series $I_k = c_0^{-1} S_2^{(k)} - \sum_a
     \alpha_a V^{(k)}_a$, with $I_k[C] = A_{C,k} = \sum_{m \mid C} GV_{C/m,k}/m^2$ after
     subtracting lower curves (exponent 2 and no factor of $\deg C$, against 3 and $\deg C$
     for threefolds). $\exp(C\cdot\alpha)$ does not depend on $k$, so each curve still costs
     one exponential, emitted with a vector of scalars. A two-negative cone can be skipped only
     if its form vanishes for every $k$.
   - This suggests one design: a list of $M$ bilinear forms giving an $M$-component $I$
     ($M = 1$ with $K = \sum_t w_t\kappa_t$ for threefolds, $M$ = number of reference surfaces
     for n-folds), plus a multicover rule. **It has to be decided before `table` is written**,
     like the key width, because it changes the value layout of every $I$ table. Memory for $I$
     scales with $M$.
   - The n-fold formulas above are derived from the reference code (`compute_inst`,
     `invert_series`), not taken from the C, and are unverified until they pass differential
     tests against the reference on `examples/fourfold.yaml` and the `fourfold` bench model.

   **`min_points`: dropped for the new method (decided).** It never enumerates the semigroup,
   so there is no count to grow towards. A request with `min_points` and the new method is an
   error; the reference method keeps supporting it.

   **`target_points`: in the main plan (decided), for both methods.** #100 implements it
   ("lightcone") and shows it is straightforward: the set $D(P) = \{e : P - e \in
   \text{cone}\}$ of curves below a target is closed under taking summands, so truncating to it
   is valid for the same reason degree truncation and face restriction are, and membership is
   $H e \le H P$ with $H$ the Mori facets (several targets: the union).
   - `modular`: filter at the C's three points (candidates, curves applied in the extraction,
     output). Filtering the terms inside each exponential as well would save more but costs a
     membership test in the scatter loop; unmeasured, left as a tuning option. The truncation
     test still goes through one function, so that option stays cheap to add.
   - `reference`: keeps `Semigroup::with_target_points`, which implements the same set for the
     generated semigroup; under saturation the two agree.
   - Interface: the existing `target_points` field and argument keep their name and meaning;
     `max_deg` becomes optional when targets are given (default: the largest target degree).
   - Needs the facets, so with `modular` it depends on cone data (phase 6) or native facets
     (phase 7a). Not offloaded to a future GPU backend (the C keeps it on the CPU).

   **Phases of the ambient toric variety, including vex fans (#97): in the main plan, core and
   fan preprocessing (decided).**
   - *Core*: an optional `strata` input, a list of (three divisor indices, rational class
     $d_S \in \mathbb{Q}^{h^{1,1}}$), one per vex 3-cone. `modular` adds the $m < 0$ branch to
     the period step, the extra $K_T$ per stratum, and the larger factorial tables. Small and
     local, the same shape as in the C. Threefold hypersurfaces only; other combinations with
     `strata` are an error. The certifier (6.5) rejects `strata` until its bounds are worked
     out, as the C's does.
   - *Also in `reference` (decided, 6.11)*: the same `strata` input. In the reference pipeline
     the term enters each instanton series directly: with $P_t$ the series holding
     $U_0\, d_{S,t}$ at every curve $C$ on a vex stratum, $\mathrm{inst}_t \mathrel{+}= c_0^{-1}
     P_t$. Contracting with the grading vector gives back the C's $S_2$ term $U_0\,(w\cdot
     d_S)$, and since the invariants do not depend on the grading, the per-$t$ split is forced.
     This is my derivation from the C; the review checked it against `gv.c` and the reference's
     $\beta = c_0^{-1} c_2$, and it is still to be validated against `refs_vex`. What else the
     reference needs, matching the C:
     - a bucket for curves with three negative intersections next to `group_by_neg_int`'s
       0/1/2 (`src/fundamental_period.rs:18`), which today drops them;
     - *every* curve with $m < 0$ kept out of `compute_c_1neg` / `compute_c_2neg` (some have
       two negatives and land in `neg2` today, where `as u32` wraps), contributing only through
       the stratum term when $S$ matches and nothing otherwise;
     - $m < 0$ with fewer than two negative intersections is an error (`gv.c` does the same).
   - *Derivation from a fan: in the main plan (phase 7c, decided).* Deriving `strata`, the
     default Mori rays and grading, and the `kappa` check from the fan, i.e. the content of
     `cgv_phase.py`, natively and without normaliz. Its only normaliz uses are two small
     convex hulls, which the double description of phase 7a covers. Callers can still supply
     `strata` directly.
   - In `modular`, a stratum whose rational class has a denominator divisible by one of the
     primes makes that prime unusable: skip to the next prime instead of aborting as the C
     does. (With ~62-bit primes this essentially never happens, but it must not be fatal.)
   - Validation: the C oracle and the 16 references in `cgv/tests/refs_vex` (6.6: sync `cgv/`).
   - The formula is taken from the C and the cited papers; we have not derived it ourselves.

4. **Key width.** The key is a lossless encoding, not a hash: coordinate $t$ occupies a signed
   bit field, $\text{key} = \sum_t v_t\, 2^{bt}$ in two's complement, which is what makes
   $\text{key}(a+b) = \text{key}(a)+\text{key}(b)$ and lets the curve be unpacked again (the
   period step, the gcd in the multicover and the output all unpack). A separate hash of the
   key picks the table slot. Consequently a coordinate that does not fit its field is a hard
   failure, not a rare collision, and the C checks a bound before starting: every cone point
   of degree ≤ `max_deg` has $\lvert x_t\rvert \le$ `max_deg` $\cdot \max_i \lvert g_{it}\rvert
   / \deg g_i$ over the generators, and it aborts if that exceeds the field. The C gives every
   coordinate the same $b = \lfloor 128/h^{1,1} \rfloor$ bits (capped at 40): 11 bits at
   $h^{1,1} = 11$, 6 at 20, 4 at 30.

   **Direction (decided): a `Key` trait, with the representation chosen per run.** The tables
   are generic over a small trait (add, compare, hash, unpack, build from coordinates), which
   Rust monomorphizes, so there is no runtime cost and the entry point dispatches once, the way
   it will for the lane count. Planned implementations:
   - *Packed integer with per-coordinate field widths*, `u64` or `u128`. Each field is sized
     from its own bound instead of a uniform width; the bound is already computed per
     coordinate. Used whenever the fields fit.
   - *Inline array* `[i16; N]` (possibly `[i8; N]`) for inputs that do not fit. No carries
     between fields and unpacking is free, at the price of larger table entries.
   - Not a heap-allocated `Vec` per key: an allocation and a pointer chase per key, as in the
     current `HashSet<DVector<i32>>`.

   Bits needed with per-coordinate fields, from the C's bound at the largest benchmark degrees:
   `h11_8` at 20: 51; `h11_9` at 18: 56; `h11_10` at 10: 55; `h11_11` at 16: 70. So three of
   the four fit a 64-bit key, which makes entries smaller than the C's in a stage limited by
   memory latency. The largest coordinate bound in those runs is 60, so `i8` would also hold
   them, with little headroom; the bound grows linearly with `max_deg`, so sparse gradings run
   to degrees in the hundreds probably need `i16` (not computed for `h11_9_plike`).

   Unmeasured, to settle in phase 5: the cost of 32-byte keys against 16 (estimate: 10–30% on
   extraction), whether `u64` keys give a real gain, and whether an `i8` variant is worth
   having next to `i16`. Phase 1 implements the trait and the packed `u128`; the others follow
   once there is a benchmark to compare them on.
5. **Exactness. Direction (decided): a tunable number of check primes, default 1, minimum 1.** The primes
   in a run split into those needed to hold the largest $\lvert GV\rvert$ and extra ones that
   only verify the lift: the result is accepted when the lift from the needed primes agrees
   with every check prime's residue, for every curve. The parameter (`check_primes` or similar)
   sets how many extra ones there are:
   - *1 (default)*: what the C does. A wrong lift passes only if it happens to agree modulo an
     unrelated ~62-bit prime.
   - *0*: **not allowed**, unless phase 5 shows that dropping the check lane is significantly
     faster. Without a check the number of primes rests entirely on the probe's size estimate,
     and a low estimate gives a wrong result (reduced modulo the product of the primes) with
     nothing to detect it. The measurement below suggests the saving is around 20%, which does
     not justify that.
   - *n > 1*: each extra prime lowers the chance of an undetected wrong lift by a further
     factor of about $2^{-62}$. This is very strong evidence, but it is still not a proof: no
     number of check primes bounds $\lvert GV\rvert$.

   A proof needs the certifier, which is **planned** (phase 8). It is a separate switch, not
   a value of this parameter. How it works, following the C's `CGV_MAJ` build and
   `cgv/tools/certify.py`:
   - The modular result is exact modulo $M$, the product of the primes, so the only possible
     error is $\lvert GV\rvert \ge M/2$. A proof therefore needs an upper bound on every
     $\lvert GV\rvert$.
   - The same pipeline is run once more over non-negative reals rounded upward, with
     subtraction replaced by addition and negation by the identity, and division applied only
     to exact positive quantities. Every value is then an upper bound on the absolute value of
     the true one.
   - To stop the bounds compounding, the candidate GVs from the modular run are fed in: what a
     curve passes on to higher degrees is the value the candidates imply, not its bound. By
     induction on degree, if every bound is below $M/2$, every candidate is the true value.
   - If some bound is too large it reports how many primes are needed, and the modular run is
     repeated with at least that many.
   - It proves the lift was not truncated. It does not check the input or the code, and it
     assumes the true invariants are integers.
   - Cost per the C's README: about one extra CPU run, with more memory.

   Structurally it needs the hooks in 4.1.5 (alpha along charge rows, a per-point source
   override in the extraction, negation in the coefficient trait); those are designed in from
   phases 2–3 so phase 8 does not force a refactor. It rejects `strata`, as the C's does.

   To investigate when we get there: which number type provides the rigorous upward rounding.
   The C uses `long double` with the hardware rounding mode. `rug::Float` has directed rounding
   but `rug` becomes optional (6.10), so either the certifier is a `rug`-only feature or it
   uses another route (a pure-Rust type with directed rounding, or `f64` with explicit outward
   rounding of each operation). The C also re-expresses $\alpha$ along the rows of the charge
   matrix to keep the bounds tight, and extends to n-folds and complete intersections only by
   analogy, which needs checking.

   Implementation notes: lanes per pass = needed + check, capped at 4 as in the C, with further
   passes when more are required; a failed check triggers another pass exactly as now. The
   parameter is one more argument through every dispatch site (`compute_gvgw_strings`,
   `io::Input`, CLI, Python).

   Cost of a lane, measured with the C (`-l`, 8 threads, max_deg 18, best of two): `h11_9`
   1.10 / 1.33 / 2.26 s and `h11_8` 1.34 / 1.68 / 2.93 s for 2 / 3 / 4 lanes, with peak memory
   0.13 / 0.16 / 0.17 GB and 0.11 / 0.14 / 0.15 GB. So a third lane costs about 20–25% in time
   and memory. The fourth costs about 40% more once the C's slow big-number CRT path (0.4–0.6 s
   of those runs, taken above three primes) is set aside; `rug::Integer` should remove most of
   that. Going from two lanes to one could not be measured (the C's normal build has at least
   two), so the ~20% for dropping the check lane is an extrapolation. It also shows what each
   *additional* check prime costs.
6. **The PR's glue on this branch. Done: dropped, `cgv/` kept as an oracle.** `build.rs`, the
   `cgv` cargo feature and `cc` build dependency, `src/cgv.rs`, `_cgv_executable` in
   `src/python.rs`, `python/cygv/_cgv_run.py`, `backend="cgv"` in `python/cygv/hkty.py`,
   `python/tests/test_cgv.py` and `.github/workflows/cgv-windows.yml` are removed; those files
   are back to their `origin/main` state. The only differences from `main` outside `cgv/` are
   three excludes that keep the directory out of the linters and out of the published crate
   (`.pre-commit-config.yaml`, `[tool.ruff]`, `exclude` in `Cargo.toml`). In the working tree,
   not committed.

   `cgv/` is built with its own Makefile (`make -C cgv cgv`) and driven through
   `cgv/tools/cgv_run.py`, which needs `numpy`, and `normaliz` for the cone path. It is removed
   in phase 9. Synced to the version in #100 (d794626, phase 0d), so the oracle covers
   lightcones and vex phases and carries `tests/refs_vex`. Its input format is "cgv 2" and it
   always needs normaliz (development machines only). After the sync, `tests/regress.py` (38
   frozen cygv outputs) and `tests/vex_regress.py` (16 phase references, with lightcone
   checks) both pass here. What replaces it as a check after phase 9 is section 5.2:
   differential tests against the reference pipeline, plus a selection of `cgv/tests/refs` and
   `refs_vex` converted into Rust fixtures before `cgv/` goes.
7. **Naming, default and layout (decided).**
   - *Names*: the existing pipeline is `reference`, the new one is `modular`. Not "legacy": the
     old pipeline is kept on purpose as the faithful HKTY implementation, gets the same new
     features (6.11), and remains the only route for `min_points` and anything `modular`
     cannot represent.
   - *Default*: a method option with three values. `auto` (default) uses `modular` when the
     request is supported and `reference` otherwise; `modular` forces it and makes an
     unsupported request an error; `reference` forces the old pipeline. Existing calls keep
     working and most get the speed-up.
   - *Layout*: the existing stage modules move under `src/reference/`. This breaks the public
     Rust paths (`cygv::semigroup`, `cygv::polynomial`, ...), which is accepted: nobody uses
     the crate directly yet and it is pre-1.0. The Python API and the CLI are unaffected.
   - Consequences to put in the release notes: results can change on incomplete generator
     lists (the default now saturates, 6.2), and the Rust module paths moved.

   Still open: what `auto` does when `prec` is given. `modular` is exact, so either fall back
   to `reference` (reproduces today's floating-point output exactly) or run `modular` and
   format the exact result as a float (faster). Also whether the crate root keeps re-exporting
   `Semigroup`, `Polynomial`, `PolynomialProperties` and `PolynomialCoeff`, or those are
   reached through `cygv::reference`.
8. **Features of the C that are not ported now (decided): future options, not ruled out.** The
   GPU extraction stage, the distributed mode (one run split by primes across machines), the
   low-memory switch and the C's tuning knobs are not part of this work, but the Rust code is
   structured so that adding any of them later is an addition, not a refactor. The resulting
   constraints are in section 4.1 and apply from phase 1. PR #95 (the GPU variant through
   `build.rs`, `device=` in Python) is the reference for the GPU interface when it comes.

9. **Dropping floating-point computation (planned).** Remove `prec` and the `rug::Float`
   variants everywhere, including the reference pipeline, so every result is an exact integer
   (GV) or rational (GW).
   - What it removes: half of the dispatch (the four `*_float_*` wrappers and the `Float` arm
     of `compute_gvgw_strings`), `prec` in `io::Input`, the CLI and Python, the `gv-float` /
     `gw-float` bench variants, and the Python dependency on `mpmath`. It also settles the open
     question in 6.7 about `auto` with `prec`.
   - The reference pipeline stays generic over `PolynomialCoeff`. That genericity is what 6.10
     uses to swap the exact number type, so it is not collapsed to one concrete type. What can
     go is what only floats needed: `zero_cutoff` as a tolerance, the rounding check in
     `invert_series`, `RoundMut`, `AbsMut`, the float comparisons in the trait bounds.
   - What is lost: bounded-precision arithmetic as a way to cap memory and time when exact
     rationals grow large. `modular` makes that unnecessary where it applies, but requests
     that stay on `reference` (`min_points`, `target_points`, anything `modular` rejects) have
     no float option any more.
   - To check *before phase 0c starts*: whether CYTools or other callers pass `prec`. And a
     deprecation path for Python: for one release `compute_gv` / `compute_gw` keep accepting
     `prec`, emit a `DeprecationWarning` and return the exact values, then the keyword goes.
   - Done as its own commit right after phase 0b (phase 0c), before the `Method` argument is
     added, so the dispatch is rewritten once.

10. **`rug` becomes an optional feature (planned).** It was there for speed, and speed now
   comes from `modular`. The default build uses a pure-Rust big integer / rational crate, so
   it needs no GMP or MPFR and no C toolchain.
   - *Shape*: the pure-Rust crate is always a dependency. `modular` uses it directly, and only
     at the edges (the CRT lift and GW rationals), so its choice barely affects `modular`'s
     speed. The reference pipeline uses it as its default coefficient type through
     `PolynomialCoeff`; the `rug` feature adds the `rug::Rational` implementation and selects
     it, for anyone who wants the reference pipeline at today's speed.
   - *Work in the reference pipeline*: `PolynomialCoeff` is written against `rug`'s in-place
     traits (`rug::Assign`, the `*Assign` operators with primitive right-hand sides,
     `RecipMut`). The pure-Rust type needs a thin wrapper implementing those, and `Assign`
     needs to become our own trait or be re-exported conditionally. The typed entry points
     return `rug::Integer` / `rug::Rational` today and need backend-neutral result types
     (the string dispatch used by the CLI and Python is unaffected).
   - *What it buys*: no from-source GMP/MPFR build, which is most of cold CI time and the
     reason for the `GMP_MPFR_SYS_CACHE` handling, the `m4`/`make` requirement and the `CC`
     consistency rule. Possibly more on Windows, where the MSYS2/MINGW64 route and its
     workarounds (`--no-build-isolation`, `PYO3_USE_RAW_DYLIB=0`) may no longer be needed for
     the default build; I believe GMP is why that route exists but have not confirmed it. The
     memory bench no longer needs the GMP allocator hooks without `rug`.
   - *What it costs*: the reference pipeline gets slower in the default build. By how much is
     unmeasured and depends on the crate; it matters for the requests that only `reference`
     can serve. For that reason the switch waits until `modular` and `auto` are in place.
   - *Certifier* (6.5, planned): the plan assumed `rug::Float` with directed rounding. Without
     `rug` it needs another route to rigorous upward rounding, or it becomes a `rug`-only
     feature.

   Still open:
   - *Which crate.* Candidates found on crates.io: `num-bigint` + `num-rational` (MIT/Apache,
     by far the most used), `malachite` (algorithms derived from GMP/FLINT, so likely the
     fastest; LGPL-3.0-only, whose fit with this crate's GPL-3.0-or-later should be confirmed),
     `dashu` (MIT/Apache). To be chosen by running the reference benches with each.
   - *What the wheels ship*: built with or without `rug`.
   - Whether the reference benches keep a `rug` variant as the performance baseline.

11. **Feature parity in the reference pipeline (decided).** Since the reference pipeline stays,
   it is the main thing `modular` is compared against, in CI, long after the C oracle is gone.
   So every new capability lands on both sides unless it is inherently tied to the modular
   method:
   - *Vex phases*: `strata` in the reference period/instanton step (6.3), phase 4b. Done before
     the modular version, so phase 6c is checked against two independent implementations.
   - *Target points*: already in the reference (`Semigroup::with_target_points`).
   - *Saturation and the Mori cone by any description* (6.2, 6.1): the reference needs the
     lattice points of the Mori cone up to `max_deg`, i.e. the saturated semigroup, from rays
     or inequalities. The 7b machinery provides them, directly (polytope enumeration) or from
     the Mori cone's Hilbert basis (Pottier); phase 7d. If 7b is dropped, the reference
     saturates only when that basis is supplied.
   - *Fan to strata* (7c) and the cone computations (7a, 7b) are shared modules used by both.
   - *Not mirrored*: check primes and the certifier (the reference is exact), the lane probe,
     the future GPU and distributed modes. `min_points` and floats go the other way: dropped
     from `modular` and from both respectively.
   With this, `auto` mostly chooses between two implementations of the same feature set, and
   `reference` is the fallback only for `min_points` and inputs `modular` cannot represent.

## 7. Notes and risks

- Porting detail: the C computes dot products and $K_{ab}$ in `int` (`gv.c`, `dot`). Use `i64`
  (or check bounds) in Rust.
- The plan was reviewed by a separate agent on 2026-10-10; its findings are folded into 4.1.5,
  5.1, 5.2, 6.2, 6.3, 6.5, 6.9, the phase list and this section.

- **Latent bug in the reference pipeline, worth fixing on `main` independently of this work.**
  `compute_c_1neg` and `compute_c_2neg` (`src/fundamental_period.rs`) convert $q_0 \cdot C$ with
  `as u32`. In a vex phase that value can be negative, so it wraps to a huge factorial argument
  instead of failing. The same holds for any inconsistent input. Phase 4b gives it defined
  behaviour; a minimal fix (an error) could also go to `main` on its own.

- **Faces of the Mori cone must keep working.** Passing only the generators of a
  lower-dimensional face computes the invariants along that face, because a face $F$ satisfies
  $a + b \in F \Rightarrow a, b \in F$, so dropping every monomial outside $F$ commutes with
  products, reciprocals and exponentials. The new algorithm inherits this as long as the
  candidates are the lattice points of $K_T \cap F$. Measured on a face of the
  $\mathbb{P}(1,1,2,2,2)[8]$ model and three facets of `h11_9` (max_deg 10): the reference and
  the C enumeration path agree everywhere; the C *cone* path does not (it returns the full-cone
  result on the toy model and fails with "CRT did not stabilize" on two of the three facets).
  The likely cause is in `cgv_run.py`, which keeps normaliz's support hyperplanes for the face
  but not the equations of its linear span, so the $K_T$ are built in the ambient space; not
  verified further. Requirements for us: the cone description used in phases 6 and 7 carries
  equations as well as inequalities (a face is not full-dimensional), saturation is taken
  within the face's span, and face inputs are part of the differential tests.

- The grading-vector contraction assumes the $\mathrm{inst}_t$ come from one prepotential. The
  PR notes that the $h^{1,1}=2$ smoke-test input of `python/tests/test_hkty.py` is not a
  consistent geometry and gives grading-dependent results with the new method, so the
  existing toy fixtures cannot all be reused for differential tests.
- The C uses `int` slot indices and `h11, ndiv <= 64` stack arrays; worth lifting rather than
  porting.
- Peak memory is per-thread working tables plus $I$; the C's README quotes 3–6 GB at degree
  26–28 for $h^{1,1}$ = 10–11. `benches/memory.rs` only counts GMP and Rust allocations, which
  is fine here since everything will be Rust-allocated.
- What happens to PR #94 itself (close, or merge first and replace later) is your call and
  affects how the author is credited.
