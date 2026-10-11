# Plan: a Rust implementation of the fast GV algorithm from PR #94

Status: **investigation done, phases 0 (0–0d), 1, 2, 3, 3b, 4, 4b, 6–6c, 5, 5b, 7a, 7b and 7d done.** Progress is
tracked by the check boxes in section 5. Where sections 1–5 mention `rug` for the new code, read
"big integers": phase 5b replaced `rug` with `malachite` everywhere (6.10). Branch
`gv-modular-rust` is PR #94 (`natemacfadden:cgv-cc-cpu`, 11 commits) rebased onto `origin/main`
(6fe293a, v0.3.0), plus our commits. The PR branch itself (`cgv-cc-cpu`) is untouched. The follow-up PRs #95, #97 and #100
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
- [x] **0b. Move the existing pipeline under `src/reference/`** (6.7), as a commit of its own
  with no functional change. Done before any new code so that the diff stays a pure move.
  - [x] `git mv` the modules and fix the `use` paths in `src`, `tests/`, `benches/` and the
    CLI.
  - [x] Split `src/hkty.rs` into `reference.rs` (`run_hkty`) and `dispatch.rs` (the string
    dispatch and the eight typed wrappers).
  - [x] Update CLAUDE.md.
  - The crate root still re-exports `Semigroup`, `Polynomial`, `PolynomialProperties`,
    `PolynomialCoeff`, `run_hkty` and everything in `dispatch`; whether it keeps doing so is the
    open question of 6.7.
- [x] **0c. Drop floating-point computation** (6.9), as its own commit.
  - [x] Check whether CYTools or other callers pass `prec` (before starting). CYTools
    (`origin/main`, bf598ab) does not, and no public code on GitHub does.
  - [x] Remove `prec` and the `Float` variants from the Rust code, `io::Input` and the CLI.
    `zero_cutoff`, `RoundMut`, `AbsMut` and the float trait bounds went too; `IsInteger`
    replaces the rounding check, and `PolynomialCoeff` now requires `Default` (zero).
    Outputs checked identical to the previous commit on both examples, `h11_9` and `h11_11`.
  - [x] Remove the `gv-float` / `gw-float` bench variants (the remaining ids keep their
    `-rational` suffix so that saved baselines still match).
  - [x] Python: `prec` removed outright (no deprecation release, since nobody passes it),
    `mpmath` dropped.
- [x] **0d. Sync the oracle** (6.6): `cgv/` is the #100 version, without the integration glue.
- [x] **1. Primitives.**
  - [x] `field`, with unit tests of its operations against the big-number crate. `Field`
    (Montgomery, as in the C), `Lanes<NL>` and `LaneField<NL>`, the prime sequence (the C's,
    so residues compare lane for lane), and the `Ring` trait of 4.1.5, which the stages take
    as a context because a residue does not know its modulus. Table construction is not in
    the trait yet; it comes with the period tables in phase 2.
  - [x] `key`, with pack/unpack round trips and overflow detection. The `Key` and
    `KeyLayout` traits, `PackedLayout` (`u128`, per-coordinate widths of at most 126 bits in
    total, which leaves a sentinel key free), and `coordinate_bounds`, the C's bound in exact
    integer arithmetic and per coordinate.
  - [x] `table`: `SlotMap`, `SparsePoly`, `DegreeHeap`, `LevelTable`. Slots and terms carry
    a run-time number of values, which settles the value layout question of 6.3. The
    per-degree workspace built from level tables belongs to the extraction (phase 3).
- [x] **2. Period and instanton, enumeration path only.**
  - [x] Candidates by full enumeration: `generated_semigroup` (the C's old `reduce_generators`
    and `enumerate_semigroup` fused into one pass in order of degree, dropping redundant
    generators on the way) and `candidates_by_enumeration` (at most two negatives).
    Single-threaded on `h11_9` in release: 1.4 s at degree 14 and 15.7 s at degree 18 for
    1.6 M and 10 M points, against 2.9 s and 21.0 s for the C's old enumeration (section 2).
  - [x] `period`: `Tables` (built from `Ring::inverses`, the table construction hook of
    4.1.5, whose default takes one inverse) and $c_0$, $c_{1,a}$, $S_2^{(k)}$ per candidate.
    Written for several anticanonical rows and $M$ bilinear forms from the start (`Geometry`
    in `modular.rs`), since that was cheaper than retrofitting it in 3b. $m < 0$ is an error
    with fewer than two negatives and contributes nothing otherwise, until 6c. 0.1 s at
    degree 14.
  - [x] `instanton`: `mul_into` / `mul_scalar`, `recip`, $\alpha$ as one series with
    $h^{1,1}$ values per term (that is the union of the supports the extraction takes), $V$
    and $I$ with $M$ values per term. 25 ms at degree 14. (Phase 5 takes
    $\sum_a \alpha_a V_a$ by pairs of columns of $\alpha$ and halves it by symmetry.)
  - [x] $I$, `c0` and `c1` tested against the reference stages mod p, and $S_2$ against the
    reference's `c2` contracted with the forms (`modular/reference_tests.rs`), on two lanes,
    at every element of the semigroup: the quintic, P(1,1,2,2,2)[8], the toy threefold
    (one and two negatives), `h11_9` at degree 5, and the fourfold of `examples/` both as
    the hypersurface it describes and as the CICY of the bench model (14 forms, two
    anticanonical rows). Checked to fail when a sign or the 1/2 in $V$ is broken.
- [x] **3. Extraction, multicover, CRT**, single-threaded and simple (`apply_curve` only), for
  threefold hypersurfaces.
  - [x] Provisional choice of the pure-Rust big-number crate, confined to `crt` and
    `multicover` (6.10, phase 5b): `num-bigint` (MIT/Apache, already in `Cargo.lock` through
    `nalgebra`; `malachite` since phase 5b). Only `crt` needs it: the multicover runs on
    residues, and the stability check works on Garner's mixed-radix digits in word-sized
    arithmetic, so big integers are only built for the final values.
  - [x] `extraction`: `extract` takes the layout, $\alpha$, $I$, the inverses and `max_deg`
    and returns the coefficients $n_C$ (what is left of $I$ at each curve once the lower
    curves are subtracted), with as many values per curve as $I$. Plain function for now;
    the trait of 4.1.1 can come with the second implementation (phase 5). The division by
    the degree belongs to `multicover`, so the extraction does not depend on the CY kind.
  - [x] `multicover`, with the threefold rule (divide by the degree, exponent 3) and the
    n-fold one (exponent 2), by gcd level, adding the multiples of curves with nonzero GV
    that have no coefficient. `Geometry` now records the CY kind.
  - [x] `crt` (`Residues`: per-prime columns, merging passes that find different curves) and
    `compute_gv` in `modular.rs`, the pass loop, with `Options { check_primes, lanes,
    min_primes, max_passes }`. The test knob is `lanes: 1` with `min_primes: 5` (two lanes and
    six primes since phase 5).
  - [x] End-to-end differential tests against `run_hkty` (`modular/reference_tests.rs`):
    the quintic at degree 12 (two passes with the default options), P(1,1,2,2,2)[8] at 12,
    the four `benches/data` models at degree 5–6, several passes with one lane (two since
    phase 5) and with three lanes and two check primes. Checked to fail when the exponent or the division by
    the degree is broken in the multicover, or the degree weight in the Euler recurrence.
    Also compared by hand with the C oracle (normaliz cone data, 8 threads): identical GVs
    for `h11_9` at degree 18 (77458), `h11_8` at 16 (42300) and `h11_10` at 12 (15008).
    Single-threaded release timings on `h11_9`: at degree 14, enumeration 1.7 s, then
    period 0.10 s, instanton 0.04 s, extraction 0.15 s, multicover 3 ms; at degree 18,
    enumeration 17 s, then 0.5 / 0.3 / 2.9 / 0.02 s.
- [x] **3b. Complete intersections and n-folds** (6.3). Done before tuning so that phase 5
  optimises the final data layout.
  - [x] Multi-row `q0`. Done in phase 2: `period` is the only stage that reads it.
  - [x] The $M$-component $I$, consumed by the extraction (phase 3).
  - [x] The exponent-2 multicover (phase 3).
  - [x] Check which geometry `examples/fourfold.yaml` and the `fourfold` bench model each
    are (5.2). They are the same data (charges, intersection numbers, generators and grading)
    except that the bench model has the nef partition `[0..3], [4..11]`. Read as the
    hypersurface it says it is, the example gives the reference non-integer GV invariants
    at degree 6, so its intersection numbers are presumably those of the complete
    intersection, and the example should get the nef partition and a corrected comment (its
    CLI test uses `min_points: 100`, where the reference happens to succeed).
  - [x] Differential tests against the reference on the fourfold CICY, end to end at degree
    10 (183 invariants, 30 of them on curves with gcd above one). The stage-level ones exist
    since phase 2.
  - [x] Fix `examples/fourfold.yaml` as above, and its expected output in `tests/cli.rs` if
    that changes. It has the nef partition now, and its comment says what it is (as a
    hypersurface the same data would even be a fivefold: 12 rays and 6 charges). Its output
    changed (56 invariants instead of 79), and the CLI test pins the first one and the count.
    The stage-level test still also runs the data without the nef partition, which covers an
    n-fold with a single row.
- [x] **4. Front end.**
  - [x] The `Method { Auto, Modular, Reference }` argument through `dispatch`, `io::Input`,
    the CLI and Python. `auto` falls back to `reference` for requests `modular` does not
    support; an explicit `modular` returns an error for them. `Method` is in `src/lib.rs`;
    every entry point takes it last and returns `Result<_, ComputeError>` (only `modular`'s
    errors get there; `run_hkty` still panics). Unsupported: `min_points`, `target_points`
    (until 6b), and no truncation at all; `auto` also falls back when the curve classes do not fit the
    keys (`KeyError::TooWide` / `BoundTooLarge`). `method` in the YAML (checked eagerly for
    `modular`), `--method` / `-m` in the CLI, `method=` in Python (`None` means the Rust
    default). Python exceptions from the subprocess now keep their type (`ValueError`).
  - [x] GW. Derived from the lifted GV invariants (`multicover::gw_from_gv`, $\sum_{m \mid C}
    \mathrm{GV}_{C/m}/m^e$ in `num-rational`, also provisional; `malachite` since 5b), since GW
    invariants are not integers; the reference's GW output is exactly those, nonzero only. Differential tests
    through `compute_gvgw_strings` (so the formatting of fractions is compared too) on the
    quintic, P11222, the fourfold CICY and `h11_8`.
  - [x] `n_threads` exposed in the Python API.
  - [x] Error messages for non-integer results that name the likely cause (an inconsistent
    input, or a cone that is not a face) instead of only "CRT did not stabilize". One text
    (`NON_INTEGER_CAUSES` in `src/lib.rs`) for `modular::Error::NotStable` and the reference's
    `NonIntegerGVError`.
  - [x] `auto` made the default, once the differential tests of phases 3 and 3b pass.
    Measured first, release CLI, `modular` on one thread against `reference` on eight, outputs
    byte-identical in every case: `h11_8` deg 12 / 14: 0.45 / 1.6 s against 11 / 62 s;
    `h11_9` deg 10 / 14: 0.17 / 2.1 s against 4.0 / 147 s; `h11_10` deg 8 / 10: 0.23 / 1.1 s
    against 5.3 / 22 s; `h11_11` deg 10 / 12: 0.76 / 3.7 s against 6.1 / 32 s; fourfold
    deg 15: 0.12 s against 0.74 s; GW on `h11_9` deg 12: 0.62 s against 24 s. Consequence for
    the release notes: on inconsistent input with `max_deg` (such as the toy threefold the
    README used) the default now gives different numbers, since neither pipeline's output
    means anything there. The README examples now use P11222 instead. The benchmarks still
    measure `reference`, explicitly.
  - Moved to phase 6: `generators` optional in `io::Input` when `mori_hyperplanes` is given,
    since that field arrives with the cone data.
- [x] **4b. Vex phases in the reference pipeline** (6.3, 6.11).
  - [x] The `strata` input through the front end. `Stratum { divisors: [usize; 3], class:
    Vec<(i64, i64)> }` in `src/lib.rs`, checked by `check_strata` (threefold hypersurfaces only,
    distinct in-range divisors, one entry per Kähler parameter, positive denominators, no two
    strata on the same divisors) in `dispatch` (`ComputeError::InvalidStrata`), `io` and
    `compute_omega`. A new argument of `run_hkty` and of every entry point, after `intnums`;
    `strata:` in the YAML (`{divisors, class}` mappings, entries integers or `"p/q"`), `strata=`
    in Python (`(divisors, class)` pairs, entries anything `Fraction` takes). `modular` does not
    support strata until 6c, so `auto` falls back to the reference for them.
  - [x] The three-negative bucket and the extra term in the reference instanton step:
    `group_by_neg_int` returns `CurveGroups`, whose `vex` list holds the $m < 0$ curves on a
    stratum; `compute_c_vex` gives $U_0\, d_{S,t}$ into `FundamentalPeriod::vex[t]`, and
    `compute_instanton_data` adds $c_0^{-1}$ `vex[t]` to instanton correction $t$. The
    derivation of section 6.3 (the per-$t$ split) holds: see the tests below.
  - [x] $m < 0$ kept out of the 1/2-negative routines (this also fixes the `as u32` wrap of
    section 7), and the error for $m < 0$ with fewer than two negatives
    (`FundamentalPeriodError::NegativeAnticanonicalDegree`).
  - [x] Tests against `cgv/tests/refs_vex` and the C oracle, development only (scratch scripts,
    since the inputs need `cgv_phase.py` and normaliz): the reference, given the Hilbert basis
    of the Mori cone as generators and the strata from `cgv_phase.prepare`, reproduces all 16
    references at their stored degrees (20 for `liam_h2_*` and the pairs, 14 for `liam_h3_*`;
    for the pairs also the FRST side under the mapped grading), and the C oracle at higher
    degrees (`liam_h2_fan1` at 30, `liam_h3_fan7/8/9` at 22/20/20, `pair_h4_1` at 24). Without
    the strata every vex case fails with non-integer invariants. Committed: unit tests of the
    grouping and of $U_0$, `tests/vex.rs` (`liam_h2_fan1` with a few invariants that are those
    of the familiar fiber class, and `pair_h4_1`, vex against FRST through `T` at degree 10),
    and YAML and Python tests. These use the geometry of two references (charges, intersection
    numbers, derived generators and strata), not their stored tables; whether that needs the
    PR author's confirmation too is open (phase 9).
- [x] **6. Cone fast path through optional input.** Runs before phase 5, and fixes the
  interface that phase 7 plugs into.
  - [x] Cone data (6.1) accepted by the Rust API and YAML. The Mori cone may be described by
    spanning vectors or by inequalities, plus equations for faces (section 7). `ConeData` in
    the new `src/cone.rs` (the shared module of 4): `hyperplanes` and `equations` (columns,
    $h \cdot C \ge 0$ and $e \cdot C = 0$) and `cones`, a `CandidateCone` (divisors $T$,
    Hilbert basis as columns) for every $T$ of one or two divisors. One more positional
    argument, `cone_data: Option<ConeData>`, after `strata` at every dispatch site; YAML and
    Python fields `mori_hyperplanes`, `mori_equations`, `cones` (mappings / pairs of
    `divisors` and `hilbert_basis`). The reference ignores it. Spanning vectors stay the
    `generators`; nothing is derived from them natively until 7a. Checked by
    `check_cone_data`: a cone for every $T$ (a forgotten one would silently drop curves), no
    duplicates, every Hilbert basis element in its $K_T$ and every generator in the Mori cone.
    Completeness of a basis cannot be checked: wrong data silently drops curves (a face test
    with deliberately wrong bases gave no invariants at all). Elements of degree above
    `max_deg` may be left out, which keeps fixtures small. Cones of three divisors are rejected
    until 6c.
  - [x] `generators` optional in `io::Input` when `mori_hyperplanes` is given, with validation
    reworked accordingly (moved from phase 4). In practice: when the cones are given (the
    hyperplanes alone are not enough to run anything before 7b). Without generators the Rust
    matrix has no columns, and a request that needs the reference (method `reference`, or one
    `modular` does not support) is `ComputeError::MissingGenerators`, or `MissingField` in YAML.
    `auto` does not fall back on key overflow then. `h11` comes from the grading vector now.
  - [x] `candidates` enumerates each $K_T$ from its Hilbert basis, with the containment and
    vanishing-form skips. The candidate source is an interface, so phase 7b can plug in a
    native one (Hilbert bases or direct enumeration) without touching the rest.
    `CandidateSource { Generators, Cones }` gives the coordinate bounds (from the basis
    elements up to `max_deg`, as the C) and the candidates; `compute_gv`/`compute_gw` take it.
    Each used cone is enumerated with `generated_semigroup` and the points are merged. The
    vanishing-form test is $Q_s K^{(k)} Q_u = 0$ for every form $k$, so it covers n-folds and
    complete intersections too. Negative intersections for the containment skip only count
    elements of degree up to `max_deg`, which is tighter and still sound.
  - [x] ~~The Python wrapper may produce the cone data with `normaliz`.~~ Decided against
    (6.1).
  - [x] Tests: candidate sets with and without the skips, face inputs (5.2). Unit tests against
    the enumeration with Hilbert bases computed by brute force (P(1,1,2,2,2)[8], a toy
    $h^{1,1} = 3$ geometry under two gradings, threefold and n-fold, where the vanishing-form
    skip drops points), and the same check on `h11_9` with Normaliz's bases
    (`tests/data/h11_9_cones.yaml`, truncated at degree 10, 62 KB). `tests/cones.rs`: `h11_9`
    at degree 8, all methods and with or without generators; the saturation test above; a facet
    with its equation (`tests/data/h11_9_facet.yaml`, 34 KB) against the reference with the
    facet's generators; errors. `cone.rs`, YAML, CLI (`examples/cones.yaml`) and Python tests.
    The fixtures come from a scratch script around Normaliz 3.11.0 (the way `cgv_run.py` does
    it, but keeping the equations of a face, whose omission is the C driver's face bug of
    section 7).
  - Measured (release, one thread, outputs identical to the enumeration path): `h11_8` deg 12 /
    14 / 16: 0.09 / 0.33 / 1.17 s against 0.43 / 1.51 / 4.42 s; `h11_9` deg 10 / 14 / 16: 0.06
    / 0.30 / 0.96 s against 0.15 / 1.89 / 5.93 s; `h11_10` deg 8 / 10: 0.12 / 0.22 s against
    0.19 / 0.97 s; `h11_11` deg 10 / 12: 0.05 / 0.11 s against 0.66 / 3.62 s. The C with cone
    data took 0.15 s for `h11_9` at 14 on 8 threads (section 2). Normaliz itself took 2.5–7 s
    per geometry for all the cones.
- [x] **6b. Target points for `modular`** (6.3).
  - [x] The facet test at the three filter points (listing candidates, choosing the curves to
    apply in the extraction, the output). `modular::lightcone::Lightcone`, built by
    `CandidateSource::candidates` and applied in `candidates`, `extraction::extract` and after
    the multicover of each pass (so the CRT never sees the garbage values outside). Two tests:
    with cone data, $h \cdot C \le h \cdot P$ over the Mori cone's hyperplanes (targets off the
    span of a face have nothing below them); without, membership of $C$ and $P - C$ in the
    generated semigroup, which is then enumerated up to the targets' degree, i.e. exactly the
    reference's `with_target_points` set, so it needs no facets and no phase 7a.
    `CandidateSource::Cones` now holds the whole `ConeData` for the hyperplanes. GW is filtered
    too: `gw_from_gv` also writes the multiples of a curve, which leave the lightcones with
    incomplete sums. Without zero among the candidates (every target off the cone) the result
    is empty instead of failing to invert $c_0$.
  - [x] Optional `max_deg`, defaulting to the largest degree of the targets. `max_deg` and
    `target_points` may now be combined everywhere (YAML, CLI, Python), the degree lowering the
    targets' largest; only `min_points` still excludes the others. `compute_gv`/`compute_gw`
    take `targets` after `max_deg` and lower it themselves; the dispatch passes no degree as
    `u32::MAX`. The reference applies `max_deg` to its results, since its semigroup only takes
    the targets. Targets are checked up front (`check_target_points`: dimensions, positive
    degree, `ComputeError::InvalidTargetPoints`), since `run_hkty` panics on a nonpositive one.
  - [x] Differential tests against the full computation and the C's lightcone mode. Unit tests
    of both lightcone tests; `reference_tests::check_targets` (GV and GW, both methods through
    `compute_gvgw_strings`, with and without a lower `max_deg`, against the full run restricted
    to the lightcones) on the quintic, P11222, `h11_9` and the fourfold CICY; `tests/cones.rs`
    on `h11_9` with cones, with and without generators, and on the facet with targets on and
    off it; YAML and Python tests. The reference's `with_target_points` is slow on thousands of
    generators (87 s in debug on `h11_9`, against 1.2 s with those of degree up to the
    targets'), so the tests drop the generators above it. Against the C's lightcone mode
    (scratch script, Normaliz cone data, one thread): identical on `h11_9` (degree 10, one to
    three targets, and with `max_deg` 8), `h11_8` (16, three targets: 881 invariants, 0.52 s
    against the C's 0.46 s; with `max_deg` 13) and `h11_11` (14 and 12). The toy threefold's
    target-point smoke tests in Python now ask for the reference: `auto` routes them to
    `modular`, where that inconsistent geometry has no integer invariants.
- [x] **6c. Vex phases in `modular`** (6.3).
  - [x] The $m < 0$ period branch. `Geometry::with_strata` (threefold hypersurfaces only,
    `Error::InvalidStrata`) keeps each stratum's sorted divisors and the terms of $w \cdot d_S$;
    the period evaluator reduces them modulo the primes and adds $U_0\,(w \cdot d_S)$ to $S_2$
    at a curve of anticanonical degree $m < 0$ whose three negative intersections are a
    stratum's, as the C does. That is the reference's $c_0^{-1}$ `vex[t]` contracted with the
    grading, since $S_2$ only enters $I$ as $c_0^{-1} S_2$. The "larger factorial tables" of
    6.3 came for free: the tables were already sized by the intersections of the candidates.
  - [x] The extra $K_T$ per stratum. `ConeData` takes a three-divisor cone for every stratum,
    and `check_cone_data` (which now takes the strata) requires one per stratum and rejects any
    other. The containment skip needs no change: the stratum's cone contains those of its
    subsets. The candidates keep the curves on a stratum (`Geometry::carries_period`) on both
    paths; from the cones, points of a stratum's cone with three negatives that carry nothing
    are dropped. `modular_unsupported` lost its strata argument, so `auto` routes vex phases to
    `modular`.
  - [x] Tests against the reference (phase 4b), the C oracle and `cgv/tests/refs_vex`.
    `reference_tests::test_vex` on `liam_h2_fan1` and the vex side of `pair_h4_1`: the stages
    curve by curve (the reference's $S_2$ plus its stratum terms contracted with $w$, which
    catches a flipped sign of $U_0$), GV and GW end to end, target points, and cone data built
    by brute force (with the facets, now computed by a test helper) against the enumeration and
    the reference, also below a target on a stratum. `tests/vex.rs` compares both methods; a
    new example, `examples/vex.yaml` (cone data with the stratum's cone, no generators), is
    checked against the reference by `tests/cli.rs`, as is the error without that cone; YAML and
    Python tests. Development only (scratch scripts, the inputs need `cgv_phase.py` and
    Normaliz): all 16 references of `refs_vex` at their stored degrees, through the generators
    (1 s for all) and through Normaliz cone data with the stratum cones, plus the lightcones of
    three stored curves per reference as `vex_regress.py` does; the C oracle at higher degrees,
    identical on `liam_h2_fan1` at 50, `liam_h3_fan4/6/8` at 30, `pair_h4_1` at 30 and
    `pair_h5_1/2` at 28 (up to 5156 invariants, 0.1–3.2 s on one thread against the C's
    0.3–1.7 s on eight); the cone path against the generators path at 26–40 on all of them.
- [x] **5. Parallelism and tuning.** Target parity with the C *with cone data* (section 2).
  Measured on a 4-core machine with hyperthreading (8 threads), old enough that the tuning is
  provisional and worth redoing on newer hardware.
  - [x] Rayon over candidates / product chunks / curves per layer. Each run builds a private
    pool of `Options::threads` workers (`n_threads` from the dispatch). Candidates: the cones in
    parallel, the union sorted. Period: chunks of candidates, appended in order. Instanton:
    products in chunks of similar cost, each into a map of its own, added to a `ShardedMap`;
    $\sum_a \alpha_a V_a$ as products of pairs of columns $a \le b$ with some $K_{ab} \ne 0$,
    each pair of terms once (each $\alpha_a$ has few terms: 1.9 → 0.58 s on one thread for
    `h11_10` at 16, against the C's 0.79 s). Extraction: the curves of each layer in parallel,
    one worker per thread with its own level tables, with the C's super-layers (degrees
    $[d, d + \delta)$). Software prefetch in the scatter (x86-64), as the C does: about 10 %.
    The multicover and the lift stay serial (under 0.15 s everywhere measured).
  - [x] Sharded $I$ with per-shard locks: $4 \times$ threads shards by the high bits of the key
    hash, per-shard buffers of 1024 terms per worker, the rest folded per shard in parallel at
    the end of each layer, as the C does.
  - [x] The level-synchronous mode for few large curves (`extraction::Tuning`): each degree is
    finalized for every curve, then the scatters of every (curve, target level) run as tasks,
    largest first. It only pays for layers of fewer than 2 curves per thread (the C uses 4, but
    here 25 curves on 8 threads went faster curve by curve), with $\alpha$ having at least
    20000 terms up to the layer's top degree. It gives about 2× on such layers; the layers just
    above them, with a few dozen curves, still scale only 2–3×, on both sides.
  - [x] Lane probe (`choose_lanes`), the C driver's rule: two primes per pass up to degree 12,
    above that a probe at 3/4 of the degree extrapolates the largest bits per degree, plus the
    check primes. `Options::lanes` is now `None` by default, and 2–4 primes: one per pass was
    never efficient (it always needs a second pass for the check) and cost compile time, so
    `NL = 1` is gone and the test knob is two lanes with six primes. On `h11_8` at 20 this turns
    two passes of two primes into one of three (44 → 27 s on one thread).
  - [x] Bench scenarios: a `Pipeline` per scenario (reference, modular, modular with cone data;
    the ids of the modular ones end in `/modular` or `/modular-cones`), cone data for the four
    realistic models in `benches/data/cones/` (Normaliz, truncated at the highest degree of the
    modular scenarios, 570 KB, excluded from the crate). Quick: the fourfold at 15 (GV and GW),
    `h11_9` at 12 and `h11_9_plike` at 400 from the generators, the realistic models at 12–16
    with cones; heavy: `h11_9` at 14 and `h11_9_plike` at 500 from the generators, the
    realistic models at 16–20 with cones.
  - [x] Stress test (5.2): `tests/stress.rs`, ignored by default, 40 runs per case at 2, 3, 4
    and 8 threads with both extraction modes against one thread (`h11_9` with cones and from
    its generators, the fourfold), 24 s in release. `reference_tests::test_threads` forces both
    modes on 2 and 3 threads against the reference in every test run.
  - [x] Compile time and binary size measured (release, the crate alone, the CLI binary): before
    phase 5 51 s and 2.33 MB; with the parallel code and lanes 1–4 74 s and 3.19 MB; with
    lanes 2 only 49 s and 2.60 MB; as committed (lanes 2–4) 66 s and 3.00 MB. Each lane count
    costs about 8 s and 200 KB, the parallel code (rayon closures, monomorphized per lane
    count) about a third on top of what it replaced.
  - Results, through the CLI (YAML in and out) against the C (`cgv`, its own lane probe, plain
    text out), outputs identical in every case:

    | Case | cygv, 8 threads | C, 8 threads | cygv, 1 thread | C, 1 thread |
    |------|-----------------|--------------|----------------|-------------|
    | `h11_8`, deg 18 | 2.2 s | 2.0 s | 6.2 s | 5.7 s |
    | `h11_8`, deg 20 | 9.4 s | 8.0 s | 27.5 s | 24.0 s |
    | `h11_9`, deg 18 | 1.5 s | 1.2 s | 3.8 s | 3.3 s |
    | `h11_9`, deg 20 | 4.9 s | 4.4 s | 14.4 s | 13.3 s |
    | `h11_10`, deg 14 | 2.1 s | 1.7 s | 4.9 s | 4.5 s |
    | `h11_10`, deg 16 | 10.2 s | 9.5 s | 30.3 s | 27.3 s |
    | `h11_11`, deg 16 | 0.6 s | 0.4 s | 1.3 s | 1.0 s |
    | `h11_11`, deg 18 | 1.8 s | 1.4 s | 4.9 s | 3.9 s |

    So within 5–25 %. The extraction does exactly the C's work (7.23e7 scatter steps on
    `h11_11` at 18) about 15 % slower per step on one thread, and at parity layer by layer on
    eight; neither `target-cpu=native` (the C builds with `-march=native`), unchecked indexing
    in the level tables nor other prefetch distances or hints changed that measurably on this
    machine (`perf` is not available on it). The rest is the CLI's YAML. Also checked with the
    parallel code: the C's lightcone mode, the 16 `refs_vex` references with their lightcones,
    and the cone path against the generators path on the vex pairs at 30. Not done: the stats
    struct of 4.1.7, so the timings above come from temporary instrumentation; and the key
    sizes of 6.4 (whether `u64` keys, which three of the four models fit, pay off in a stage
    bound by memory latency), which belongs with the tuning on newer hardware.
- [x] **5b. Replace `rug` with a pure-Rust crate** (6.10). Planned as making `rug` optional;
  the benchmarks made it a replacement.
  - [x] The crate, chosen from benchmarks: the reference pipeline on `rug`, `num-rational`,
    `malachite` and `dashu` through one generic wrapper, every output identical (temporary
    harness, not committed). Single-threaded:

    | Scenario | `rug` | `num-rational` | `malachite` | `dashu` |
    |---|---|---|---|---|
    | `threefold`, deg 30 | 0.85 s | 16.6 s | 0.94 s | 1.20 s |
    | `h11_9`, deg 10 | 8.3 s | 15.6 s | 5.9 s | 7.9 s |
    | `h11_9_plike`, deg 200 | 12.8 s | 65.2 s | 9.3 s | 14.8 s |
    | `h11_10`, deg 8 | 6.8 s | 8.2 s | 5.4 s | 6.4 s |
    | `h11_11`, deg 10 | 7.5 s | 8.5 s | 6.6 s | 7.2 s |

    On eight threads `h11_11` at 12 took 31.7 / 29.8 / 31.7 s (`rug` / `malachite` /
    `dashu`), `h11_9` at 12 23.0 / 16.8 / 18.4 s with peak RSS 955 / 871 / 856 MB.
    `num-rational` is out (its `gcd` is binary, and every operation normalises), and
    `malachite` beats `dashu` almost everywhere. `malachite` stores one-limb numbers inline, so
    small values do not allocate, which is where it gains on `rug`; it loses when the numbers
    get huge (the toy `threefold`). The license fits (6.10).
  - [x] `PolynomialCoeff` implemented for it: `reference::Rational`, a newtype over
    `malachite::Rational`, with `Assign` now our own trait.
  - [x] Backend-neutral result types: with one crate left, the typed entry points return
    `malachite::Integer` / `malachite::Rational`, and `modular` uses the same types (`crt`,
    `gw_from_gv`), so `num-bigint` and `num-rational` are gone too. `malachite` is re-exported.
  - [x] ~~The `rug` feature~~: dropped (decided with you), since `malachite` is as fast.
  - [x] CI workflows and CLAUDE.md simplified: no GMP cache or `CC`, Windows on the stock
    MSVC toolchain (Rust and Python CI with `setup-python`, the wheels with cibuildwheel like
    the other platforms, 32-bit Windows skipped), so no MSYS2, `--no-build-isolation` or
    `PYO3_USE_RAW_DYLIB`; `Swatinem/rust-cache` with a `dev` / `release` key instead of the
    hand-rolled caches. Not run yet: CI only triggers on `main` and PRs to it.

    Before and after, release, eight threads, `cargo bench --bench memory` (peak RSS, time):

    | Scenario | `rug` | `malachite` |
    |---|---|---|
    | `fourfold`, deg 15, GW | 44 MB, 1.30 s | 44 MB, 0.64 s |
    | `h11_9`, deg 10 | 130 MB, 5.3 s | 133 MB, 2.7 s |
    | `h11_9_plike`, deg 240 | 45 MB, 34.4 s | 46 MB, 28.1 s |
    | `h11_8`, deg 14 | 514 MB, 60.6 s | 469 MB, 33.3 s |
    | `h11_9`, deg 14 | 665 MB, 120.8 s | 629 MB, 75.3 s |
    | `h11_10`, deg 10 | 549 MB, 19.6 s | 488 MB, 15.2 s |
    | `h11_11`, deg 12 | 927 MB, 28.3 s | 844 MB, 26.9 s |

    The toy `threefold` is the exception: criterion gives 0.48 / 0.55 s at degree 30 and 5.66
    / 6.74 s at degree 50 (GV; GW 4.71 / 5.63 s), i.e. 15–20 % slower. (The memory bench
    makes that look like 2x: at degree 50 `malachite` makes 151 million allocations to
    `rug`'s 10 million, and every one goes through the bench's shared atomic counters.) The
    modular pipeline is unaffected, as expected (`h11_9` at 14: 1.61 / 1.46 s). Builds from
    clean, locally: dev 96 s with GMP built from source (what an uncached CI job paid) against
    49 s; release 159 s against 163 s, since `malachite` is large and the release profile is
    one LTO unit; the CLI binary grows from 3.0 to 4.5 MB.

- [ ] **7a. Double description (committed, not a prototype).** Can start any time. Support
  hyperplanes from spanning vectors, in any dimension, with exact integers. Needed regardless
  of how 7b turns out: it gives the Mori cone's facets (target points, `mori_rays` input, the
  default grading) and the hull facets that phase 7c needs.
  - [x] Decide: our own implementation or `howzat`. Our own (decided with you), from this
    evaluation.
    `howzat` 0.4 (with `dashu`, not its default `rug`):
    - It is AGPL-3.0-only, as are its dependencies `calculo` and `hullabaloo-core`. GPL-3.0
      code may be combined with it, but the combination would put AGPL's network clause on
      every downstream service built on `cygv`.
    - Correct and fast on the full-dimensional Mori cones: on the generators of the four
      `benches/data` files it gives Normaliz's facets in 0.06–0.95 s.
    - It drops facets on a lower-dimensional cone. On the generators of
      `tests/data/h11_9_facet.yaml` (an eight-dimensional cone in $\mathbb{Z}^9$), fed the
      cdd way (the origin as a point, the generators as rays, as its own tests do), it returns
      the equation and 7 of the 9 facets, whatever the order of the generators. Normaliz's 9
      are all facets (each is tight on generators of rank 7). Fed as rays only, it returns
      the 9 facets but no equation.

    Our own, as a scratch prototype: the dual cone in exact `i128` with overflow checks
    (none overflowed), the lineality handled as cdd does, Motzkin's iteration with the
    combinatorial adjacency test, about 100 lines. It matches Normaliz on all five inputs,
    the facet's equation included (its facet vectors differ from Normaliz's by multiples of
    the equation, so the output needs a canonical reduction modulo the equations). The order
    of the generators decides its speed: as given, 0.01–1.2 s (12–22 s on `h11_11` shuffled);
    sorted by the 1-norm of the generators, 4–68 ms, i.e. 3–50x faster than `howzat`.
  - [x] Implementation and tests. `src/cone/hull.rs`: `cone::facets` (generators to
    `Inequalities`, i.e. hyperplanes and equations, the `ConeData` convention) and
    `cone::rays` (inequalities to `Generators`, extreme rays and a lineality basis), both
    through the dual cone of a set of vectors, public. Exact `i128` with overflow checks; an
    entry of the result that does not fit 32 bits, or an intermediate one 128, is the error
    `TooLarge` rather than a fallback to `malachite` (outputs go into `i32` matrices anyway).
    Canonical output: zero, repeated and redundant inputs are dropped; the equations (the
    lineality basis) are in reduced row echelon form, each row primitive with a positive
    pivot; each hyperplane (ray) is moved along them to vanish at their pivots and made
    primitive; both sorted. The input is sorted by 1-norm, then lexicographically, so the
    result and the run time do not depend on its order: 10 / 30 / 130 / 10 ms on the
    generators of `h11_8` / `h11_9` / `h11_10` / `h11_11` (release, shuffled or not; slower
    than the prototype's best on `h11_10`, 68 ms, whose ties were in file order). Cones with
    lines and lower-dimensional ones are handled: the adjacency test is the combinatorial one
    alone, without the rank precheck, whose count would need the true dimension. Tests: by
    hand (orthant, square, a face, a ray, half-plane, plane, $\{0\}$, faces through
    `rays`), random cones in dimensions 2–5 against a brute force (normals of every $n - 1$
    generators; a third of the cones on $x_n = x_1$, against the brute force of their
    projection), the round trip `facets(rays(facets(G))) == facets(G)`, overflow, and
    `tests/cones.rs::test_facets`: Normaliz's hyperplanes on the four bench geometries and
    on the facet of `h11_9` (Normaliz's hyperplanes reduced along the equation the same way).
    Not wired into the pipelines or the inputs yet: the consumers are 7b (the $K_T$ by
    inequalities, and their rays for coordinate bounds), 7c (the hull facets) and an optional
    `mori_hyperplanes` derived from `generators`.
- [x] **7b. Candidate lattice points without normaliz: time-boxed comparison, then integrate
  or drop** (5.1). Both approaches matched Normaliz; B integrated (decided with you).
  - [x] Direct enumeration of polytope lattice points (A), as a prototype, dropped. One
    polytope per set $T$ of divisors holding exactly the points whose negative intersections
    are with $T$ (so the polytopes are disjoint): the Mori facets, $Q_r \cdot C \ge 0$ for
    $r \notin T$, $Q_r \cdot C \le -1$ for $r \in T$ and the degree bound, over $T$ empty,
    single, the pairs whose form does not vanish, and the strata. Vertices by `cone::rays`
    of the homogenized cone, then the projections onto the first $k$ coordinates, so the
    range of each coordinate given the earlier ones is exact. The enumeration was never the
    cost (milliseconds, almost no empty branches); the projections were. Computing them as
    `cone::facets` of the projected vertices blew up the double description (9–37 s on
    `h11_10` / `h11_11`; reversing the insertion order helped, to 1.5–3.4 s, cdd's cutoff
    rules did not). Fourier–Motzkin from the polytope's facets down, keeping after each step
    only the inequalities that the projected vertices show to be facets (tight on vertices of
    rank one less than the projection's), keeps every step at its true facet count (at most
    15): 21–250 ms per geometry.
  - [x] A correctly specified Pottier algorithm (B): `src/cone/hilbert.rs`, described in
    CLAUDE.md. The reduction ($z$ dropped if an element $w$ of the same sign has
    $|\lambda(w)| \le |\lambda(z)|$ and $z - w$ in the old cone) is unsound when the old cone
    has a lineality space: started from $\pm e_i$ as 5.1 proposed, it collapsed `h11_8`'s
    basis to 2 elements. Started from a unimodular cone containing the Mori cone (a lattice
    basis $m, a_i + K m$ with $m$ the sum of the facets, $K = 4$), its first cut ran over 15
    minutes with more than 45000 elements and growing. Started from the simplicial cone of
    the first $n$ independent Mori facets (fundamental parallelepiped of 1 point on three
    geometries and 16 on `h11_10`), it is immediate: the Mori cone's basis in 1–4 ms (at most
    a few dozen elements throughout), all the $K_T$ in 14–147 ms.
  - [x] Comparison on the benchmark geometries (5.1), release, on the protocol's inputs:
    candidate sets identical to Normaliz's cone data on the four geometries at their Quick
    and Heavy modular degrees (`h11_8` 16 / 20, `h11_9` 16 / 20, `h11_10` 12 / 16, `h11_11`
    16 / 18; the Huge tier has no modular scenarios), on the facet of `h11_9` (A only; B's
    face support came with the integration), on the vex phases h2 and h4 against cone data by
    brute force, and on `h11_9` under its sparse grading at 300 / 400 / 500 against the
    semigroup enumeration, up to the pairs with vanishing forms that both skip. Without
    `[1,0,0,0,0,0,1,0,0]` the facets are the same, so the candidates stay saturated, where
    the semigroup enumeration loses two. Times (A; B's Mori basis + cones; the enumeration of
    the bases into candidates, the same for B and Normaliz's data): `h11_8` 19 / 25 ms,
    1 + 21 ms, 6 / 15 ms; `h11_9` 58 / 94 ms, 1 + 14 / 1 + 23 ms, 11 / 68 ms; `h11_10`
    130 / 186 ms, 4 + 26 / 4 + 37 ms, 10 / 107 ms; `h11_11` 222 / 227 ms, 3 + 122 / 3 + 147
    ms, 10 / 31 ms. Normaliz took 2.5–7 s per geometry. Memory: neither has an intermediate
    beyond a few dozen inequalities or basis elements; the candidates dominate and are the
    same either way. B was chosen as the hybrid it is (compact bases, then the existing,
    tested cone path, and the reference saturated by the Mori cone's basis); its unknowns at
    larger $h^{1,1}$ are the parallelepiped (capped at $2^{22}$ points, `ConeError::TooLarge`)
    and the size of the Mori cone's basis, A's the projections, one polytope per $T$.
  - [x] The winner integrated. `cone::cone_data(hyperplanes, equations, q, strata, grading,
    max_deg)` gives complete `ConeData` (every $T$, vanishing pairs included), checked against
    Normaliz's bases on the four bench geometries and the facet (`tests/cones.rs`). The
    dispatch describes the Mori cone by the cone data's inequalities or `cone::facets` of the
    generators, and without cone data computes it, truncated at `max_deg` lowered to the
    targets' degree, and runs `CandidateSource::Cones`; the target lightcones are then the
    hyperplane ones too. `CandidateSource::Generators` stays for library callers. Errors:
    `ComputeError::MoriCone(ConeError)` (a line in the cone, a grading that is not positive on
    it, overflow). No opt-out and no size threshold yet (6.2).
    Before and after, `cargo bench --bench memory`, eight threads (peak RSS, time; the result
    counts are identical): modular without cone data, which used to enumerate the semigroup,
    `h11_9` deg 12 105 MB, 0.48 s against 27 MB, 0.11 s; deg 14 223 MB, 1.74 s against 38 MB,
    0.17 s; `h11_9_plike` deg 400 69 MB, 0.52 s against 53 MB, 0.18 s; deg 500 329 MB, 3.55 s
    against 69 MB, 0.76 s. The reference, now on the Mori cone's Hilbert basis instead of
    thousands of generators, `h11_9` deg 10 3.5 s against 1.9 s and `h11_10` deg 8 5.0 s
    against 1.6 s, the same peak; the rest, and modular with supplied cone data, unchanged.
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
- [x] **7d. Saturation for the reference pipeline** (6.11). Done with 7b; `auto` has one
  saturation semantics (6.2).
  - [x] The lattice points of the Mori cone up to `max_deg`, which is exactly the saturated
    semigroup the reference enumerates, from the 7b machinery: the reference runs on
    `cone::hilbert_basis` of the Mori cone's inequalities instead of the generators. Except
    without any of `max_deg`, `min_points` and `target_points`, where the generators are the
    curve classes (`Semigroup::from_data`) and are passed as given.
  - [x] `mori_hyperplanes` accepted by both methods: with cone data the generators are
    optional for the reference too (YAML, CLI, Python), except for such an untruncated
    request. `tests/cones.rs::test_saturation`: without the Hilbert basis element of GV 20,
    both methods now give the full result, with or without cone data.
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
  `compute_c_2neg` on `main`, independently of this work (section 7). Fixed on this branch by
  phase 4b.

Follow-ups after the main plan, not tracked here: GPU and the other C features (6.8).

### 5.1 Phase 7 in detail: cone computations without normaliz

Goal: produce the cone data of 6.1 inside the crate, so the fast path does not depend on an
external program. No suitable crate was found on crates.io (no Hilbert basis crate, no normaliz
or 4ti2 bindings; `howzat` does double description only, and was evaluated under phase 7a).
Outcome: both candidates below worked; B (Pottier) was integrated, see phase 7b.

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
- Unsaturated input (`h11_9` without `[1,0,0,0,0,0,1,0,0]`, see 6.2; the hole in `h11_11.yaml`
  does not change any invariant): the candidates must come out saturated.

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
   look for published tables to use as fixtures. (`examples/fourfold.yaml` used to describe
   itself as a hypersurface; it is the complete intersection of the `fourfold` bench model,
   and has its nef partition since phase 3b.) And before shipping `refs_vex` as fixtures,
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
   just rays) or inequalities (a new `mori_hyperplanes` field). **Decided (2026-10-10): the
   Python wrapper does not call `normaliz`**, not even as an opt-in; it accepts cone data the
   caller computed, and waits for phase 7b.
2. **Saturation semantics. Direction: saturate by default.** The reference uses the semigroup
   the generators *generate* (non-negative integer combinations); the correct object is *all
   lattice points of the cone they span*, $\sigma \cap \mathbb{Z}^n$. The two differ whenever
   the input misses a Hilbert basis element of the cone, e.g. generators (1,0), (1,2) miss (1,1).
   This is a common problem for cygv users today: CYTools passes the rays plus a finite sample
   of lattice points, which can miss an element. `benches/data/h11_11.yaml` does (its 2032
   generators miss `[2,0,0,2,1,0,0,-1,-1,0,0]`, degree 6; not an artefact of the file's
   trimming at degree 16), and the reference then silently computes over a semigroup with holes.
   (Found in phase 6: that point meets three divisors negatively and has GV 0, and the results
   with and without it agree up to degree 12, so this hole is invisible. A missing Hilbert basis
   element with nonzero GV is not: dropping `[1,0,0,0,0,0,1,0,0]` (GV 20) from `h11_9` makes
   the reference panic on a non-integer and the enumeration path fail with `NotStable`, while
   the cones give the full result; `tests/cones.rs` tests the latter two.)

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
     need saturated inputs, or the reference fed the completed generator set. (`h11_11` turned
     out not to show the difference; see above for the test that does.)
   - Reporting which Hilbert basis elements the input missed (the C's old note) would need the
     Mori cone's basis again. Dropped, as in #100; it can come back as a diagnostic if wanted.

   **Done (phases 7b and 7d)**: both methods saturate, from the Mori cone's inequalities
   (the cone data's, or the facets of the generators), with the exception of requests without
   any truncation, whose generators are the curve classes as given. No warnings are needed,
   since `auto` never mixes semantics.
   Still open: what "too large" means (a dimension threshold, a time budget, or both; nothing
   on the benchmark geometries comes close, so no threshold yet), whether skipping for size
   should warn or error, and the name of the opt-out; none exists yet.
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
     scales with $M$. (Settled in phase 1: the slots and terms of the tables in `table` hold
     a number of values fixed when the table is made.)
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

   Unmeasured, to settle in phase 5 (still open after it): the cost of 32-byte keys against 16 (estimate: 10–30% on
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
   but `rug` is gone (6.10), so the certifier needs another route: a pure-Rust type with
   directed rounding (`malachite`'s `Float`, behind its `floats` feature, has rounding modes),
   or `f64` with explicit outward rounding of each operation. The C also re-expresses $\alpha$ along the rows of the charge
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

   Still open: whether the crate root keeps re-exporting `Semigroup`, `Polynomial`,
   `PolynomialProperties` and `PolynomialCoeff`, or those are reached through
   `cygv::reference`. (The question of what `auto` does with `prec` went away with phase 0c.)
8. **Features of the C that are not ported now (decided): future options, not ruled out.** The
   GPU extraction stage, the distributed mode (one run split by primes across machines), the
   low-memory switch and the C's tuning knobs are not part of this work, but the Rust code is
   structured so that adding any of them later is an addition, not a refactor. The resulting
   constraints are in section 4.1 and apply from phase 1. PR #95 (the GPU variant through
   `build.rs`, `device=` in Python) is the reference for the GPU interface when it comes.

9. **Dropping floating-point computation (done, phase 0c).** Remove `prec` and the `rug::Float`
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
   - Checked before phase 0c: neither CYTools nor any other public caller passes `prec`, so it
     was removed outright instead of going through a deprecation release.
   - Done as its own commit right after phase 0b (phase 0c), before the `Method` argument is
     added, so the dispatch is rewritten once.

10. **`rug` replaced by `malachite` (done, phase 5b).** The plan was to make `rug` an optional
   feature, with a pure-Rust crate by default, since speed now comes from `modular`. The
   benchmarks of phase 5b changed that: `malachite` is as fast as `rug` on the reference
   pipeline, faster on most scenarios, so `rug` was dropped altogether rather than kept as a
   feature.
   - *Shape*: `malachite` (with only its integers and rationals) is the only big-number crate.
     `modular` uses `malachite::Integer` in the CRT lift and `malachite::Rational` for GW. The
     reference pipeline stays generic over `PolynomialCoeff`, whose `Assign` is now our own
     trait, and its one implementation is `reference::Rational`, a newtype over
     `malachite::Rational` (the orphan rule forbids implementing `AddAssign<i32>` and the like on
     the foreign type). The typed entry points return `malachite::Integer` and
     `malachite::Rational`, and the crate re-exports `malachite`.
   - *License*: `malachite` is LGPL-3.0-only. LGPL-3.0 is GPL-3.0 plus additional permissions,
     so GPL-3.0-or-later code may use it, and the combination (the binary, the wheels) is
     distributed under GPL-3.0; the source of the whole is available, which covers the LGPL's
     relinking requirement. `rug` and GMP were LGPL-3.0-or-later, so nothing changes for anyone
     who uses the crate under GPL-3.0; the only thing given up is distributing a combined work
     under a hypothetical GPL-4.0.
   - *MSRV*: `malachite` 0.13 needs Rust 1.90, which is now `rust-version`.
   - *What it bought*: no from-source GMP build, so no C toolchain, `m4` or `make`, no
     `GMP_MPFR_SYS_CACHE` handling and no `CC` consistency rule. On Windows the MSYS2/MINGW64
     route existed for GMP (gmp-mpfr-sys only builds with MinGW), so CI and the wheels now use
     the stock MSVC toolchain and cibuildwheel, and `--no-build-isolation` and
     `PYO3_USE_RAW_DYLIB=0` went with it. With no custom shell left, the cargo caches are
     `Swatinem/rust-cache` again. The memory bench no longer needs GMP allocator hooks.
   - *Certifier* (6.5, planned): it assumed `rug::Float` with directed rounding; see there.

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
  behaviour (done); a minimal fix (an error) could also go to `main` on its own.

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
  porting. Phase 1 uses `u32` slot numbers, i.e. about four billion slots, which is more than
  fits in memory.
- Peak memory is per-thread working tables plus $I$; the C's README quotes 3–6 GB at degree
  26–28 for $h^{1,1}$ = 10–11. `benches/memory.rs` counts Rust allocations, which since phase
  5b are all of them.
- What happens to PR #94 itself (close, or merge first and replace later) is your call and
  affects how the author is credited.
