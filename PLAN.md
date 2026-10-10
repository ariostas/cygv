# Plan: a Rust implementation of the fast GV algorithm from PR #94

Status: **investigation done, phases 0 and 0d done, nothing implemented yet.** Where sections 1–5
mention `rug` for the new code, read "big integers": decision 6.10 makes `rug` optional. Branch
`gv-modular-rust` is PR #94 (`natemacfadden:cgv-cc-cpu`, 11 commits) rebased onto `origin/main`
(6fe293a, v0.3.0), plus our commits. The PR branch itself (`cgv-cc-cpu`) is untouched. The
follow-up PRs #95, #97 and #100 are covered in section 1.1.

Section 6 lists the decisions that need an answer before implementation starts.

## 1. What PR #94 does

It adds `cgv/`, a standalone C program (`gv.c`, 2226 lines, plus a 88-line driver), and wires it
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
6.3 (target points decided; phases added), phases 6–6c and 7a–7c, and section 7.

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

src/reference.rs            run_hkty (today's src/hkty.rs minus the dispatch)
src/reference/semigroup.rs, polynomial.rs (+ polynomial/), factorial.rs,
              fundamental_period.rs, instanton.rs, series_inversion.rs, misc.rs
                            moved as they are, with their error modules and tests

src/modular.rs              entry point: run(...), options and stats structs, pass loop, lane probe
src/modular/field.rs        Montgomery arithmetic, prime search, Lanes<const NL: usize>
src/modular/key.rs          Key trait; packed u64/u128 with per-coordinate fields, inline arrays; bound check, hash
src/modular/table.rs        slot map, sparse polynomial, per-degree level tables, heap
src/modular/candidates.rs   generator reduction, semigroup enumeration, cone union, saturation check
src/modular/hilbert.rs      (phase 7, if it passes) support hyperplanes, truncated Hilbert bases of the K_T
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
5. **Stages are generic over the coefficient, not hard-wired to residues.** The certifier
   (6.5) is the same pipeline over upward-rounded non-negative reals; the C gets it by swapping
   one type at compile time. So the stages are written against a small trait (zero, from
   integer, add, sub, mul, fused multiply-add, exact inverse, is-zero) with the residue lanes
   as the main implementer. Only the lane implementer needs to satisfy constraint 2.
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

0. **Branch hygiene.** Done (6.6): the PR's integration glue is removed and `cgv/*.c` stays as
   a development oracle and performance baseline.
0b. **Move the existing pipeline under `src/reference/`** (6.7), as a commit of its own with
   no functional change: `git mv`, fix the `use` paths in `src`, `tests/`, `benches/` and the
   CLI, split `src/hkty.rs` into `reference.rs` and `dispatch.rs`, update CLAUDE.md. Done
   before any new code so that the diff stays a pure move.
0c. **Drop floating-point computation** (6.9), as its own commit: `prec`, the `Float` variants,
   the float bench variants, `mpmath`.
0d. **Sync the oracle.** Done (6.6): `cgv/` is the #100 version, without the integration glue.
1. **Primitives**: `field`, `key`, `table`, with unit tests (field ops against `rug`, pack/unpack
   round trips and overflow detection).
2. **Period and instanton, enumeration path only**: candidates by full enumeration, `period`,
   `instanton`. Test $I$, `c0`, `c1` against the reference stages mod p.
3. **Extraction, multicover, CRT**, single-threaded and simple (`apply_curve` only), for
   threefold hypersurfaces. End-to-end differential tests against `run_hkty` on `examples/` and
   `benches/data`.
3b. **Complete intersections and n-folds** (decision 6.3): multi-row `q0`, the $M$-component
   $I$, the exponent-2 multicover. Differential tests against the reference on the fourfold
   models. Done before tuning so that phase 5 optimises the final data layout.
4. **Front end**: the `Method { Auto, Modular, Reference }` argument through `dispatch`,
   `io::Input`, the CLI and Python, plus GW. `auto` falls back to `reference` for requests
   `modular` does not support; an explicit `modular` returns an error for them. `auto` only
   becomes the default once the differential tests of phases 3 and 3b pass.
5. **Parallelism and tuning**: rayon over candidates / product chunks / curves per layer, sharded
   $I$ with per-shard locks, the level-synchronous mode for few large curves, lane probe. Add
   bench scenarios; target parity with the C on the table in section 2.
5b. **Make `rug` optional** (6.10): choose the pure-Rust crate from benchmarks, implement
   `PolynomialCoeff` for it, backend-neutral result types, the `rug` feature, then simplify
   the CI workflows and CLAUDE.md accordingly. After phase 4, so that `auto` already routes
   most requests to `modular` when the reference pipeline's default build gets slower.
   `modular` is written against the pure-Rust crate from phase 3 on, so a provisional choice
   is needed there; it is confined to `crt` and `multicover`.
6. **Cone fast path through optional input.** `candidates` takes cone data (defined in 6.1) and
   enumerates each $K_T$ from its Hilbert basis, with the containment and vanishing-form skips.
   The Mori cone may be described by spanning vectors or by inequalities. The data comes from
   the caller at this point: the Rust API and YAML accept it, and the Python wrapper may produce
   it with `normaliz`. This fixes the interface that phase 7 plugs into.
6b. **Target points for `modular`** (6.3): the facet test at the three filter points, optional
   `max_deg`, differential tests against the full computation and the C's lightcone mode.
6c. **Vex phases, core** (6.3): the `strata` input, the $m < 0$ period branch, the extra $K_T$;
   tests against the C oracle and `cgv/tests/refs_vex`. Also make `reference` reject negative
   anticanonical degree.
7. **Native cone computations.** Independent of phases 1–5 and can run alongside them. See
   section 5.1.
   - *7a, double description (committed, not a prototype)*: support hyperplanes from spanning
     vectors, in any dimension, with exact integers. Needed regardless of how 7b turns out: it
     gives the Mori cone's facets (target points, `mori_rays` input, the default grading) and
     the hull facets that phase 7c needs. Our own implementation or `howzat`, decided then.
   - *7b, Hilbert bases: time-boxed prototype, then integrate or drop.* The $K_T$ bases.
7c. **Phases from a fan** (6.3): derive `strata`, the default Mori rays and grading, and
   optionally check `kappa`, from the fan's maximal cones and `q`, without normaliz. Port of
   `cgv_phase.py`: rays as the integer Gale dual of `q`; walls (3-faces shared by two maximal
   cones) and their classes by Cox–Little–Schenck 6.4.4; vex cones as the fan's 2- and 3-faces
   whose rays share no facet of the hull of the rays (facets from 7a); the toric intersection
   ring for the `kappa` check. Exact rational linear algebra on small matrices with the
   big-number crate of 6.10. Estimated 400–600 lines with tests. Tested piece by piece against
   `cgv_phase.py` (normaliz needed on the development machine only) and end to end against
   `cgv/tests/refs_vex`. A directly supplied `strata` input stays accepted. Rejects inputs
   outside the paper's assumptions (a fine fan on the boundary points of a reflexive
   4-polytope not interior to facets), e.g. vex 2-cones (MacFadden–Sheridan Prop. 5). Note
   that the vex path is only fully normaliz-free if 7b succeeds: each stratum adds a $K_T$
   whose Hilbert basis is still needed.
8. **Certifier** (6.5): the majorant run, fed with the candidate GVs, behind its own switch;
   tested against `cgv_maj` on the benchmark geometries. Optional extras in the same phase:
   generator-reduction back-port to `Semigroup`.
9. **Wrap up**: remove `cgv/`, update README / CLAUDE.md / docs, credit the PR author.

Follow-ups after the main plan: GPU and the other C features (6.8).

Phases 1–5 need no decision except 6.5–6.7 and deliver the ~70× column of section 2.

### 5.1 Phase 7 in detail: cone computations without normaliz

Goal: produce the cone data of 6.1 inside the crate, so the fast path does not depend on an
external program. No suitable crate was found on crates.io (no Hilbert basis crate, no normaliz
or 4ti2 bindings; `howzat` does double description only and has not been evaluated).

The problem is much narrower than the general one normaliz solves:

- The Mori cone's own Hilbert basis is not needed (6.1). Only its support hyperplanes are, from
  whatever description the input gives.
- Each $K_T$ (including one per vex 3-cone) is that cone cut by a few halfspaces $Q_r \cdot C \ge 0$. Cutting a cone with a
  known Hilbert basis by one halfspace is Pottier's algorithm (normaliz's "dual mode"): split
  the basis by the sign of $Q_r \cdot C$, keep adding irreducible sums of a positive and a
  negative element until nothing new appears, keep the non-negative ones.
- Only elements of degree ≤ `max_deg` are needed, and sums only increase degree, so the
  completion can be truncated. The result then depends on the grading vector and `max_deg`,
  unlike normaliz's output.
- The cones are nested ($K_\emptyset \subset K_{\{r\}} \subset K_{\{r,s\}}$), so intermediate
  bases can be shared.

Pieces:

1. *Support hyperplanes of the Mori cone*, needed for the membership test inside the
   irreducibility check. This is phase 7a, which is committed independently of the prototype
   (also needed for target points and phase 7c). The caller may still supply them.
2. *Truncated halfspace cutting* as above, with the reduction test.
3. *Driver* over all $T$ with $\lvert T\rvert \le 2$, reusing intermediate bases and skipping
   two-divisor cones with vanishing form before computing them.

Prototype protocol:

- Run on the four `benches/data` geometries (and `h11_9` under its sparse grading) at their
  Heavy and Huge degrees.
- Correctness oracle: `normaliz` (installed locally) on the same cones, compared after
  truncating its output to `max_deg`; plus identical candidate sets and GVs end to end.
- Unsaturated input (`h11_11.yaml`): the $K_T$ bases must come out saturated whatever lattice
  points the input lists, per decision 6.2.

Exit criteria:

- *Integrate* if every benchmark geometry takes at most a few seconds (normaliz: ~1.7 s for the
  91 cones of `h11_9`) and matches the oracle. It becomes `src/modular/hilbert.rs`, feeding the
  phase 6 interface; the optional-input path stays for callers that already have the data.
- *Drop* if intermediate bases blow up or the time is not competitive with the enumeration
  fallback. The optional-input path of phase 6 remains the only fast path, and the findings are
  recorded here.

Main risk: Pottier's algorithm can produce large intermediate bases, and its behaviour at
$h^{1,1}$ = 10–11 is unknown until measured.

## 6. Decisions needed

1. **Hilbert bases (the main one).** The last ~15–20× needs the Hilbert basis of the Mori cone
   and of every $K_T$ (91 cones for `h11_9`), plus the Mori cone's support hyperplanes. The PR
   gets them from the `normaliz` binary in Python.

   "Cone data", following #100: (i) the Mori cone's support hyperplanes (facets), needed for
   target points and the coordinate bound; (ii) the Hilbert basis of each $K_T$: every $T$ with
   $\lvert T\rvert \le 2$ (13 + 78 = 91 cones for `h11_9`, of which 19 are enumerated after the
   skips) plus one per vex 3-cone (6.3); (iii) the divisor set $T$ of each cone, which enables
   the vanishing-form skip. Independent of the grading vector and `max_deg` when computed
   exactly. The Mori cone's own Hilbert basis is *not* needed: the $K_T$ already contain every
   candidate, and the basis was only ever used for the saturation check, which #100 deletes.
   The data is trusted: a wrong basis for some $K_T$ silently drops candidate curves.

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
   - This raises the stakes of phase 7: a native Hilbert basis computation is no longer only a
     speed-up, it is what makes the default correct without an external program.
   - Unlike the C after #100, we keep the enumeration of the *generated* semigroup, but only as
     the path for the opt-out cases above (and as the phase 2 stepping stone). It needs lattice
     points as input, so it is unavailable when the cone is given by rays or inequalities only.
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
   - *Not in `reference`*: a request with `strata` under `auto` goes to `modular`, and
     `reference` rejects it. Separately, `reference` must reject negative anticanonical degree
     instead of wrapping it (section 7).
   - *Derivation from a fan: in the main plan (phase 7c, decided).* Deriving `strata`, the
     default Mori rays and grading, and the `kappa` check from the fan, i.e. the content of
     `cgv_phase.py`, natively and without normaliz. Its only normaliz uses are two small
     convex hulls, which the double description of phase 7a covers. Callers can still supply
     `strata` directly.
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
   checks) both pass here. Still open: whether to keep some of `cgv/tests/refs` (38 frozen outputs, 6.2 MB)
   as test data after that, or rely on the reference implementation and `benches/data`.
7. **Naming, default and layout (decided).**
   - *Names*: the existing pipeline is `reference`, the new one is `modular`. Not "legacy": the
     old pipeline is kept on purpose as the faithful HKTY implementation and remains the only
     route for `min_points`, `target_points` (for now) and anything `modular` cannot handle.
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
   - To check before it ships: whether CYTools or other callers pass `prec`.
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

## 7. Notes and risks

- **Latent bug in the reference pipeline, worth fixing on `main` independently of this work.**
  `compute_c_1neg` and `compute_c_2neg` (`src/fundamental_period.rs`) convert $q_0 \cdot C$ with
  `as u32`. In a vex phase that value can be negative, so it wraps to a huge factorial argument
  instead of failing. The same holds for any inconsistent input. It should be an error.

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
