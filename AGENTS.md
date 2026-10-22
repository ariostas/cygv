# AGENTS.md

This file provides guidance to agents when working with code in this repository.

## Project

`cygv` computes Gopakumar-Vafa (GV) and Gromov-Witten (GW) invariants of Calabi-Yau manifolds
realized as hypersurfaces or complete intersections in toric varieties, via the HKTY procedure.
It is a Rust reimplementation (with improvements) of the algorithm in
"[Computational Mirror Symmetry](https://arxiv.org/abs/2303.00757)".

It ships as both a Rust crate (`cygv` on crates.io) and a Python package (`cygv` on PyPI) built
with maturin/PyO3. Licensed GPLv3+.

## Commands

```bash
# Rust
cargo build
cargo test                            # unit tests in src/**, integration tests in tests/
cargo test test_threefold             # single test by name
cargo test --test full_hkty           # single integration test file
cargo doc --open
cargo run -- --file examples/threefold.yaml   # the CLI
cargo build --no-default-features     # library only, i.e. without the cli feature

# Benchmarks (see the Benchmarks section below)
cargo bench                           # both targets, quick scenarios only
cargo bench --bench bench             # wall clock (criterion)
cargo bench --bench memory            # allocations and peak memory
cargo bench -- fourfold               # only the scenarios whose id contains "fourfold"
CYGV_BENCH_HEAVY=1 cargo bench        # add the high-degree scenarios
CYGV_BENCH_HUGE=1 cargo bench --bench memory  # add the largest ones (~45 min, ~12 GB peak)
CYGV_BENCH_THREADS=1 cargo bench      # pin the worker thread count

# Python (builds the Rust extension via maturin)
uv sync --extra tests                 # creates .venv and builds the extension
uv run pytest
uv run pytest python/tests/test_hkty.py::test_threefold

# Lint (runs rustfmt, cargo-check, clippy, ruff, mypy, codespell)
prek -a --quiet
```

Every dependency is pure Rust, so a stock Rust toolchain builds everything, on Windows the
default MSVC one. The big integers and rationals come from `malachite` (LGPL-3.0-only, which
GPL-3.0-or-later code may link: the combination is distributed under GPL-3.0). It replaced
`rug`, whose GMP had to be compiled from source with a C toolchain, `m4` and `make` (and on
Windows through MSYS2/MINGW64), and whose build and its caching dominated CI; see `PLAN.md`
phase 5b for the comparison that chose it.

Both CI workflows cache `~/.cargo/registry` and the dependencies in `target/` with
`Swatinem/rust-cache`, which keys on the job id, the toolchain and the lockfile. Both workflows
name their job `build-and-test`, so each passes its own `key` (`dev` for the Rust workflow, which
builds the dev profile, `release` for the Python one); sharing an entry would leave whichever job
saved second restoring a `target/` built with the wrong profile.

`[tool.uv] cache-keys` in `pyproject.toml` lists `src/**/*.rs`; without it uv would not rebuild the
editable extension when only Rust sources change, and `uv run pytest` would test a stale build.
`uv.lock` is gitignored — it pins only the dev environment and does not affect consumers.

`Cargo.lock`, by contrast, is tracked. The usual "libraries do not commit a lockfile" rule assumes
consumers re-resolve from `Cargo.toml`, which holds for the crates.io half of this project but not
for the PyPI half: those users install wheels that CI compiled, so the lockfile decides what they
actually get. It also keeps the benchmarks honest — comparing two revisions is only meaningful if
the dependency set is held fixed — and keeps CI failures attributable to commits rather than to a
dependency that drifted underneath them. The cost is that CI no longer continuously exercises the
newest semver-compatible dependencies; the monthly grouped cargo Dependabot PR is what covers that,
so it should be treated as a real test run rather than rubber-stamped.

Every workflow that compiles Rust enforces the lockfile, so a manifest change committed without its
lockfile update fails instead of silently re-resolving. The Rust workflow passes `--locked` to
cargo directly; the Python and deploy workflows reach cargo through maturin and set
`MATURIN_PEP517_ARGS=--locked` instead. Two wrinkles there: cibuildwheel needs it repeated inside
`CIBW_ENVIRONMENT` to reach the container builds, and maturin's `build_sdist` never reads the
variable, so the sdist job is unaffected by it.

PR titles must follow conventional commits (enforced by `.github/workflows/pr_title.yml`).

## Architecture

The modules split into the pipeline and the front end shared by everything that calls it.
There are two pipelines. `src/reference.rs` and `src/reference/` hold a faithful implementation
of the HKTY procedure in exact rationals (named *reference* because the faster, differently
structured `modular` pipeline sits beside it; see `PLAN.md`). `src/problem.rs` holds the public
input (`Problem`, `Options`) and `src/dispatch.rs` the public entry points, which pick a
pipeline by the `Method` enum of `src/lib.rs`, and `src/io.rs`, the CLI and `src/python.rs` sit
on top of those (see "One problem, one dispatch" below). `src/cone.rs` holds the *cone data*
(`ConeData`: the Mori cone's inequalities and equations, and the Hilbert basis of each cone
$K_T$, the curve classes in the Mori cone that only meet the divisors in $T$ negatively, for
every $T$ of one or two divisors and, in a vex phase, the three divisors of every stratum) and
its validation (`check_cone_data`), and the native cone computations of `PLAN.md` phase 7,
shared by both pipelines. Callers may still supply the Hilbert bases, computed with an
external program such as Normaliz (`Options::hilbert_bases`, with the cone given by
`MoriCone::Inequalities` or generators), and they are used as given. `src/cone/hull.rs` (phase 7a):
`cone::facets` and `cone::rays` convert between generators and inequalities by the double
description method, in `i128` with overflow checks (`TooLarge`), with a canonical output
(equations in reduced row echelon form, each hyperplane reduced to vanish at their pivots).
It sorts its input by 1-norm first, which is what makes it fast on Mori cones of thousands of
generators. `src/cone/hilbert.rs` (phase 7b): `cone::hilbert_basis` and `cone::cone_data`
compute Hilbert bases by Pottier's algorithm (Normaliz's dual algorithm): cut a pointed cone
with a known basis by one inequality at a time, pairing elements of opposite signs and
dropping the sums an element of the same sign lies below. Its reduction is only sound for a
pointed cone, so it starts from the simplicial cone of the first independent facets (its
basis from the lattice points of the fundamental parallelepiped, at most
`MAX_PARALLELEPIPED`), not from $\pm e_i$, and a face is handled in a lattice basis of its
span. `cone_data` cuts the Mori cone's basis by the rows $Q_r$ down a tree (each row cut or
left out, at most two left out, three for a stratum), dropping sums above `max_deg`, which
is sound because every cone is then inside the Mori cone. As in #102's `ed93fb4`, only the
rows negative on some element of that basis enter the tree ($K_T = K_{T \cap N}$), and each
distinct cone is given once, labeled by a $T$ with the same negative set $M$ on its basis
that is not a pair when there is one, so that the candidates' skip of vanishing pair forms
keeps exactly the union of the full set (the Mori cone alone, as $K_{\{0\}}$, when nothing is
negative): a face of the Mori cone at $h^{1,1} = 491$ goes from minutes to under a second.
User-supplied Hilbert bases must still have every $T$. Errors are `ConeError`. Phase 7b
also tried enumerating the lattice points of each $K_T$ directly (`PLAN.md`); both matched
Normaliz, and the Hilbert bases won because they feed the existing cone path unchanged.
`src/phase.rs` (phase 7c, a port of `cgv_phase.py` of the C program of PR #94): `phase::Fan`
takes the maximal cones of the fan of the ambient toric fourfold and the charges, with the
rays as the integer Gale dual of the charges (`kernel_basis` of `src/cone/hilbert.rs`), and
checks the paper's assumptions as far as the fan shows them (simplicial, complete, a reflexive polytope of rays
none of which is interior to a facet, regular through a pointed Mori cone, no vex 2-cones;
that every boundary point not interior to a facet is a ray is not checked). `Fan::strata`
gives the strata (the 3-cones whose rays share no facet of the polytope, from
`cone::facets`, with their wall classes, Cox–Little–Schenck 6.4.4), `Fan::mori_cone` the
Mori cone (the facets of the cone of the wall classes) and the default grading vector (the
sum of its facets), `Fan::check_grading` checks a given grading vector on the wall classes,
and `Fan::check_intnums` compares the intersection numbers with the toric intersection ring
(`Ring`, recursive by linear equivalence, over the neighbors of each divisor) as sparse
trilinear forms on $h^{1,1}$ independent divisors. The wall classes go through a left
inverse of the charges computed once (`left_inverse`, free for CYTools' identity block).
Everything but `mori_cone`, a double description in dimension $h^{1,1}$, takes about a
second at $h^{1,1} = 491$, so the dispatch only calls it for a default it needs. Exact
`malachite::Rational` arithmetic; errors are `FanError`.
`tests/phase.rs` compares it with what `cgv_phase.prepare` derived for the 20 fans of
`tests/data/refs_vex.yaml`.

`src/modular.rs` and `src/modular/` hold the faster pipeline, which is still being extended and
tuned. It computes modulo word-sized primes and lifts the results with the CRT, and shares no
code with `reference`. Its primitives are `field` (Montgomery arithmetic,
`Lanes<NL>` residues, and the `Ring` trait the stages are written against, which takes the
primes as a context), `bound` (`BoundRing`, the other `Ring`: upward-rounded nonnegative
bounds for the certifier of `PLAN.md` phase 8, an `f32` mantissa with an `i32` exponent, 8
bytes, moved up by `next_up` after each operation, subtraction being addition), `key` (curve classes packed
into keys of 128, 256 or 512 bits that add like the classes, `PackedLayout<W>` over a
`Word`, `u128` or `Wide<N>` of `N` words, with field widths sized per run from coordinate
bounds; `Setup` takes the narrowest word that fits, `Tuning::min_key_bits` forces a wider
one for the tests, `with_word!` dispatches, and the `Lift` trait runs wide keys with four
lanes only, which keeps the build from compiling the pipeline for every lane count and
width; 64-bit keys were measured and gave nothing, see `PLAN.md` 5.3), `lattice` (after the
C's `setup_lattice`: when the vectors the candidates are sums of span a lattice of rank
below $h^{1,1}$, as on a face of the Mori cone, the run keys the curves by their coordinates
in an LLL-reduced basis of it, which takes a face of rank 24 at $h^{1,1} = 491$ from ~300
bits to ~120; `Lattice::spanned` goes through the vectors shortest first, solving for the
coordinates in floating point and checking them exactly, adding a vector independent modulo
a 61-bit prime, joining one of fractional coordinates exactly, and reducing by Cohen's
integral LLL in malachite integers; `Setup` moves the source, the targets, the grading
vector and the geometry into it (`Reduced`, `Geometry::in_lattice`), `pass` restricts
$\alpha$ to it after the instanton (`restrict_series`) unless the extraction goes along the
rows, and the output goes back; `Tuning::lattice` turns it off), `table` (open-addressing tables and sparse series over those keys, each entry
holding a run-time number of values; the index of a `SlotMap` has 8-byte buckets of a tag
and a slot) and `level` (the extraction's level tables, `PLAN.md` phase 8b, after the C of
PR #102: they key the offsets of an exponential by the linear hash $H = \varphi(m) M$, with
$\varphi$ packing a box of the offsets, less the coordinate the degree fixes, into a `u64`
when it fits, a `u128` otherwise, and for wide keys into a word as wide as theirs
(`Key::Level`), so that a scatter's target is one addition, and
`OffsetKeys` turns a hash back into a key only for the emitted terms; the tree's nodes keep
their terms in `CompactTable`s, of any capacity and filled to 13/16, which take two thirds of
the memory of a `LevelTable` and are as fast). A *pass* (`pass`, generic over the `Ring`) runs the
stages modulo `NL` primes: `candidates` (from a `CandidateSource`: either
the sums of the Hilbert basis elements of each $K_T$ of the cone data, skipping cones contained
in another and two-divisor cones whose form $Q_s K^{(k)} Q_u$ vanishes, as the C does, or the
whole semigroup of the generators up to `max_deg`, enumerated in order of degree and filtered
to the classes with at most two negative intersections; the cones are much faster and give
every lattice point of the Mori cone, the enumeration only those the generators generate, so
the dispatch always passes cones, computed when the caller gave none, and `Generators` is
left to library callers and the tests; in
a vex phase both also keep the curves of negative anticanonical degree on a stratum,
`Geometry::carries_period`), `period` (tables of factorials and harmonic numbers, then $c_0$,
$c_1$ and the second derivatives already contracted with the bilinear forms, $S_2$, to which
a curve on a stratum adds $U_0\,(w \cdot d_S)$, the reference's $c_0^{-1}$ `vex[t]` term
contracted with the grading; the strata reach `Geometry` through `with_strata`), `instanton` (sparse
products and reciprocal, $\alpha$, and the instanton series $I$), `extraction` (the
coefficients of $I$ in the basis $z^C \exp(C\cdot\alpha)$, curve by curve in order of degree,
each exponential built by the Euler recurrence in per-degree level tables, the emitted terms
subtracted through shards whose buffers are hashed first and prefetched while they are
folded in; along the rows of the charges, $\exp(C \cdot \alpha) = \prod_r E_r^{Q_r \cdot C}$
with $E_r = \exp \beta_r$, and `extraction/tree.rs`, `PLAN.md` phase 8e, shares the powers of
the rows of the most terms between the curves: a tree of the curves by their exponents on
those rows, by default at most two of them, whose nodes are multiplied by the power of their
row once for all of their curves, relaxed in the degree, while each curve builds the
exponential of the other rows as before; the pass takes $\beta$ along the rows whenever the
tree is used, which `extraction::Tuning::tree_depth` and `uses_tree` decide: by default for a
series of one component, while an n-fold's several components go curve by curve, the tree
taking them one at a time; the shards evict the final terms after every layer, most of them
zero, and two options of the tuning trade time for memory, which the low-memory mode turns
on: `sinks`, which sends the terms of the degrees no power reaches beyond straight to the
root, into tables split by hash that the threads fill through buffers, and `leaf_window`,
which has the leaves take a window of degrees at a time, every curve building its
exponential again for each later window) and `multicover` (the GV
invariants from those, by gcd of the coordinates, solved in place in order of degree;
exponent 3 and a division by the degree for threefolds, exponent 2 for n-folds). `compute_gv` in `src/modular.rs` repeats passes with the
next primes until `crt::Residues` gives a stable lift (the lift from all but the last
`check_primes` primes must agree with those), skipping a pass whose primes cannot invert a
denominator of the strata (`Error::NotInvertible`), with `Options` for the lane count (1 to 4 primes
per pass; by default `choose_lanes` predicts it from a probe run at 3/4 of the degree, as the
C's driver does), a minimum number of primes, the thread count, the extraction's tuning,
`certify` and `low_memory` (`PLAN.md` phase 10: one prime per pass, the sinks and leaf windows
of one degree, for about a third of the memory in two and a half times the time; the
command line and the Python module also call `return_freed_memory`, which sets glibc's
`M_MMAP_THRESHOLD` so that large freed blocks go back to the system, a third less resident
memory); big integers (`malachite::Integer`) are only built by `crt`, for the final values.
With `certify` (the certifier, `PLAN.md` phase 8, after the C's `cgv_maj` and `certify.py`),
`lifted` runs `bounds` on the lift and `primes_needed` compares every bound with half the
product of the primes, adding passes until it is below: `bounds` is a `pass` over `BoundRing`
whose extraction `source` makes every curve pass on what the candidates imply,
$\deg(C) \sum_{n \mid C} \lvert\mathrm{GV}_{C/n}\rvert / n^3$ ($n^2$ and no degree for n-folds), instead of
its own bound, and whose exponentials go along the rows of the charges ($c_1$ along the rows
from `period`, $\beta$ from `instanton`, `rows` in `extract`), so a term vanishes with
$Q_r \cdot C$; the result matches `cgv_maj`'s curve set and is at or below its bounds. Vex
phases are refused (`Error::CertifyStrata`), as by the C. `bound_gv` exposes the bounds.
`compute_gw` derives the GW invariants from the lifted GV invariants
(`multicover::gw_from_gv`, in `malachite::Rational`), since they are not integers. Both truncate by
`max_deg` and optionally by target points: `lightcone` holds the union of their backward
lightcones (the $C$ with $P - C$ in the Mori cone), tested by the Mori cone's hyperplanes with
cone data and by membership in the generated semigroup otherwise (which is then enumerated up
to the targets' degree, matching the reference), and applied at the C's three points: the
candidates, the curves the extraction subtracts, and each pass's output (GW too, since the
multiples of a curve leave the lightcones). `Geometry` in `src/modular.rs` holds the integer data
the stages read: the CY kind, the charges, one anticanonical row per part of the nef
partition, and $M$ bilinear forms, one for a threefold ($\sum_t w_t \kappa_{tab}$) and one per
reference surface for an n-fold, which give $I$ $M$ components. In a lattice the charges have
two sides: `q`/`q0` on the curves, in its coordinates (`dim`), which the candidates, the
period's intersections and the extraction's rows read, and `divisor_q`/`divisor_q0` in the
$h^{1,1}$ coordinates of the forms, which the period's derivatives, $c_1$ and the instanton
keep. `src/modular/reference_tests.rs`
compares $c_0$, $c_1$, $S_2$ and $I$ with the reference stages modulo the primes, curve by curve
(a mismatch there pins a bug to a stage), and the GV invariants of `compute_gv` with those of
`run_hkty` and the GW invariants of the two methods, which only agree on consistent geometries
(not the toy threefold of `examples/` and the tests, where the methods give different results).

The computation is a four-stage pipeline, orchestrated by `reference::run_hkty`
(`src/reference.rs`), whose stages are submodules of `reference`:

1. **`semigroup`** — builds the truncated affine semigroup $S_\sigma = \sigma \cap \mathbb{Z}^n$ of
   curve classes in the Mori cone. Four constructors select which points get computed:
   `from_data`, `with_max_degree`, `with_min_elements`, `with_target_points`. Elements are stored
   as **columns** of a `DMatrix<i32>`, sorted by degree under the grading vector; downstream code
   relies on that sort order and on element 0 being the identity.
2. **`fundamental_period`** — computes $\omega$, $\omega^{-1}$, and first/second derivatives of its
   coefficients (`FundamentalPeriod { c0, c1, c2, c0_inv, vex }`). Curves are bucketed by how many
   negative intersections they have with the GLSM charge matrix `q` (`group_by_neg_int`), and each
   bucket has its own coefficient routine (`compute_c_0neg` / `_1neg` / `_2neg`). Curves of
   negative anticanonical degree $m$ stay out of those: in a vex phase of the toric variety
   (`Stratum` in `src/lib.rs`, the `strata` input; threefold hypersurfaces only) one that meets
   exactly the three divisors of a stratum negatively gets the stratum term of `compute_c_vex`
   in `vex[t]`, every other one contributes nothing, and $m < 0$ with fewer than two negative
   intersections is an error. The formula comes from the C program of PR #94.
3. **`instanton`** — derives the alpha, beta, F, instanton, and exp(alpha) polynomials from the
   fundamental period, returning `InstantonData { inst, expalpha }`. The stratum terms enter
   instanton correction $t$ as $c_0^{-1}$ `vex[t]`.
4. **`series_inversion`** — `invert_series` walks degrees in increasing order, extracting GV/GW
   invariants level by level; it keeps a sliding window of `previous_qn` levels whose size scales
   with $h^{1,1}$.

`misc::process_int_nums` preprocesses the intersection-number dict into the canonicalized dict, the
set of relevant index pairs, and the index count. Threefolds and higher-dimensional CYs canonicalize
indices differently (the first index is a reference surface in the n-fold case).

### Key cross-cutting patterns

- **One problem, one dispatch.** The public input is `src/problem.rs` (decided in
  `INTERFACE.md`, phase 9): a `Problem` of a `Geometry` (`q` with one row per divisor,
  `intnums`, `nefpart` as divisor indices, and `ambient`: `Plain`, `Fan(cones)` or
  `Strata(list)`), a `MoriCone` (`Generators`, `Inequalities { hyperplanes, equations }` or
  `FromFan`), a `Selection` (`MaxDeg`, `TargetPoints { points, max_deg }`, `MinPoints`,
  `CurveClasses`) and an optional `grading_vector`, plus `Options` (`method`, `n_threads`,
  `certify`, `low_memory`, `hilbert_bases`, and `tuning: modular::Tuning`). Every field has the name of the
  Python keyword and YAML key that set it, and every choice is an enum, so conflicting input
  cannot be built. `compute_gv`, `compute_gw` and the runtime-kind `compute` (returning
  `Computed::Gv` or `Computed::Gw`, for Python and the CLI) in `src/dispatch.rs` (a private
  module, re-exported at the root) take `&Problem` and `&Options` and return `Invariants<T>`:
  the `dimension`, the `grading_vector` used, the nonzero `entries` (`Invariant { curve,
  surface: Option<usize>, degree, value }`, sorted by degree, curve, surface) and a `Report`
  (`#[non_exhaustive]`: the method that ran, `certified`, `grading_vector_derived`,
  `elapsed`). `Problem::validate` (also run by `compute`, through `Problem::check`) checks
  everything that is cheap before computing — `check_nefpart`, the dimension (with a hint
  when `q` looks transposed), `check_intnums`, `check_strata`, the fan (`Fan::new`, its
  intersection ring, its strata), the shapes of the cone and the selection, a given grading
  vector against the fan and the generators, the target degrees, the Hilbert bases when the
  cone is given by inequalities, `method: modular` with a selection it does not support, and
  certify in a vex phase — and returns `Error::Input { field, reason }` naming the field
  (`Error::Unsupported`, `Error::MoriCone`, `Error::Modular` otherwise). `compute` then
  resolves the cone's inequalities (given, `cone::facets` of the generators, or
  `Fan::mori_cone`) and the default grading vector (the primitive sum of those inequalities),
  and picks the pipeline: `Auto` runs `modular` for `MaxDeg` and `TargetPoints` and the
  reference otherwise; curve classes too wide for the modular keys are `Error::Unsupported`
  naming `method: reference`, never a silent fallback. `modular` takes its candidates from
  `Options::hilbert_bases` (checked by `check_cone_data`) or from `cone::cone_data` of the
  inequalities, truncated at the maximum degree lowered to the targets'; the reference runs on
  `cone::hilbert_basis` of the inequalities, or for `CurveClasses` on those classes as given
  (checked to lie in the cone, zero added, duplicates dropped, `Semigroup::from_data`), and
  ignores the Hilbert bases. With `TargetPoints` and `max_deg` the reference's results are
  truncated by degree, since its semigroup only takes the targets. `run_hkty<T>` is generic over
  the coefficient type, which is only ever `reference::Rational` (`src/reference/rational.rs`)
  today, takes the `InvariantKind` and `CYKind` enums at runtime, and still panics on invalid
  input, which is why the checks come first. Adding an option means a field of `Options` (or
  `Problem`), its key in `src/io.rs` and its keyword in `src/python.rs` and `hkty.py`. The
  computation is always exact: floating-point arithmetic (the old `prec` option) was removed.
- **`PolynomialCoeff<T>`** (`src/reference/polynomial/coefficient.rs`) is the compound trait the
  pipeline is written against, so that the exact number type can be swapped. It is built almost
  entirely from `*Assign` traits (`Assign` itself is ours): coefficient arithmetic is done in
  place to avoid allocating bignums. When adding a numeric operation, add it as an in-place trait
  (see `RecipMut`, `IsInteger`) rather than a by-value one. `Default` must give zero.
  `reference::Rational` implements it as a newtype over `malachite::Rational`, since the orphan
  rule forbids implementing `AddAssign<i32>` and the like on the foreign type directly. The typed
  entry points return `malachite::Integer` and `malachite::Rational`, and the crate re-exports
  `malachite` so callers can name them.
- **`Polynomial<T>`** (`src/reference/polynomial.rs`) is sparse over the fixed monomial set of the
  semigroup: `coeffs: HashMap<monomial index, T>` plus a sorted `nonzero: Vec<usize>` index list.
  Both must be kept consistent — code that inserts into `coeffs` directly (e.g. when assembling
  from worker threads) must rebuild and re-sort `nonzero` afterwards, and `clean_up` drops the
  terms that cancelled to zero. `PolynomialProperties` carries the semigroup reference, the
  monomial→index map, and a `zero` coefficient that new ones are cloned from.
- **Threading** uses rayon throughout stages 2–4. `run_hkty` builds a private `ThreadPool` with
  `n_threads` workers and `install`s the pipeline in it, so the stage functions take no thread
  count and simply run in whatever pool they are called from (the global one, when called
  directly as the unit tests do). Deliberately not the global pool: that would ignore `n_threads`
  and share workers with any rayon use of a library consumer. Results are mostly `collect`ed and
  then assembled on the calling thread — collecting a `Vec` from an indexed iterator keeps its
  order, so `alpha[t]` and friends stay put. Errors are propagated by collecting into a `Result`,
  which stops the remaining items early. The exception is `compute_omega`, which produces far
  more, smaller results than the other stages: collecting all of its coefficients before filing
  them into their polynomials raised its peak memory by half, so `file_while` streams them over
  a bounded `mpsc` channel to a filing thread instead. That thread is a plain `std` thread outside
  the pool on purpose; a pool worker blocked on the channel could starve the producers, and
  deadlocks outright with `n_threads = 1`. The bound makes the producers wait whenever the filing
  thread falls behind (a busy machine is enough), which keeps the peak at its floor instead of
  swinging by tens of MiB from run to run. Where a worker needs bignum scratch space it gets it from
  `for_each_init`, which builds it once per rayon job instead of once per item (see `CScratch` in
  `src/reference/fundamental_period.rs`).
- **Threading in `modular`** follows the same pattern: `modular::run` builds a private pool with
  `Options::threads` workers (the dispatch passes `n_threads`, which must be positive) and installs
  candidates, passes and lift in it. The candidates enumerate the cones in parallel, the period
  splits the candidates into chunks, and the products of `instanton` and the extraction use
  shards: a `table::ShardedMap` (or the extraction's own shards, which also carry a degree heap)
  splits a table by the high bits of the key hash, each shard behind a `Mutex`, and workers add
  to it in batches (a chunk's own map, or per-shard buffers of `FLUSH` terms), so each lock is
  taken rarely. The extraction's tree builds the curves' exponentials in parallel by curve,
  sums up the tree by parent, subtracts by table and multiplies by target table (its products
  into one table are one task, but for the tables of the root and its children, which are in
  one part per thread and take the products in pieces of their sources through buffers).
  Lists of tasks sorted largest first go through rayon with `with_max_len(1)`: left to itself
  it runs such a list in about twice as many runs as threads, each in order, and the first run
  takes the largest tasks one after the other. The curve-by-curve
  extraction applies the curves of a layer (degrees $[d, d + \delta)$, $\delta$
  the lowest degree of $\alpha$) in parallel, each worker indexed by
  `rayon::current_thread_index` with its own level tables; layers of few large curves (the first
  ones) instead go level-synchronously, with the scatter of each curve into each of its levels as
  a task (`extraction::Tuning`). Results do not depend on the thread count, since the arithmetic
  is exact and every table is sorted when taken. The certifier's bounds round, so their sums
  depend on their order, and for a `Ring` that is not `EXACT` the stages fix it: the
  instanton products use chunks that do not depend on the threads and merge them in order
  (`ShardedMap::add_in_order`), the exponentials sort their sources by hash, the tree sorts
  each layer by key before numbering its nodes, adds the curves' terms in their order (kept
  a group of about a million terms at a time, `KEEP_TERMS`) and folds the subtractions in the
  order of pieces of their tables, a batch at a time (`EMIT_PIECE`, `EMIT_BATCH`; each depth
  groups the bounds' sums differently, so the bounds depend on it), and the curve-by-curve
  extraction sorts each layer and applies it in groups of curves (`GROUP_TERMS`,
  `GROUP_CURVES`), whose terms record their curve and are kept until the group is done, when
  each shard takes its terms in the order of the curves (`Shard::fold_ordered`;
  `modular::tests::test_bounds_threads`). That costs a certified run 3–13% of its time and
  4–14% of its memory (measured again in `PLAN.md` phase 10, and kept); `tests/stress.rs` (ignored, run with
  `cargo test --release --test stress -- --ignored` before releases and after changes to the
  parallel code) checks that over repeated runs, and `reference_tests::test_threads` forces
  trees of several depths and both modes of the curve-by-curve extraction on a few threads in
  every test run. The tuning constants were measured on a
  4-core machine with hyperthreading, and checked on a 16-thread Zen 5 (`PLAN.md` 5.3, the
  measurement session).

### Command line interface

`src/bin/cygv/main.rs` is a thin clap front end: it reads YAML from a file or stdin, parses it with
`io::Input::load_all`, applies the overrides (`--threads`, `--method`, `--max-deg` through
`Input::set_max_deg`, `--invariants`, `--certify`/`--no-certify`, `--low-memory`) and validates every document
before computing any, and writes one YAML document of results per input document to a file or
stdout, flushing after each one. All of the parsing and formatting lives in `src/io.rs`, whose
`Input` holds a `Problem`, the `InvariantKind` and the `Options`, with the field names of the
Python API (each inner list is a *column*), in a flat document. `Input::from_yaml` rejects
unknown fields (with what replaced the removed `cones` and `is_threefold`), parses the one-of
groups (`fan`/`strata`; `generators`/`mori_hyperplanes`(+`mori_equations`)/nothing with a fan;
`max_deg`/`target_points`(+`max_deg`)/`min_points`/`curve_classes`) into their enums
(`IoError::Conflicting`, `IoError::MissingOneOf`), and calls `Problem::validate`, mapping
`Error::Input { field, .. }` to `IoError::InvalidField`; it does not check anything itself. The
output carries `dimension`, the `grading_vector` and a `report` line (method, certified,
grading_vector_derived; no timing, so outputs stay comparable), and `surface` for n-folds.
Example inputs live in `examples/*.yaml` and are exercised by `tests/cli.rs`.
`tests/data/h11_9_*.yaml` holds the inequalities and Hilbert bases computed with Normaliz for
`benches/data/h11_9.yaml` (the whole Mori cone, and a facet with its generators), which
`tests/cones.rs` and the modular pipeline's comparison tests append to that file in place of
its generators. `tests/data/refs.yaml` and `refs_vex.yaml` are the frozen fixtures (`PLAN.md`
5.2, layer 4): CLI inputs whose GV invariants neither pipeline computed (cygv 0.2.3 for
threefolds up to $h^{1,1} = 10$, the C program of PR #94 for vex phases and FRSTs, with the
FRST side of each pair mapped from its vex phase), with the expected results in the output
format in the matching `*_gv.yaml`, compared by `tests/fixtures.rs` (the reference method too,
on the cheap ones). The files say where each geometry comes from, and adding a fixture is a
document in each. `tests/data/` is excluded from the published crate.

The whole CLI sits behind the default `cli` feature: `clap` and `yaml-rust2` are optional
dependencies, `[[bin]]` declares `required-features = ["cli"]`, and both `io` (in `src/lib.rs`) and
`tests/cli.rs` are `#[cfg]`-gated on it. The wheel does not need any of it, so `[tool.maturin]` sets
`no-default-features = true` — which means the Python build never exercises `src/io.rs`, and
`cargo build --no-default-features` (run in CI) is what keeps the gate honest.

### Python layer

`src/python.rs` exposes a single `_compute(kind, args)` PyO3 function (module built with `gil_used =
false`) taking the keywords as one dict of plain lists and ints, which it parses into a `Problem`
and `Options` (the one-of groups as in the YAML layer) and passes to `compute`; it returns the
entries with their values as **strings**, the dimension, the grading vector and the report, and
the Python side parses the values into `int` or `fractions.Fraction`. A new field is one key there
and one keyword in `hkty.py`.

`python/cygv/hkty.py` wraps it: `Problem` is a dataclass of the geometry and the Mori cone (`q`,
`intnums` as a dict or `[i, j, k, value]` entries, `nefpart`, `fan`, `strata`, `generators`,
`mori_hyperplanes`, `mori_equations`, `grading_vector`) whose `gv` and `gw` methods take the
selection and the options as keywords; `compute_gv(q, intnums, **kwargs)` and `compute_gw` split
their keywords between the two. The result, `Invariants`, is a `list` of `(key, value)` pairs (the
key a curve tuple, or `(curve, surface)` for an n-fold) with `kind`, `dimension`,
`grading_vector` and a `Report` dataclass attached. The call runs in a `multiprocessing.Process`
so Ctrl-C interrupts the computation without killing the interpreter (the Rust side installs a
ctrlc handler that calls `exit(1)`). Exceptions raised in the subprocess keep their type: invalid
or unsupported input is a `ValueError`, a failed computation (`NotStable`, cones too large) a
`cygv.ComputationError` (a `RuntimeError` created in `src/python.rs`), and panics a
`RuntimeError`. The subprocess's stderr is redirected to a temp file and surfaced in the raised
`RuntimeError` when the process dies.

Note the transpose convention: the Rust API takes column-major matrices, while the Python API
takes row-oriented lists (`to_matrix` in `src/python.rs` treats each inner list as a column).

## Benchmarks

Both bench targets are driven by the same scenario list in `benches/common/mod.rs`: a model
crossed with one or both `Variant`s (GV or GW), at a couple of maximum degrees. A scenario id looks
like `fourfold/deg15/gw-rational` (the suffix predates the removal of the floating-point variants
and is kept so that saved criterion baselines still match), and any argument after `--` that is
not a flag filters on it. Adding a model or a degree means adding a row to `scenarios()`; nothing else has to
change. Each row also names a `Pipeline`: the reference (`Method::Reference`, explicitly, since
the default is `Auto`), or the modular one without or with the supplied cone data (the
inequalities as the cone, and the Hilbert bases as `Options::hilbert_bases`), or in its
low-memory mode without it. The modular
ids end in `/modular`, `/modular-cones` or `/modular-lowmem` (the reference ones have no suffix, so that older
baselines still match), so `cargo bench -- modular` runs only those. The modular scenarios go
to much higher degrees, mostly with the cone data. Those without it compute the cone data from
the generators (the facets, then the Hilbert bases); before phase 7b they enumerated the
semigroup of the generators, which took most of their time, so their numbers are not
comparable across that change.

There are two kinds of model. `threefold` ($h^{1,1} = 2$ hypersurface) and `fourfold` ($h^{1,1} = 6$
CICY) are small, hardcoded, and run both variants. `h11_8`, `h11_9`, `h11_10` and `h11_11`
(named after their $h^{1,1}$) are realistic threefold hypersurfaces loaded from
`benches/data/*.yaml`, which are ordinary CLI inputs (so the benches need the `cli` feature, and
each file also runs with `cygv --file`). Their generators are thousands of Mori cone lattice points,
the way CYTools passes them, so they also exercise semigroup construction at scale; that is what
catches something like an accidentally cubic `find_generators`. `h11_9_plike` is `h11_9` under a
sparse grading vector with large entries, which spreads the curve classes over hundreds of degree
levels and stresses the series inversion rather than the bignum arithmetic. The realistic models
only run `gv-rational`. Their data comes from an external project's test set and was cross-checked
against an independent HKTY implementation; the files record where each geometry comes from, and
some drop generators above the highest degree any scenario reaches. `benches/data/cones/` holds
their cone data, computed with Normaliz and truncated at the highest degree of the modular
scenarios; like `tests/data/`, it is excluded from the published crate. `h11_9_plike` has none,
since the truncation is by the other grading.

- `benches/bench.rs` measures wall clock with criterion. Each scenario carries a `sample_size` and
  an `expected_secs` that only sizes criterion's measurement window;
  criterion will suggest raising the target time, which is expected — the window is deliberately
  set just under one run per sample so slow cases are not sampled twice over.
- `benches/memory.rs` measures allocations. `cargo bench --bench memory` prints peak live bytes,
  total bytes allocated, allocation count, peak RSS and wall clock per scenario.

All of the memory goes through Rust's global allocator, the big numbers included, so
`benches/memory.rs` counts it with a `GlobalAlloc` wrapper. Its counters are shared atomics, so
its wall clock is inflated on allocation-heavy runs (twice criterion's on `threefold` at degree 50);
take timings from criterion. Each scenario runs in a child process
(`--run <id>`, spawned by the parent), so peaks do not leak between scenarios and `VmHWM` reflects
one scenario only.

Scenarios come in three `Tier`s. The high-degree ones are gated behind `CYGV_BENCH_HEAVY` because
criterion has to run each of them ten or more times; the memory target runs everything once, so
enabling them there is cheap. The `Huge` ones, behind `CYGV_BENCH_HUGE`, take minutes and several
GB each even for a single run, so only the memory target runs them; they are what to check when a
change is meant to reduce peak memory.
`CYGV_BENCH_THREADS` pins the worker count (default: one thread per core), which cuts the noise
when comparing two revisions.

## Docs

`src/lib.rs` includes `README.md` as the crate-level doc, and doc comments use `$...$` LaTeX
rendered by KaTeX via `docs-header.html` (wired up through `[package.metadata.docs.rs]`).
