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
structured `modular` pipeline sits beside it; see `PLAN.md`). `src/dispatch.rs` holds the public
entry points, which pick a pipeline by the `Method` enum of `src/lib.rs`, and `src/io.rs`, the
CLI and `src/python.rs` sit on top of those. `src/cone.rs` holds the optional *cone data*
(`ConeData`: the Mori cone's inequalities and equations, and the Hilbert basis of each cone
$K_T$, the curve classes in the Mori cone that only meet the divisors in $T$ negatively, for
every $T$ of one or two divisors and, in a vex phase, the three divisors of every stratum) and
its validation (`check_cone_data`), and the native cone computations of `PLAN.md` phase 7,
shared by both pipelines. Callers may still supply cone data, computed with an external
program such as Normaliz, and it is used as given. `src/cone/hull.rs` (phase 7a):
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
is sound because every cone is then inside the Mori cone. Errors are `ConeError`. Phase 7b
also tried enumerating the lattice points of each $K_T$ directly (`PLAN.md`); both matched
Normaliz, and the Hilbert bases won because they feed the existing cone path unchanged.
`src/phase.rs` (phase 7c, a port of `cgv/tools/cgv_phase.py`): `phase::Fan` takes the maximal
cones of the fan of the ambient toric fourfold and the charges, with the rays as the integer
Gale dual of the charges (`kernel_basis` of `src/cone/hilbert.rs`), and checks the paper's
assumptions as far as the fan shows them (simplicial, complete, a reflexive polytope of rays
none of which is interior to a facet, regular through a pointed Mori cone, no vex 2-cones;
that every boundary point not interior to a facet is a ray is not checked). `Fan::phase`
gives the Mori cone (the extreme rays of the cone of the wall classes, Cox–Little–Schenck
6.4.4), the default grading vector (the sum of its facets) and the strata (the 3-cones whose
rays share no facet of the polytope, from `cone::facets`, with their wall classes), and
`Fan::check_intnums` compares the intersection numbers with the toric intersection ring
(`Ring`, recursive by linear equivalence) as trilinear forms on $h^{1,1}$ independent
divisors. Exact `malachite::Rational` arithmetic on small systems; errors are `FanError`.
`tests/phase.rs` compares it with `cgv_phase.prepare` on the 20 fans of `cgv/tests/refs_vex`.

`src/modular.rs` and `src/modular/` hold the faster pipeline, which is still being extended and
tuned. It computes modulo word-sized primes and lifts the results with the CRT, and shares no
code with `reference`. Its primitives are `field` (Montgomery arithmetic,
`Lanes<NL>` residues, and the `Ring` trait the stages are written against, which takes the
primes as a context), `key` (curve classes packed into fixed-size keys that add like the
classes, with a layout sized per run from coordinate bounds) and `table` (open-addressing
tables and sparse series over those keys, each entry holding a run-time number of values).
A *pass* runs the stages modulo `NL` primes: `candidates` (from a `CandidateSource`: either
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
each exponential built by the Euler recurrence in per-degree tables) and `multicover` (the GV
invariants from those, by gcd of the coordinates; exponent 3 and a division by the degree for
threefolds, exponent 2 for n-folds). `compute_gv` in `src/modular.rs` repeats passes with the
next primes until `crt::Residues` gives a stable lift (the lift from all but the last
`check_primes` primes must agree with those), with `Options` for the lane count (2 to 4 primes
per pass; by default `choose_lanes` predicts it from a probe run at 3/4 of the degree, as the
C's driver does), a minimum number of primes, the thread count and the extraction's tuning;
big integers (`malachite::Integer`) are only built by `crt`, for the final values.
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
reference surface for an n-fold, which give $I$ $M$ components. `src/modular/reference_tests.rs`
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

- **One dispatch, several entry points.** Every entry point — the four
  `compute_g{v,w}_rat_{threefold,nfold}` wrappers in `src/dispatch.rs`, and
  `compute_gvgw_strings` in the same file, which is the runtime dispatch used by both
  `src/python.rs` and the CLI — goes through the private `dispatch::compute`, which picks the
  pipeline from the `Method`: `Auto` (the default) runs `modular` when the request has
  `max_deg` or `target_points` and no `min_points` (see `modular_unsupported`), and falls back
  to `reference` otherwise or when the curve classes do not fit `modular`'s keys; `Modular`
  makes an unsupported request a `ComputeError::Unsupported`; `Reference` always runs
  `run_hkty`. The entry points return `Result<_, ComputeError>`, but only `modular`'s errors
  and invalid strata, cone data and target points (checked up front by `check_strata`,
  `check_cone_data` and `check_target_points`) get there: `run_hkty` still panics on invalid
  input. `max_deg` may come with `target_points` (the reference then truncates its results by
  degree, since its semigroup only takes the targets). Both give only the nonzero invariants.
  Both pipelines use every lattice point of the Mori cone (decision 6.2, phases 7b and 7d):
  the Mori cone is described by the cone data's inequalities, or else by `cone::facets` of
  the generators; `modular` takes its candidates from the cone data, or from
  `cone::cone_data` of those inequalities (truncated at `max_deg`, lowered to the targets'
  largest degree); the reference runs on `cone::hilbert_basis` of them instead of the
  generators. The exception is a request with none of `max_deg`, `min_points` and
  `target_points`, where the reference's semigroup is the generators as given
  (`Semigroup::from_data`), so they are passed unchanged. With cone data the generators may
  have no columns, except for such a request (`ComputeError::MissingGenerators`). Failures of
  the cone computations are `ComputeError::MoriCone`. The optional `fan` (threefold
  hypersurfaces only) is turned into a `phase::Phase` first: its strata replace the `strata`
  argument, which must then be empty, its Mori cone spans the cone when the generators have
  no columns, and its grading vector is used when `grading_vector` is empty
  (`ComputeError::MissingGradingVector` without a fan); the intersection numbers are checked
  against it, and its errors are `ComputeError::InvalidFan`.
  `run_hkty<T>` is generic over the coefficient type, which is only ever `reference::Rational`
  (`src/reference/rational.rs`) today, and takes the `InvariantKind` and `CYKind` enums at runtime. Adding an option means touching
  all of the dispatch sites. The computation is always exact: floating-point arithmetic (the old
  `prec` option) was removed.
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
  `Options::threads` workers (the dispatch passes `n_threads`, zero meaning one) and installs
  candidates, passes and lift in it. The candidates enumerate the cones in parallel, the period
  splits the candidates into chunks, and the products of `instanton` and the extraction use
  shards: a `table::ShardedMap` (or the extraction's own shards, which also carry a degree heap)
  splits a table by the high bits of the key hash, each shard behind a `Mutex`, and workers add
  to it in batches (a chunk's own map, or per-shard buffers of `FLUSH` terms), so each lock is
  taken rarely. The extraction applies the curves of a layer (degrees $[d, d + \delta)$, $\delta$
  the lowest degree of $\alpha$) in parallel, each worker indexed by
  `rayon::current_thread_index` with its own level tables; layers of few large curves (the first
  ones) instead go level-synchronously, with the scatter of each curve into each of its levels as
  a task (`extraction::Tuning`). Results do not depend on the thread count, since the arithmetic
  is exact and every table is sorted when taken; `tests/stress.rs` (ignored, run with
  `cargo test --release --test stress -- --ignored` before releases and after changes to the
  parallel code) checks that over repeated runs, and `reference_tests::test_threads` forces both
  extraction modes on a few threads in every test run. The tuning constants were measured on a
  4-core machine with hyperthreading and are provisional.

### Command line interface

`src/bin/cygv/main.rs` is a thin clap front end: it reads YAML from a file or stdin, parses it with
`io::Input::load_all`, and writes one YAML document of results per input document to a file or
stdout, flushing after each one. `--threads` and `--method` override the `n_threads` and
`method` fields of every document. All of the parsing, validation, and formatting lives in
`src/io.rs`, whose `Input` struct mirrors the arguments of `run_hkty` and whose field conventions
match the Python API (each inner list is a *column*). `Input::from_yaml` rejects unknown fields and
validates shapes eagerly — including a `process_int_nums` dry run, `method: modular` without
`max_deg` or `target_points`, the cone data, the target points, the fan (from which it derives the
grading vector when `grading_vector` is left out, since the output needs it), and leaving out
`generators` without cone data or a fan, or without any of `max_deg`, `min_points` and
`target_points` — because
`run_hkty` itself still unwraps its errors. `Input::compute` sorts the results by degree and then lexicographically by
curve class, since the order in which curve classes of equal degree come out of `run_hkty` is not
deterministic. Example inputs live in `examples/*.yaml` and are exercised by `tests/cli.rs`.
`tests/data/*.yaml` holds cone data computed with Normaliz for `benches/data/h11_9.yaml` (the
whole Mori cone, and a facet with its generators), which `tests/cones.rs` and the modular
pipeline's comparison tests append to that file; it is excluded from the published crate.

The whole CLI sits behind the default `cli` feature: `clap` and `yaml-rust2` are optional
dependencies, `[[bin]]` declares `required-features = ["cli"]`, and both `io` (in `src/lib.rs`) and
`tests/cli.rs` are `#[cfg]`-gated on it. The wheel does not need any of it, so `[tool.maturin]` sets
`no-default-features = true` — which means the Python build never exercises `src/io.rs`, and
`cargo build --no-default-features` (run in CI) is what keeps the gate honest.

### Python layer

`src/python.rs` exposes a single `_compute_gvgw` PyO3 function (module built with `gil_used =
false`) that forwards to `dispatch::compute_gvgw_strings`; results are returned as **strings** and
parsed on the Python side into `int` or `fractions.Fraction` depending on `find_gv`.

`python/cygv/hkty.py` wraps it: `compute_gv` / `compute_gw` normalize inputs to numpy arrays, infer
`is_threefold` from the shape of `q` and the nef partition, and run the Rust call in a
`multiprocessing.Process` so Ctrl-C interrupts the computation without killing the interpreter
(the Rust side installs a ctrlc handler that calls `exit(1)`). They take `method` (a name, or
`None` for the Rust default) and `n_threads`, and `generators` and `grading_vector` may be
`None` (with cone data or a `fan`; `None` becomes an empty list). Exceptions raised in the subprocess keep their type
(an invalid `method` or a `ComputeError` is a `ValueError`), except panics, which become a
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
the default is `Auto`), or the modular one without or with the supplied cone data. The modular
ids end in `/modular` or `/modular-cones` (the reference ones have no suffix, so that older
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
