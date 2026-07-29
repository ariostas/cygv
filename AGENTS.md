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
cargo test                            # unit tests in src/**, integration test in tests/full_hkty.rs
cargo test test_threefold             # single test by name
cargo test --test full_hkty           # single integration test file
cargo doc --open

# Benchmarks (see the Benchmarks section below)
cargo bench                           # both targets, quick scenarios only
cargo bench --bench bench             # wall clock (criterion)
cargo bench --bench memory            # allocations and peak memory
cargo bench -- fourfold               # only the scenarios whose id contains "fourfold"
CYGV_BENCH_HEAVY=1 cargo bench        # add the high-degree scenarios
CYGV_BENCH_THREADS=1 cargo bench      # pin the worker thread count

# Python (builds the Rust extension via maturin)
uv sync --extra tests                 # creates .venv and builds the extension
uv run pytest
uv run pytest python/tests/test_hkty.py::test_threefold

# Lint (runs rustfmt, cargo-check, clippy, ruff, mypy, codespell)
prek -a --quiet
```

`rug` depends on GMP/MPFR, which are compiled from source, so the first build is slow and needs a C
toolchain plus `m4` and `make` (`configure: error: No usable m4 in $PATH` means `m4` is missing).
CI sets `CC=clang` on all platforms; locally the default compiler is usually fine, and forcing
`CC=clang` fails outright if clang is not installed. Windows builds go through MSYS2/MINGW64.

That source build dominates CI (~13 min on Windows, where it used to be more than half the job),
so both CI workflows cache it. `gmp-mpfr-sys` keeps the built libraries in the directory named by
`GMP_MPFR_SYS_CACHE`, under a subdirectory keyed by version, target triple, `CC` and `CFLAGS`, and
copies them back instead of rebuilding — about 3.5 MB, and a ~240 s build becomes ~8 s. The
workflows point that variable at `$RUNNER_TEMP/gmp-mpfr-sys` and restore it with `actions/cache`.
Two consequences worth remembering: `CC` must be set identically for every cargo invocation in a
job (a step that drops it looks under a different `CC-*` subdirectory and rebuilds from source),
and the cache does not depend on the cargo profile, so the whole Python matrix and the Rust
workflow share one key per OS — which is also why the key hashes `Cargo.lock` rather than anything
build-specific. The Rust workflow additionally uses `Swatinem/rust-cache` for the
registry and `target/`, but not on Windows — that action invokes `cargo` outside the msys2 shell,
where the MINGW64 toolchain is not on `PATH`.

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

PR titles must follow conventional commits (enforced by `.github/workflows/pr_title.yml`).

## Architecture

The computation is a four-stage pipeline, orchestrated by `hkty::run_hkty` (`src/hkty.rs`):

1. **`semigroup`** — builds the truncated affine semigroup $S_\sigma = \sigma \cap \mathbb{Z}^n$ of
   curve classes in the Mori cone. Four constructors select which points get computed:
   `from_data`, `with_max_degree`, `with_min_elements`, `with_target_points`. Elements are stored
   as **columns** of a `DMatrix<i32>`, sorted by degree under the grading vector; downstream code
   relies on that sort order and on element 0 being the identity.
2. **`fundamental_period`** — computes $\omega$, $\omega^{-1}$, and first/second derivatives of its
   coefficients (`FundamentalPeriod { c0, c1, c2, c0_inv }`). Curves are bucketed by how many
   negative intersections they have with the GLSM charge matrix `q` (`group_by_neg_int`), and each
   bucket has its own coefficient routine (`compute_c_0neg` / `_1neg` / `_2neg`).
3. **`instanton`** — derives the alpha, beta, F, instanton, and exp(alpha) polynomials from the
   fundamental period, returning `InstantonData { inst, expalpha }`.
4. **`series_inversion`** — `invert_series` walks degrees in increasing order, extracting GV/GW
   invariants level by level; it keeps a sliding window of `previous_qn` levels whose size scales
   with $h^{1,1}$.

`misc::process_int_nums` preprocesses the intersection-number dict into the canonicalized dict, the
set of relevant index pairs, and the index count. Threefolds and higher-dimensional CYs canonicalize
indices differently (the first index is a reference surface in the n-fold case).

### Key cross-cutting patterns

- **Const generics select variants at compile time.** `run_hkty<T, FIND_GV, IS_THREEFOLD>` is
  monomorphized over `T ∈ {Rational, Float}` and the two bools. Every entry point — the eight
  `compute_g{v,w}_{rat,float}_{threefold,nfold}` wrappers in `src/hkty.rs` and the `match (find_gv,
  is_threefold)` in `src/python.rs` — is just a dispatch into that one generic. Adding an option
  here means touching all of those dispatch sites.
- **`PolynomialCoeff<T>`** (`src/polynomial/coefficient.rs`) is the compound trait that lets the
  same code run over exact `rug::Rational` and arbitrary-precision `rug::Float`. It is built almost
  entirely from `*Assign` traits: coefficient arithmetic is done in place to avoid allocating
  bignums. When adding a numeric operation, add it as an in-place trait (see `RecipMut`,
  `RoundMut`, `AbsMut`) rather than a by-value one.
- **`NumberPool<T>`** (`src/pool.rs`) recycles allocated bignums. Almost every function threads an
  `&mut NumberPool<T>`, and `Polynomial` has pool-aware `drop`/`clear`/`clone`/`move_into` that
  return coefficients to the pool instead of freeing them. Do not use `std::mem::drop` or `Clone`
  on polynomials where the pool-aware version exists. `run_hkty` builds
  `all_pools: (main_pool, Vec<per-thread pools>)` and passes it down by `&mut`.
- **`Polynomial<T>`** (`src/polynomial.rs`) is sparse over the fixed monomial set of the semigroup:
  `coeffs: HashMap<monomial index, T>` plus a sorted `nonzero: Vec<usize>` index list. Both must be
  kept consistent — code that inserts into `coeffs` directly (e.g. when assembling from worker
  threads) must rebuild and re-sort `nonzero` afterwards. `PolynomialProperties` carries the
  semigroup reference, the monomial→index map, and the `zero_cutoff` used by `clean_up` to drop
  numerically-zero terms.
- **Threading** uses one uniform pattern throughout stages 2–4: `Arc<Mutex<slice::Iter>>` as a work
  queue, `thread::scope` to spawn one worker per per-thread `NumberPool`, an `mpsc` channel back to
  the main thread which assembles results, and `drop(tx)` to terminate the receive loop. No rayon.

### Python layer

`src/python.rs` exposes a single `_compute_gvgw` PyO3 function (module built with `gil_used =
false`); results are returned as **strings** and parsed on the Python side into `int`,
`fractions.Fraction`, or `mpmath.mpf` depending on `find_gv`/`prec`.

`python/cygv/hkty.py` wraps it: `compute_gv` / `compute_gw` normalize inputs to numpy arrays, infer
`is_threefold` from the shape of `q` and the nef partition, and run the Rust call in a
`multiprocessing.Process` so Ctrl-C interrupts the computation without killing the interpreter
(the Rust side installs a ctrlc handler that calls `exit(1)`). The subprocess's stderr is redirected
to a temp file and surfaced in the raised `RuntimeError` when the process dies.

Note the transpose convention: the Rust API takes column-major matrices, while the Python API
takes row-oriented lists (`to_matrix` in `src/python.rs` treats each inner list as a column).

## Benchmarks

Both bench targets are driven by the same scenario list in `benches/common/mod.rs`: a model
(`threefold`, $h^{1,1} = 2$ hypersurface; `fourfold`, $h^{1,1} = 6$ CICY) crossed with the four
`Variant`s (GV/GW × rational/float), at a couple of maximum degrees. A scenario id looks like
`fourfold/deg15/gw-float`, and any argument after `--` that is not a flag filters on it. Adding a
model or a degree means adding a row to `scenarios()`; nothing else has to change.

- `benches/bench.rs` measures wall clock with criterion. Each scenario carries a `sample_size` and
  an `expected_secs` (measured on a 6-core machine) that only sizes criterion's measurement window;
  criterion will suggest raising the target time, which is expected — the window is deliberately
  set just under one run per sample so slow cases are not sampled twice over.
- `benches/memory.rs` measures allocations. `cargo bench --bench memory` prints peak live bytes,
  total bytes allocated, allocation count, peak RSS and wall clock per scenario.

Memory needs care because the bulk of it lives in GMP/MPFR numbers, which C code allocates and
Rust's global allocator therefore never sees. `benches/memory.rs` installs its own allocation
functions through `gmp_mpfr_sys::gmp::set_memory_functions` (MPFR routes through the same hooks
when built against GMP, which is how `gmp-mpfr-sys` builds it) plus a `GlobalAlloc` wrapper, so
both sides feed the same counters. GMP passes the block size back on free and realloc, so no
bookkeeping headers are needed. The hooks must be installed before any GMP object exists, which is
why it happens first thing in `main`. Each scenario then runs in a child process (`--run <id>`,
spawned by the parent), so peaks do not leak between scenarios and `VmHWM` reflects one scenario
only. The `gmp-mpfr-sys` dev-dependency must stay compatible with the version `rug` links against,
and must request exactly the features `rug` does (`default-features = false`, `mpc` and `mpfr`).
Cargo fingerprints feature *names*, not what they expand to, so enabling `default` here — which
resolves to the very same two features — splits it into a second build unit and compiles GMP and
MPFR from source a second time, doubling every cold CI build.

The high-degree scenarios are gated behind `CYGV_BENCH_HEAVY` because criterion has to run each of
them ten or more times; the memory target runs everything once, so enabling them there is cheap.
`CYGV_BENCH_THREADS` pins the worker count (default: one thread per core), which cuts the noise
when comparing two revisions.

## Docs

`src/lib.rs` includes `README.md` as the crate-level doc, and doc comments use `$...$` LaTeX
rendered by KaTeX via `docs-header.html` (wired up through `[package.metadata.docs.rs]`).
