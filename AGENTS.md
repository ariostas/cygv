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

`rug` depends on GMP, which is compiled from source, so the first build is slow and needs a C
toolchain plus `m4` and `make` (`configure: error: No usable m4 in $PATH` means `m4` is missing).
CI sets `CC=clang` on all platforms; locally the default compiler is usually fine, and forcing
`CC=clang` fails outright if clang is not installed. Windows builds go through MSYS2/MINGW64.

That source build dominates CI, so both CI workflows cache it. (The timings below date from when
MPFR was built too, before `rug`'s float feature was dropped, so GMP alone takes less. Back then
it took ~13 min on Windows, where it used to be more than half the job.) `gmp-mpfr-sys` keeps the built libraries in the directory named by
`GMP_MPFR_SYS_CACHE`, under a subdirectory keyed by version, target triple, `CC` and `CFLAGS`, and
copies them back instead of rebuilding — about 3.5 MB, and a ~240 s build becomes ~8 s. The
workflows point that variable at `$RUNNER_TEMP/gmp-mpfr-sys` and restore it with `actions/cache`.
Two consequences worth remembering: `CC` must be set identically for every cargo invocation in a
job (a step that drops it looks under a different `CC-*` subdirectory and rebuilds from source),
and the cache does not depend on the cargo profile, so the whole Python matrix and the Rust
workflow share one key per OS — which is also why the key hashes `Cargo.lock` rather than anything
build-specific.

Both workflows separately cache `~/.cargo/registry` and `target/` with a hand-rolled
`actions/cache` step rather than `Swatinem/rust-cache`, because that action runs `cargo` from node
instead of through the msys2 shell and would therefore key the Windows cache on the runner's
preinstalled MSVC toolchain rather than the MINGW64 one doing the build. Two things that action
handled for us have to be done explicitly. The toolchain version goes in the key (and in the
restore-key prefix): a compiler upgrade invalidates every artifact in `target/`, and because
`actions/cache` skips saving whenever the primary key hits exactly, a stale entry would otherwise
never be replaced. And the two workflows use different key prefixes — `cargo-dev-` and
`cargo-release-` — because the Rust workflow builds the dev profile and the Python workflow builds
release; sharing a key would leave whichever job saved second restoring a `target/` with the wrong
profile in it. Nothing prunes these entries, so they are considerably larger than the GMP one.

The Windows Python job additionally installs with `--no-build-isolation`. maturin — the only entry
in `build-system.requires` — has no MINGW64 wheel on PyPI, so pip's isolated build environment
compiled it from source on every run, five minutes of work to duplicate the
`mingw-w64-x86_64-python-maturin` that pacman had already installed. `deploy.yml` skips isolation
on Windows for the same reason. The catch is that pip no longer enforces `build-system.requires`
there, so the version pacman ships has to keep satisfying `>=1.5,<2.0` on its own.

Both Windows jobs also set `PYO3_USE_RAW_DYLIB=0`. pyo3 0.29 started linking the Python DLL with
rustc's `raw-dylib`, which writes the DLL name straight into the import table rather than taking it
from an import library, and it takes that name from the interpreter configuration. For abi3 builds
against pyo3 0.29 or newer, maturin stops letting pyo3 probe the interpreter and writes a
`PYO3_CONFIG_FILE` of its own, which carries neither a `lib_name` nor any hint that the interpreter
is a MinGW one — so pyo3 falls back to the MSVC name `python3.dll`, while MSYS2 ships
`libpython3.dll`. Nothing fails while building, because raw-dylib needs no import library; the
extension simply fails to load, with `ImportError: DLL load failed while importing cygv: The
specified module could not be found`. Opting out restores what pyo3 0.28 did: link `-lpython3`, let
the MinGW linker resolve it to `/mingw64/lib/libpython3.dll.a`, and import `libpython3.dll`. The
variable is only read by pyo3's build script on Windows, so setting it job-wide is harmless
elsewhere. It can go once maturin records the MinGW library name in the config file it generates.

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
CLI and `src/python.rs` sit on top of those.

`src/modular.rs` and `src/modular/` hold the faster pipeline, which is still being extended and
tuned. It computes modulo word-sized primes and lifts the results with the CRT, and shares no
code with `reference`. Its primitives are `field` (Montgomery arithmetic,
`Lanes<NL>` residues, and the `Ring` trait the stages are written against, which takes the
primes as a context), `key` (curve classes packed into fixed-size keys that add like the
classes, with a layout sized per run from coordinate bounds) and `table` (open-addressing
tables and sparse series over those keys, each entry holding a run-time number of values).
A *pass* runs the stages modulo `NL` primes: `candidates` (for now the whole semigroup up to
`max_deg`, enumerated in order of degree, filtered to the classes with at most two negative
intersections), `period` (tables of factorials and harmonic numbers, then $c_0$, $c_1$ and the
second derivatives already contracted with the bilinear forms, $S_2$), `instanton` (sparse
products and reciprocal, $\alpha$, and the instanton series $I$), `extraction` (the
coefficients of $I$ in the basis $z^C \exp(C\cdot\alpha)$, curve by curve in order of degree,
each exponential built by the Euler recurrence in per-degree tables) and `multicover` (the GV
invariants from those, by gcd of the coordinates; exponent 3 and a division by the degree for
threefolds, exponent 2 for n-folds). `compute_gv` in `src/modular.rs` repeats passes with the
next primes until `crt::Residues` gives a stable lift (the lift from all but the last
`check_primes` primes must agree with those), with `Options` for the lane count and a minimum
number of primes; big integers come from `num-bigint`, a provisional choice confined to `crt`.
`compute_gw` derives the GW invariants from the lifted GV invariants
(`multicover::gw_from_gv`, in `num-rational`), since they are not integers. It only supports
truncation by `max_deg`, and everything is single-threaded so far. `Geometry` in `src/modular.rs` holds the integer data
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

- **One dispatch, several entry points.** Every entry point — the four
  `compute_g{v,w}_rat_{threefold,nfold}` wrappers in `src/dispatch.rs`, and
  `compute_gvgw_strings` in the same file, which is the runtime dispatch used by both
  `src/python.rs` and the CLI — goes through the private `dispatch::compute`, which picks the
  pipeline from the `Method`: `Auto` (the default) runs `modular` when the request has
  `max_deg` and no `min_points` or `target_points` (see `modular_unsupported`), and falls back
  to `reference` otherwise or when the curve classes do not fit `modular`'s keys; `Modular`
  makes an unsupported request a `ComputeError::Unsupported`; `Reference` always runs
  `run_hkty`. The entry points return `Result<_, ComputeError>`, but only `modular`'s errors
  get there: `run_hkty` still panics on invalid input. Both give only the nonzero invariants.
  `run_hkty<T>` is generic over the coefficient type, which is only ever `rug::Rational` today,
  and takes the `InvariantKind` and `CYKind` enums at runtime. Adding an option means touching
  all of the dispatch sites. The computation is always exact: floating-point arithmetic (the old
  `prec` option) was removed.
- **`PolynomialCoeff<T>`** (`src/reference/polynomial/coefficient.rs`) is the compound trait the
  pipeline is written against, so that the exact number type can be swapped (`PLAN.md` plans a
  pure-Rust one next to `rug::Rational`). It is built almost entirely from `*Assign` traits:
  coefficient arithmetic is done in place to avoid allocating bignums. When adding a numeric
  operation, add it as an in-place trait (see `RecipMut`, `IsInteger`) rather than a by-value
  one. `Default` must give zero.
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

### Command line interface

`src/bin/cygv/main.rs` is a thin clap front end: it reads YAML from a file or stdin, parses it with
`io::Input::load_all`, and writes one YAML document of results per input document to a file or
stdout, flushing after each one. `--threads` and `--method` override the `n_threads` and
`method` fields of every document. All of the parsing, validation, and formatting lives in
`src/io.rs`, whose `Input` struct mirrors the arguments of `run_hkty` and whose field conventions
match the Python API (each inner list is a *column*). `Input::from_yaml` rejects unknown fields and
validates shapes eagerly — including a `process_int_nums` dry run, and `method: modular` without
`max_deg` — because `run_hkty` itself still unwraps its errors. `Input::compute` sorts the results by degree and then lexicographically by
curve class, since the order in which curve classes of equal degree come out of `run_hkty` is not
deterministic. Example inputs live in `examples/*.yaml` and are exercised by `tests/cli.rs`.

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
`None` for the Rust default) and `n_threads`. Exceptions raised in the subprocess keep their type
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
change. Every scenario runs the reference pipeline (`Method::Reference`, explicitly, since the
default is `Auto`); the modular pipeline has no benchmarks yet.

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
some drop generators above the highest degree any scenario reaches.

- `benches/bench.rs` measures wall clock with criterion. Each scenario carries a `sample_size` and
  an `expected_secs` that only sizes criterion's measurement window;
  criterion will suggest raising the target time, which is expected — the window is deliberately
  set just under one run per sample so slow cases are not sampled twice over.
- `benches/memory.rs` measures allocations. `cargo bench --bench memory` prints peak live bytes,
  total bytes allocated, allocation count, peak RSS and wall clock per scenario.

Memory needs care because the bulk of it lives in GMP numbers, which C code allocates and Rust's
global allocator therefore never sees. `benches/memory.rs` installs its own allocation functions
through `gmp_mpfr_sys::gmp::set_memory_functions` plus a `GlobalAlloc` wrapper, so
both sides feed the same counters. GMP passes the block size back on free and realloc, so no
bookkeeping headers are needed. The hooks must be installed before any GMP object exists, which is
why it happens first thing in `main`. Each scenario then runs in a child process (`--run <id>`,
spawned by the parent), so peaks do not leak between scenarios and `VmHWM` reflects one scenario
only. The `gmp-mpfr-sys` dev-dependency must stay compatible with the version `rug` links against,
and must request exactly the features `rug` does: `default-features = false` and none at all,
since `rug` is used with only its `rational` and `std` features. Cargo fingerprints feature
*names*, not what they expand to, so a different list here splits it into a second build unit and
compiles GMP from source a second time, doubling every cold CI build. Turning `rug`'s `float`
feature back on would bring MPFR back, and this list would then have to name `mpfr` too.

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
