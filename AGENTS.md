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

`[tool.uv] cache-keys` in `pyproject.toml` lists `src/**/*.rs`; without it uv would not rebuild the
editable extension when only Rust sources change, and `uv run pytest` would test a stale build.
`uv.lock` is gitignored — it pins only the dev environment and does not affect consumers.

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

## Docs

`src/lib.rs` includes `README.md` as the crate-level doc, and doc comments use `$...$` LaTeX
rendered by KaTeX via `docs-header.html` (wired up through `[package.metadata.docs.rs]`).
