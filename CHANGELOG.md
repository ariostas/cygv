# Changelog

All notable changes to this project are documented here. The project follows
[semantic versioning](https://semver.org/); before 1.0, a minor version can break the API.

## Unreleased

A new, much faster pipeline, a new interface around one problem type, checks of the input,
and support for vex phases. The interface changed in every layer; the tables below say what to
write instead of what 0.3.0 took.

### Added

- The modular pipeline (`method="modular"`, the default through `"auto"`), which computes
  modulo word-sized primes and lifts the results by the CRT. It is much faster than the
  original pipeline, now called the reference pipeline (`method="reference"`), and reaches
  $h^{1,1}$ in the hundreds on faces of the Mori cone. `certify` proves its lift from residues.
- Vex phases of the ambient toric variety, given by their `fan` or their `strata` (threefold
  hypersurfaces only). The fan also gives the Mori cone and the default grading vector, and
  the intersection numbers are checked against it.
- The Mori cone can be given by its inequalities (`mori_hyperplanes`, with `mori_equations`
  for a face), and the Hilbert bases of its cones $K_T$ can be given (`hilbert_bases`) or are
  computed natively. Every lattice point of the cone is used, even when the generators do not
  generate them all.
- `curve_classes` selects exactly the listed classes; `target_points` works with both methods.
- Consistency checks of the input: the modular pipeline checks its result with a second
  contraction of the intersection numbers (`verify`, `"quick"` by default), the reference
  pipeline checks every instanton correction, and cheap checks run before computing (the cone
  against the classes that meet every divisor nonnegatively, a vex phase without its fan or
  strata, a sample of given Hilbert bases, an optional `dimension`, degrees above $2^{20}$).
- A report with each result: the method that ran, how its lift is known to be right
  (`lift`: exact, stable or certified), which consistency check ran (`consistency`), whether
  the grading vector was derived, and the time.
- Progress output (`verbose` in Python, `-v`/`-vv` on the command line, the `log` facade in
  Rust), ctrl+c that interrupts the computation in Python, and `low_memory`.

### Changed

- The big integers are `malachite`'s instead of `rug`'s, so no C toolchain is needed to build.
- The Python computation runs in the calling process, on its own threads with the GIL
  released, instead of in a `multiprocessing` subprocess.
- GW invariants are always exact (`Fraction`s); floating-point arithmetic (`prec`) is gone, and
  with it the `mpmath` dependency.

### Migrating from 0.3.0

#### Python

| 0.3.0 | Now |
| --- | --- |
| `compute_gv(generators, grading_vector, q, intnums, max_deg=...)` | `compute_gv(q, intnums, generators=generators, grading_vector=grading_vector, max_deg=...)`; everything after `intnums` is a keyword, and the grading vector can be left out |
| no `max_deg`, `min_points` or `target_points`: compute at the generators | `curve_classes=generators` (zero is ignored) |
| `nefpart=...` (by position too) | `nefpart=...`, by keyword only, also in `Problem(q, intnums, nefpart=...)` |
| `prec=...` | gone: the computation is always exact |
| a `list` of `(key, value)` pairs | an `Invariants` sequence of the same pairs, with `invariants[key]`, `key in invariants`, `invariants.get(key)`, `dict(invariants)`, the numpy arrays `curves`, `degrees`, `values` (and `surfaces` for n-folds), and `dimension`, `grading_vector` and `report` |
| GW values as `Fraction`s (or `mpmath` floats with `prec`) | `Fraction`s |
| a `RuntimeError` for every failure | `ValueError` for invalid input, `cygv.ComputationError` for a failed computation (e.g. inconsistent input) |
| | `Problem(q, intnums, ...)` with `.gv(...)` and `.gw(...)` holds the geometry for several computations |

#### Command line (YAML documents)

| 0.3.0 | Now |
| --- | --- |
| no `max_deg`, `min_points` or `target_points`: compute at the generators | `curve_classes: <the generators>` |
| `is_threefold` (input and output) | gone from the input (the dimension is derived from `q` and `nefpart`; `dimension` checks it); the output has `dimension` |
| `prec` | gone |
| `grading_vector` required | optional: the primitive sum of the inequalities of the Mori cone by default |
| | new keys: `fan`, `strata`, `mori_hyperplanes`, `mori_equations`, `hilbert_bases`, `curve_classes`, `dimension`, `method`, `certify`, `verify`, `low_memory`; new flags `--method`, `--max-deg`, `--invariants`, `--certify`, `--verify`, `--low-memory`, `-v`; the output has a `report` line |

#### Rust

| 0.3.0 | Now |
| --- | --- |
| `cygv::run_hkty`, `compute_gvgw_strings` and the `compute_{gv,gw}_{rat,float}_{threefold,nfold}` functions | `cygv::compute_gv(&problem, &options)`, `compute_gw` or `compute` (a `Problem` of a `Geometry`, a `MoriCone` and a `Selection`, built with `Problem::new` and `Geometry::new`, and `Options`), returning `Invariants` with a `Report`; `cygv::modular::compute_gv` for the low-level modular pipeline |
| the stage modules `factorial`, `fundamental_period`, `instanton`, `misc`, `polynomial`, `semigroup`, `series_inversion` and `hkty` | private, inside `cygv::reference`; only their error types are public |
| `Semigroup`, `PolynomialCoeff`, `PrettyPrintPolynomial` | gone from the public API |
| `rug` integers and rationals | `malachite::Integer` and `malachite::Rational`, re-exported as `cygv::malachite` |
| `io::Input` with flat fields | `io::Input { name, problem, kind, options }` |
| | an `Interrupt` in `Options::interrupt` stops a computation from another thread |
| struct literals of the input types | `Problem`, `Geometry`, `Options` and the tuning types are `#[non_exhaustive]`: build them with their constructors or `Default` and set their fields |
