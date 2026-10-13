# The user-facing interface of `cygv`: review and redesign proposal

Status: design document for review, written against branch `gv-modular-rust` at `fbc6ea6`
(40 commits ahead of `main`, which is the released `v0.3.0`). Nothing here is implemented.
File and line references are to that revision.

The question asked was whether the interface "got a little convoluted" over the modular,
cone, vex-phase and certifier work, and how to make it as clear as possible for users. The
short answer is that the *computation* is now in good shape and most of the recent options
are each justified, but they were added one at a time to a flat argument list, so the
interface now exposes implementation structure (two pipelines, two Mori cone descriptions
plus an acceleration structure, three truncation modes with different pipeline support,
threefold/n-fold variants, fans vs strata) rather than the three things a user actually has:
a Calabi-Yau, a cone of curve classes, and a choice of which classes to compute.

Sections: 1 what exists today; 2 assessment (good, convoluted, per surface);
3 the redesign (concepts, then Rust, Python, YAML/CLI); 4 compatibility and migration;
5 prioritized list. Appendix A is the option interaction matrix, Appendix B the probe
transcript the examples come from.

---

## 1. The interface today

### 1.1 Surfaces

| Surface | Entry points | Where |
|---|---|---|
| Rust | `compute_gv_rat_threefold`, `compute_gw_rat_threefold`, `compute_gv_rat_nfold`, `compute_gw_rat_nfold` (14 positional args each), `compute_gvgw_strings` (16 args, string results), plus `reference::run_hkty` (12 args, panics on bad input), `modular::compute_gv/gw/bound_gv`, `cone::*`, `phase::Fan` | `src/dispatch.rs:551-709`, `src/dispatch.rs:491-549`, `src/reference.rs:24-40`, `src/modular.rs:479-539`, `src/lib.rs:229-258` |
| Python | `compute_gv`, `compute_gw` (4 positional + 12 keyword args) | `python/cygv/hkty.py:224-325`, `src/python.rs:27-51` |
| YAML/CLI | 19 fields per document, `--threads`, `--method`, `--certify` | `src/io.rs:107-127`, `src/bin/cygv/main.rs:81-100` |

### 1.2 Inputs, as a user sees them

Every argument, with the rule a user has to know. The "rule" column is the convolution.

| Argument | Meaning | Rules a user must know today |
|---|---|---|
| `q` | GLSM charges | Python/YAML: one vector per Kähler parameter (h11 vectors of ndiv entries). Rust: `DMatrix` with one *row* per divisor (`src/python.rs:11-25` transposes). The same matrix is described as having divisors as "rows of q" (`src/io.rs:150`, `src/lib.rs:75`) and as "columns of q" (`examples/fan.yaml:16`). |
| `intnums` | Intersection numbers | Python: `dict` only (a list of 4-tuples is a `TypeError` from PyO3, Appendix B); YAML: list of `[i,j,k,v]` *or* mapping `"i,j,k": v`. One ordering per triple; for n-folds the first index is a reference surface and only the last two are unordered (`src/lib.rs:154-185`). |
| `nefpart` | Nef partition | Empty = hypersurface. Typed `Sized | None` in Python (`hkty.py:232`). Decides the dimension together with `q`'s shape. |
| `is_threefold` | CY kind | Inferred from shapes in Python (`hkty.py:148-151`), overridable only in YAML (`src/io.rs:310-316`); the Rust API has four separate functions *and* a `CYKind` argument to the string entry point. |
| `generators` | Mori cone | Three meanings depending on other arguments: with any truncation, "any vectors spanning the cone" (saturated: every lattice point is used, `dispatch.rs:323-333`); without truncation, "exactly the curve classes to compute" (`dispatch.rs:316-321`, `untruncated_classes`); may be `None`/no columns only with cone data or a fan *and* a truncation (`dispatch.rs:256`, `:319`, `src/io.rs:370-374`). |
| `grading_vector` | Degree | Must be given unless a `fan` is (`ComputeError::MissingGradingVector`, `dispatch.rs:230`), although the Mori cone's facets would give one in every case (`PLAN.md` 5.3 follow-up). |
| `max_deg`, `min_points`, `target_points` | Which classes | Three modes. `max_deg`+`target_points` combine; `min_points` with either is rejected by the CLI (`src/io.rs:317-319`) but silently accepted by Python/Rust, where `run_hkty` then ignores `min_points` (`src/reference.rs:41-45`; Appendix B). `min_points` and "no truncation" force the reference pipeline (`dispatch.rs:124-136`). |
| `method` | Pipeline | `auto` / `modular` / `reference`. `auto` also falls back to the reference, silently, when the curve classes do not fit the modular keys (`dispatch.rs:310-313`), which is exactly the high-h11 case the reference cannot finish (`PLAN.md` 5.3, first item). |
| `n_threads` | Threads | `None` = one per core, `0` = one (`dispatch.rs:426-427`); `modular::Options::threads` is `Option<usize>` with a different doc (`src/modular.rs:437-439`). |
| `certify` | Proof | Only acts on the modular pipeline; silently ignored with `method="reference"` or when `auto` falls back (Appendix B). Refused in vex phases (`Error::CertifyStrata`). |
| `strata` | Vex phase | Threefold hypersurfaces only; Python `(divisors, class)` pairs with `Fraction`s; YAML mappings with `"-1/2"` strings; Rust `Vec<(i64, i64)>` numerator/denominator pairs (`src/lib.rs:73-80`) although the crate re-exports `malachite::Rational`. Cannot be given with `fan`. |
| `fan` | Toric variety | Threefold hypersurfaces only. Gives strata, a Mori cone (used only if `generators` is empty and there is no cone data) and a grading vector (used only if the given one is empty), and checks `intnums` and a given grading vector (`dispatch.rs:197-229`). |
| `mori_hyperplanes`, `mori_equations`, `cones` | "Cone data" | Three keyword arguments for one Rust struct (`ConeData`, `src/cone.rs:44-54`), all-or-nothing (`src/python.rs:93-113`). Bundles two different things: a *description of the Mori cone* by inequalities (an alternative to `generators`) and the *Hilbert bases of the K_T* (a precomputed acceleration structure that is otherwise computed internally, `dispatch.rs:271-294`). Trusted, not verified beyond membership. |
| `invariants` (YAML) / function name (Python, Rust) | GV or GW | One per call/document. |
| `name` (YAML) | Label | Echoed. |

### 1.3 Outputs

- Rust: `Vec<(DVector<i32>, Integer)>` for threefolds and `Vec<Labeled<Integer>>` with
  `Labeled<T> = ((DVector<i32>, usize), T)` for n-folds (`dispatch.rs:49`); strings from
  `compute_gvgw_strings`. No degree, no indication of which pipeline ran, how many primes,
  whether the result is certified.
- Python: list of `(tuple, int)` or `((tuple, surface), int)` pairs, "in no particular
  order" (`README.md:100`); `Fraction`s for GW. Order differs between pipelines (the
  modular one is sorted by degree, the reference is not), so equal inputs can print
  differently.
- YAML: sorted; echoes `invariants`, `is_threefold`, `grading_vector`; GW as quoted
  fractions (`src/io.rs:489-521`). Nothing about the method, primes or certification.

### 1.4 How the main consumer calls it

CYTools (`cytools/calabiyau.py:1991-2090`) calls `cygv.compute_gv(generators=...,
grading_vector=..., q=..., intnums=..., max_deg=..., min_points=..., target_points=...)`,
*all by keyword*. It builds `generators` as the Mori cone rays stacked with
`100 * h11` lattice points (`:2049-2050`), always computes a grading vector itself
(`:2054`), and converts the result with `dict(invariants)` or by iterating pairs
(`:2073-2076`). Three consequences on this branch:

1. `cy.compute_gvs()` with no truncation sends thousands of points as the exact class list
   to the *reference* pipeline (`dispatch.rs:316`): the slow path, with holes if the sample
   missed a class. `min_points=` does the same. Only `max_deg=`/`target_points=` reach the
   modular pipeline. CYTools does not know this.
2. The stacked lattice points are now pointless (saturation uses every lattice point of the
   cone the rays span) and only cost a `cone::facets` call over thousands of columns.
3. Because the call is by keyword, any redesign that keeps those seven keyword names (or
   aliases them) leaves CYTools working unchanged.

PR #102's Python interface (`cgv/tools/cgv_run.py:316-345`, `cgv/README.md`) is the
contrast: `compute_gvs(cy_or_dict, max_deg, grading_vec=None, lightcone=None, ...)`, where
the dict is `mori_rays` *or* `mori_hyperplanes`, `q`, `intnums`, optional `grading_vector`;
the grading vector defaults to the primitive sum of the facets (`cgv_run.py:177-189`), and
`cy_input(cy)` (`:270-283`) reads a CYTools object in four lines. `compute_gv_fan(cones, q,
intnums, max_deg, mori_rays=None, grading_vector=None, lightcone=None)` is the fan form.
It is small because the C does one thing (threefold hypersurfaces, GV, modular); but its
shape, a geometry plus one cone description plus one truncation with everything else
defaulted, is the right one and carries over.

---

## 2. Assessment

### 2.1 What is good and should be kept

- **One dispatch.** Every surface funnels into `dispatch::compute` (`dispatch.rs:170-355`),
  validation happens there for both pipelines, and `Method::Auto` means users do not have to
  know the pipelines exist for the common case. Keep that shape; the redesign only changes
  what goes in and comes out.
- **Saturation by default** (`PLAN.md` 6.2): "generators" means any vectors spanning the
  cone, every lattice point is used. This removed the most common silent-wrong-answer in
  the released version. Keep it, and make the one exception (the explicit class list)
  explicit rather than a mode of `generators`.
- **Up-front validation with reasons** (`check_nefpart`, `check_intnums`, `check_strata`,
  `check_cone_data`, `check_target_points`, the fan checks) and the `NON_INTEGER_CAUSES`
  hint (`src/lib.rs:189-192`). The messages are good; where they attach is the problem.
- **The fan as input.** Deriving strata, Mori cone and grading from the maximal cones and
  checking `intnums` against the intersection ring (`src/phase.rs`) is exactly the kind of
  default that removes options. It should become the *preferred* way to describe a toric
  hypersurface, not a special case bolted on the side.
- **Python conventions**: row-oriented lists / numpy array-likes, `Fraction`s for GW,
  subprocess isolation for ctrl-C, `ValueError` for bad input, panics as `RuntimeError`.
- **YAML**: unknown fields rejected, eager shape checks, both `intnums` spellings, stream
  of documents, flushed per document, examples exercised by `tests/cli.rs`.
- **`modular::Options` as a struct** (`src/modular.rs:423-446`), decision 4.1.6: the
  right pattern, applied so far only inside `modular`.
- **Exact results, `malachite` re-exported**, `Integer`/`Rational` as the typed results.

### 2.2 What is convoluted

#### Cross-cutting

**C1. Three overlapping ways to say what the Mori cone is, plus one that is not a cone
description but looks like one.** `generators` (spanning vectors), `mori_hyperplanes` +
`mori_equations` (inequalities), `fan` (wall classes), and `cones` (Hilbert bases of the
K_T, which only accelerate the modular candidates). The precedence is in code, not in
the user's head: cone data wins over the fan, which wins over generators
(`dispatch.rs:251-262`); generators may be omitted only with the first two *and* a
truncation; the fan's cone is computed only when needed (`:212-218`). PR #102 has one
field (`mori_rays` or `mori_hyperplanes`), and the fan form derives it. The `ConeData`
bundle is the root of the three-keyword Python API and of the "needs both
mori_hyperplanes and cones" error (Appendix B).

**C2. `generators` means two things.** With a truncation it is a cone description; without
one it is the exact class list, zero added, duplicates dropped, positive degree required
(`dispatch.rs:357-383`). The README needs a paragraph for it (`README.md:85-88`), the
error for a missing one needs a clause for it (`dispatch.rs:102-106`: "without cone data
or a fan, or without max_deg, min_points and target_points"), and the CLI help needs the
same sentence twice (`main.rs:62-63`, `:75-79`).

**C3. Truncation is three options with a support matrix.** `max_deg`, `min_points`,
`target_points`; `max_deg` combines with `target_points` but `min_points` combines with
nothing, which the CLI enforces (`src/io.rs:317`) and Python/Rust do not (Appendix B:
`min_points=5, max_deg=3` returns degree-3 results, `min_points` ignored by
`reference.rs:41-45`). Which pipeline can run depends on the choice
(`modular_unsupported`, `dispatch.rs:124-136`), which leaks into `method`'s documentation
on every surface.

**C4. `method` and `certify` interact invisibly.** `certify=True` with `method="reference"`
is accepted and ignored (Appendix B); with `auto` and `min_points` the same happens
silently. `auto` falls back to the reference on `KeyError::TooWide | BoundTooLarge`
(`dispatch.rs:310-313`), with no signal, for inputs where the reference "cannot finish"
(`PLAN.md` 5.3). The user asked for a fast, proven result and gets a slow, unproven one
with no message.

**C5. Threefold vs n-fold is encoded twice.** A runtime `CYKind` *and* separate
`_threefold`/`_nfold` functions in Rust; inferred and un-overridable in Python, inferred
but overridable in YAML; the result shape changes (`(class, v)` vs `((class, surface), v)`,
`README.md:101-103`) instead of carrying a surface field that is `None`/absent for
threefolds.

**C6. Orientation is explained in three incompatible frames.** Rust `q` has one row per
divisor (`phase.rs:145-146`), Python/YAML `q` is one vector per Kähler parameter
(`README.md:27`), and divisor indices are "rows of q" in `io.rs:150`/`lib.rs:75` but
"columns of q" in `examples/fan.yaml:16`. A transposed `q` fails with "the dimension of
the CY must be at least three, but it is -5" (Appendix B), which does not name the likely
cause.

**C7. Nothing comes back but numbers.** No degree (the CLI computes it for sorting,
`io.rs:428`, Python users recompute it), no method used, primes, passes, certification
status, candidate count, timings. Decision 4.1.7 ("reporting is data") was never built
(`PLAN.md` 5.3). The only way to learn that `auto` fell back is the wall clock.

**C8. Defaults that require expertise.** The grading vector is mandatory without a fan
although both pipelines already compute the facets it can be derived from
(`dispatch.rs:257-261`; `Fan::mori_cone` already applies the rule, `phase.rs:394`).
CYTools and PR #102 both compute one themselves (`calabiyau.py:2052`,
`cgv_run.py:177-189`), which is a sign the library should.

**C9. Documentation load.** `AGENTS.md`'s "Key cross-cutting patterns" bullet on the
dispatch is ~60 lines of rules; the README's Python section (`README.md:39-98`) is five
dense paragraphs of "X may be None unless Y and Z"; the CLI `--help` is 68 lines
(`main.rs:12-79`). Every one of those sentences is a rule the redesign should make
unnecessary rather than document better.

#### Rust (`src/lib.rs`, `src/dispatch.rs`, public modules)

**R1. Five near-identical entry points of 14-16 positional arguments**
(`dispatch.rs:551-709`), each `#[allow(clippy::too_many_arguments)]`, none with a doc
comment, differing only in two enum values that the string variant takes as arguments.
Adding one option touched five signatures plus `io::Input` plus `python.rs` plus
`hkty.py` twice; that is the mechanism by which the interface grew flat.

**R2. The `_rat_` infix is vestigial** (the float variants went with decision 6.9) and the
names say nothing about kind of input.

**R3. Sentinels instead of types.** Empty `RowDVector` means "no grading vector"
(`dispatch.rs:230`), zero columns means "no generators" (`:253`), `Option<DMatrix>` with
zero columns is also "no targets" (`:240`), `fan: Option<Vec<Vec<usize>>>` is a raw nested
vector, `Stratum::class` is `Vec<(i64, i64)>` fractions (`lib.rs:79`) although
`malachite::Rational` is re-exported for exactly this purpose (`lib.rs:252`), and
`Labeled<T>` is a nested tuple (`dispatch.rs:49`).

**R4. Two option vocabularies.** The dispatch takes `n_threads: Option<u32>`, `method`,
`certify` as positional arguments and builds a `modular::Options` from them
(`dispatch.rs:425-430`); `Options` itself (lanes, check primes, min primes, extraction
tuning) is unreachable from the entry points, the CLI and Python. Decision 4.1.6 says
options travel in one struct through the dispatch; they do not yet.

**R5. Validation is duplicated.** `io::Input::from_yaml` re-runs `check_nefpart`,
`check_intnums`, `check_strata`, `check_cone_data`, `check_target_points`, the fan checks
and the `Method::Modular` support check (`src/io.rs:317-396`) that `dispatch::compute`
runs again (`dispatch.rs:188-268`), because `ComputeError` carries `&'static str` reasons
without the field name, and the YAML layer wants field-tagged errors. `check_cone_data`
already returns `(field, reason)` for this reason (`cone.rs:84-89`), the others do not.

**R6. The public surface is the whole crate.** `pub mod` for all eight `reference`
submodules and twelve `modular` submodules (`reference.rs:3-10`, `modular.rs:25-36`),
root re-exports of `Polynomial`, `PolynomialProperties`, `PrettyPrintPolynomial`,
`Semigroup`, `PolynomialCoeff` (`lib.rs:229-239`; decision 7 left this open), and
`run_hkty` at the root, documented as `// TODO: Need to handle errors properly`
(`reference.rs:22`) and panicking on invalid input while sitting next to entry points that
promise not to. A `cargo doc` reader cannot tell the API from the internals.

**R7. Strings as the runtime-typed result** (`compute_gvgw_strings`, `dispatch.rs:491`).
It exists because the result type depends on `InvariantKind` at runtime; a two-variant
enum of typed results would serve Python and the CLI equally and keep the numbers typed to
the FFI boundary.

#### Python (`python/cygv/hkty.py`, `src/python.rs`)

**P1. Positional arguments that may be `None`.** `compute_gv(None, None, q, intnums,
max_deg=3, fan=fan)` (`README.md:69`) is the documented way to use a fan. The first two
positionals exist for compatibility with 0.3.0's `(generators, grading_vector, q,
intnums)`; everything that made them optional is unreleased.

**P2. The `ConeData` trio** `mori_hyperplanes=`, `mori_equations=`, `cones=` with
all-or-nothing coupling (`python.rs:93-113`), and `cones` as `(divisors, basis)` pairs
while `strata` is `(divisors, class)` pairs: two list-of-pairs conventions with different
second elements.

**P3. Inconsistent acceptance.** `intnums` must be a `dict` (list → `TypeError` from PyO3,
not `ValueError`); `target_points` may be 1-D (one point, `hkty.py:165`) but `generators`
may not (`hkty.py:195`); `strata` classes accept strings, `Fraction`s and ints (test at
`test_hkty.py:237`); `method` is a case-insensitive string (`lib.rs:222-226`) with no
enum or constants exposed.

**P4. All failures are `ValueError`** (`python.rs:141`), including "did not stabilize"
(`Error::NotStable`), which is a computation failure, and unknown-method, which is a
programming error. Callers cannot distinguish "fix your input" from "the geometry is
inconsistent" from "give it more primes".

**P5. Unordered results, shape depends on the kind** (C5, C7). CYTools has to compute
degrees itself; `test_hkty.py` sorts everything before comparing.

**P6. Everything is marshalled twice**: `hkty.py` regularizes to lists, pickles them across
a `multiprocessing` pipe, and `python.rs` rebuilds matrices column by column. Fine for
correctness, but each new argument is added in six places (`hkty.py:21-38`, `:79-96`,
`:103-123`, `:224-240`, `:276-292`, `python.rs:30-50`).

#### YAML / CLI (`src/io.rs`, `src/bin/cygv/main.rs`)

**Y1. Nineteen flat fields** (`io.rs:107-127`) with the exclusivity rules in prose
(`main.rs:60-79`). The one-of groups (how the Mori cone is given; how the classes are
chosen) are invisible in the schema.

**Y2. `is_threefold` in and out** (`io.rs:118`, `:500`): an implementation detail; the
dimension, or nothing, is what a reader wants.

**Y3. Output lacks provenance**: no method, primes, certification, candidate count
(C7). `surface_index` appears only for n-folds (`io.rs:512-514`).

**Y4. CLI overrides exist for `--threads`, `--method`, `--certify` but not for the two
things one actually varies between runs, the maximum degree and GV vs GW**
(`main.rs:84-100`). `--certify` is one-way (`input.certify |= args.certify`,
`main.rs:150`).

**Y5. Divergent validation.** `ConflictingTruncations` exists only here (`io/error.rs:17`);
the fan-derived grading is computed here *and* in the dispatch (`io.rs:320-355`,
`dispatch.rs:207-228`), the second a double-description in dimension h11.

---

## 3. Proposed redesign

### 3.1 Principles

1. **Three nouns, one verb.** A *geometry* (what the CY is), a *cone* (where the curve
   classes live), a *selection* (which classes to compute), and `compute`. Options are a
   fourth, separate noun with a complete default.
2. **One way to say each thing; one-of groups are enums.** The Mori cone is *one* enum
   value. The selection is *one* enum value. Mutual exclusivity becomes impossible to get
   wrong rather than documented.
3. **Derive, don't require.** Grading vector, CY kind, Mori cone (from a fan), strata (from
   a fan), Hilbert bases: all derived unless given. A given value is *checked* against the
   derived one where that is cheap.
4. **Accelerators are not inputs.** Precomputed Hilbert bases are an optional hint that
   never changes the answer, so they live in the options, not next to the cone.
5. **Report what happened.** Results carry degree, kind, grading vector, method used,
   primes, certification, counts, timings.
6. **The same names on every surface.** Rust field = Python keyword = YAML key.
7. **Validate once, with the field name.** One `validate` on the problem struct returns
   errors that name the field; the YAML layer stops re-checking.

### 3.2 The conceptual model

```
Problem
├── geometry: Geometry            # the Calabi-Yau
│   ├── charges (q)               # one row per toric divisor, one column per Kähler parameter
│   ├── intnums                   # {(i,j,k): v}; n-fold: i is a reference surface
│   ├── nef_partition             # [] = hypersurface
│   └── ambient: Ambient          # Plain | Fan(maximal cones) | Strata(list)
├── cone: MoriCone                # Generators(vectors) | Inequalities{hyperplanes, equations} | FromFan
├── grading_vector: Option<_>     # default: primitive sum of the cone's facet normals
└── selection: Selection          # UpToDegree(d) | Below{targets, max_deg} | MinClasses(n) | Classes(list)

Options
├── method: Method                # Auto | Modular | Reference
├── threads: Option<usize>
├── certify: bool
├── hilbert_bases: Option<Vec<CandidateCone>>   # precomputed K_T (Normaliz); never changes the answer
└── modular: modular::Tuning      # lanes, check_primes, min_primes, extraction tuning

Invariants<T>
├── kind, grading_vector
├── entries: Vec<{ curve, surface: Option<usize>, degree, value: T }>   # sorted by degree, then class
└── report: Report                # method used, primes, passes, certified, counts, timings, warnings
```

Why these choices:

- `Ambient` makes "fan XOR strata" a type. `Plain` is today's default (no strata, FRST
  assumed). `Fan(cones)` derives strata (possibly none), can derive the cone and the grading,
  and checks `intnums`. `Strata(list)` is for callers without a fan (today's `strata`). Two
  variants cannot be given at once, so the "the strata are derived from the fan, so they
  cannot be given with it" error (`dispatch.rs:202-206`) disappears.
- `MoriCone::FromFan` is only valid with `Ambient::Fan`; `validate` says so. The three
  variants cover today's `generators`, `mori_hyperplanes`+`mori_equations`, and the fan's
  cone. The precedence rules of `dispatch.rs:251-262` go away because there is one value.
- `Selection::Classes(list)` is today's "no truncation: generators as given", named for what
  it is. `generators` then always means "vectors spanning the cone" (C2 resolved).
  `MinClasses(n)` is `min_points`; it stays reference-only and says so in its doc, and
  `Method::Auto` reports the fallback instead of hiding it.
- `hilbert_bases` in `Options`, not in `MoriCone`: it is the only piece of today's cone
  data that is an accelerator. Its validation (`check_cone_data`'s per-cone checks,
  `cone.rs:127-170`) stays; the hyperplane part of that check moves to `MoriCone`.
- `grading_vector` optional everywhere, with the default PR #102 and `Fan::mori_cone`
  already use (primitive sum of facet normals: positive on a pointed cone), which is also
  what the 5.3 follow-up asks for. When the user gives one it is checked for positivity on
  the facets' rays (cheap: `cone::rays` of the facets, or the generators themselves).
- `Invariants.entries[i].surface` is `None` for threefolds instead of a different tuple
  shape (C5). `report` is decision 4.1.7.

### 3.3 Rust API

```rust
// src/lib.rs — the whole public surface, re-exported from `cygv::`

pub struct Geometry {
    /// One row per toric divisor, one column per Kähler parameter (h11).
    pub charges: DMatrix<i32>,
    /// Keyed by [i, j, k]; one ordering per triple; for an n-fold i is a reference surface.
    pub intnums: HashMap<[usize; 3], i32>,
    /// Disjoint parts covering every divisor; empty for a hypersurface.
    pub nef_partition: Vec<Vec<usize>>,
    pub ambient: Ambient,
}

pub enum Ambient {
    /// No information beyond the charges (a standard phase is assumed).
    Plain,
    /// The maximal cones of the fan, as divisor indices. Threefold hypersurfaces only.
    /// Gives the strata, and defaults for the cone and the grading vector.
    Fan(Vec<Vec<usize>>),
    /// The strata of a vex phase, for callers that have them but not the fan.
    Strata(Vec<Stratum>),
}

pub struct Stratum {
    pub divisors: [usize; 3],
    pub class: Vec<malachite::Rational>,
}

pub enum MoriCone {
    /// Any vectors spanning the cone (rays, or rays plus lattice points), one per column.
    Generators(DMatrix<i32>),
    /// h . C >= 0 for every column h, and e . C = 0 for every column e (a face).
    Inequalities { hyperplanes: DMatrix<i32>, equations: DMatrix<i32> },
    /// The cone spanned by the wall classes of `Ambient::Fan`.
    FromFan,
}

pub enum Selection {
    UpToDegree(u32),
    /// The backward lightcones of the targets, up to their largest degree or `max_deg`.
    Below { targets: DMatrix<i32>, max_deg: Option<u32> },
    /// At least this many classes (reference pipeline only).
    MinClasses(u32),
    /// Exactly these classes (zero is added if missing).
    Classes(DMatrix<i32>),
}

pub struct Problem {
    pub geometry: Geometry,
    pub cone: MoriCone,
    pub selection: Selection,
    /// Default: the primitive sum of the facet normals of the cone.
    pub grading_vector: Option<RowDVector<i32>>,
}

impl Problem {
    /// Checks everything the pipelines assume and resolves the defaults (grading vector,
    /// strata from the fan, the cone's inequalities), returning the error with the field
    /// at fault. `compute_*` call it; the YAML layer calls it for its own messages.
    pub fn validate(&self) -> Result<Resolved, Error>;
    pub fn kind(&self) -> CYKind;
    pub fn h11(&self) -> usize;
}

pub struct Options {
    pub method: Method,                        // Auto
    pub threads: Option<usize>,                // None: one per core
    pub certify: bool,
    /// Precomputed Hilbert bases of the cones K_T (e.g. from Normaliz). Trusted; never
    /// changes the result, only skips computing them.
    pub hilbert_bases: Option<Vec<CandidateCone>>,
    pub modular: modular::Tuning,              // lanes, check_primes, min_primes, extraction
}
impl Default for Options { ... }

pub fn compute_gv(problem: &Problem, options: &Options) -> Result<Invariants<Integer>, Error>;
pub fn compute_gw(problem: &Problem, options: &Options) -> Result<Invariants<Rational>, Error>;
/// For callers that only know the kind at run time (Python, the CLI).
pub fn compute(problem: &Problem, kind: InvariantKind, options: &Options) -> Result<Computed, Error>;
pub enum Computed { Gv(Invariants<Integer>), Gw(Invariants<Rational>) }

pub struct Invariants<T> {
    pub kind: CYKind,
    pub grading_vector: RowDVector<i32>,
    /// Sorted by degree, then by curve class, then by surface.
    pub entries: Vec<Invariant<T>>,
    pub report: Report,
}
pub struct Invariant<T> {
    pub curve: Vec<i32>,
    /// The reference surface for an n-fold; `None` for a threefold.
    pub surface: Option<usize>,
    pub degree: i64,
    pub value: T,
}
pub struct Report {
    pub method: Method,                  // the one that ran (never Auto)
    pub certified: bool,                 // true for the reference, or a certified modular run
    pub primes: usize, pub passes: usize,
    pub candidates: usize, pub classes: usize,
    pub grading_vector_derived: bool,
    pub timings: Vec<(&'static str, Duration)>,
    pub warnings: Vec<Warning>,          // e.g. FellBackToReference(reason), CertifyIgnored
}

pub enum Error {
    /// Invalid input; names the field (as in the Rust struct and the YAML/Python keys).
    Input { field: &'static str, reason: String },
    /// The request cannot run with the requested method.
    Unsupported { method: Method, what: &'static str },
    /// The modular lift did not stabilize, the strata could not be inverted, ...
    Computation(modular::Error),
    MoriCone(cone::ConeError),
    Fan(phase::FanError),
}
```

Usage:

```rust
let problem = Problem {
    geometry: Geometry { charges, intnums, nef_partition: vec![], ambient: Ambient::Fan(cones) },
    cone: MoriCone::FromFan,
    selection: Selection::UpToDegree(10),
    grading_vector: None,
};
let gv = cygv::compute_gv(&problem, &Options::default())?;
for e in &gv.entries { println!("{:?} (deg {}) {}", e.curve, e.degree, e.value); }
eprintln!("{:?}, {} primes, certified: {}", gv.report.method, gv.report.primes, gv.report.certified);
```

Plain structs with public fields rather than a builder: the fields are few, every one has
an obvious default or is required, and `..Default::default()` covers the rest. A builder
can be added later without breaking this.

Other Rust changes that go with it:

- `dispatch.rs` keeps `compute` and the pipeline choice; its five entry points and
  `compute_gvgw_strings` go. `reference::run_hkty` stays public *under `reference`* only,
  taking a `&Resolved` instead of twelve arguments, and stops panicking on input that
  `validate` already rejects. `modular::compute_gv/gw/bound_gv` stay as the low-level API.
- `CYKind::from_is_threefold` and `Method::from_name`/`as_str` stay (the Python and YAML
  layers need them); `InvariantKind::as_str` too.
- Public module surface: `pub mod cone; pub mod phase; pub mod modular; pub mod reference;`
  with the stage submodules `pub(crate)` or behind `#[doc(hidden)]`, and the root
  re-exports of `Polynomial`, `PolynomialProperties`, `PrettyPrintPolynomial`,
  `Semigroup`, `PolynomialCoeff` dropped (reachable as `cygv::reference::...` for anyone
  who wants them; decision 7's open question).
- `ConeData` becomes `MoriCone::Inequalities` + `Options::hilbert_bases`;
  `CandidateCone` keeps its name and shape.

### 3.4 Python API

Primary form, keyword-first, same names as Rust/YAML:

```python
def compute_gv(
    q: ArrayLike,                              # one charge vector per Kähler parameter, as today
    intnums: Mapping[tuple[int,int,int], int] | Iterable[Sequence[int]],   # dict, or [i,j,k,v] rows
    *,
    # which classes (exactly one; max_deg may accompany target_points)
    max_deg: int | None = None,
    target_points: ArrayLike | None = None,
    min_points: int | None = None,
    curve_classes: ArrayLike | None = None,
    # the Mori cone (exactly one, or none with a fan)
    generators: ArrayLike | None = None,
    mori_hyperplanes: ArrayLike | None = None,
    mori_equations: ArrayLike | None = None,
    # the geometry
    nefpart: Iterable[Iterable[int]] | None = None,
    fan: Iterable[Iterable[int]] | None = None,
    strata: Iterable[tuple[Iterable[int], Iterable[int | str | Fraction]]] | None = None,
    grading_vector: ArrayLike | None = None,    # default: from the cone's facets
    # options
    method: str | None = None,                  # "auto" | "modular" | "reference"
    n_threads: int | None = None,
    certify: bool = False,
    hilbert_bases: Iterable[tuple[Iterable[int], ArrayLike]] | None = None,
) -> Invariants: ...

def compute_gw(...same...) -> Invariants: ...
```

`Invariants` is a result class that *is* a sequence of `(class, value)` pairs, so
`dict(result)`, `for cls, v in result`, `sorted(result)`, `len(result)` and the README's
examples keep working unchanged, and it also carries the report:

```python
class Invariants(Sequence[tuple[Key, Value]]):
    kind: str  # "gv" | "gw"
    dimension: int  # 3, 4, ...
    grading_vector: tuple[int, ...]
    # .method, .certified, .primes, .passes, .candidates, .timings, .warnings
    report: Report

    def dok(self) -> dict[Key, Value]: ...
    def degree(self, curve) -> int: ...
    # (curve_tuple, value) for threefolds, ((curve_tuple, surface), value) for n-folds
    def __iter__(self): ...
```

Keeping the n-fold key as `(curve, surface)` preserves `dict()` round-trips for CYTools;
whether to switch to a flat key with a `surface` attribute per entry is a separate,
smaller decision (I would keep the pair: it is a hashable key, which `dict` needs).

Compatibility shim for the 0.3.0 positional order `(generators, grading_vector, q,
intnums)`: detectable because the old second positional is a 1-D integer vector where the
new API expects `intnums`. For one release, accept it with a `DeprecationWarning` that
prints the new spelling; after that, a `TypeError` with the same message. CYTools is
unaffected either way (it passes keywords only; Section 1.4).

Errors: `ValueError` for invalid input (as now), a new `cygv.ComputationError(RuntimeError)`
for failures of a valid request (`NotStable`, `NotInvertible`, `TooLarge`), `TypeError`
for wrong Python types (already what PyO3 raises; the wrapper should convert `intnums`
lists before the FFI so that it is a `ValueError` only when the *content* is wrong).
`method` keeps its string form; expose `cygv.METHODS = ("auto", "modular", "reference")`.

Marshalling: pass one dict through the pipe and one `#[pyfunction]` taking a `&PyDict` or a
`#[derive(FromPyObject)]` struct, so a new field is added in one Python place and one Rust
place instead of six (P6).

**CYTools objects.** Two ways to accept `cy.compute_gvs(...)`'s data:

- *A duck-typed adapter in cygv*: `cygv.Problem.from_cytools(cy, grading_vector=None)`
  reading `cy.curve_basis(include_origin=False, as_matrix=True)`,
  `cy.intersection_numbers(in_basis=True, format="dok")` and
  `cy.mori_cone_cap(in_basis=True).rays()` (the three calls CYTools and PR #102 both make,
  `calabiyau.py:2061-2068`, `cgv_run.py:270-283`), and, when available, the fan with the
  divisor relabelling `cgv_phase.phase_from_polytope` does (`cgv_phase.py:300-311`). No
  import of CYTools, no dependency. Convenient; but it ties cygv to CYTools' method names
  and to the ray-labelling convention, which cygv cannot test in CI without CYTools.
- *CYTools builds the `Problem`*: cygv provides the dataclass and a clear mapping in the
  README; CYTools' `_compute_gvs_gws` becomes `cygv.compute_gv(q=..., intnums=...,
  generators=mori.rays(), max_deg=...)`.

Recommendation: ship the `Problem` dataclass and document the mapping (second option),
and offer `from_cytools` only if CYTools' maintainers want the fan-relabelling logic to
live in cygv; the FRST-only three-call version is small enough to live in either.

```python
@dataclass
class Problem:
    q: ArrayLike
    intnums: Mapping | Iterable
    nefpart: ... = None
    generators: ... = None
    mori_hyperplanes: ... = None
    mori_equations: ... = None
    fan: ... = None
    strata: ... = None
    grading_vector: ... = None
    def gv(self, *, max_deg=None, target_points=None, min_points=None, curve_classes=None, **options) -> Invariants
    def gw(self, ...same...) -> Invariants
```

### 3.5 YAML and CLI

One-of groups become nested mappings, everything else stays flat with the same names as
Python/Rust. Released field names (`generators`, `grading_vector`, `max_deg`, `min_points`,
`target_points`, `nefpart`, `invariants`, `name`, `n_threads`) stay valid at the top level,
so the 0.3.0 example documents still run.

```yaml
---
name: vex phase of P(1,1,2,2,2)[8]
q: [[1, 0, -1, 1, 0, 1], [0, 1, 1, 1, 1, 0]]
intnums: [[0, 1, 1, 4], [1, 1, 1, 2]]       # or the "i,j,k": v mapping
# nefpart: [[...], [...]]                   # complete intersections only
fan: [[0, 1, 2, 3], [0, 1, 2, 4], ...]      # or strata: [{divisors: [1, 2, 4], class: [1, -1]}]
# the Mori cone: one of generators | mori_hyperplanes (+ mori_equations) | nothing with a fan
# generators: [[0, 1], [1, -1]]
# grading_vector: [2, 1]                    # default: from the facets
# which classes: one of max_deg | target_points (+ max_deg) | min_points | curve_classes
max_deg: 10
invariants: [gv, gw]                        # one or both; default gv
options:
  method: auto
  n_threads: 4
  certify: true
  hilbert_bases:                            # optional accelerator (renamed from cones)
    - {divisors: [0], hilbert_basis: [[0, 1], [1, 1]]}
```

Output:

```yaml
---
name: vex phase of P(1,1,2,2,2)[8]
dimension: 3
grading_vector: [2, 1]
report: {method: modular, certified: false, primes: 3, passes: 2, candidates: 61, classes: 58, seconds: 0.41}
gv:
  - {curve_class: [1, -1], degree: 1, gv: 1}
  - {curve_class: [0, 1], degree: 1, gv: 640}
gw:
  - {curve_class: [1, -1], degree: 1, gw: '1'}
```

`is_threefold` becomes `dimension`; `surface_index` becomes `surface` and appears only
for n-folds, as now. A single-kind document gets a single `gv:`/`gw:` list; allowing both
in one document saves re-reading large inputs twice.

CLI: keep `--file/--output/--threads/--method/--certify`; add `--max-deg N` and
`--invariants gv|gw|gv,gw` as per-document overrides (Y4), and `--no-certify` for symmetry.
`--help` shrinks to the field list with one line each once the rules are in the schema.

The YAML layer calls `Problem::validate()` and maps `Error::Input { field, .. }` to its
`InvalidField { field, reason }`; the duplicated checks in `io.rs:317-396` go.

### 3.6 Behavioural changes that remove rules

| Today | Proposed |
|---|---|
| Grading vector required without a fan | Derived from the facets (cone data, generators or fan) when omitted; checked when given; `report.grading_vector_derived` |
| `generators` = classes when nothing truncates | `Selection::Classes` / `curve_classes:`; `generators` without any selection is an error for one release with a message, then simply an error |
| `min_points` + `max_deg` accepted by Python/Rust, rejected by CLI | `Selection` is one value; impossible |
| `certify` silently ignored on the reference / after fallback | `report.certified` is `true` for the reference (exact) and `report.warnings` has `CertifyIgnored` when `auto` fell back; `Method::Reference` + `certify` is fine (exact anyway) |
| `auto` falls back to the reference on key width silently | `report.warnings: FellBackToReference(KeyTooWide)`; the Python wrapper emits a `warnings.warn`; the CLI prints to stderr. (Whether it should be an error instead is open: the reference does finish at small h11 with large coordinates, so a warning seems right.) |
| Transposed `q` → "dimension is -5" | `validate` tries the transpose: "q has 2 charge vectors of 6 entries but 6 Kähler parameters were implied; it looks transposed" |
| `intnums` list is a `TypeError` in Python | Accepted, as in YAML |
| Results unordered | Sorted by degree, then class, on every surface |
| `fan` + `strata` error | `Ambient` enum |
| `mori_hyperplanes` without `cones` error | `MoriCone::Inequalities` alone is complete; `hilbert_bases` is optional |

---

## 4. Compatibility and migration

### 4.1 What is released and what is not

`v0.3.0` is `main` (same `Cargo.toml` version, tag on `main` only). Its surface:

- Python: `compute_gv/gw(generators, grading_vector, q, intnums, max_deg=None,
  min_points=None, target_points=None, nefpart=None, prec=None)` (`main:python/cygv/hkty.py:141-150`).
- YAML: `name, generators, grading_vector, q, intnums, nefpart, invariants, max_deg,
  min_points, target_points, is_threefold, prec, n_threads` (`main:src/io.rs:59-73`).
- Rust: `cygv::hkty::compute_g{v,w}_rat_{threefold,nfold}`, module paths `cygv::semigroup`
  etc. (`main:src/lib.rs`).

Everything else on this branch (`method`, `certify`, `strata`, `fan`, `mori_hyperplanes`,
`mori_equations`, `cones`, `n_threads` in Python, optional `generators`/`grading_vector`,
`ComputeError`, `dispatch`, `cone`, `phase`, `modular`) has never shipped and can be
renamed or reshaped **without any deprecation period**. This is the main reason to do the
redesign *now*: after the next release every one of those names is a commitment.

The next release already breaks `prec` (decision 6.9), the Rust module layout (decision 7)
and the result semantics on unsaturated input (6.2). It is a breaking bump anyway.

### 4.2 What stays, what is aliased, what breaks

| Surface | Stays as is | Deprecated alias (one release) | Breaks |
|---|---|---|---|
| Python | keyword names `generators, grading_vector, q, intnums, max_deg, min_points, target_points, nefpart`; list-of-pairs iteration of the result; `int`/`Fraction` values; `(curve, surface)` keys for n-folds | positional `(generators, grading_vector, q, intnums)` (detected, warns); `generators` with no selection → treated as `curve_classes` with a warning | `prec` (already gone); unreleased kwargs `cones`→`hilbert_bases` (never shipped, no alias needed) |
| CYTools | **everything**: it calls by keyword with released names and uses `dict()`/pair iteration (`calabiyau.py:2058-2076`) | the no-truncation default path hits the `generators`-as-classes alias and warns once | nothing; recommended follow-up in CYTools: pass `generators=mori.rays()` and prefer `max_deg` to `min_points`, and stop computing its own grading vector |
| YAML | the 13 released fields and their meanings; `name`; both `intnums` spellings | `is_threefold` (input: accepted and checked against the inferred kind; output: replaced by `dimension`); `generators` with no selection → `curve_classes` with a stderr note | `prec`; unreleased `mori_*`/`cones`/`strata`/`fan` move into their new places |
| Rust | `InvariantKind`, `CYKind`, `Method`, `Stratum` (field types change), `malachite` re-export, `modular::*` low-level API, `cone::*`, `phase::Fan` | none (pre-1.0, decision 7: break freely) | the five entry points, `compute_gvgw_strings`, root `run_hkty`, root re-exports of `Polynomial` etc., `ConeData` |

### 4.3 Order of steps

Each step leaves the tree green and is independently reviewable.

1. **Rust core types** (`Problem`, `Geometry`, `Ambient`, `MoriCone`, `Selection`,
   `Options`, `Invariants`, `Report`, `Error`) in `src/lib.rs`; `dispatch::compute` rewritten
   over them; the old five entry points reimplemented as thin wrappers for one commit so the
   tests keep running, then deleted together with `compute_gvgw_strings`. `validate` absorbs
   `check_*` and returns field-tagged errors. `Report` filled by both pipelines (the modular
   `run` already knows primes and passes; the reference reports `certified: true`).
2. **Defaults**: grading vector from facets; transposed-`q` message; sorting in `compute`.
3. **YAML** over `Problem::validate`; nested `options:`; `dimension`/`report` in the output;
   `--max-deg`, `--invariants`; examples and `tests/cli.rs` updated; `is_threefold` and
   untruncated `generators` aliases.
4. **Python**: keyword-first signature with the positional shim; `Invariants` class;
   `Problem` dataclass; single-dict marshalling; `ComputationError`; list `intnums`;
   `python/tests` extended with the shim and the report.
5. **Docs**: README rewritten around the three nouns (the current five paragraphs of rules
   become one table); `AGENTS.md` dispatch section shortened to match; `PLAN.md` 6.2 "opt-out
   of saturation" closed by `Selection::Classes`; release notes listing Section 4.2.
6. **CYTools PR** (separate repository): keywords unchanged, but pass the rays only, use
   `max_deg`, read `result.report`.

Steps 1-2 are the bulk; 3 and 4 are mostly mechanical once 1 exists.

---

## 5. Prioritized changes

Costs are rough, for one person familiar with the code. "Benefit" is for users unless
noted.

### Quick wins (each under half a day, independent of the restructuring)

| # | Change | Cost | Benefit | Where |
|---|---|---|---|---|
| Q1 | Default grading vector from the facets when none is given (already a 5.3 follow-up) | small | Removes a required argument and `MissingGradingVector`; matches PR #102 and CYTools | `dispatch.rs:230`, `:257-262` |
| Q2 | Reject `min_points` with `max_deg`/`target_points` in `dispatch` (today only the CLI does) | tiny | Removes a silent-ignore (Appendix B) | `dispatch.rs:263`, `reference.rs:41-45` |
| Q3 | Accept `intnums` as a list of `[i,j,k,v]` in Python | tiny | Consistency with YAML; CYTools' coordinate-list output format | `hkty.py`, `python.rs:36` |
| Q4 | Sort results by degree then class in `compute` for both pipelines | tiny | Deterministic output, no per-caller sorting | `dispatch.rs:465-483` |
| Q5 | Warn (Python `warnings`, CLI stderr) when `auto` falls back to the reference, and when `certify` has no effect | small | The two most surprising silent behaviours | `dispatch.rs:264-268`, `:310-313` |
| Q6 | Detect a transposed `q` and say so | small | The most likely first-time error gets a useful message | `dispatch.rs:189-192`, `io.rs:380-386` |
| Q7 | Fix the rows/columns wording to one sentence used everywhere ("a divisor index is a position inside each charge vector; in Rust, a row of `q`") | tiny | C6 | `io.rs:150`, `examples/fan.yaml:16`, `lib.rs:75`, `README.md:27` |
| Q8 | `--max-deg` and `--invariants` CLI overrides | small | Y4 | `main.rs:81-100` |
| Q9 | Rename `cones` → `hilbert_bases` and `surface_index` → `surface` while unreleased | tiny | Names say what they are | `io.rs`, `hkty.py`, `python.rs`, README, examples |
| Q10 | `Stratum::class: Vec<malachite::Rational>` | small | Uses the type the crate already re-exports for this | `lib.rs:79`, `python.rs:72-83`, `io.rs:672-677` |

Q1, Q2, Q4, Q9 and Q10 are subsumed by the restructuring but are worth doing first if the
restructuring is deferred, because Q1/Q2/Q9/Q10 change unreleased behaviour that becomes
frozen at the next release.

### Medium (one to three days each)

| # | Change | Cost | Benefit |
|---|---|---|---|
| M1 | `Problem`/`Options`/`Invariants`/`Report`/`Error` and one `compute` (Section 3.3, step 1) | 2-3 days incl. tests | Removes R1-R5, R7, C1-C5, C7; makes every later option a field. The single largest clarity gain and the prerequisite for M2-M4 |
| M2 | YAML over `validate`, nested `options`, `dimension`/`report` output, aliases (step 3) | 1 day | Y1-Y3, Y5; halves `io.rs` validation |
| M3 | Python keyword-first signature, positional shim, `Invariants` class, `Problem` dataclass, single-dict marshalling, `ComputationError` (step 4) | 1-2 days | P1-P6; CYTools keeps working; new users never see `compute_gv(None, None, q, ...)` |
| M4 | `Report` contents: primes, passes, candidates, timings, warnings from both pipelines (4.1.7) | 1 day | C7; also what the benches and a future profiler read |
| M5 | README/AGENTS.md rewrite around the model (step 5) | 1 day | C9; the test of the redesign is that the rules paragraphs disappear |

### Larger / optional

| # | Change | Cost | Benefit | Note |
|---|---|---|---|---|
| L1 | Trim the public Rust surface (`pub(crate)` stage modules, drop root re-exports, `run_hkty` under `reference` with a `Resolved` input, no panics) | 1-2 days (doc links, tests that reach into modules) | R6; a `cargo doc` page that is the API | Decision 7's open question; pre-1.0 so free to do |
| L2 | `Problem.from_cytools` adapter with fan relabelling | 1 day + CYTools coordination | One-line use from CYTools incl. vex phases | Maintenance coupling; see 3.4 for the alternative |
| L3 | Both invariants in one YAML document / one Python call (`compute(problem, kinds=...)`) sharing the GV pass | 1 day | GW is derived from GV (`multicover::gw_from_gv`), so it is nearly free; today it reruns everything | Nice-to-have |
| L4 | Certifying vex phases; faces of high h11 in a reduced lattice basis | large | Removes the two remaining "unsupported" branches (`Error::CertifyStrata`, key width fallback) | Computation work listed in `PLAN.md` 5.3, not interface work; the interface should already report them as warnings/errors with the right field |

### What I would not change

- The pipeline split and `Method::Auto` as the default: users should not have to know,
  and after M1 they only meet it in `report.method`.
- Column-major `DMatrix` inside Rust: it is the crate's convention and every stage relies
  on it; the fix for C6 is one sentence of documentation at the boundary, not a transpose.
- Python running in a subprocess, `Fraction`s, `ValueError` for bad input.
- The YAML stream-of-documents format and the two `intnums` spellings.

---

## Appendix A. Option interaction matrix today

Rows: how the Mori cone is given. Columns: the selection. Cells: what runs under
`Method::Auto`, and the gotcha.

| cone \ selection | `max_deg` | `target_points` (+`max_deg`) | `min_points` | none |
|---|---|---|---|---|
| `generators` | modular, saturated (`dispatch.rs:283-294`); reference on key overflow, silently | same | reference, saturated via Hilbert basis (`:324-333`) | reference on the generators **as given** (`:316-321`) |
| `mori_hyperplanes`+`mori_equations`+`cones` (all three or error) | modular from the given bases; generators optional | same | reference; needs no generators | **error**: `MissingGenerators` (`:319`), although the cone is fully described |
| `fan` | modular; strata, cone and grading derived; `strata` forbidden | same | reference | **error**: `MissingGenerators` |
| `fan` + `generators` | fan's cone ignored, generators used (`:253`); grading still derived | same | reference | reference on the generators as given |
| nothing | `MissingGenerators` | | | |

`certify` has an effect only in the modular cells. `grading_vector` is required in every
cell without a fan. `strata` is allowed in every cell without a fan, threefold
hypersurfaces only; `cones` must then include each stratum's cone.

## Appendix B. Probe transcript

Run with `uv run python` on this branch (`examples`' P(1,1,2,2,2)[8], `q = [[0,0,1,1,1,1],
[1,1,0,0,0,-2]]`, `intnums = {(0,0,0): 8, (0,0,1): 4}`, `gens = [[1,0],[0,1]]`):

```
no truncation, gens as given:   OK [((1, 0), 640), ((0, 1), 4)]
no truncation, None gens:       ValueError: the generators are needed, without cone data or a fan, or without max_deg, min_points and target_points
no grading vector:              ValueError: the grading vector is needed without a fan
modular + min_points:           ValueError: the modular method does not support min_points; use the auto or reference method
modular, no truncation:         ValueError: the modular method does not support computations without max_deg or target_points; ...
min_points=5 + max_deg=3:       OK (min_points silently ignored; the CLI rejects this combination)
q transposed (divisor rows):    ValueError: the dimension of the CY must be at least three, but it is -5
only mori_hyperplanes:          ValueError: the cone data needs both mori_hyperplanes and cones
generators 1-D:                 ValueError: generators must be a 2D array-like of ints
target_points 1-D:              OK (one point)
grading not positive:           ValueError: the grading vector must be positive on every nonzero point of the Mori cone
certify=True + reference:       OK (certify silently ignored)
n_threads=0:                    OK (means one thread)
intnums as list of 4-tuples:    TypeError: 'list' object is not an instance of 'dict'
fan with n-fold-shaped q:       ValueError: a fan is only supported for threefold hypersurfaces
```

CLI (`./target/release/cygv` on stdin):

```
q + intnums + generators + max_deg, no grading_vector:  Error: the required field "grading_vector" is missing
q + intnums + grading_vector + max_deg, no generators:  Error: the required field "generators" is missing
```
