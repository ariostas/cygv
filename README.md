# cygv

[![Rust CI](https://github.com/ariostas/cygv/actions/workflows/rust.yml/badge.svg)](https://github.com/ariostas/cygv/actions/workflows/rust.yml) [![Python CI](https://github.com/ariostas/cygv/actions/workflows/python.yml/badge.svg)](https://github.com/ariostas/cygv/actions/workflows/python.yml) [![Crates.io Version](https://img.shields.io/crates/v/cygv)](https://crates.io/crates/cygv) [![PyPI - Version](https://img.shields.io/pypi/v/cygv)](https://pypi.org/project/cygv/) [![PyPI - Downloads](https://img.shields.io/pypi/dm/cygv)](https://pypi.org/project/cygv/)

This project implements an efficient algorithm to perform the HKTY procedure [[1], [2]] to compute Gopakumar-Vafa (GV) and Gromov-Witten (GW) invariants of Calabi-Yau (CY) manifolds obtained as hypersurfaces or complete intersections in toric varieties. This project is based on the work presented in the paper "[Computational Mirror Symmetry]", but written in the Rust programming language and with some additional improvements.

Every computation needs three things:

- **the Calabi-Yau**: its GLSM charges `q` and its intersection numbers `intnums`, the nef
  partition `nefpart` of a complete intersection, and for a vex phase of the ambient toric
  variety either its `fan` or its `strata`;
- **its Mori cone**: vectors that span it (`generators`), or the inequalities that cut it out
  (`mori_hyperplanes`, with `mori_equations` for a face), or nothing when the `fan` gives it;
- **which curve classes**: those up to a degree (`max_deg`), those below some classes
  (`target_points`), the lowest degrees holding a number of them (`min_points`), or exactly the
  ones listed (`curve_classes`).

Everything else has a default: the degree is taken with the `grading_vector`, by default the
primitive sum of the inequalities of the Mori cone; the computation picks its method; and every
lattice point of the Mori cone is used, even when the generators do not generate them all. The
names are the same in Python, in the YAML input of the command line interface, and in Rust.

## Python package

The `cygv` package on PyPI ships prebuilt wheels, so no Rust toolchain is needed to install it.

```bash
pip install cygv
```

A `Problem` holds the Calabi-Yau and its Mori cone, and its `gv` and `gw` methods compute the
invariants of the selected curve classes. Every vector argument is a list of vectors, or any
array-like that numpy accepts.

```python
from cygv import Problem

# The CY threefold hypersurface in P(1,1,2,2,2), resolved, with h^{1,1} = 2.
problem = Problem(
    # The GLSM charges, one charge vector per Kähler parameter.
    q=[[0, 0, 1, 1, 1, 1], [1, 1, 0, 0, 0, -2]],
    # The intersection numbers, one ordering of each triple (or a list of [i, j, k, value]).
    intnums={(0, 0, 0): 8, (0, 0, 1): 4},
    # Vectors that span the Mori cone.
    generators=[[1, 0], [0, 1]],
)

gv = problem.gv(max_deg=3)
# [((0, 1), 4), ((1, 0), 640), ((1, 1), 640), ((2, 0), 10032), ((2, 1), 72224), ...]
dict(gv)[(1, 0)]
# 640
gv.grading_vector, gv.report
# ((1, 1), Report(method='modular', certified=False, grading_vector_derived=True, seconds=0.003))

problem.gw(max_deg=3)
# [((0, 1), Fraction(4, 1)), ((1, 0), Fraction(640, 1)), ((0, 2), Fraction(1, 2)), ...]
```

The results are `(curve_class, invariant)` pairs sorted by degree and then by curve class, with
GV invariants as `int`s and GW invariants as `fractions.Fraction`s. For CYs of dimension greater
than three, an invariant is labeled by a curve class and by a reference surface (the first index
of the intersection numbers), as in `(((1, 0, -4, 0, 0, 0), 2), 80)`. They also carry the
`dimension`, the `grading_vector` and a `report` of how they were computed.

`compute_gv` and `compute_gw` do the same in one call, taking the arguments of `Problem` and of
`gv` or `gw` together: `compute_gv(q, intnums, generators=..., max_deg=3)`.

`target_points` selects the curve classes $C$ with $P - C$ in the Mori cone for one of the given
classes $P$, up to the largest degree of the $P$, or up to `max_deg` if it is given too. A
complete intersection takes its nef partition as `nefpart`, lists of divisor indices (positions
in the charge vectors) that split the divisors, each into exactly one part.

A vex phase of the ambient toric variety (MacFadden–Sheridan, arXiv:2512.14817), whose curves of
negative anticanonical degree need the class of the stratum they force into the CY, needs those
strata, as `strata`: `(divisors, class)` pairs of the three divisors and the class of the
stratum, whose entries may be `Fraction`s. They can instead be derived from the fan of the
ambient toric variety, as `fan`: its maximal cones, lists of four divisors, which must form a
fine regular fan on the boundary points of a reflexive polytope that are not interior to facets
(FRST or vex). The fan also gives the Mori cone (spanned by the classes of its walls) and the
default grading vector, and `intnums` is checked against its intersection ring. Both are only
supported for threefold hypersurfaces.

```python
fan = [
    (0, 1, 2, 3),
    (0, 1, 2, 4),
    (0, 1, 3, 4),
    (0, 2, 3, 4),
    (1, 2, 3, 5),
    (1, 2, 4, 5),
    (1, 3, 4, 5),
    (2, 3, 4, 5),
]
q = [[1, 0, -1, 1, 0, 1], [0, 1, 1, 1, 1, 0]]
Problem(q, {(0, 1, 1): 4, (1, 1, 1): 2}, fan=fan).gv(max_deg=3)
# [((0, 1), 640), ((1, -1), 1), ((0, 2), 10032), ((1, 0), 20), ...]
```

Two implementations of the computation are available, selected with `method`. The `modular` one
computes modulo word-sized primes and recovers the exact integers by the Chinese remainder
theorem, which is much faster, but only supports `max_deg` and `target_points`, and curve classes
whose coordinates fit its keys of up to 512 bits. On a face of the Mori cone it takes their
coordinates in a reduced basis of the lattice they span, which keeps them small: faces of rank
up to a few dozen fit even when $h^{1,1}$ is in the hundreds. The `reference` one is a faithful
implementation of the HKTY procedure with exact rational arithmetic. The default, `auto`, uses
the modular one for `max_deg` and `target_points` and the reference one otherwise; curve
classes too large for the modular one are an error rather than a slow run of the reference.
`n_threads` sets the number of threads, which defaults to one per core.

The modular implementation accepts its result once one more prime agrees with it, which makes a
wrong result vanishingly unlikely but does not rule it out. With `certify=True` it proves the
result instead: it bounds every invariant from above, by running the computation once more over
upward-rounded reals, and takes more primes until every bound is below half their product, which
costs about as much as the computation itself. Vex phases cannot be certified yet.

For runs that do not fit in memory otherwise, `low_memory=True` makes the modular implementation
prefer memory to speed: it takes about a third of the memory, for about two and a half times
the time.

The modular implementation takes the curve classes from the Hilbert bases of the cones $K_T$ of
the classes in the Mori cone that only meet the divisors in $T$ negatively, which it computes
itself. A program such as [Normaliz](https://www.normaliz.uni-osnabrueck.de/) can compute them
instead, passed as `hilbert_bases`: `(divisors, hilbert_basis)` pairs for every set $T$ of one or
two divisors, and in a vex phase also for the three divisors of every stratum. Elements of degree
above `max_deg` can be left out. They are trusted, so a missing element silently drops curve
classes.

Invalid input raises a `ValueError`, and a computation that fails (for instance because the
invariants are not integers, which usually means inconsistent input) a `cygv.ComputationError`.
The computation runs in a subprocess, so it can be interrupted with ctrl+c without losing the
Python session.

## Rust crate

The same model, as types: a `Problem` of a `Geometry`, a `MoriCone` and a `Selection`, computed
by `compute_gv` or `compute_gw` with `Options`.

```rust
use cygv::{compute_gv, Geometry, MoriCone, Options, Problem, Selection};
use nalgebra::dmatrix;
use std::collections::HashMap;

let problem = Problem {
    geometry: Geometry {
        // One row per toric divisor, one column per Kähler parameter.
        q: dmatrix![0, 1; 0, 1; 1, 0; 1, 0; 1, 0; 1, -2],
        intnums: HashMap::from([((0, 0, 0), 8), ((0, 0, 1), 4)]),
        ..Geometry::default()
    },
    // One generator per column.
    cone: MoriCone::Generators(dmatrix![1, 0; 0, 1]),
    selection: Selection::MaxDeg(3),
    grading_vector: None,
};
let gv = compute_gv(&problem, &Options::default()).unwrap();
assert_eq!(gv.grading_vector, [1, 1]);
assert_eq!((gv.entries[1].curve.as_slice(), gv.entries[1].value.clone()), (&[1, 0][..], 640.into()));
```

## Command line interface

This project also ships a `cygv` executable, so that it can be used without Python. It reads
YAML-formatted data from a file, or from the standard input, and writes the resulting invariants
to a file, or to the standard output.

```bash
cargo install cygv
cygv --file input.yaml --output invariants.yaml
```

It is part of the default `cli` feature. Library users that do not need it can depend on this
crate with `default-features = false` to avoid pulling in the argument and YAML parsers.

The input is a stream of YAML documents, each of which specifies a single computation, with the
same names as the Python arguments.

```yaml
---
name: an example threefold
q: [[0, 0, 1, 1, 1, 1], [1, 1, 0, 0, 0, -2]]
intnums:
  - [0, 0, 0, 8]
  - [0, 0, 1, 4]
generators: [[1, 0], [0, 1]]
invariants: gv
max_deg: 3
```

The intersection numbers can equivalently be given as a mapping from indices to values, which is
easier to read when there are many of them.

```yaml
intnums:
  "0,0,0": 8
  "0,0,1": 4
```

The results are written back as one YAML document per input document, sorted by the degree of
the curve classes under the grading vector.

```yaml
---
name: an example threefold
invariants: gv
dimension: 3
grading_vector: [1, 1]
report: {method: modular, certified: false, grading_vector_derived: true}
results:
  - {curve_class: [0, 1], degree: 1, gv: 4}
  - {curve_class: [1, 0], degree: 1, gv: 640}
  - {curve_class: [1, 1], degree: 2, gv: 640}
  - {curve_class: [2, 0], degree: 2, gv: 10032}
  - {curve_class: [2, 1], degree: 3, gv: 72224}
  - {curve_class: [3, 0], degree: 3, gv: 288384}
```

Run `cygv --help` for the full list of input fields and options. `--method`, `--threads`,
`--max-deg`, `--invariants`, `--certify` and `--low-memory` set those for every document, and
the `examples` directory has complete inputs.

## License

Licensed under the GNU General Public License v3 or later (GPLv3+)
([LICENSE](https://github.com/ariostas/cygv/blob/main/LICENSE) or <https://www.gnu.org/licenses/gpl-3.0.en.html#license-text>).

## Contribution

Any contribution submitted for inclusion in the work by you, shall be
licensed as under the GPLv3+, without any additional terms or conditions.
See [CONTRIBUTING.md](https://github.com/ariostas/cygv/blob/main/CONTRIBUTING.md) for guidelines.

[1]: https://arxiv.org/abs/hep-th/9308122
[2]: https://arxiv.org/abs/hep-th/9406055
[Computational Mirror Symmetry]: https://arxiv.org/abs/2303.00757
