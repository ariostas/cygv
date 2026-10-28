# cygv

[![Rust CI](https://github.com/ariostas/cygv/actions/workflows/rust.yml/badge.svg)](https://github.com/ariostas/cygv/actions/workflows/rust.yml) [![Python CI](https://github.com/ariostas/cygv/actions/workflows/python.yml/badge.svg)](https://github.com/ariostas/cygv/actions/workflows/python.yml) [![Crates.io Version](https://img.shields.io/crates/v/cygv)](https://crates.io/crates/cygv) [![PyPI - Version](https://img.shields.io/pypi/v/cygv)](https://pypi.org/project/cygv/) [![PyPI - Downloads](https://img.shields.io/pypi/dm/cygv)](https://pypi.org/project/cygv/)

cygv computes Gopakumar-Vafa (GV) and Gromov-Witten (GW) invariants of Calabi-Yau (CY)
manifolds obtained as hypersurfaces or complete intersections in toric varieties, by the HKTY
procedure [[1], [2]]. It is a Rust implementation, with some improvements, of the algorithm of
"[Computational Mirror Symmetry]", usable from Python, from Rust, and from the command line.

Every computation needs three things, with the same names in all three:

- **the CY**: its GLSM charges `q` and intersection numbers `intnums`; for a complete
  intersection its nef partition `nefpart`, and for a vex phase of the ambient toric variety its
  `fan` or its `strata`;
- **its Mori cone**: vectors that span it (`generators`), the inequalities that cut it out
  (`mori_hyperplanes`, with `mori_equations` for a face), or nothing when the `fan` gives it;
- **the curve classes**: those up to a degree (`max_deg`), those below some classes
  (`target_points`), the lowest degrees holding a number of them (`min_points`), or exactly the
  listed ones (`curve_classes`).

The degree is taken with the `grading_vector`, by default the primitive sum of the inequalities
of the Mori cone. Every lattice point of the Mori cone is used, even when the generators do not
generate them all.

## Python

The `cygv` package on PyPI ships prebuilt wheels, so no Rust toolchain is needed.

```bash
pip install cygv
```

A `Problem` holds the CY and its Mori cone, and its `gv` and `gw` methods compute the invariants.
Vectors can be lists or any array-like that numpy accepts.

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
# Invariants(gv, [((0, 1), 4), ((1, 0), 640), ((1, 1), 640), ((2, 0), 10032), ((2, 1), 72224), ... (6 in all)])
gv[(1, 0)], gv.get((5, 5))
# (640, None)
gv.curves[:2], gv.values[:2]
# (array([[0, 1], [1, 0]], dtype=int32), array([  4, 640]))
gv.grading_vector, gv.report
# ((1, 1), Report(method='modular', certified=False, grading_vector_derived=True, seconds=0.003))

list(problem.gw(max_deg=3))
# [((0, 1), Fraction(4, 1)), ((1, 0), Fraction(640, 1)), ((0, 2), Fraction(1, 2)), ...]
```

The result is a sequence of `(curve_class, invariant)` pairs sorted by degree, with GV invariants
as `int`s and GW invariants as `Fraction`s, that can be indexed by curve class or turned into a
`dict`. The same data are numpy arrays, much faster for many invariants: `curves`, `degrees`,
`values` (`int64` when they all fit, Python objects otherwise) and, for CYs of dimension greater
than three, `surfaces`, since there an invariant is also labeled by a reference surface (the
first index of the intersection numbers), as in `(((1, 0, -4, 0, 0, 0), 2), 80)`.
`compute_gv(q, intnums, generators=..., max_deg=3)` and `compute_gw` do the same in one call.

`target_points` selects the classes $C$ with $P - C$ in the Mori cone for one of the given $P$,
up to the largest degree of the $P$ or up to `max_deg`. `nefpart` lists, for each part, the
divisors (positions in the charge vectors) it contains.

**Vex phases.** In a vex phase (MacFadden–Sheridan, arXiv:2512.14817), curves of negative
anticanonical degree need the strata they force into the CY. Give either the `strata`, as
`(divisors, class)` pairs of three divisors and the stratum's class (entries may be `Fraction`s),
or the `fan` of the ambient toric variety: its maximal cones, which must form a fine regular fan
on the boundary points of a reflexive polytope that are not interior to facets. The fan also
gives the Mori cone and the default grading vector, and `intnums` is checked against it. Both are
only supported for threefold hypersurfaces.

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
# Invariants(gv, [((0, 1), 640), ((1, -1), 1), ((0, 2), 10032), ((1, 0), 20), ((0, 3), 288384), ... (6 in all)])
```

**Options.**

| Keyword | Effect |
| --- | --- |
| `method` | `"modular"` computes modulo primes and lifts by the CRT, and is much faster; `"reference"` is a faithful implementation in exact rationals. The default, `"auto"`, picks `modular` for everything but `min_points`, which it does not support. |
| `certify` | Proves the modular result, by bounding every invariant and taking primes until the bounds fit, for about twice the time. Without it the result is accepted once one more prime agrees, which makes an error vanishingly unlikely. Not for vex phases yet. |
| `low_memory` | About a third of the memory, for about two and a half times the time. |
| `n_threads` | The number of threads, one per logical CPU by default. |
| `verbose` | Writes the progress to `sys.stderr`; `2` also writes each degree. |
| `hilbert_bases` | The Hilbert bases of the cones the modular method takes its curve classes from, computed elsewhere (by [Normaliz](https://www.normaliz.uni-osnabrueck.de/), say) instead of by cygv: `(divisors, basis)` pairs for every set of one or two divisors (and the three of each stratum). They are trusted, so a missing element silently drops classes. |

The modular method takes curve classes of up to 512 bits; on a face of the Mori cone it uses a
reduced basis of the lattice they span, so faces of rank up to a few dozen fit even when
$h^{1,1}$ is in the hundreds. Larger classes are an error, not a slow run of the reference.

**CYTools.** A CY of [CYTools](https://cy.tools) (which calls cygv in `cy.compute_gvs`) gives its
data in the form cygv takes, with the origin of the polytope left out of the charges:

```python
gv = compute_gv(
    cy.curve_basis(include_origin=False, as_matrix=True),
    cy.intersection_numbers(in_basis=True, format="dok"),
    generators=cy.mori_cone_cap(in_basis=True).rays(),
    max_deg=10,
)
```

**Errors and interrupts.** Invalid input raises a `ValueError`, and a failed computation (for
instance invariants that are not integers, which usually means inconsistent input) a
`cygv.ComputationError`. The computation runs on threads of its own with the GIL released, and
ctrl+c (or interrupting a Jupyter kernel) stops it within a fraction of a second, raising
`KeyboardInterrupt`; a second ctrl+c stops waiting. Since it runs in the Python process, running
out of memory takes the process or the kernel with it, which `low_memory=True` helps avoid.

## Rust

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

The crate re-exports `nalgebra` and `malachite`, whose matrices and big integers it takes and
returns. Progress goes through the [`log`](https://docs.rs/log) facade; `cygv::log_to_stderr`
installs a minimal logger.

## Command line

The `cygv` executable reads a stream of YAML documents, one computation each, with the same names
as the Python arguments, and writes one document of results for each.

```bash
cargo install cygv
cygv --file input.yaml --output invariants.yaml
```

```yaml
---
name: an example threefold
q: [[0, 0, 1, 1, 1, 1], [1, 1, 0, 0, 0, -2]]
intnums: [[0, 0, 0, 8], [0, 0, 1, 4]]  # or {"0,0,0": 8, "0,0,1": 4}
generators: [[1, 0], [0, 1]]
invariants: gv
max_deg: 3
```

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

`--method`, `--threads`, `--max-deg`, `--invariants`, `--certify` and `--low-memory` override
those for every document, and `-v` (or `-vv`) writes the progress to stderr. See `cygv --help`
and the `examples` directory. The executable is the default `cli` feature; library users can
leave it out with `default-features = false`.

## License

Licensed under the GNU General Public License v3 or later (GPLv3+)
([LICENSE](https://github.com/ariostas/cygv/blob/main/LICENSE) or <https://www.gnu.org/licenses/gpl-3.0.en.html#license-text>).
Any contribution submitted for inclusion in the work shall be licensed under the GPLv3+, without
any additional terms or conditions; see
[CONTRIBUTING.md](https://github.com/ariostas/cygv/blob/main/CONTRIBUTING.md).

[1]: https://arxiv.org/abs/hep-th/9308122
[2]: https://arxiv.org/abs/hep-th/9406055
[Computational Mirror Symmetry]: https://arxiv.org/abs/2303.00757
