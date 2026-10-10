# cygv

[![Rust CI](https://github.com/ariostas/cygv/actions/workflows/rust.yml/badge.svg)](https://github.com/ariostas/cygv/actions/workflows/rust.yml) [![Python CI](https://github.com/ariostas/cygv/actions/workflows/python.yml/badge.svg)](https://github.com/ariostas/cygv/actions/workflows/python.yml) [![Crates.io Version](https://img.shields.io/crates/v/cygv)](https://crates.io/crates/cygv) [![PyPI - Version](https://img.shields.io/pypi/v/cygv)](https://pypi.org/project/cygv/) [![PyPI - Downloads](https://img.shields.io/pypi/dm/cygv)](https://pypi.org/project/cygv/)

This project implements an efficient algorithm to perform the HKTY procedure [[1], [2]] to compute Gopakumar-Vafa (GV) and Gromov-Witten (GW) invariants of Calabi-Yau (CY) manifolds obtained as hypersurfaces or complete intersections in toric varieties. This project is based on the work presented in the paper "[Computational Mirror Symmetry]", but written in the Rust programming language and with some additional improvements.

## Python package

The `cygv` package on PyPI ships prebuilt wheels, so no Rust toolchain is needed to install it.

```bash
pip install cygv
```

It exposes two functions, `compute_gv` and `compute_gw`, which take the data describing the CY and
return the corresponding invariants. Every vector argument is given as a list of vectors, or as any
array-like that numpy accepts.

```python
from cygv import compute_gv, compute_gw

# The CY threefold hypersurface in P(1,1,2,2,2), resolved, with h^{1,1} = 2.
# The generators of the Mori cone.
generators = [[1, 0], [0, 1]]
# A vector that has a positive inner product with every generator.
grading_vector = [1, 1]
# The GLSM charge matrix, given as one charge vector per Kähler parameter.
q = [[0, 0, 1, 1, 1, 1], [1, 1, 0, 0, 0, -2]]
# The triple intersection numbers. Only one ordering of each triplet may be given.
intnums = {(0, 0, 0): 8, (0, 0, 1): 4}

compute_gv(generators, grading_vector, q, intnums, max_deg=3)
# [((1, 0), 640), ((0, 1), 4), ((2, 0), 10032), ((1, 1), 640), ((3, 0), 288384), ...]

compute_gw(generators, grading_vector, q, intnums, max_deg=3)
# [((0, 1), Fraction(4, 1)), ((1, 0), Fraction(640, 1)), ((0, 2), Fraction(1, 2)), ...]
```

How many curve classes are computed is controlled by `max_deg`, `min_points`, or `target_points`,
of which at most one may be given. When none of them is given, the generators are taken to be the
complete list of curve classes to use. Complete intersections are specified by passing the nef
partition as `nefpart`. A vex phase of the ambient toric variety (MacFadden–Sheridan,
arXiv:2512.14817), whose curves of negative anticanonical degree need the class of the stratum
they force into the CY, is specified by passing those strata as `strata`, a list of
`(divisors, class)` pairs: the three toric divisors (indices of the rays) and the class of the
stratum, whose entries may be `Fraction`s. Strata are only supported for threefold
hypersurfaces. The computation is exact.

Two implementations of the computation are available, selected with `method`. The `modular` one
computes modulo word-sized primes and recovers the exact integers by the Chinese remainder theorem,
which is much faster, but it only supports `max_deg`. The `reference` one is a faithful
implementation of the HKTY procedure with exact rational arithmetic, and supports everything. The
default, `auto`, uses the modular one whenever it can. `n_threads` sets the number of threads,
which defaults to one per core (the modular implementation only uses one so far).

The results are returned as a list of `(curve_class, invariant)` pairs, in no particular order. GV
invariants are Python `int`s, and GW invariants are `fractions.Fraction`s. For CYs of dimension greater than three the invariants are labeled by a curve
class and by the index of a reference surface, which is the first index of the intersection
numbers, so each entry looks like `(((1, 0, -4, 0, 0, 0), 2), 80)` instead.

The computation runs in a subprocess, so it can be interrupted with ctrl+c without losing the
Python session.

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

The input is a stream of YAML documents, each of which specifies a single CY manifold.

```yaml
---
name: an example threefold
generators: [[1, 0], [0, 1]]
grading_vector: [1, 1]
q: [[0, 0, 1, 1, 1, 1], [1, 1, 0, 0, 0, -2]]
intnums:
  - [0, 0, 0, 8]
  - [0, 0, 1, 4]
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
is_threefold: true
grading_vector: [1, 1]
results:
  - {curve_class: [0, 1], degree: 1, gv: 4}
  - {curve_class: [1, 0], degree: 1, gv: 640}
  - {curve_class: [1, 1], degree: 2, gv: 640}
  - {curve_class: [2, 0], degree: 2, gv: 10032}
  - {curve_class: [2, 1], degree: 3, gv: 72224}
  - {curve_class: [3, 0], degree: 3, gv: 288384}
```

Run `cygv --help` for the full list of input fields and options, such as `method` (which can also
be set for every document with `--method`), and see the `examples` directory for complete inputs.

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
