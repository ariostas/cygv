from __future__ import annotations

import os
import tempfile
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, fields
from fractions import Fraction
from multiprocessing import Pipe, Process
from multiprocessing.connection import Connection, wait
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from cygv.cygv import _compute

# An invariant is labeled by its curve class, and for an n-fold also by its reference surface.
Key = tuple[int, ...] | tuple[tuple[int, ...], int]

METHODS = ("auto", "modular", "reference")


@dataclass(frozen=True)
class Report:
    """How the invariants were computed."""

    #: The pipeline that ran, "modular" or "reference".
    method: str
    #: Whether the results are proven exact (the reference pipeline is exact; the modular one
    #: with ``certify=True``).
    certified: bool
    #: Whether the grading vector is the default rather than given.
    grading_vector_derived: bool
    #: The wall-clock time of the computation, in seconds.
    seconds: float


class Invariants(list):  # type: ignore[type-arg]
    """The nonzero invariants, as ``(key, value)`` pairs sorted by degree and then by curve
    class, where the key is the curve class (a tuple), or for an n-fold the pair of the curve
    class and the reference surface. ``dict(invariants)`` maps the keys to the values.

    Besides the pairs it holds what is needed to read them: ``kind`` ("gv" or "gw"),
    ``dimension``, ``grading_vector`` and the ``report``.
    """

    kind: str
    dimension: int
    grading_vector: tuple[int, ...]
    report: Report

    def __init__(
        self,
        pairs: Iterable[tuple[Key, int | Fraction]],
        kind: str,
        dimension: int,
        grading_vector: tuple[int, ...],
        report: Report,
    ) -> None:
        super().__init__(pairs)
        self.kind = kind
        self.dimension = dimension
        self.grading_vector = grading_vector
        self.report = report

    def degree(self, curve: Sequence[int]) -> int:
        """The degree of a curve class under the grading vector."""
        return sum(int(c) * g for c, g in zip(curve, self.grading_vector, strict=True))

    def __reduce__(self) -> Any:
        return (
            Invariants,
            (list(self), self.kind, self.dimension, self.grading_vector, self.report),
        )


def _run_subprocess(
    conn: Connection, stderr_path: Path, kind: str, args: dict[str, Any]
) -> None:
    with stderr_path.open("wb") as f:
        os.dup2(f.fileno(), 2)
    try:
        conn.send(_compute(kind, args))
    # Errors that the Rust side reports, such as invalid input, keep their type.
    except Exception as e:  # noqa: BLE001
        conn.send(e)
    # The catch has to be this broad: a panic on the Rust side surfaces as a
    # PanicException, which derives from BaseException rather than Exception, and
    # it still needs to be reported back to the parent process.
    except BaseException as e:  # noqa: BLE001
        conn.send(RuntimeError(str(e)))
    conn.close()


# The computation runs in a subprocess, so that ctrl+c interrupts it without exiting the main
# Python process.
def _compute_in_subprocess(kind: str, args: dict[str, Any]) -> Any:
    parent_conn, child_conn = Pipe(duplex=False)
    with tempfile.NamedTemporaryFile(delete=False) as f:
        stderr_path = Path(f.name)
    process = Process(
        target=_run_subprocess, args=(child_conn, stderr_path, kind, args)
    )
    process.start()
    child_conn.close()

    ready = wait([parent_conn, process.sentinel])
    if parent_conn in ready:
        result = parent_conn.recv()
        process.join()
        stderr_path.unlink()
        if isinstance(result, BaseException):
            raise result
        return result

    process.join()
    with stderr_path.open("rb") as f:
        stderr_msg = f.read().decode(errors="replace").strip()
    stderr_path.unlink()
    msg = f"Computation failed (exit code {process.exitcode})"
    if stderr_msg:
        msg += f":\n{stderr_msg}"
    raise RuntimeError(msg)


def _vectors(vectors: ArrayLike | None, name: str) -> list[Any] | None:
    """A list of vectors, which may be empty, as a list of lists."""
    if vectors is None:
        return None
    array = np.array(vectors, dtype=int)
    if array.size == 0:
        return []
    if array.ndim != 2:
        msg = f"{name} must be a 2D array-like of ints"
        raise ValueError(msg)
    res: list[Any] = array.tolist()
    return res


def _points(points: ArrayLike | None, name: str) -> list[Any] | None:
    """One point or a list of points as a list of lists."""
    if points is None:
        return None
    array = np.array(points, dtype=int)
    if array.ndim == 1 and array.size > 0:
        array = array.reshape(1, -1)
    return _vectors(array, name)


def _vector(vector: ArrayLike | None, name: str) -> list[int] | None:
    if vector is None:
        return None
    array = np.array(vector, dtype=int)
    if array.ndim != 1:
        msg = f"{name} must be a 1D array-like of ints"
        raise ValueError(msg)
    res: list[int] = array.tolist()
    return res


def _index_lists(lists: Iterable[ArrayLike] | None) -> list[list[int]] | None:
    if lists is None:
        return None
    return [[int(i) for i in np.asarray(x).ravel()] for x in lists]


def _intnums(intnums: Mapping[Any, int] | Iterable[Sequence[int]]) -> dict[Any, int]:
    """The intersection numbers as a dict from index triples, from a dict or from a list of
    ``[i, j, k, value]`` entries."""
    if isinstance(intnums, Mapping):
        items = intnums.items()
    else:
        entries = [tuple(int(x) for x in e) for e in intnums]
        if any(len(e) != 4 for e in entries):
            msg = "intnums must be a dict or a list of [i, j, k, value] entries"
            raise ValueError(msg)
        items = ((e[:3], e[3]) for e in entries)  # type: ignore[assignment]
    res: dict[Any, int] = {}
    for key, value in items:
        indices = tuple(int(i) for i in key)
        if len(indices) != 3 or any(i < 0 for i in indices):
            msg = (
                "the indices of the intersection numbers must be three nonnegative ints"
            )
            raise ValueError(msg)
        if indices in res:
            msg = f"the intersection number {indices} is given twice"
            raise ValueError(msg)
        res[indices] = int(value)
    return res


def _strata(
    strata: Iterable[tuple[ArrayLike, Iterable[Any]]] | None,
) -> list[Any] | None:
    """The strata as (divisors, class) pairs, with the class as (num, den) pairs."""
    if strata is None:
        return None
    res = []
    for divisors, cls in strata:
        fracs = [Fraction(x) for x in cls]
        res.append(
            (
                [int(d) for d in np.asarray(divisors).ravel()],
                [(f.numerator, f.denominator) for f in fracs],
            )
        )
    return res


def _hilbert_bases(
    bases: Iterable[tuple[ArrayLike, ArrayLike]] | None,
) -> list[Any] | None:
    if bases is None:
        return None
    return [
        (
            [int(d) for d in np.asarray(divisors).ravel()],
            _vectors(basis, "the Hilbert bases"),
        )
        for divisors, basis in bases
    ]


@dataclass
class Problem:
    """A Calabi-Yau hypersurface or complete intersection in a toric variety, and its Mori
    cone, whose GV or GW invariants :meth:`gv` and :meth:`gw` compute.

    The vectors are given as lists (or array-likes) of ints, one vector per inner list:

    - ``q``: the GLSM charges, one charge vector per Kähler parameter. A divisor index is a
      position in these vectors.
    - ``intnums``: the intersection numbers, as a dict ``{(i, j, k): value}`` or a list of
      ``[i, j, k, value]`` entries, each triple in one ordering. For an n-fold, ``i`` is a
      reference surface.
    - ``nefpart``: the nef partition of a complete intersection, as lists of divisor
      indices; ``None`` for a hypersurface.
    - ``fan``: the maximal cones of the fan of the ambient toric variety, as lists of divisor
      indices (threefold hypersurfaces only). It gives the strata of a vex phase, the Mori
      cone and the default grading vector, and ``intnums`` is checked against it.
    - ``strata``: instead of the fan, the strata of a vex phase, as ``(divisors, class)``
      pairs, with three divisors and a class whose entries may be fractions.
    - The Mori cone, one of: ``generators``, vectors that span it (every lattice point of the
      cone they span is used); ``mori_hyperplanes``, the inequalities ``h . C >= 0`` that cut
      it out, with ``mori_equations``, the equations of its span if it is a face; or neither,
      with a fan.
    - ``grading_vector``: the degree is the inner product with it, which must be positive on
      the Mori cone. By default, the primitive sum of the inequalities of the Mori cone.
    """

    q: ArrayLike
    intnums: Mapping[Any, int] | Iterable[Sequence[int]]
    nefpart: Iterable[ArrayLike] | None = None
    fan: Iterable[ArrayLike] | None = None
    strata: Iterable[tuple[ArrayLike, Iterable[Any]]] | None = None
    generators: ArrayLike | None = None
    mori_hyperplanes: ArrayLike | None = None
    mori_equations: ArrayLike | None = None
    grading_vector: ArrayLike | None = None

    def gv(
        self,
        *,
        max_deg: int | None = None,
        target_points: ArrayLike | None = None,
        min_points: int | None = None,
        curve_classes: ArrayLike | None = None,
        method: str | None = None,
        n_threads: int | None = None,
        certify: bool = False,
        hilbert_bases: Iterable[tuple[ArrayLike, ArrayLike]] | None = None,
    ) -> Invariants:
        """The GV invariants, as integers, of the curve classes selected by one of:

        - ``max_deg``: those of degree at most this;
        - ``target_points``: those below the given classes ``P`` (the ``C`` with ``P - C`` in
          the Mori cone), up to their largest degree, or ``max_deg`` if it is also given and
          lower;
        - ``min_points``: the lowest degrees that hold at least this many classes;
        - ``curve_classes``: exactly the classes listed (zero is added if it is missing).

        ``method`` is "auto" (the default), "modular" or "reference": auto runs the much
        faster modular pipeline for ``max_deg`` and ``target_points``, and the reference one
        otherwise. ``certify`` proves the results exact. ``n_threads`` defaults to one per
        core. ``hilbert_bases`` are the Hilbert bases of the cones ``K_T`` from a program such
        as Normaliz, as ``(divisors, basis)`` pairs, which the modular pipeline then uses
        instead of computing them.

        Invalid input raises a ``ValueError``, and a computation that fails a
        :class:`ComputationError`.
        """
        return self._compute(
            "gv",
            max_deg=max_deg,
            target_points=target_points,
            min_points=min_points,
            curve_classes=curve_classes,
            method=method,
            n_threads=n_threads,
            certify=certify,
            hilbert_bases=hilbert_bases,
        )

    def gw(
        self,
        *,
        max_deg: int | None = None,
        target_points: ArrayLike | None = None,
        min_points: int | None = None,
        curve_classes: ArrayLike | None = None,
        method: str | None = None,
        n_threads: int | None = None,
        certify: bool = False,
        hilbert_bases: Iterable[tuple[ArrayLike, ArrayLike]] | None = None,
    ) -> Invariants:
        """The GW invariants, as fractions, of the curve classes selected as for :meth:`gv`,
        which describes the other arguments too."""
        return self._compute(
            "gw",
            max_deg=max_deg,
            target_points=target_points,
            min_points=min_points,
            curve_classes=curve_classes,
            method=method,
            n_threads=n_threads,
            certify=certify,
            hilbert_bases=hilbert_bases,
        )

    def _compute(
        self,
        kind: str,
        *,
        max_deg: int | None,
        target_points: ArrayLike | None,
        min_points: int | None,
        curve_classes: ArrayLike | None,
        method: str | None,
        n_threads: int | None,
        certify: bool,
        hilbert_bases: Iterable[tuple[ArrayLike, ArrayLike]] | None,
    ) -> Invariants:
        q = _vectors(self.q, "q")
        if not q:
            msg = "q must not be empty"
            raise ValueError(msg)
        args = {
            "q": q,
            "intnums": _intnums(self.intnums),
            "nefpart": _index_lists(self.nefpart),
            "fan": _index_lists(self.fan),
            "strata": _strata(self.strata),
            "generators": _vectors(self.generators, "generators"),
            "mori_hyperplanes": _vectors(self.mori_hyperplanes, "mori_hyperplanes"),
            "mori_equations": _vectors(self.mori_equations, "mori_equations"),
            "grading_vector": _vector(self.grading_vector, "grading_vector"),
            "max_deg": max_deg,
            "target_points": _points(target_points, "target_points"),
            "min_points": min_points,
            "curve_classes": _vectors(curve_classes, "curve_classes"),
            "method": method,
            "n_threads": n_threads,
            "certify": certify,
            "hilbert_bases": _hilbert_bases(hilbert_bases),
        }
        entries, dimension, grading_vector, report = _compute_in_subprocess(kind, args)
        value = int if kind == "gv" else Fraction
        if dimension == 3:
            pairs = [(tuple(c), value(v)) for ((c, _), v) in entries]
        else:
            pairs = [((tuple(c), s), value(v)) for ((c, s), v) in entries]
        return Invariants(
            pairs, kind, dimension, tuple(grading_vector), Report(*report)
        )


_PROBLEM_FIELDS = tuple(f.name for f in fields(Problem))


def _split(kwargs: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """The keywords of :class:`Problem`, and the others."""
    problem = {k: v for k, v in kwargs.items() if k in _PROBLEM_FIELDS}
    rest = {k: v for k, v in kwargs.items() if k not in _PROBLEM_FIELDS}
    return problem, rest


def compute_gv(q: ArrayLike, intnums: Any, **kwargs: Any) -> Invariants:
    """The GV invariants: ``Problem(q, intnums, ...).gv(...)``, with the keywords of
    :class:`Problem` and of :meth:`Problem.gv` together."""
    problem, rest = _split(kwargs)
    return Problem(q, intnums, **problem).gv(**rest)


def compute_gw(q: ArrayLike, intnums: Any, **kwargs: Any) -> Invariants:
    """The GW invariants: ``Problem(q, intnums, ...).gw(...)``, with the keywords of
    :class:`Problem` and of :meth:`Problem.gw` together."""
    problem, rest = _split(kwargs)
    return Problem(q, intnums, **problem).gw(**rest)
