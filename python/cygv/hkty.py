from __future__ import annotations

import operator
import os
import sys
import tempfile
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, fields
from fractions import Fraction
from multiprocessing import Pipe, Process
from multiprocessing.connection import Connection
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


def _forward(f: Any) -> None:
    """Writes what the subprocess wrote to its stderr since the last call to ours."""
    data = f.read()
    if data:
        sys.stderr.write(data.decode(errors="replace"))
        sys.stderr.flush()


# The computation runs in a subprocess, so that ctrl+c interrupts it without exiting the main
# Python process. Its stderr goes to a file, which is forwarded to ours as it grows if
# `verbose`, and otherwise only shown if it dies.
def _compute_in_subprocess(
    kind: str, args: dict[str, Any], *, verbose: bool = False
) -> Any:
    parent_conn, child_conn = Pipe(duplex=False)
    with tempfile.NamedTemporaryFile(delete=False) as f:
        stderr_path = Path(f.name)
    process = Process(
        target=_run_subprocess, args=(child_conn, stderr_path, kind, args)
    )
    try:
        process.start()
        child_conn.close()
        # The result, or the end of the pipe if the child died without sending it.
        try:
            if verbose:
                with stderr_path.open("rb") as f:
                    while not parent_conn.poll(0.2):
                        _forward(f)
                    _forward(f)
            result = parent_conn.recv()
        except EOFError:
            process.join()
            with stderr_path.open("rb") as f:
                stderr_msg = f.read().decode(errors="replace").strip()
            msg = f"Computation failed (exit code {process.exitcode})"
            if stderr_msg:
                msg += f":\n{stderr_msg}"
            raise RuntimeError(msg) from None
        process.join()
    # Interrupted, ctrl+c or otherwise: the computation stops with the caller.
    except BaseException:
        if process.pid is not None:
            process.terminate()
            process.join()
        raise
    finally:
        parent_conn.close()
        stderr_path.unlink(missing_ok=True)
    if isinstance(result, BaseException):
        raise result
    return result


def _int(x: Any, name: str) -> int:
    """An integer, from an int or an integral float."""
    try:
        return operator.index(x)
    except TypeError:
        pass
    if isinstance(x, float | np.floating) and float(x).is_integer():
        return int(x)
    msg = f"{name} must hold integers, not {x!r}"
    raise ValueError(msg)


def _count(x: Any, name: str) -> int | None:
    """A nonnegative integer below 2**32, or None."""
    if x is None:
        return None
    v = _int(x, name)
    if not 0 <= v < 2**32:
        msg = f"{name} must be a nonnegative int below 2**32, not {v}"
        raise ValueError(msg)
    return v


def _int_array(x: ArrayLike, name: str) -> Any:
    """An array of integers that fit in 32 bits."""
    try:
        array = np.asarray(x)
    except ValueError:
        msg = f"{name} must be a rectangular array of ints, with vectors of the same length"
        raise ValueError(msg) from None
    if array.dtype.kind == "f":
        if not np.all(np.isfinite(array) & (array == np.round(array))):
            msg = f"{name} must hold integers"
            raise ValueError(msg)
        array = array.astype(np.int64)
    elif array.dtype.kind not in "iub":
        values = [_int(v, name) for v in array.ravel()]
        if any(not -(2**31) <= v < 2**31 for v in values):
            msg = f"{name} must hold integers that fit in 32 bits"
            raise ValueError(msg)
        array = np.array(values, dtype=np.int64).reshape(array.shape)
    array = array.astype(np.int64)
    if array.size and (array.min() < -(2**31) or array.max() >= 2**31):
        msg = f"{name} must hold integers that fit in 32 bits"
        raise ValueError(msg)
    return array


def _vectors(vectors: ArrayLike | None, name: str) -> list[Any] | None:
    """A list of vectors, which may be empty, as a list of lists."""
    if vectors is None:
        return None
    array = _int_array(vectors, name)
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
    array = _int_array(points, name)
    if array.ndim == 1 and array.size > 0:
        array = array.reshape(1, -1)
    return _vectors(array, name)


def _vector(vector: ArrayLike | None, name: str) -> list[int] | None:
    if vector is None:
        return None
    array = _int_array(vector, name)
    if array.ndim != 1:
        msg = f"{name} must be a 1D array-like of ints"
        raise ValueError(msg)
    res: list[int] = array.tolist()
    return res


def _index_lists(
    lists: Iterable[ArrayLike] | None, name: str
) -> list[list[int]] | None:
    if lists is None:
        return None
    return [_int_array(x, name).ravel().tolist() for x in lists]


def _intnums(intnums: Mapping[Any, int] | Iterable[Sequence[int]]) -> dict[Any, int]:
    """The intersection numbers as a dict from index triples, from a dict or from a list of
    ``[i, j, k, value]`` entries."""
    if isinstance(intnums, Mapping):
        items = intnums.items()
    else:
        entries = [tuple(_int(x, "intnums") for x in e) for e in intnums]
        if any(len(e) != 4 for e in entries):
            msg = "intnums must be a dict or a list of [i, j, k, value] entries"
            raise ValueError(msg)
        items = ((e[:3], e[3]) for e in entries)  # type: ignore[assignment]
    res: dict[Any, int] = {}
    for key, value in items:
        indices = tuple(_int(i, "intnums") for i in key)
        if len(indices) != 3 or any(i < 0 for i in indices):
            msg = (
                "the indices of the intersection numbers must be three nonnegative ints"
            )
            raise ValueError(msg)
        if indices in res:
            msg = f"the intersection number {indices} is given twice"
            raise ValueError(msg)
        v = _int(value, "intnums")
        if not -(2**31) <= v < 2**31:
            msg = "the intersection numbers must fit in 32 bits"
            raise ValueError(msg)
        res[indices] = v
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
                _int_array(divisors, "strata").ravel().tolist(),
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
            _int_array(divisors, "hilbert_bases").ravel().tolist(),
            _vectors(basis, "hilbert_bases"),
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
        low_memory: bool = False,
        hilbert_bases: Iterable[tuple[ArrayLike, ArrayLike]] | None = None,
        verbose: bool | int = False,
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
        otherwise. ``certify`` proves the results exact. ``low_memory`` prefers memory to
        speed, for the runs that do not fit otherwise: the modular pipeline then takes about a
        third of the memory, for about two and a half times the time. ``n_threads`` defaults
        to one per core. ``hilbert_bases`` are the Hilbert bases of the cones ``K_T`` from a program such
        as Normaliz, as ``(divisors, basis)`` pairs, which the modular pipeline then uses
        instead of computing them. ``verbose`` writes the progress of the computation to
        ``sys.stderr``: with ``True`` (or 1) its stages and passes, with 2 also the degrees of
        each pass.

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
            low_memory=low_memory,
            hilbert_bases=hilbert_bases,
            verbose=verbose,
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
        low_memory: bool = False,
        hilbert_bases: Iterable[tuple[ArrayLike, ArrayLike]] | None = None,
        verbose: bool | int = False,
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
            low_memory=low_memory,
            hilbert_bases=hilbert_bases,
            verbose=verbose,
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
        low_memory: bool,
        hilbert_bases: Iterable[tuple[ArrayLike, ArrayLike]] | None,
        verbose: bool | int,
    ) -> Invariants:
        q = _vectors(self.q, "q")
        if not q:
            msg = "q must not be empty"
            raise ValueError(msg)
        args = {
            "q": q,
            "intnums": _intnums(self.intnums),
            "nefpart": _index_lists(self.nefpart, "nefpart"),
            "fan": _index_lists(self.fan, "fan"),
            "strata": _strata(self.strata),
            "generators": _vectors(self.generators, "generators"),
            "mori_hyperplanes": _vectors(self.mori_hyperplanes, "mori_hyperplanes"),
            # No equations is the same as none given.
            "mori_equations": _vectors(self.mori_equations, "mori_equations") or None,
            "grading_vector": _vector(self.grading_vector, "grading_vector"),
            "max_deg": _count(max_deg, "max_deg"),
            "target_points": _points(target_points, "target_points"),
            "min_points": _count(min_points, "min_points"),
            "curve_classes": _vectors(curve_classes, "curve_classes"),
            "method": method,
            "n_threads": _count(n_threads, "n_threads"),
            "certify": certify,
            "low_memory": low_memory,
            "hilbert_bases": _hilbert_bases(hilbert_bases),
            "verbose": _count(verbose, "verbose"),
        }
        entries, dimension, grading_vector, report = _compute_in_subprocess(
            kind, args, verbose=bool(verbose)
        )
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
