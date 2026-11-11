from __future__ import annotations

import operator
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import KW_ONLY, dataclass, fields
from fractions import Fraction
from functools import cached_property
from typing import Any, cast, overload

import numpy as np
from numpy.typing import ArrayLike, NDArray

from cygv.cygv import _compute

# An invariant is labeled by its curve class, and for an n-fold also by its reference surface.
Key = tuple[int, ...] | tuple[tuple[int, ...], int]

METHODS = ("auto", "modular", "reference")
VERIFY = ("off", "quick", "full")


@dataclass(frozen=True)
class Report:
    """How the invariants were computed."""

    #: The pipeline that ran, "modular" or "reference".
    method: str
    #: How the values are known to be those the input determines, which is about the
    #: arithmetic, not about the input: "exact" (the reference pipeline), "stable" (the
    #: modular one, whose lift from residues one more prime agreed with) or "certified" (the
    #: modular one with ``certify=True``, which bounds every invariant).
    lift: str
    #: Which check of the consistency of the input ran: "quick" or "full" (see ``verify`` in
    #: :meth:`Problem.gv`; the reference pipeline checks fully), or "unchecked" (turned off,
    #: or for n-folds and ``h11 = 1``, where it does not apply).
    consistency: str
    #: Whether the grading vector is the default rather than given.
    grading_vector_derived: bool
    #: The wall-clock time of the computation, in seconds.
    seconds: float


# The pairs are typed loosely, as by the list this replaced: a key is a curve tuple or a pair
# of one and a surface, a value an int or a Fraction.
class Invariants(Sequence[tuple[Any, Any]]):
    """The nonzero invariants, as a sequence of ``(key, value)`` pairs sorted by degree and
    then by curve class, where the key is the curve class (a tuple), or for an n-fold the pair
    of the curve class and the reference surface. ``dict(invariants)`` maps the keys to the
    values.

    ``invariants[key]`` and ``invariants.get(key)`` look up the value of a key (with a lookup
    table built on the first use), and an int index gives the pair at that position.

    The same data are available as numpy arrays, which is much faster for many invariants:
    ``curves``, one row per curve class; ``surfaces``, the reference surfaces of an n-fold
    (``None`` for a threefold); ``degrees``; and ``values``, of ``int64`` when they all fit,
    and of Python ints otherwise (``Fraction``s for GW invariants). Besides those it holds
    what is needed to read them: ``kind`` ("gv" or "gw"), ``dimension``, ``grading_vector``
    and the ``report``.
    """

    kind: str
    dimension: int
    grading_vector: tuple[int, ...]
    report: Report
    #: The curve classes, one per row.
    curves: NDArray[np.int32]
    #: The reference surfaces of an n-fold, ``None`` for a threefold.
    surfaces: NDArray[np.int64] | None
    #: The degrees of the curve classes under the grading vector.
    degrees: NDArray[np.int64]

    def __init__(
        self,
        kind: str,
        curves: NDArray[np.int32],
        surfaces: NDArray[np.int64] | None,
        degrees: NDArray[np.int64],
        values: NDArray[np.int64] | list[int],
        denominators: NDArray[np.int64] | list[int] | None,
        dimension: int,
        grading_vector: tuple[int, ...],
        report: Report,
    ) -> None:
        self.kind = kind
        self.curves = curves
        self.surfaces = surfaces
        self.degrees = degrees
        # The values, or the numerators of GW invariants, and their denominators.
        self._values = values
        self._denominators = denominators
        self.dimension = dimension
        self.grading_vector = grading_vector
        self.report = report
        self._index: dict[Key, int] | None = None
        for array in (curves, surfaces, degrees, values, denominators):
            if isinstance(array, np.ndarray):
                array.flags.writeable = False

    @cached_property
    def values(self) -> NDArray[Any]:
        """The invariants, in an array of ``int64`` when they are integers that all fit, and
        of Python objects otherwise."""
        if self._denominators is None and isinstance(self._values, np.ndarray):
            return self._values
        values = np.empty(len(self), dtype=object)
        values[:] = self._value_list()
        values.flags.writeable = False
        return values

    def _value_list(self) -> list[int] | list[Fraction]:
        """The invariants as Python ints or fractions."""
        values = _to_list(self._values)
        if self._denominators is None:
            return values
        denominators = _to_list(self._denominators)
        return [Fraction(n, d) for n, d in zip(values, denominators, strict=True)]

    # Not `keys`, which would have `dict()` take the invariants for a mapping.
    def _keys(self) -> list[Key]:
        """The keys of the invariants: the curve classes as tuples, or for an n-fold the pairs
        of the curve class and the reference surface."""
        curves = map(tuple, self.curves.tolist())
        if self.surfaces is None:
            return list(curves)
        return list(zip(curves, self.surfaces.tolist(), strict=True))

    def __len__(self) -> int:
        return len(self.degrees)

    def __iter__(self) -> Iterator[tuple[Any, Any]]:
        return iter(zip(self._keys(), self._value_list(), strict=True))

    @overload
    def __getitem__(self, i: int) -> tuple[Any, Any]: ...

    @overload
    def __getitem__(self, i: slice) -> list[tuple[Any, Any]]: ...

    @overload
    def __getitem__(self, i: Sequence[Any] | NDArray[Any]) -> Any: ...

    def __getitem__(
        self, i: int | slice | Sequence[Any] | NDArray[Any]
    ) -> tuple[Any, Any] | list[tuple[Any, Any]] | Any:
        # A key can be any sequence of ints (a list, an array) as well as a tuple.
        if not isinstance(i, int | np.integer | slice):
            if self._index is None:
                self._index = {k: j for j, k in enumerate(self._keys())}
            try:
                k = self._key(tuple(i))
            except (TypeError, ValueError):
                raise KeyError(i) from None
            j = self._index.get(k)
            if j is None:
                raise KeyError(i)
            return self[j][1]
        if isinstance(i, slice):
            return [self[j] for j in range(*i.indices(len(self)))]
        i = range(len(self))[i]
        curve = tuple(self.curves[i].tolist())
        key: Key = curve if self.surfaces is None else (curve, int(self.surfaces[i]))
        value: int | Fraction = int(self._values[i])
        if self._denominators is not None:
            value = Fraction(value, int(self._denominators[i]))
        return key, value

    def __contains__(self, key: object) -> bool:
        """Whether ``key`` is the key of a nonzero invariant here, as ``invariants[key]``
        looks it up (a curve class, or for an n-fold a ``(curve, surface)`` pair, from any
        sequences of ints), or one of the ``(key, value)`` pairs that iterating gives."""
        if isinstance(key, int | np.integer | slice | str):
            return False
        if isinstance(key, tuple) and self._is_pair(key):
            return any(p == key for p in self)
        try:
            self[cast(Any, key)]
        except KeyError:
            return False
        return True

    def _is_pair(self, x: tuple[Any, ...]) -> bool:
        """Whether ``x`` is a ``(key, value)`` pair rather than a key: its first entry is
        a key, which is a sequence, of ints for a threefold and a pair for an n-fold."""
        if len(x) != 2 or isinstance(x[0], int | np.integer | str):
            return False
        if self.surfaces is None:
            return True
        first = x[0]
        return (
            isinstance(first, Sequence)
            and len(first) == 2
            and not isinstance(first[0], int | np.integer)
        )

    def get(self, key: Any, default: Any = None) -> Any:
        """The value of a key (a curve class, or for an n-fold the pair of a curve class and a
        reference surface), or ``default`` if it has no nonzero invariant here."""
        try:
            return self[tuple(key)]
        except (KeyError, TypeError):
            return default

    def _key(self, key: tuple[Any, ...]) -> Key:
        """A key as :meth:`_keys` gives them, from any sequences of ints."""
        if self.surfaces is None:
            return tuple(_int(x, "a curve class") for x in key)
        if len(key) != 2:
            msg = "the keys of an n-fold's invariants are (curve, surface) pairs"
            raise KeyError(msg)
        return (tuple(int(x) for x in key[0]), int(key[1]))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Sequence):
            return NotImplemented
        return len(self) == len(other) and list(self) == list(other)

    # Unhashable, as the list this replaced, since it compares by its pairs.
    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        head = ", ".join(repr(p) for p in self[:5])
        more = f", ... ({len(self)} in all)" if len(self) > 5 else ""
        return f"Invariants({self.kind}, [{head}{more}])"

    def degree(self, curve: Sequence[int]) -> int:
        """The degree of a curve class under the grading vector."""
        if len(curve) != len(self.grading_vector):
            msg = f"the curve class must have {len(self.grading_vector)} entries"
            raise ValueError(msg)
        return sum(int(c) * g for c, g in zip(curve, self.grading_vector, strict=True))

    def __reduce__(self) -> Any:
        return (
            Invariants,
            (
                self.kind,
                self.curves,
                self.surfaces,
                self.degrees,
                self._values,
                self._denominators,
                self.dimension,
                self.grading_vector,
                self.report,
            ),
        )


def _to_list(values: NDArray[np.int64] | list[int]) -> list[int]:
    """Python ints from an array or a list of them."""
    if isinstance(values, np.ndarray):
        res: list[int] = values.tolist()
        return res
    return values


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
    elif array.dtype.kind not in "iub":
        values = [_int(v, name) for v in array.ravel()]
        if any(not -(2**31) <= v < 2**31 for v in values):
            msg = f"{name} must hold integers that fit in 32 bits"
            raise ValueError(msg)
        array = np.array(values, dtype=np.int64).reshape(array.shape)
    # Checked before the cast, which would wrap large unsigned and float values.
    if array.size and (array.min() < -(2**31) or array.max() >= 2**31):
        msg = f"{name} must hold integers that fit in 32 bits"
        raise ValueError(msg)
    return array.astype(np.int64)


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


def _indices(x: ArrayLike, name: str) -> list[int]:
    """A list of indices, which must not be negative."""
    res: list[int] = _int_array(x, name).ravel().tolist()
    if any(i < 0 for i in res):
        msg = f"{name} must hold nonnegative indices"
        raise ValueError(msg)
    return res


def _index_lists(
    lists: Iterable[ArrayLike] | None, name: str
) -> list[list[int]] | None:
    if lists is None:
        return None
    return [_indices(x, name) for x in lists]


def _flag(x: Any, name: str) -> bool:
    """A bool, from a bool (numpy's too) or the int 0 or 1."""
    if isinstance(x, bool | np.bool_):
        return bool(x)
    try:
        v = operator.index(x)
    except TypeError:
        v = None
    if v not in (0, 1):
        msg = f"{name} must be True or False, not {x!r}"
        raise ValueError(msg)
    return bool(v)


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
        if any(not (abs(f.numerator) < 2**63 and f.denominator < 2**63) for f in fracs):
            msg = "the classes of the strata must be fractions of integers that fit in 64 bits"
            raise ValueError(msg)
        res.append(
            (
                _indices(divisors, "strata"),
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
            _indices(divisors, "hilbert_bases"),
            _vectors(basis, "hilbert_bases"),
        )
        for divisors, basis in bases
    ]


def _method(method: Any) -> str | None:
    if method is not None and not isinstance(method, str):
        msg = f'method must be "auto", "modular" or "reference", not {method!r}'
        raise ValueError(msg)
    return method


def _verify(verify: Any) -> str:
    """The consistency check, from its name or a bool."""
    if isinstance(verify, bool | np.bool_):
        return "quick" if verify else "off"
    if not isinstance(verify, str) or verify.lower() not in VERIFY:
        msg = f'verify must be "off", "quick" or "full" (or a bool), not {verify!r}'
        raise ValueError(msg)
    return verify.lower()


def _verbosity(verbose: Any) -> int:
    """The verbosity, 0 to 2, from a bool or an int (larger ones are 2)."""
    if isinstance(verbose, bool | np.bool_):
        return int(verbose)
    v = _count(verbose, "verbose")
    return min(v or 0, 2)


# Not compared by value, since the fields may be arrays.
@dataclass(eq=False)
class Problem:
    """A Calabi-Yau hypersurface or complete intersection in a toric variety, and its Mori
    cone, whose GV or GW invariants :meth:`gv` and :meth:`gw` compute.

    ``q`` and ``intnums`` can be given by position, the others only by keyword. The vectors
    are given as lists (or array-likes) of ints, one vector per inner list:

    - ``q``: the GLSM charges, one charge vector per Kähler parameter. A divisor index is a
      position in these vectors.
    - ``intnums``: the intersection numbers, as a dict ``{(i, j, k): value}`` or a list of
      ``[i, j, k, value]`` entries, each triple in one ordering. For an n-fold, ``i`` is a
      reference surface.
    - ``nefpart``: the nef partition of a complete intersection, as lists of divisor
      indices; ``None`` for a hypersurface.
    - ``dimension``: optionally, the dimension of the CY, which is then checked against the
      one that ``q`` and ``nefpart`` imply (a complete intersection given without its nef
      partition would be taken for a hypersurface of higher dimension).
    - ``fan``: the maximal cones of the fan of the ambient toric variety, as lists of divisor
      indices (threefold hypersurfaces only). It gives the strata of a vex phase, the Mori
      cone and the default grading vector, and ``intnums`` is checked against it.
    - ``strata``: instead of the fan, the strata of a vex phase, as ``(divisors, class)``
      pairs, with three divisors and a class whose entries may be fractions.
    - The Mori cone, one of: ``generators``, vectors that span it (every lattice point of the
      cone they span is used); ``mori_hyperplanes``, the inequalities ``h . C >= 0`` that cut
      it out, with ``mori_equations``, the equations of its span if it is a face; or neither,
      with a fan.
    - ``hilbert_bases``: optionally, the Hilbert bases of the cones ``K_T`` of the Mori cone
      from a program such as Normaliz, as ``(divisors, basis)`` pairs, which the modular
      pipeline then uses instead of computing them. They are trusted, so a missing element
      silently drops classes.
    - ``grading_vector``: the degree is the inner product with it, which must be positive on
      the Mori cone. By default, the primitive sum of the inequalities of the Mori cone.
    """

    q: ArrayLike
    intnums: Mapping[Any, int] | Iterable[Sequence[int]]
    # The rest by keyword only: given by position, the generators would be taken for the nef
    # partition.
    _: KW_ONLY
    nefpart: Iterable[ArrayLike] | None = None
    dimension: int | None = None
    fan: Iterable[ArrayLike] | None = None
    strata: Iterable[tuple[ArrayLike, Iterable[Any]]] | None = None
    generators: ArrayLike | None = None
    mori_hyperplanes: ArrayLike | None = None
    mori_equations: ArrayLike | None = None
    hilbert_bases: Iterable[tuple[ArrayLike, ArrayLike]] | None = None
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
        verify: str | bool = "quick",
        low_memory: bool = False,
        verbose: bool | int = False,
    ) -> Invariants:
        """The GV invariants, as integers, of the curve classes selected by one of:

        - ``max_deg``: those of degree at most this;
        - ``target_points``: those below the given classes ``P`` (the ``C`` with ``P - C`` in
          the Mori cone), up to their largest degree, or ``max_deg`` if it is also given and
          lower;
        - ``min_points``: the lowest degrees that hold at least this many classes;
        - ``curve_classes``: exactly the classes listed (computed as ``target_points``, whose
          lightcones their invariants depend on, keeping only these).

        ``method`` is "auto" (the default), "modular" or "reference": auto runs the much
        faster modular pipeline for every selection but ``min_points``, and the reference one
        for that.

        ``verify`` checks that the input is consistent, for a threefold with ``h11 >= 2``: the
        modular pipeline only uses the intersection numbers contracted with the grading
        vector, so wrong ones (or a Mori cone too small, or a vex phase without its fan or
        strata) can give wrong integers, which a second pass with another contraction then
        detects, raising a :class:`ComputationError`. "quick" (the default, or ``True``) checks
        up to about half the maximum degree, for a few percent of the time, "full" up to the
        maximum degree, for about half a pass more, and "off" (or ``False``) not at all. It
        cannot detect intersection numbers all scaled by a common factor, nor anything for
        ``h11 = 1`` or an n-fold. The reference pipeline always checks.

        ``certify`` proves that the lift from residues is right, by bounding every invariant
        and taking primes until the bounds fit, for about twice the time; it says nothing
        about the input, and without it a wrong lift is already vanishingly unlikely. It is
        not supported for vex phases yet. ``low_memory`` prefers memory to
        speed, for the runs that do not fit otherwise: the modular pipeline then takes about a
        third of the memory, for about two and a half times the time. ``n_threads`` defaults
        to one per logical CPU. ``verbose`` writes the progress of the computation to
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
            verify=verify,
            low_memory=low_memory,
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
        verify: str | bool = "quick",
        low_memory: bool = False,
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
            verify=verify,
            low_memory=low_memory,
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
        verify: str | bool,
        low_memory: bool,
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
            "dimension": _count(self.dimension, "dimension"),
            "fan": _index_lists(self.fan, "fan"),
            "strata": _strata(self.strata),
            "generators": _vectors(self.generators, "generators"),
            "mori_hyperplanes": _vectors(self.mori_hyperplanes, "mori_hyperplanes"),
            # No equations is the same as none given.
            "mori_equations": _vectors(self.mori_equations, "mori_equations") or None,
            "hilbert_bases": _hilbert_bases(self.hilbert_bases),
            "grading_vector": _vector(self.grading_vector, "grading_vector"),
            "max_deg": _count(max_deg, "max_deg"),
            "target_points": _points(target_points, "target_points"),
            "min_points": _count(min_points, "min_points"),
            "curve_classes": _points(curve_classes, "curve_classes"),
            "method": _method(method),
            "n_threads": _count(n_threads, "n_threads"),
            "certify": _flag(certify, "certify"),
            "verify": _verify(verify),
            "low_memory": _flag(low_memory, "low_memory"),
            "verbose": _verbosity(verbose),
        }
        (
            curves,
            surfaces,
            degrees,
            values,
            denominators,
            dimension,
            grading_vector,
            report,
        ) = _compute(kind, args)
        if kind == "gw" and denominators is None:
            denominators = []
        return Invariants(
            kind,
            curves,
            surfaces,
            degrees,
            values,
            denominators,
            dimension,
            tuple(grading_vector),
            Report(*report),
        )


# The keywords of cygv 0.3 that are gone, with what to write instead.
_REMOVED = {
    "prec": "the computation is always exact now; leave prec out",
    "cones": "it is called hilbert_bases now, and the Mori cone it came with is given by "
    "mori_hyperplanes alone",
    "is_threefold": "the dimension is derived from q and nefpart; leave it out, or give "
    "dimension to have it checked",
    "find_gv": "call compute_gv or compute_gw",
}

# The order of the positional arguments of cygv 0.3.
_OLD_ORDER = (
    "cygv 0.3 took (generators, grading_vector, q, intnums, max_deg, min_points, "
    "target_points, nefpart, prec) by position; now q and intnums come first and everything "
    "else is a keyword, e.g. compute_gv(q, intnums, generators=generators, "
    "grading_vector=grading_vector, max_deg=10)"
)


def _compute_checked(
    kind: str,
    q: ArrayLike,
    intnums: Any,
    args: tuple[Any, ...],
    keywords: dict[str, Any],
    removed: dict[str, Any],
) -> Invariants:
    """The invariants for :func:`compute_gv` and :func:`compute_gw`, with errors that say
    what to write instead of the arguments of cygv 0.3."""
    name = f"compute_{kind}"
    if args:
        msg = f"{name}() takes 2 positional arguments but {2 + len(args)} were given; {_OLD_ORDER}"
        raise TypeError(msg)
    for key in removed:
        hint = _REMOVED.get(key)
        msg = f"{name}() got an unexpected keyword argument {key!r}"
        raise TypeError(f"{msg} ({hint})" if hint else msg)
    problem = {k: keywords.pop(k) for k in _PROBLEM_FIELDS if k not in ("q", "intnums")}
    selection = ("max_deg", "target_points", "min_points", "curve_classes")
    if (
        all(keywords[k] is None for k in selection)
        and problem["generators"] is not None
    ):
        msg = (
            "which curve classes: give max_deg, target_points, min_points or curve_classes "
            "(cygv 0.3 computed the invariants at the generators when none of the first "
            "three was given, which is curve_classes=generators now)"
        )
        raise ValueError(msg)
    p = Problem(q, intnums, **problem)
    return p.gv(**keywords) if kind == "gv" else p.gw(**keywords)


_PROBLEM_FIELDS = tuple(f.name for f in fields(Problem))


def compute_gv(
    q: ArrayLike,
    intnums: Mapping[Any, int] | Iterable[Sequence[int]],
    *args: Any,
    nefpart: Iterable[ArrayLike] | None = None,
    dimension: int | None = None,
    fan: Iterable[ArrayLike] | None = None,
    strata: Iterable[tuple[ArrayLike, Iterable[Any]]] | None = None,
    generators: ArrayLike | None = None,
    mori_hyperplanes: ArrayLike | None = None,
    mori_equations: ArrayLike | None = None,
    hilbert_bases: Iterable[tuple[ArrayLike, ArrayLike]] | None = None,
    grading_vector: ArrayLike | None = None,
    max_deg: int | None = None,
    target_points: ArrayLike | None = None,
    min_points: int | None = None,
    curve_classes: ArrayLike | None = None,
    method: str | None = None,
    n_threads: int | None = None,
    certify: bool = False,
    verify: str | bool = "quick",
    low_memory: bool = False,
    verbose: bool | int = False,
    **removed: Any,
) -> Invariants:
    """The GV invariants: ``Problem(q, intnums, ...).gv(...)``, with the keywords of
    :class:`Problem` and of :meth:`Problem.gv` together, which describe them."""
    keywords = {
        k: v
        for k, v in locals().items()
        if k not in ("q", "intnums", "args", "removed")
    }
    return _compute_checked("gv", q, intnums, args, keywords, removed)


def compute_gw(
    q: ArrayLike,
    intnums: Mapping[Any, int] | Iterable[Sequence[int]],
    *args: Any,
    nefpart: Iterable[ArrayLike] | None = None,
    dimension: int | None = None,
    fan: Iterable[ArrayLike] | None = None,
    strata: Iterable[tuple[ArrayLike, Iterable[Any]]] | None = None,
    generators: ArrayLike | None = None,
    mori_hyperplanes: ArrayLike | None = None,
    mori_equations: ArrayLike | None = None,
    hilbert_bases: Iterable[tuple[ArrayLike, ArrayLike]] | None = None,
    grading_vector: ArrayLike | None = None,
    max_deg: int | None = None,
    target_points: ArrayLike | None = None,
    min_points: int | None = None,
    curve_classes: ArrayLike | None = None,
    method: str | None = None,
    n_threads: int | None = None,
    certify: bool = False,
    verify: str | bool = "quick",
    low_memory: bool = False,
    verbose: bool | int = False,
    **removed: Any,
) -> Invariants:
    """The GW invariants: ``Problem(q, intnums, ...).gw(...)``, with the keywords of
    :class:`Problem` and of :meth:`Problem.gw` together, which describe them."""
    keywords = {
        k: v
        for k, v in locals().items()
        if k not in ("q", "intnums", "args", "removed")
    }
    return _compute_checked("gw", q, intnums, args, keywords, removed)
