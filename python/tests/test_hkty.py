from __future__ import annotations

import _thread
import inspect
import pickle
import threading
import time
from dataclasses import replace
from fractions import Fraction
from typing import Any

import numpy as np
import pytest

from cygv import (
    METHODS,
    ComputationError,
    Invariants,
    Problem,
    compute_gv,
    compute_gw,
)


def test_threefold() -> None:
    generators = [[0, -1], [1, 2]]
    grading_vector = [3, -1]
    q = [[1, 1, 1, 0, 1, 2], [0, 0, -1, 1, 1, -1]]
    intnums = {(0, 0, 0): 2, (0, 0, 1): 1, (0, 1, 1): -1, (1, 1, 1): -5}
    problem = Problem(q, intnums, generators=generators, grading_vector=grading_vector)
    gv = problem.gv(min_points=100)
    assert isinstance(gv, Invariants)
    assert gv.kind == "gv"
    assert gv.dimension == 3
    assert gv.grading_vector == (3, -1)
    assert gv.report.method == "reference"
    assert gv.report.lift == "exact"
    assert gv.report.consistency == "full"
    assert not gv.report.grading_vector_derived
    # Sorted by degree, then by curve class.
    keys = [(gv.degree(c), c) for c, _ in gv]
    assert keys == sorted(keys)
    problem.gw(min_points=100)
    for target_points in ([10, 10], [[10, 10], [12, 10]]):
        gv = problem.gv(target_points=target_points, method="reference")
        assert problem.gv(target_points=target_points) == gv
        gw = problem.gw(target_points=target_points, method="reference")
        assert problem.gw(target_points=target_points) == gw

    # With 5 for -5 the geometry is inconsistent: the instanton corrections disagree, which
    # the reference method finds, and the modular one gives invariants that are not integers.
    wrong = replace(problem, intnums={**intnums, (1, 1, 1): 5})
    for kind in (wrong.gv, wrong.gw):
        with pytest.raises(ComputationError, match="do not agree"):
            kind(max_deg=10, method="reference")
    with pytest.raises(ComputationError, match="not an integer"):
        wrong.gv(max_deg=10)

    generators = [
        [-1, 0, 0, 0, 0, 1, -1],
        [0, 0, 0, 0, 0, -2, 1],
        [0, 0, 0, 0, 0, -1, 1],
        [0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 1, -1, 0, 0],
        [0, 0, 1, -1, 1, 0, 0],
        [0, 0, 1, 0, 0, 1, 0],
        [1, 0, 0, 0, 1, 0, -1],
        [1, 1, 0, 0, 0, -2, 0],
    ]
    grading_vector = [-2, 3, 2, 5, 4, 0, 1]
    q = [
        [-1, 1, 0, 0, -1, 0, 0, 0, -1, 0, 0],
        [1, 0, 1, 0, 3, 0, -2, 0, -1, 0, 0],
        [1, 0, 0, 1, -1, 0, 1, 0, 0, 0, 0],
        [2, 0, 0, 0, 2, 1, -2, 0, 1, 0, 0],
        [2, 0, 0, 0, 3, 0, -3, 1, 1, 0, 0],
        [0, 0, 0, 0, 1, 0, -1, 0, -1, 1, 0],
        [1, 0, 0, 0, 2, 0, -2, 0, 0, 0, 1],
    ]
    intnums = {
        (0, 5, 6): 1,
        (1, 2, 3): 2,
        (1, 3, 4): 2,
        (1, 4, 6): 1,
        (1, 5, 6): 1,
        (0, 0, 2): -2,
        (0, 0, 6): -1,
        (0, 5, 5): -2,
        (0, 6, 6): -1,
        (1, 2, 2): 4,
        (1, 3, 3): -2,
        (1, 4, 4): -2,
        (1, 5, 5): -2,
        (1, 6, 6): -2,
        (2, 2, 2): 24,
        (2, 2, 3): 8,
        (2, 3, 3): 2,
        (3, 3, 3): -12,
        (3, 3, 4): 10,
        (3, 4, 4): -8,
        (4, 4, 4): 6,
        (4, 4, 6): 1,
        (4, 6, 6): -3,
        (5, 5, 5): 8,
        (5, 5, 6): -2,
        (6, 6, 6): 7,
    }
    # For now, these are just smoke tests
    problem = Problem(q, intnums, generators=generators, grading_vector=grading_vector)
    problem.gv(min_points=100)
    problem.gw(min_points=100)


def test_fourfold() -> None:
    generators = [
        [0, 0, 0, -2, 3, 0],
        [1, 0, -4, 0, 0, 0],
        [0, 0, 1, 1, 0, 0],
        [0, 1, 0, 0, 1, -2],
        [0, -2, 0, 1, 0, 2],
        [0, 0, 0, 0, -2, 1],
    ]
    grading_vector = [73, -30, 18, -17, -11, -21]
    q = [
        [1, 1, 1, 1, 0, 0, 8, 12, 4, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, -4, -6, -3, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, -2, -3, -2, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, -2, -2, -1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, -3, -4, -2, 0, 0, 1],
    ]
    intnums = {
        (0, 0, 2): 1,
        (0, 1, 1): -8,
        (0, 1, 3): 4,
        (0, 1, 5): 8,
        (0, 2, 2): -4,
        (0, 3, 3): -8,
        (0, 4, 4): -16,
        (0, 4, 5): 8,
        (0, 5, 5): -16,
        (1, 0, 1): -8,
        (1, 0, 3): 4,
        (1, 0, 5): 8,
        (1, 1, 3): 16,
        (1, 1, 5): -32,
        (1, 3, 3): -16,
        (1, 5, 5): 64,
        (2, 0, 0): 1,
        (2, 0, 2): -4,
        (2, 2, 2): 16,
        (3, 0, 1): 4,
        (3, 0, 3): -8,
        (3, 1, 1): 16,
        (3, 1, 3): -16,
        (4, 0, 4): -16,
        (4, 0, 5): 8,
        (4, 4, 4): -128,
        (4, 4, 5): 32,
        (5, 0, 1): 8,
        (5, 0, 4): 8,
        (5, 0, 5): -16,
        (5, 1, 1): -32,
        (5, 1, 5): 64,
        (5, 4, 4): 32,
        (5, 5, 5): -128,
        (6, 0, 0): -8,
        (6, 0, 3): 16,
        (6, 0, 5): -32,
        (6, 1, 1): -128,
        (6, 1, 3): 64,
        (6, 1, 5): 128,
        (6, 3, 3): -64,
        (6, 5, 5): -256,
        (7, 0, 0): 4,
        (7, 0, 1): 16,
        (7, 0, 3): -16,
        (7, 1, 1): 64,
        (7, 1, 3): -64,
        (7, 3, 3): 64,
        (8, 0, 0): 8,
        (8, 0, 1): -32,
        (8, 0, 5): 64,
        (8, 1, 1): 128,
        (8, 1, 5): -256,
        (8, 5, 5): 512,
        (9, 0, 0): -4,
        (9, 0, 2): 16,
        (9, 2, 2): -64,
        (10, 0, 0): -8,
        (10, 0, 1): -16,
        (10, 1, 1): -64,
        (10, 1, 3): 64,
        (10, 3, 3): -128,
        (11, 0, 0): -16,
        (11, 0, 4): -128,
        (11, 0, 5): 32,
        (11, 4, 4): -768,
        (11, 4, 5): 128,
        (12, 0, 0): 8,
        (12, 0, 4): 32,
        (12, 4, 4): 128,
        (13, 0, 0): -16,
        (13, 0, 1): 64,
        (13, 0, 5): -128,
        (13, 1, 1): -256,
        (13, 1, 5): 512,
        (13, 5, 5): -1024,
    }
    # For now, these are just smoke tests
    problem = Problem(
        q,
        intnums,
        nefpart=[[0, 1, 2, 3], [4, 5, 6, 7, 8, 9, 10, 11]],
        generators=generators,
        grading_vector=grading_vector,
    )
    gv = problem.gv(min_points=100)
    assert gv.dimension == 4
    # The keys of an n-fold carry the reference surface.
    (curve, surface), value = gv[0]
    assert isinstance(curve, tuple)
    assert isinstance(surface, int)
    assert isinstance(value, int)
    problem.gw(min_points=100)


QUINTIC = Problem([[1, 1, 1, 1, 1]], {(0, 0, 0): 5}, generators=[[1]])


def test_methods() -> None:
    # The quintic, which both methods support with max_deg.
    expected = QUINTIC.gv(max_deg=5, method="reference", n_threads=2)
    assert expected[:2] == [((1,), 2875), ((2,), 609250)]
    assert expected.grading_vector == (1,)
    assert expected.report.grading_vector_derived
    for method in (None, *METHODS, "Modular"):
        assert QUINTIC.gv(max_deg=5, method=method) == expected
    assert QUINTIC.gv(max_deg=5).report.method == "modular"
    assert QUINTIC.gv(max_deg=5).report.lift == "stable"
    # The quintic has h11 = 1, where the consistency check does not apply.
    assert QUINTIC.gv(max_deg=5).report.consistency == "unchecked"
    certified = QUINTIC.gv(max_deg=5, certify=True)
    assert certified == expected
    assert certified.report.lift == "certified"
    assert QUINTIC.gv(max_deg=5, low_memory=True) == expected
    expected = QUINTIC.gw(max_deg=5, method="reference")
    assert expected[1] == ((2,), Fraction(609250 * 8 + 2875, 8))
    assert QUINTIC.gw(max_deg=5, method="modular") == expected
    assert dict(expected)[(2,)] == Fraction(609250 * 8 + 2875, 8)

    # The results survive pickling, report included.
    copy = pickle.loads(pickle.dumps(expected))
    assert copy == expected
    assert copy.report == expected.report

    # Target points, below which both methods compute the same invariants, optionally
    # with a lower maximum degree.
    expected = QUINTIC.gv(target_points=[4], method="reference")
    assert expected == QUINTIC.gv(max_deg=4, method="reference")
    assert QUINTIC.gv(target_points=[4], method="modular") == expected
    assert QUINTIC.gv(target_points=[4], max_deg=3) == expected[:3]

    # The modular method does not support min_points, which auto runs with the reference.
    expected = QUINTIC.gv(min_points=3, method="reference")
    assert expected
    assert QUINTIC.gv(min_points=3, method="auto") == expected
    assert QUINTIC.gv(curve_classes=[[1], [2]]) == expected
    with pytest.raises(ValueError, match="min_points"):
        QUINTIC.gv(min_points=3, method="modular")
    with pytest.raises(ValueError, match="unknown method"):
        QUINTIC.gv(max_deg=3, method="fast")

    # Exactly one selection.
    with pytest.raises(ValueError, match="which curve classes"):
        QUINTIC.gv()
    with pytest.raises(ValueError, match="only one of"):
        QUINTIC.gv(max_deg=3, min_points=3)


def test_curve_classes() -> None:
    # The invariants of a class depend on the classes below it, which need not be listed.
    p11222 = Problem(
        [[0, 0, 1, 1, 1, 1], [1, 1, 0, 0, 0, -2]],
        {(0, 0, 0): 8, (0, 0, 1): 4},
        generators=[[1, 0], [0, 1]],
    )
    expected = dict(p11222.gv(max_deg=3))
    assert expected[(1, 1)] == 640
    for method in ["auto", "reference"]:
        res = p11222.gv(curve_classes=[[2, 1], [1, 1]], method=method)
        assert res.report.method == ("modular" if method == "auto" else "reference")
        assert list(res) == [((1, 1), 640), ((2, 1), expected[(2, 1)])]
    with pytest.raises(ValueError, match="outside the Mori cone"):
        p11222.gv(curve_classes=[[2, -1]])


def test_compute_functions() -> None:
    # compute_gv and compute_gw take the keywords of the problem and of the computation.
    expected = QUINTIC.gv(max_deg=4)
    q, intnums = [[1, 1, 1, 1, 1]], {(0, 0, 0): 5}
    assert compute_gv(q, intnums, generators=[[1]], max_deg=4) == expected
    assert compute_gv(q=q, intnums=intnums, generators=[[1]], max_deg=4) == expected
    assert compute_gw(q, intnums, generators=[[1]], max_deg=4) == QUINTIC.gw(max_deg=4)
    # The intersection numbers may also be [i, j, k, value] entries.
    assert compute_gv(q, [[0, 0, 0, 5]], generators=[[1]], max_deg=4) == expected
    with pytest.raises(TypeError):
        compute_gv(q, intnums, generators=[[1]], max_deg=4, bad_keyword=1)
    # The positional order of earlier versions, (generators, grading_vector, q, intnums),
    # is rejected, with what to write instead, and so are their removed keywords.
    with pytest.raises(TypeError, match=r"cygv 0\.3 took"):
        compute_gv([[1]], [1], q, intnums, max_deg=4)  # type: ignore[list-item]
    with pytest.raises(TypeError, match="always exact"):
        compute_gv(q, intnums, generators=[[1]], max_deg=4, prec=100)
    with pytest.raises(TypeError, match="hilbert_bases"):
        compute_gw(q, intnums, generators=[[1]], max_deg=4, cones=[])
    # Without a selection, 0.3 computed at the generators.
    with pytest.raises(ValueError, match="curve_classes=generators"):
        compute_gv(q, intnums, generators=[[1]], max_deg=None, min_points=None)
    # The keywords are explicit, for help() and completion.
    params = inspect.signature(compute_gv).parameters
    assert {"generators", "max_deg", "verify", "certify", "fan"} <= set(params)
    # Problem's fields after intnums are keywords only: by position, generators would have
    # been taken for the nef partition.
    with pytest.raises(TypeError):
        Problem(q, intnums, [[1]])  # type: ignore[call-arg]


def test_verbose(capsys: pytest.CaptureFixture[str]) -> None:
    # The progress goes to stderr, and the results are those of a quiet run.
    expected = QUINTIC.gv(max_deg=4)
    assert not capsys.readouterr().err
    assert QUINTIC.gv(max_deg=4, verbose=True) == expected
    assert "pass over primes" in capsys.readouterr().err
    with pytest.raises(ValueError, match="verbose"):
        QUINTIC.gv(max_deg=4, verbose=-1)
    assert QUINTIC.gv(max_deg=4, verbose=300) == expected
    capsys.readouterr()


def test_invalid_input() -> None:
    # Invalid input is a ValueError with either method, not a crash of the computation.
    bad_intnums = {(0, 0, 0): 5, (0, 0, 1): 1}
    for method in ("reference", "modular"):
        with pytest.raises(ValueError, match="every toric divisor"):
            compute_gv(
                [[1, 1, 1, 1, 1]],
                {(0, 0, 0): 5},
                generators=[[1]],
                nefpart=[[0, 1, 2, 3]],
                max_deg=3,
                method=method,
            )
        with pytest.raises(ValueError, match="below h11"):
            compute_gv(
                [[1, 1, 1, 1, 1]],
                bad_intnums,
                generators=[[1]],
                max_deg=3,
                method=method,
            )
    # A transposed q says so.
    with pytest.raises(ValueError, match="transposed"):
        compute_gv(
            [[1], [1], [1], [1], [1]], {(0, 0, 0): 5}, generators=[[1]], max_deg=3
        )
    # So does the origin of the polytope among the divisors, as CYTools includes by default.
    with pytest.raises(ValueError, match="include_origin=False"):
        compute_gv([[-5, 1, 1, 1, 1, 1]], {(0, 0, 0): 5}, generators=[[1]], max_deg=3)
    with pytest.raises(ValueError, match="Mori cone is needed"):
        compute_gv([[1, 1, 1, 1, 1]], {(0, 0, 0): 5}, max_deg=3)
    # Options of the wrong type, and negative indices, are ValueErrors too.
    for kwargs in (
        {"certify": "yes"},
        {"low_memory": 2},
        {"method": 1},
        {"nefpart": [[0, -1]]},
    ):
        with pytest.raises(ValueError, match=next(iter(kwargs))):
            compute_gv(
                [[1, 1, 1, 1, 1]], {(0, 0, 0): 5}, generators=[[1]], max_deg=3, **kwargs
            )
    assert QUINTIC.gv(max_deg=2, certify=1) == QUINTIC.gv(max_deg=2)  # type: ignore[arg-type]
    # Problems compare by identity, since their fields may be arrays.
    q = np.array([[1, 1, 1, 1, 1]])
    assert Problem(q, {(0, 0, 0): 5}) != Problem(q, {(0, 0, 0): 5})
    with pytest.raises(ValueError, match="not both"):
        compute_gv(
            [[1, 1, 1, 1, 1]],
            {(0, 0, 0): 5},
            generators=[[1]],
            mori_hyperplanes=[[1]],
            max_deg=3,
        )
    # No equations are none.
    assert compute_gv(
        [[1, 1, 1, 1, 1]],
        {(0, 0, 0): 5},
        generators=[[1]],
        mori_equations=[],
        max_deg=1,
    ) == [((1,), 2875)]


def test_integers() -> None:
    def quintic(q: Any = None, intnums: Any = None, **kwargs: Any) -> Any:
        args = {"generators": [[1]], "max_deg": 2, **kwargs}
        return compute_gv(
            [[1, 1, 1, 1, 1]] if q is None else q,
            {(0, 0, 0): 5} if intnums is None else intnums,
            **args,
        )

    expected = quintic()
    # Integral floats are integers, and anything else is an error, not truncated.
    assert quintic(q=[[1.0, 1, 1, 1, 1]], max_deg=2.0) == expected
    for q in ([[1.5, 1, 1, 1, 1]], [["1", 1, 1, 1, 1]]):
        with pytest.raises(ValueError, match="q must hold integers"):
            quintic(q=q)
    with pytest.raises(ValueError, match="intnums must hold integers"):
        quintic(intnums={(0, 0, 0): 5.5})
    with pytest.raises(ValueError, match="rectangular"):
        quintic(q=[[1, 1, 1, 1, 1], [1, 1]])
    for big in (
        [[2**40, 1, 1, 1, 1]],
        np.array([[2**64 - 1, 1, 1, 1, 1]], dtype=np.uint64),
        np.array([[1e30, 1, 1, 1, 1]]),
    ):
        with pytest.raises(ValueError, match="32 bits"):
            quintic(q=big)
    with pytest.raises(ValueError, match="32 bits"):
        quintic(intnums={(0, 0, 0): 2**40})
    for max_deg in (-1, 2.5, 2**40):
        with pytest.raises(ValueError, match="max_deg"):
            quintic(max_deg=max_deg)
    with pytest.raises(ValueError, match="n_threads"):
        quintic(n_threads=-2)
    # Target points outside the Mori cone.
    with pytest.raises(ValueError, match="outside the Mori cone"):
        compute_gv(
            [[0, 0, 1, 1, 1, 1], [1, 1, 0, 0, 0, -2]],
            {(0, 0, 0): 8, (0, 0, 1): 4},
            generators=[[1, 0], [0, 1]],
            target_points=[[2, -1]],
        )


def test_arrays() -> None:
    # The invariants are numpy arrays too, of int64 while they fit and of ints after.
    gv = QUINTIC.gv(max_deg=12)
    assert gv.curves.dtype == np.int32
    assert gv.curves.shape == (12, 1)
    assert gv.degrees.tolist() == list(range(1, 13))
    assert gv.surfaces is None
    assert gv.values.dtype == object
    assert gv.values.tolist() == [v for _, v in gv]
    assert [k for k, _ in gv] == [(d,) for d in range(1, 13)]
    small = QUINTIC.gv(max_deg=5)
    assert small.values.dtype == np.int64
    assert small.values.tolist() == [v for _, v in gv[:5]]
    # They cannot be changed behind the pairs' back.
    with pytest.raises(ValueError, match="read-only"):
        small.curves[0, 0] = 2
    gw = QUINTIC.gw(max_deg=3)
    assert gw.values.tolist() == [v for _, v in gw]
    assert gw.values is gw.values
    assert gw[-1] == list(gw)[-1]

    # The values by key.
    assert gv[(2,)] == gv.get([2]) == dict(gv)[(2,)]
    assert gv.get((13,)) is None
    with pytest.raises(KeyError):
        gv[(13,)]
    with pytest.raises(ValueError, match="1 entries"):
        gv.degree((1, 0))


def test_interrupt() -> None:
    # Ctrl+c stops the computation, which runs with the GIL released: the timer that sends
    # it is a Python thread.
    timer = threading.Timer(0.5, _thread.interrupt_main)
    timer.start()
    start = time.monotonic()
    with pytest.raises(KeyboardInterrupt):
        # About 20 s and a few MB, since the invariants need many primes.
        QUINTIC.gv(max_deg=200, n_threads=2)
    timer.join()
    assert time.monotonic() - start < 5
    # The next computation is not affected.
    assert QUINTIC.gv(max_deg=2)[0] == ((1,), 2875)


def test_vex_phase() -> None:
    # A vex phase with h11 = 2, a flop of the resolved P(1,1,2,2,2)[8], which needs its
    # stratum: three divisors and a class, whose entries may be fractions.
    q = [[1, 0, -1, 1, 0, 1], [0, 1, 1, 1, 1, 0]]
    intnums = {(0, 1, 1): 4, (1, 1, 1): 2}
    generators = [[0, 1], [1, -1]]
    strata = [((1, 2, 4), (1, Fraction(-1)))]

    gv = dict(Problem(q, intnums, strata=strata, generators=generators).gv(max_deg=4))
    assert gv[(0, 1)] == 640
    assert gv[(1, -1)] == 1
    assert gv[(1, 0)] == 20
    as_strings = [([1, 2, 4], ["1", "-2/2"])]
    problem = Problem(q, intnums, strata=as_strings, generators=generators)
    assert dict(problem.gv(max_deg=4, method="reference")) == gv

    # Both methods support strata, also with Hilbert bases, which then need the cone of the
    # stratum.
    assert dict(problem.gv(max_deg=4, method="modular")) == gv
    sets: list[tuple[int, ...]] = [(s,) for s in range(6)]
    sets += [(s, t) for s in range(6) for t in range(s + 1, 6)]
    bases = [(T, [[0, 1], [1, 0]] if 2 in T else [[0, 1], [1, 1]]) for T in sets]
    bases.append(((1, 2, 4), [[0, 1], [1, -1]]))
    by_inequalities = Problem(
        q,
        intnums,
        strata=strata,
        mori_hyperplanes=[[1, 0], [1, 1]],
        grading_vector=[2, 1],
    )
    assert dict(by_inequalities.gv(max_deg=4)) == gv
    assert dict(replace(by_inequalities, hilbert_bases=bases).gv(max_deg=4)) == gv
    with pytest.raises(ValueError, match="stratum"):
        replace(by_inequalities, hilbert_bases=bases[:-1]).gv(max_deg=4)

    # Without the stratum, the curve of negative anticanonical degree gives it away.
    with pytest.raises(ValueError, match="negative anticanonical degree"):
        Problem(q, intnums, generators=generators).gv(max_deg=6)

    # The fan gives the stratum, the Mori cone and the grading vector.
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
    from_fan = Problem(q, intnums, fan=fan).gv(max_deg=4)
    assert dict(from_fan) == gv
    assert from_fan.report.grading_vector_derived
    with pytest.raises(ValueError, match="not both"):
        Problem(q, intnums, fan=fan, strata=strata).gv(max_deg=4)
    with pytest.raises(ValueError, match="intersection ring"):
        Problem(q, {(0, 1, 1): 4, (1, 1, 1): 3}, fan=fan).gv(max_deg=4)
    with pytest.raises(ValueError, match="vex phase"):
        Problem(q, intnums, fan=fan).gv(max_deg=4, certify=True)


def test_inequalities() -> None:
    # P(1,1,2,2,2)[8], whose Mori cone is spanned by (1, 0) and (0, 1). Only the last
    # divisor can be negative, on the classes (x, y) with x < 2y.
    q = [[0, 0, 1, 1, 1, 1], [1, 1, 0, 0, 0, -2]]
    intnums = {(0, 0, 0): 8, (0, 0, 1): 4}
    sets = [(s,) for s in range(6)] + [
        (s, t) for s in range(6) for t in range(s + 1, 6)
    ]
    bases = [(T, [[1, 0], [0, 1]] if 5 in T else [[1, 0], [2, 1]]) for T in sets]

    expected = Problem(q, intnums, generators=[[1, 0], [0, 1]]).gv(max_deg=6)
    assert ((0, 1), 4) in expected
    problem = Problem(q, intnums, mori_hyperplanes=[[1, 0], [0, 1]])
    assert problem.gv(max_deg=6) == expected
    with_bases = replace(problem, hilbert_bases=bases)
    assert with_bases.gv(max_deg=6) == expected
    assert problem.gv(max_deg=6, method="reference") == expected
    res = replace(problem, hilbert_bases=dict(bases).items()).gw(max_deg=6)
    assert ((0, 2), Fraction(1, 2)) in res

    # On the face spanned by (1, 0), which meets no divisor negatively, so it is every K_T.
    face = Problem(
        q,
        intnums,
        mori_hyperplanes=[[1, 0]],
        mori_equations=[[0, 1]],
        hilbert_bases=[(T, [[1, 0]]) for T in sets],
        grading_vector=[1, 1],
    )
    assert face.gv(max_deg=6) == [(c, gv) for (c, gv) in expected if c[1] == 0]

    # Below the targets (2, 2) and (4, 0), i.e. the classes (x, y) with x, y <= 2, or
    # with x <= 4 and y = 0.
    res = with_bases.gv(target_points=[[2, 2], [4, 0]])
    assert res == [
        (c, gv)
        for (c, gv) in expected
        if (c[0] <= 2 and c[1] <= 2) or c in ((3, 0), (4, 0))
    ]

    # The Hilbert bases are checked.
    with pytest.raises(ValueError, match="every set"):
        replace(problem, hilbert_bases=bases[1:]).gv(max_deg=6)
    with pytest.raises(ValueError, match="needs mori_hyperplanes"):
        Problem(q, intnums, mori_equations=[[0, 1]]).gv(max_deg=6)


def test_consistency() -> None:
    # A spurious intersection number gives integers at degree 4 without the check, wrong ones,
    # and an error with it, in both methods.
    q = [[0, 0, 1, 1, 1, 1], [1, 1, 0, 0, 0, -2]]
    good = Problem(q, {(0, 0, 0): 8, (0, 0, 1): 4}, generators=[[1, 0], [0, 1]])
    wrong = replace(good, intnums={(0, 0, 0): 8, (0, 0, 1): 4, (0, 1, 1): 1})
    assert wrong.gv(max_deg=4, verify="off")[(1, 0)] == 896
    assert good.gv(max_deg=4)[(1, 0)] == 640
    for method in ("modular", "reference"):
        with pytest.raises(ComputationError, match="do not agree"):
            wrong.gv(max_deg=4, method=method)
    with pytest.raises(ComputationError, match="do not agree"):
        wrong.gv(max_deg=4, verify=True, certify=True)
    cases: list[tuple[str | bool, str]] = [
        ("full", "full"),
        (False, "unchecked"),
        ("quick", "quick"),
    ]
    for verify, consistency in cases:
        assert good.gv(max_deg=3, verify=verify).report.consistency == consistency
    with pytest.raises(ValueError, match="verify"):
        good.gv(max_deg=3, verify="sometimes")
    # The dimension, when given, is checked.
    with pytest.raises(ValueError, match="dimension"):
        replace(good, dimension=4).gv(max_deg=3)
    assert replace(good, dimension=3).gv(max_deg=3) == good.gv(max_deg=3)


def test_lookup() -> None:
    gv = QUINTIC.gv(max_deg=3)
    # `in` and indexing look up the same keys, from any sequence of ints.
    assert (1,) in gv  # type: ignore[comparison-overlap]
    assert [1] in gv  # type: ignore[comparison-overlap]
    assert np.array([2]) in gv
    assert (4,) not in gv  # type: ignore[comparison-overlap]
    # The (key, value) pairs that iterating gives are in it too.
    assert ((1,), 2875) in gv
    assert ((1,), 1) not in gv
    assert gv[[1]] == gv[(1,)] == gv[np.array([1])] == 2875
    assert gv.get([4]) is None
    with pytest.raises(KeyError):
        gv[(4,)]
    with pytest.raises(KeyError):
        gv[("a",)]
