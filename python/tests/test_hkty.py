from __future__ import annotations

import multiprocessing
import os
import pickle
import sys
from fractions import Fraction
from typing import Any

import pytest

from cygv import (
    METHODS,
    ComputationError,
    Invariants,
    Problem,
    compute_gv,
    compute_gw,
    hkty,
)


def test_threefold() -> None:
    generators = [[0, -1], [1, 2]]
    grading_vector = [3, -1]
    q = [[1, 1, 1, 0, 1, 2], [0, 0, -1, 1, 1, -1]]
    intnums = {(0, 0, 0): 2, (0, 0, 1): 1, (0, 1, 1): -1, (1, 1, 1): 5}
    # For now, these are just smoke tests. This is not a consistent geometry, on which the
    # modular method finds invariants that are not integers, so the reference method is
    # requested explicitly where the modular one would run.
    problem = Problem(q, intnums, generators=generators, grading_vector=grading_vector)
    gv = problem.gv(min_points=100)
    assert isinstance(gv, Invariants)
    assert gv.kind == "gv"
    assert gv.dimension == 3
    assert gv.grading_vector == (3, -1)
    assert gv.report.method == "reference"
    assert gv.report.certified
    assert not gv.report.grading_vector_derived
    # Sorted by degree, then by curve class.
    keys = [(gv.degree(c), c) for c, _ in gv]
    assert keys == sorted(keys)
    problem.gw(min_points=100)
    for target_points in ([10, 10], [[10, 10], [12, 10]]):
        problem.gv(target_points=target_points, method="reference")
        problem.gw(target_points=target_points, method="reference")

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
    assert not QUINTIC.gv(max_deg=5).report.certified
    certified = QUINTIC.gv(max_deg=5, certify=True)
    assert certified == expected
    assert certified.report.certified
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

    # The modular method does not support min_points or curve_classes, which auto runs
    # with the reference.
    expected = QUINTIC.gv(min_points=3, method="reference")
    assert expected
    assert QUINTIC.gv(min_points=3, method="auto") == expected
    # Three classes, zero included.
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
    # is rejected.
    with pytest.raises(TypeError):
        compute_gv([[1]], [1], q, intnums, max_deg=4)  # type: ignore[call-arg]


def test_verbose(capsys: pytest.CaptureFixture[str]) -> None:
    # The progress goes to stderr, and the results are those of a quiet run.
    expected = QUINTIC.gv(max_deg=4)
    assert not capsys.readouterr().err
    assert QUINTIC.gv(max_deg=4, verbose=True) == expected
    assert "pass over primes" in capsys.readouterr().err
    with pytest.raises(ValueError, match="verbose"):
        QUINTIC.gv(max_deg=4, verbose=-1)


def test_invalid_input() -> None:
    # Invalid input is a ValueError with either method, not a crash of the computation.
    bad_intnums = {(0, 0, 0): 5, (0, 0, 1): 1}
    for method in ("reference", "modular"):
        with pytest.raises(ValueError, match="every row"):
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
    with pytest.raises(ValueError, match="Mori cone is needed"):
        compute_gv([[1, 1, 1, 1, 1]], {(0, 0, 0): 5}, max_deg=3)
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
    with pytest.raises(ValueError, match="32 bits"):
        quintic(q=[[2**40, 1, 1, 1, 1]])
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


def test_child_dies(monkeypatch: pytest.MonkeyPatch) -> None:
    # Only a forked child sees the patch.
    if multiprocessing.get_start_method() != "fork":
        pytest.skip("needs the fork start method")

    def die(_kind: str, _args: Any) -> None:
        sys.stderr.write("the end\n")
        sys.stderr.flush()
        os._exit(3)

    monkeypatch.setattr(hkty, "_compute", die)
    with pytest.raises(RuntimeError, match="exit code 3"):
        compute_gv([[1, 1, 1, 1, 1]], {(0, 0, 0): 5}, generators=[[1]], max_deg=1)


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
    assert dict(by_inequalities.gv(max_deg=4, hilbert_bases=bases)) == gv
    with pytest.raises(ValueError, match="stratum"):
        by_inequalities.gv(max_deg=4, hilbert_bases=bases[:-1])

    # Without the stratum the invariants are not integers.
    with pytest.raises(ComputationError, match="not integers"):
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
    assert problem.gv(max_deg=6, hilbert_bases=bases) == expected
    assert problem.gv(max_deg=6, method="reference") == expected
    res = problem.gw(max_deg=6, hilbert_bases=dict(bases).items())
    assert ((0, 2), Fraction(1, 2)) in res

    # On the face spanned by (1, 0), which meets no divisor negatively, so it is every K_T.
    face = Problem(
        q,
        intnums,
        mori_hyperplanes=[[1, 0]],
        mori_equations=[[0, 1]],
        grading_vector=[1, 1],
    )
    face_bases = [(T, [[1, 0]]) for T in sets]
    assert face.gv(max_deg=6, hilbert_bases=face_bases) == [
        (c, gv) for (c, gv) in expected if c[1] == 0
    ]

    # Below the targets (2, 2) and (4, 0), i.e. the classes (x, y) with x, y <= 2, or
    # with x <= 4 and y = 0.
    res = problem.gv(target_points=[[2, 2], [4, 0]], hilbert_bases=bases)
    assert res == [
        (c, gv)
        for (c, gv) in expected
        if (c[0] <= 2 and c[1] <= 2) or c in ((3, 0), (4, 0))
    ]

    # The Hilbert bases are checked.
    with pytest.raises(ValueError, match="every set"):
        problem.gv(max_deg=6, hilbert_bases=bases[1:])
    with pytest.raises(ValueError, match="needs mori_hyperplanes"):
        Problem(q, intnums, mori_equations=[[0, 1]]).gv(max_deg=6)
