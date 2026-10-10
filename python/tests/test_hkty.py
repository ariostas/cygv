from __future__ import annotations

from fractions import Fraction

import pytest

from cygv import compute_gv, compute_gw


def test_threefold() -> None:
    generators = [[0, -1], [1, 2]]
    grading_vector = [3, -1]
    q = [[1, 1, 1, 0, 1, 2], [0, 0, -1, 1, 1, -1]]
    intnums = {(0, 0, 0): 2, (0, 0, 1): 1, (0, 1, 1): -1, (1, 1, 1): 5}
    # For now, these are just smoke tests. This is not a consistent geometry, on which the
    # modular method finds invariants that are not integers, so the reference method is
    # requested explicitly where the modular one would run.
    args = (generators, grading_vector, q, intnums)
    compute_gv(*args, min_points=100)
    compute_gw(*args, min_points=100)
    compute_gv(*args, target_points=[10, 10], method="reference")
    compute_gv(*args, target_points=[[10, 10], [12, 10]], method="reference")
    compute_gw(*args, target_points=[10, 10], method="reference")
    compute_gw(*args, target_points=[[10, 10], [12, 10]], method="reference")

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
    compute_gv(generators, grading_vector, q, intnums, min_points=100)
    compute_gw(generators, grading_vector, q, intnums, min_points=100)


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
    compute_gv(generators, grading_vector, q, intnums, min_points=100)
    compute_gw(generators, grading_vector, q, intnums, min_points=100)


def test_methods() -> None:
    # The quintic, which both methods support with max_deg.
    generators = [[1]]
    grading_vector = [1]
    q = [[1, 1, 1, 1, 1]]
    intnums = {(0, 0, 0): 5}
    args = (generators, grading_vector, q, intnums)

    expected = sorted(compute_gv(*args, max_deg=5, method="reference", n_threads=2))
    assert expected[:2] == [((1,), 2875), ((2,), 609250)]
    for method in (None, "auto", "modular", "Modular"):
        assert sorted(compute_gv(*args, max_deg=5, method=method)) == expected
    expected = sorted(compute_gw(*args, max_deg=5, method="reference"))
    assert expected[1] == ((2,), Fraction(609250 * 8 + 2875, 8))
    assert sorted(compute_gw(*args, max_deg=5, method="modular")) == expected

    # Target points, below which both methods compute the same invariants, optionally
    # with a lower maximum degree.
    expected = sorted(compute_gv(*args, target_points=[4], method="reference"))
    assert expected == sorted(compute_gv(*args, max_deg=4, method="reference"))
    assert sorted(compute_gv(*args, target_points=[4], method="modular")) == expected
    assert sorted(compute_gv(*args, target_points=[4], max_deg=3)) == expected[:3]

    # The modular method does not support min_points, and auto falls back to the reference.
    expected = sorted(compute_gv(*args, min_points=3, method="reference"))
    assert expected
    assert sorted(compute_gv(*args, min_points=3, method="auto")) == expected
    with pytest.raises(ValueError, match="min_points"):
        compute_gv(*args, min_points=3, method="modular")
    with pytest.raises(ValueError, match="unknown method"):
        compute_gv(*args, max_deg=3, method="fast")


def test_vex_phase() -> None:
    # A vex phase with h11 = 2, a flop of the resolved P(1,1,2,2,2)[8], which needs its
    # stratum: three divisors and a class, whose entries may be fractions.
    generators = [[0, 1], [1, -1]]
    grading_vector = [2, 1]
    q = [[1, 0, -1, 1, 0, 1], [0, 1, 1, 1, 1, 0]]
    intnums = {(0, 1, 1): 4, (1, 1, 1): 2}
    args = (generators, grading_vector, q, intnums)
    strata = [((1, 2, 4), (1, Fraction(-1)))]

    gv = dict(compute_gv(*args, max_deg=4, strata=strata))
    assert gv[(0, 1)] == 640
    assert gv[(1, -1)] == 1
    assert gv[(1, 0)] == 20
    as_strings = [([1, 2, 4], ["1", "-2/2"])]
    assert (
        dict(compute_gv(*args, max_deg=4, strata=as_strings, method="reference")) == gv
    )

    # Both methods support strata, also with cone data, which then needs the cone of the
    # stratum.
    assert dict(compute_gv(*args, max_deg=4, strata=strata, method="modular")) == gv
    sets: list[tuple[int, ...]] = [(s,) for s in range(6)]
    sets += [(s, t) for s in range(6) for t in range(s + 1, 6)]
    cones = [(T, [[0, 1], [1, 0]] if 2 in T else [[0, 1], [1, 1]]) for T in sets]
    cones.append(((1, 2, 4), [[0, 1], [1, -1]]))
    hyperplanes = [[1, 0], [1, 1]]
    cone_args = (None, grading_vector, q, intnums)
    res = compute_gv(
        *cone_args,
        max_deg=4,
        strata=strata,
        mori_hyperplanes=hyperplanes,
        cones=cones,
    )
    assert dict(res) == gv
    with pytest.raises(ValueError, match="stratum"):
        compute_gv(
            *cone_args,
            max_deg=4,
            strata=strata,
            mori_hyperplanes=hyperplanes,
            cones=cones[:-1],
        )

    # Without the stratum the invariants are not integers.
    with pytest.raises(ValueError, match="not integers"):
        compute_gv(*args, max_deg=6)


def test_cone_data() -> None:
    # P(1,1,2,2,2)[8], whose Mori cone is spanned by (1, 0) and (0, 1). Only the last
    # divisor can be negative, on the classes (x, y) with x < 2y.
    grading_vector = [1, 1]
    q = [[0, 0, 1, 1, 1, 1], [1, 1, 0, 0, 0, -2]]
    intnums = {(0, 0, 0): 8, (0, 0, 1): 4}
    hyperplanes = [[1, 0], [0, 1]]
    sets = [(s,) for s in range(6)] + [
        (s, t) for s in range(6) for t in range(s + 1, 6)
    ]
    cones = [(T, [[1, 0], [0, 1]] if 5 in T else [[1, 0], [2, 1]]) for T in sets]

    expected = sorted(
        compute_gv([[1, 0], [0, 1]], grading_vector, q, intnums, max_deg=6)
    )
    assert ((0, 1), 4) in expected
    for generators in (None, [[1, 0], [0, 1]]):
        res = compute_gv(
            generators,
            grading_vector,
            q,
            intnums,
            max_deg=6,
            mori_hyperplanes=hyperplanes,
            cones=cones,
        )
        assert sorted(res) == expected
    res = compute_gw(
        None,
        grading_vector,
        q,
        intnums,
        max_deg=6,
        mori_hyperplanes=hyperplanes,
        cones=dict(cones).items(),
    )
    assert ((0, 2), Fraction(1, 2)) in res

    # On the face spanned by (1, 0), which meets no divisor negatively, so it is every K_T.
    face_cones = [(T, [[1, 0]]) for T in sets]
    res = compute_gv(
        None,
        grading_vector,
        q,
        intnums,
        max_deg=6,
        mori_hyperplanes=[[1, 0]],
        mori_equations=[[0, 1]],
        cones=face_cones,
    )
    assert sorted(res) == [(c, gv) for (c, gv) in expected if c[1] == 0]

    # Below the targets (2, 2) and (4, 0), i.e. the classes (x, y) with x, y <= 2, or
    # with x <= 4 and y = 0.
    res = compute_gv(
        None,
        grading_vector,
        q,
        intnums,
        target_points=[[2, 2], [4, 0]],
        mori_hyperplanes=hyperplanes,
        cones=cones,
    )
    assert sorted(res) == [
        (c, gv)
        for (c, gv) in expected
        if (c[0] <= 2 and c[1] <= 2) or c in ((3, 0), (4, 0))
    ]

    # Both parts are needed, the generators are needed by the reference method, and the
    # cone data is checked.
    with pytest.raises(ValueError, match="mori_hyperplanes and cones"):
        compute_gv(None, grading_vector, q, intnums, max_deg=6, cones=cones)
    with pytest.raises(ValueError, match="generators"):
        compute_gv(
            None,
            grading_vector,
            q,
            intnums,
            max_deg=6,
            mori_hyperplanes=hyperplanes,
            cones=cones,
            method="reference",
        )
    with pytest.raises(ValueError, match="every set"):
        compute_gv(
            None,
            grading_vector,
            q,
            intnums,
            max_deg=6,
            mori_hyperplanes=hyperplanes,
            cones=cones[1:],
        )
