"""cygv-style interface to cgv: same signatures as cygv.compute_gv / cygv.compute_gw.

    from cygv_compat import compute_gv
    gvs = compute_gv(generators, grading_vector, q, intnums, max_deg=20)   # [(curve, GV), ...]

Options cgv does not support raise NotImplementedError.
"""
from cgv_run import compute_gvs


def compute_gv(generators, grading_vector, q, intnums, max_deg=None, min_points=None,
               target_points=None, nefpart=None, prec=None):
    if nefpart is not None:
        raise NotImplementedError("cgv: complete intersections (nefpart) are not supported")
    if min_points is not None:
        raise NotImplementedError("cgv: min_points is not supported; pass max_deg")
    if target_points is not None:
        raise NotImplementedError("cgv: target_points is not supported")
    if max_deg is None:
        raise NotImplementedError("cgv: max_deg is required")
    # prec is accepted and ignored: cgv is exact (modular arithmetic + CRT)
    d = dict(generators=[[int(x) for x in g] for g in generators],
             grading_vector=[int(x) for x in grading_vector],
             q=[[int(x) for x in r] for r in q],
             intnums=[[int(i), int(j), int(k), int(v)] for (i, j, k), v in intnums.items()])
    return list(compute_gvs(d, int(max_deg)).items())


def compute_gw(generators, grading_vector, q, intnums, max_deg=None, min_points=None,
               target_points=None, nefpart=None, prec=None):
    raise NotImplementedError("cgv: Gromov-Witten invariants are not supported")
