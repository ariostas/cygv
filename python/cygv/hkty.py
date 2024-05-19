from cygv.cygv import _compute_gvgw

from fractions import Fraction

import mpmath as mp
import numpy as np


def _is_threefold(q, nefpart) -> bool:
    ambient_dim = q.shape[1] - q.shape[0]
    cy_codim = 1 if nefpart is None or len(nefpart) == 0 else len(nefpart)
    return (ambient_dim - cy_codim) == 3


def compute_gv(
    generators,
    grading_vector,
    q,
    intnums,
    max_deg=None,
    min_points=None,
    nefpart=None,
    prec=None,
):
    generators = np.array(generators, dtype=int)
    grading_vector = np.array(grading_vector, dtype=int)
    q = np.array(q, dtype=int)
    intnums = np.array(intnums, dtype=int)
    is_threefold = _is_threefold(q, nefpart)
    res_tmp = _compute_gvgw(
        generators,
        grading_vector,
        q,
        intnums,
        True,
        is_threefold,
        max_deg,
        min_points,
        nefpart,
        prec,
    )
    if is_threefold:
        res = [(tuple(v), int(gv)) for ((v, _), gv) in res_tmp]
    else:
        res = [((tuple(v), c), int(gv)) for ((v, c), gv) in res_tmp]
    return res


def compute_gw(
    generators,
    grading_vector,
    q,
    intnums,
    max_deg=None,
    min_points=None,
    nefpart=None,
    prec=None,
):
    if prec is not None:
        mp.prec = prec
    generators = np.array(generators, dtype=int)
    grading_vector = np.array(grading_vector, dtype=int)
    q = np.array(q, dtype=int)
    intnums = np.array(intnums, dtype=int)
    is_threefold = _is_threefold(q, nefpart)
    res_tmp = _compute_gvgw(
        generators,
        grading_vector,
        q,
        intnums,
        False,
        is_threefold,
        max_deg,
        min_points,
        nefpart,
        prec,
    )
    if is_threefold:
        res = [
            (tuple(v), (Fraction(gw) if prec is None else mp.mpf(gw)))
            for ((v, _), gw) in res_tmp
        ]
    else:
        res = [
            ((tuple(v), c), (Fraction(gw) if prec is None else mp.mpf(gw)))
            for ((v, c), gw) in res_tmp
        ]
    return res
