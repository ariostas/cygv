from __future__ import annotations

from collections.abc import Sized
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from typing import Any

import mpmath as mp
import numpy as np
from numpy.typing import ArrayLike

from cygv.hkty_worker import _compute_gvgw_queue


# We wrap the raw `_compute_gvgw` function so that we can use ctrl+c
# to interrupt the computation without fully exiting the main python thread.
def _wrapped_compute_gvgw(
    generators: ArrayLike,
    grading_vector: ArrayLike,
    q: ArrayLike,
    intnums: dict[tuple[int, int, int], int],
    find_gv: bool,
    is_threefold: bool,
    max_deg: int | None = None,
    min_points: int | None = None,
    nefpart: Sized | None = None,
    prec: int | None = None,
) -> Any:
    with ProcessPoolExecutor() as pool:
        future = pool.submit(
            _compute_gvgw_queue,
            generators,
            grading_vector,
            q,
            intnums,
            find_gv,
            is_threefold,
            max_deg,
            min_points,
            nefpart,
            prec,
        )
    result = future.result()
    if isinstance(result, Exception):
        raise (result)
    return result


def _is_threefold(q: ArrayLike, nefpart: Sized | None) -> bool:
    ambient_dim = len(q[0]) - len(q)
    cy_codim = 1 if nefpart is None or len(nefpart) == 0 else len(nefpart)
    return (ambient_dim - cy_codim) == 3


def compute_gv(
    generators: ArrayLike,
    grading_vector: ArrayLike,
    q: ArrayLike,
    intnums: dict[tuple[int, int, int], int],
    max_deg: int | None = None,
    min_points: int | None = None,
    nefpart: Sized | None = None,
    prec: int | None = None,
) -> list[Any]:
    generators = np.array(generators, dtype=int)
    grading_vector = np.array(grading_vector, dtype=int)
    q = np.array(q, dtype=int)
    is_threefold = _is_threefold(q, nefpart)
    res_tmp = _wrapped_compute_gvgw(
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
    generators: ArrayLike,
    grading_vector: ArrayLike,
    q: ArrayLike,
    intnums: dict[tuple[int, int, int], int],
    max_deg: int | None = None,
    min_points: int | None = None,
    nefpart: Sized | None = None,
    prec: int | None = None,
) -> list[Any]:
    if prec is not None:
        mp.prec = prec
    generators = np.array(generators, dtype=int)
    grading_vector = np.array(grading_vector, dtype=int)
    q = np.array(q, dtype=int)
    is_threefold = _is_threefold(q, nefpart)
    res_tmp = _wrapped_compute_gvgw(
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
