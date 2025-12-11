from __future__ import annotations

from collections.abc import Sized
from typing import Any

from numpy.typing import ArrayLike

from cygv.cygv import _compute_gvgw


def _compute_gvgw_queue(
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
    result = None
    try:
        result = _compute_gvgw(
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
    except Exception as e:
        result = e
    return result
