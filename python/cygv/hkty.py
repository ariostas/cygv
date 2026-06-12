from __future__ import annotations

import os
import tempfile
from collections.abc import Sized
from fractions import Fraction
from multiprocessing import Pipe, Process
from multiprocessing.connection import Connection, wait
from pathlib import Path
from typing import Any

import mpmath as mp
import numpy as np
from numpy.typing import ArrayLike

from cygv.cygv import _compute_gvgw


def _compute_gvgw_subprocess(
    conn: Connection,
    stderr_path: Path,
    generators: ArrayLike,
    grading_vector: ArrayLike,
    q: ArrayLike,
    intnums: dict[tuple[int, int, int], int],
    find_gv: bool,
    is_threefold: bool,
    max_deg: int | None = None,
    min_points: int | None = None,
    target_points: ArrayLike | None = None,
    nefpart: Sized | None = None,
    prec: int | None = None,
) -> None:
    with stderr_path.open("wb") as f:
        os.dup2(f.fileno(), 2)
    try:
        conn.send(
            _compute_gvgw(
                generators,
                grading_vector,
                q,
                intnums,
                find_gv,
                is_threefold,
                max_deg,
                min_points,
                target_points,
                nefpart,
                None,
                1000,
                prec,
            )
        )
    except BaseException as e:
        conn.send(RuntimeError(str(e)))
    conn.close()


# We wrap the raw `_compute_gvgw` function so that we can use ctrl+c
# to interrupt the computation without fully exiting the main python process.
def _wrapped_compute_gvgw(
    generators: ArrayLike,
    grading_vector: ArrayLike,
    q: ArrayLike,
    intnums: dict[tuple[int, int, int], int],
    find_gv: bool,
    is_threefold: bool,
    max_deg: int | None = None,
    min_points: int | None = None,
    target_points: ArrayLike | None = None,
    nefpart: Sized | None = None,
    prec: int | None = None,
) -> Any:
    parent_conn, child_conn = Pipe(duplex=False)
    with tempfile.NamedTemporaryFile(delete=False) as f:
        stderr_path = Path(f.name)
    process = Process(
        target=_compute_gvgw_subprocess,
        args=(
            child_conn,
            stderr_path,
            generators,
            grading_vector,
            q,
            intnums,
            find_gv,
            is_threefold,
            max_deg,
            min_points,
            target_points,
            nefpart,
            prec,
        ),
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


def _is_threefold(q: ArrayLike, nefpart: Sized | None) -> bool:
    ambient_dim = len(q[0]) - len(q)
    cy_codim = 1 if nefpart is None or len(nefpart) == 0 else len(nefpart)
    return (ambient_dim - cy_codim) == 3


def _regularize_target_points(
    target_points: ArrayLike | None,
) -> np.ndarray[Any, Any] | None:
    if target_points is None:
        return None
    target_points = np.array(target_points, dtype=int)
    if target_points.size == 0:
        return None
    if target_points.ndim > 2:
        raise ValueError("target_points must be at most 2-dimensional")
    if target_points.ndim == 1:
        target_points = target_points.reshape(1, -1)
    return target_points


def compute_gv(
    generators: ArrayLike,
    grading_vector: ArrayLike,
    q: ArrayLike,
    intnums: dict[tuple[int, int, int], int],
    max_deg: int | None = None,
    min_points: int | None = None,
    target_points: ArrayLike | None = None,
    nefpart: Sized | None = None,
    prec: int | None = None,
) -> list[Any]:
    generators = np.array(generators, dtype=int)
    grading_vector = np.array(grading_vector, dtype=int)
    q = np.array(q, dtype=int)
    target_points = _regularize_target_points(target_points)
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
        target_points,
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
    target_points: ArrayLike | None = None,
    nefpart: Sized | None = None,
    prec: int | None = None,
) -> list[Any]:
    if prec is not None:
        mp.mp.prec = prec
    generators = np.array(generators, dtype=int)
    grading_vector = np.array(grading_vector, dtype=int)
    q = np.array(q, dtype=int)
    target_points = _regularize_target_points(target_points)
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
        target_points,
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
