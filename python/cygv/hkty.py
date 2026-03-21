from __future__ import annotations

import os
from collections.abc import Sized
from fractions import Fraction
from multiprocessing import Pipe, Process
from multiprocessing.connection import Connection, wait
from typing import Any

import mpmath as mp
import numpy as np
from numpy.typing import ArrayLike

from cygv.cygv import _compute_gvgw


def _compute_gvgw_subprocess(
    conn: Connection,
    stderr_fd: int,
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
) -> None:
    os.dup2(stderr_fd, 2)
    os.close(stderr_fd)
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
                nefpart,
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
    nefpart: Sized | None = None,
    prec: int | None = None,
) -> Any:
    parent_conn, child_conn = Pipe(duplex=False)
    stderr_r_fd, stderr_w_fd = os.pipe()
    process = Process(
        target=_compute_gvgw_subprocess,
        args=(
            child_conn,
            stderr_w_fd,
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
        ),
    )
    process.start()
    child_conn.close()
    os.close(stderr_w_fd)

    ready = wait([parent_conn, process.sentinel])
    if parent_conn in ready:
        result = parent_conn.recv()
        process.join()
        os.close(stderr_r_fd)
        if isinstance(result, BaseException):
            raise result
        return result

    process.join()
    stderr_msg = os.read(stderr_r_fd, 65536).decode(errors="replace").strip()
    os.close(stderr_r_fd)
    msg = f"Computation failed (exit code {process.exitcode})"
    if stderr_msg:
        msg += f":\n{stderr_msg}"
    raise RuntimeError(msg)


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
        mp.mp.prec = prec
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
