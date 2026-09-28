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
                prec,
            )
        )
    # The catch has to be this broad: a panic on the Rust side surfaces as a
    # PanicException, which derives from BaseException rather than Exception, and
    # it still needs to be reported back to the parent process.
    except BaseException as e:  # noqa: BLE001
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
        msg = "target_points must be a 1D or 2D array-like of ints"
        raise ValueError(msg)
    if target_points.ndim == 1:
        target_points = target_points.reshape(1, -1)
    return target_points


def _cgv_gvs(
    generators: ArrayLike,
    grading_vector: ArrayLike,
    q: ArrayLike,
    intnums: dict[tuple[int, int, int], int],
    max_deg: int | None,
    min_points: int | None,
    target_points: ArrayLike | None,
    nefpart: Sized | None,
) -> dict[tuple[int, ...], int]:
    """GV invariants with the bundled cgv program (cgv/ in the repository), through
    the same code as its standalone Python interface (cgv/tools/cgv_run.py)."""
    if nefpart is not None and len(nefpart) > 0:
        msg = "backend='cgv' supports hypersurfaces only (no nefpart)"
        raise NotImplementedError(msg)
    if not _is_threefold(q, nefpart):
        msg = "backend='cgv' supports threefolds only"
        raise NotImplementedError(msg)
    if min_points is not None or target_points is not None or max_deg is None:
        msg = "backend='cgv' needs max_deg (min_points and target_points are not supported)"
        raise NotImplementedError(msg)
    from cygv import _cgv_run  # noqa: PLC0415
    from cygv.cygv import _cgv_executable  # noqa: PLC0415

    d = {
        "generators": [[int(x) for x in g] for g in np.array(generators, dtype=int)],
        "grading_vector": [int(x) for x in np.array(grading_vector, dtype=int)],
        "q": [[int(x) for x in r] for r in np.array(q, dtype=int)],
        "intnums": [[int(i), int(j), int(k), int(v)] for (i, j, k), v in intnums.items()],
    }
    _cgv_run.BIN = _cgv_executable()
    threads = os.cpu_count() or 1
    try:
        gvs, _, _ = _cgv_run.run_cgv(d, int(max_deg), threads)
    except FileNotFoundError:
        # no normaliz for the cone data: cgv enumerates the cone itself (same result, slower)
        gvs, _, _ = _cgv_run.run_cgv(d, int(max_deg), threads, cones=False)
    return gvs


def _gw_from_gv(
    gvs: dict[tuple[int, ...], int], grading_vector: ArrayLike, max_deg: int
) -> dict[tuple[int, ...], Fraction]:
    """Genus-0 GW invariants of a threefold from its GV invariants:
    GW(b) = sum over k | b of GV(b/k) / k^3."""
    w = np.array(grading_vector, dtype=int)
    gw: dict[tuple[int, ...], Fraction] = {}
    for c, v in gvs.items():
        deg = int(np.array(c, dtype=int) @ w)
        k = 1
        while deg * k <= max_deg:
            b = tuple(k * x for x in c)
            gw[b] = gw.get(b, Fraction(0)) + Fraction(v, k**3)
            k += 1
    return {b: x for b, x in gw.items() if x != 0}


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
    backend: str = "cygv",
) -> list[Any]:
    if backend == "cgv":
        return list(
            _cgv_gvs(generators, grading_vector, q, intnums, max_deg, min_points, target_points, nefpart).items()
        )
    if backend != "cygv":
        msg = f"unknown backend {backend!r} (use 'cygv' or 'cgv')"
        raise ValueError(msg)
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
    backend: str = "cygv",
) -> list[Any]:
    if prec is not None:
        mp.mp.prec = prec
    if backend == "cgv":
        gvs = _cgv_gvs(generators, grading_vector, q, intnums, max_deg, min_points, target_points, nefpart)
        gws = _gw_from_gv(gvs, grading_vector, int(max_deg))  # type: ignore[arg-type]
        return [(b, (x if prec is None else mp.mpf(x.numerator) / x.denominator)) for b, x in gws.items()]
    if backend != "cygv":
        msg = f"unknown backend {backend!r} (use 'cygv' or 'cgv')"
        raise ValueError(msg)
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
