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
    device: str,
) -> dict[tuple[int, ...], int]:
    """GV invariants with the bundled cgv program (cgv/ in the repository), through
    the same code as its standalone Python interface (cgv/tools/cgv_run.py).

    device: "cpu", "gpu" (GPU 0), "gpu:N", or "auto" (a suitable GPU if this cygv was
    built with cgv's GPU variant, else the CPU); see cgv_run.pick_device."""
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
    from cygv.cygv import _cgv_executable, _cgv_gpu_executable  # noqa: PLC0415

    d = {
        "generators": [[int(x) for x in g] for g in np.array(generators, dtype=int)],
        "grading_vector": [int(x) for x in np.array(grading_vector, dtype=int)],
        "q": [[int(x) for x in r] for r in np.array(q, dtype=int)],
        "intnums": [
            [int(i), int(j), int(k), int(v)] for (i, j, k), v in intnums.items()
        ],
    }
    gpu = _cgv_gpu_executable()
    gpu_bin, hip = gpu if gpu is not None else (None, False)
    extra = _cgv_run.pick_device(device, gpu_bin, hip=hip)  # type: ignore[no-untyped-call]
    if extra:
        if gpu_bin is None:
            msg = "this cygv was built without cgv's GPU variant (build from source with nvcc or hipcc)"
            raise ValueError(msg)
        _cgv_run.BIN = gpu_bin
    else:
        _cgv_run.BIN = _cgv_executable()
    threads = os.cpu_count() or 1
    try:
        out = _cgv_run.run_cgv(d, int(max_deg), threads, extra=extra)  # type: ignore[no-untyped-call]
    except FileNotFoundError:
        # no normaliz for the cone data: cgv enumerates the cone itself (same result, slower)
        out = _cgv_run.run_cgv(d, int(max_deg), threads, extra=extra, cones=False)  # type: ignore[no-untyped-call]
    gvs: dict[tuple[int, ...], int] = out[0]
    return gvs


def _cgv_phase_gvs(
    cones: ArrayLike,
    q: ArrayLike,
    intnums: dict[tuple[int, int, int], int],
    max_deg: int,
    generators: ArrayLike | None,
    grading_vector: ArrayLike | None,
    saturate: bool,
    device: str,
) -> tuple[dict[tuple[int, ...], int], list[int]]:
    """GV invariants in the phase given by a fan, with the bundled cgv program, through the same
    code as its standalone Python interface (cgv/tools/cgv_phase.py). Returns (GVs, grading)."""
    from cygv import _cgv_phase, _cgv_run  # noqa: PLC0415
    from cygv.cygv import _cgv_executable, _cgv_gpu_executable  # noqa: PLC0415

    try:
        d = _cgv_phase.prepare(  # type: ignore[no-untyped-call]
            [[int(x) for x in c] for c in cones],
            [[int(x) for x in r] for r in np.array(q, dtype=int)],
            {(int(i), int(j), int(k)): int(v) for (i, j, k), v in intnums.items()},
            None
            if generators is None
            else [[int(x) for x in g] for g in np.array(generators, dtype=int)],
            None
            if grading_vector is None
            else [int(x) for x in np.array(grading_vector, dtype=int)],
            saturate,
        )
    except FileNotFoundError as e:
        msg = "compute_gv_phase needs normaliz (e.g. conda install -c conda-forge normaliz)"
        raise RuntimeError(msg) from e
    gpu = _cgv_gpu_executable()
    gpu_bin, hip = gpu if gpu is not None else (None, False)
    extra = _cgv_run.pick_device(device, gpu_bin, hip=hip)  # type: ignore[no-untyped-call]
    if extra and gpu_bin is None:
        msg = "this cygv was built without cgv's GPU variant (build from source with nvcc or hipcc)"
        raise ValueError(msg)
    binary = gpu_bin if extra else _cgv_executable()
    gvs, _ = _cgv_phase.run(
        d, int(max_deg), os.cpu_count() or 1, binary=binary, extra=extra
    )  # type: ignore[no-untyped-call]
    grading: list[int] = d["grading_vector"]
    return gvs, grading


def compute_gv_phase(
    cones: ArrayLike,
    q: ArrayLike,
    intnums: dict[tuple[int, int, int], int],
    max_deg: int,
    generators: ArrayLike | None = None,
    grading_vector: ArrayLike | None = None,
    saturate: bool = True,
    device: str = "auto",
) -> list[Any]:
    """GV invariants of a CY threefold hypersurface in a given phase of the ambient toric variety:
    an FRST or a vex fan (a fine regular fan that does not refine the face fan). Uses cgv.

    cones: maximal cones of the fan, as tuples of column indices of q. q, intnums: as for compute_gv.
    generators: Mori cone generators, or a subset (lightcone GVs); default: the fan's wall curves.
    saturate: True uses every lattice point of the cone they span (normaliz Hilbert basis); False uses
    exactly their semigroup. grading_vector: default an interior point of the dual cone. The result is
    in the same format as compute_gv's. Needs normaliz. See cgv/README.md, "Any phase"."""
    gvs, _ = _cgv_phase_gvs(
        cones, q, intnums, max_deg, generators, grading_vector, saturate, device
    )
    return list(gvs.items())


def compute_gw_phase(
    cones: ArrayLike,
    q: ArrayLike,
    intnums: dict[tuple[int, int, int], int],
    max_deg: int,
    generators: ArrayLike | None = None,
    grading_vector: ArrayLike | None = None,
    saturate: bool = True,
    device: str = "auto",
    prec: int | None = None,
) -> list[Any]:
    """Genus-0 GW invariants in a given phase (see compute_gv_phase), from its GV invariants."""
    if prec is not None:
        mp.mp.prec = prec
    gvs, grading = _cgv_phase_gvs(
        cones, q, intnums, max_deg, generators, grading_vector, saturate, device
    )
    gws = _gw_from_gv(gvs, grading, max_deg)
    return [
        (b, (x if prec is None else mp.mpf(x.numerator) / x.denominator))
        for b, x in gws.items()
    ]


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
    device: str = "auto",
) -> list[Any]:
    if backend == "cgv":
        return list(
            _cgv_gvs(
                generators,
                grading_vector,
                q,
                intnums,
                max_deg,
                min_points,
                target_points,
                nefpart,
                device,
            ).items()
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
    device: str = "auto",
) -> list[Any]:
    if prec is not None:
        mp.mp.prec = prec
    if backend == "cgv":
        gvs = _cgv_gvs(
            generators,
            grading_vector,
            q,
            intnums,
            max_deg,
            min_points,
            target_points,
            nefpart,
            device,
        )
        assert max_deg is not None  # checked by _cgv_gvs
        gws = _gw_from_gv(gvs, grading_vector, max_deg)
        return [
            (b, (x if prec is None else mp.mpf(x.numerator) / x.denominator))
            for b, x in gws.items()
        ]
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
