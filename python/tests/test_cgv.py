from __future__ import annotations

import gzip
import json
from fractions import Fraction
from math import gcd
from pathlib import Path
from typing import Any

import pytest

from cygv import compute_gv, compute_gw

# cgv/tools/cgv_run.py (bundled unchanged) leaves a few files for the garbage collector to close
pytestmark = [
    pytest.mark.filterwarnings("ignore::ResourceWarning"),
    pytest.mark.filterwarnings("ignore::pytest.PytestUnraisableExceptionWarning"),
]

REPO = Path(__file__).resolve().parents[2]
REFS = REPO / "cgv" / "tests" / "refs"


def test_bundled_cgv_run_is_a_copy() -> None:
    """python/cygv/_cgv_run.py must stay identical to cgv/tools/cgv_run.py."""
    tools = REPO / "cgv" / "tools" / "cgv_run.py"
    if not tools.exists():
        pytest.skip("not a repository checkout")
    bundled = REPO / "python" / "cygv" / "_cgv_run.py"
    assert bundled.read_bytes() == tools.read_bytes()


def test_cgv_backend_matches() -> None:
    """backend='cgv' gives cygv's GV and GW invariants exactly."""
    for name in [
        "quintic_D10",
        "h10_72_3998_0_default_K5",
        "h10_92_226_0_plike_K7",
        "h10_52_8839_0_default_K8",
    ]:
        path = REFS / f"{name}.json.gz"
        if not path.exists():
            pytest.skip("not a repository checkout")
        with gzip.open(path) as f:
            ref = json.load(f)
        d = ref["input"]
        kw: dict[str, Any] = {
            "generators": d["generators"],
            "grading_vector": d["grading_vector"],
            "q": d["q"],
            "intnums": {(i, j, k): x for i, j, k, x in d["intnums"]},
            "max_deg": ref["max_deg"],
        }
        assert dict(compute_gv(**kw, backend="cgv")) == dict(compute_gv(**kw)), name
        assert dict(compute_gw(**kw, backend="cgv")) == dict(compute_gw(**kw)), name


def test_cgv_backend_two_parameter_model() -> None:
    """P(1,1,2,2,2)[8] (h11 = 2), with only the two Mori cone generators.

    Not the h11 = 2 input of test_hkty.py: those smoke-test numbers are not a consistent
    geometry, and cgv, which combines the h11 instanton series weighted by the grading vector,
    then depends on the grading vector. On real geometries both backends agree.
    """
    kw: dict[str, Any] = {
        "generators": [[1, 0], [0, 1]],
        "grading_vector": [1, 1],
        "q": [[1, 0, 0, 0, 1, -2], [0, 1, 1, 1, 0, 1]],
        "intnums": {(0, 1, 1): 4, (1, 1, 1): 8},
        "max_deg": 10,
    }
    gv = dict(compute_gv(**kw, backend="cgv"))
    assert gv[(1, 0)] == 4
    assert gv[(0, 1)] == 640  # the classic values of this model
    assert gv == dict(compute_gv(**kw))
    assert dict(compute_gw(**kw, backend="cgv")) == dict(compute_gw(**kw))


def test_cgv_backend_devices(monkeypatch: pytest.MonkeyPatch) -> None:
    """device='cpu' always; device='gpu' if this cygv was built with cgv's GPU variant."""
    from cygv.cygv import _cgv_gpu_executable  # noqa: PLC0415

    kw: dict[str, Any] = {
        "generators": [[1, 0], [0, 1]],
        "grading_vector": [1, 1],
        "q": [[1, 0, 0, 0, 1, -2], [0, 1, 1, 1, 0, 1]],
        "intnums": {(0, 1, 1): 4, (1, 1, 1): 8},
        "max_deg": 12,
    }
    ref = dict(compute_gv(**kw))
    assert dict(compute_gv(**kw, backend="cgv", device="cpu")) == ref
    if _cgv_gpu_executable() is None:
        with pytest.raises(ValueError, match="without cgv's GPU variant"):
            compute_gv(**kw, backend="cgv", device="gpu")
    else:
        monkeypatch.setenv(
            "CGV_GPU_MIN", "0"
        )  # else cgv keeps a problem this small on the CPU
        assert dict(compute_gv(**kw, backend="cgv", device="gpu")) == ref


def test_cgv_backend_rejects_unsupported() -> None:
    kw: dict[str, Any] = {
        "generators": [[0, -1], [1, 2]],
        "grading_vector": [3, -1],
        "q": [[1, 1, 1, 0, 1, 2], [0, 0, -1, 1, 1, -1]],
        "intnums": {(0, 0, 0): 2, (0, 0, 1): 1, (0, 1, 1): -1, (1, 1, 1): 5},
    }
    with pytest.raises(NotImplementedError):
        compute_gv(**kw, min_points=100, backend="cgv")
    with pytest.raises(ValueError, match="unknown backend"):
        compute_gv(**kw, max_deg=10, backend="nope")


def test_bundled_cgv_phase_is_a_copy() -> None:
    """python/cygv/_cgv_phase.py must stay identical to cgv/tools/cgv_phase.py."""
    tools = REPO / "cgv" / "tools" / "cgv_phase.py"
    if not tools.exists():
        pytest.skip("not a repository checkout")
    bundled = REPO / "python" / "cygv" / "_cgv_phase.py"
    assert bundled.read_bytes() == tools.read_bytes()


def _phase_ref(name: str) -> dict[str, Any]:
    import shutil  # noqa: PLC0415

    path = REPO / "cgv" / "tests" / "refs_vex" / f"{name}.json.gz"
    if not path.exists():
        pytest.skip("not a repository checkout")
    if shutil.which("normaliz") is None:
        pytest.skip("compute_gv_phase needs normaliz")
    with gzip.open(path) as f:
        ref: dict[str, Any] = json.load(f)
    return ref


def test_phase_frst_and_vex() -> None:
    """compute_gv_phase in FRST and vex phases (cgv/tests/refs_vex), GV and GW."""
    from cygv import compute_gv_phase, compute_gw_phase  # noqa: PLC0415

    for name in ["liam_h2_fan0", "liam_h2_fan1", "liam_h3_fan4", "liam_h3_fan9"]:
        r = _phase_ref(name)
        kappa = {(i, j, k): v for i, j, k, v in r["kappa"]}
        kw: dict[str, Any] = {"grading_vector": r["grading"]}
        gv = dict(compute_gv_phase(r["cones"], r["q"], kappa, r["max_deg"], **kw))
        assert gv == {tuple(k): v for k, v in r["gvs"]}, name
        gw = dict(compute_gw_phase(r["cones"], r["q"], kappa, r["max_deg"], **kw))
        assert all(isinstance(x, Fraction) for x in gw.values())
        primitive = [c for c in gv if gcd(*c) == 1]
        assert primitive, name
        assert all(gw[c] == gv[c] for c in primitive), (
            name
        )  # GW = GV on primitive classes


def test_phase_vex_equals_frst_of_another_polytope() -> None:
    """A vex phase and an FRST of a different polytope with the same CY give the same GVs."""
    import numpy as np  # noqa: PLC0415

    from cygv import compute_gv_phase  # noqa: PLC0415

    r = _phase_ref("pair_h4_1")
    T = np.array(r["T"])
    w_vex = np.array(r["grading_vex"])
    w_frst = np.round(np.linalg.inv(T) @ w_vex).astype(int)
    sides = {}
    for side, w in (("vex", w_vex), ("frst", w_frst)):
        s = r[side]
        kappa = {(i, j, k): v for i, j, k, v in s["kappa"]}
        sides[side] = dict(
            compute_gv_phase(
                s["cones"], s["q"], kappa, r["max_deg"], grading_vector=w.tolist()
            )
        )
    mapped = {
        tuple(int(x) for x in T.T @ np.array(n)): v for n, v in sides["vex"].items()
    }
    assert mapped == sides["frst"]
    assert len(mapped) > 100
