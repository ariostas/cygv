from __future__ import annotations

import gzip
import json
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
