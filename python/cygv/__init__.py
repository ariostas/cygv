from __future__ import annotations

from importlib.metadata import version

from cygv.cygv import ComputationError
from cygv.hkty import (
    METHODS,
    Invariants,
    Problem,
    Report,
    compute_gv,
    compute_gw,
)

__version__ = version("cygv")

__all__ = [
    "METHODS",
    "ComputationError",
    "Invariants",
    "Problem",
    "Report",
    "compute_gv",
    "compute_gw",
]
