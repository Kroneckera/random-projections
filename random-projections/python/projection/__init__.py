"""
Random‑Projections · Convex‑polygon geometry toolkit
====================================================

This package wraps a high‑performance C++17 backend (compiled via pybind11)
and provides a pure‑Python fallback for environments where the native
extension cannot be built.

Main symbols
------------
Polygon               — convex polygon container
Region                — enum {INTERIOR, BOUNDARY}
ProjectionCalculator  — exact / Monte‑Carlo projection utilities

Helpers available in **both** backends:
    is_convex(points)
Fallback‑only helpers:
    create_polygon(points)              (validate & return Polygon)
    create_regular_polygon(n, radius=1) (convenience generator)
"""

from __future__ import annotations

import importlib
import os
import warnings

__version__ = "0.1.0"

# ───────────────────────────  Backend selection  ────────────────────────────
_FORCE_FALLBACK = os.environ.get("PROJECTION_FORCE_FALLBACK", "0").lower() \
                  in {"1", "true", "yes"}

_USING_CPP_BACKEND: bool

if not _FORCE_FALLBACK:
    try:
        # 1. Import the compiled extension.  Raises ModuleNotFoundError if the
        #    C++ module is missing (not built) or ImportError if a dependency
        #    such as libstdc++ cannot be found.
        _core = importlib.import_module("projection._core")
        # 2. Python façade built directly on top of the C++ objects
        from .api import Polygon, Region, ProjectionCalculator
        _USING_CPP_BACKEND = True
    except (ModuleNotFoundError, ImportError) as exc:
        warnings.warn(
            "C++ backend unavailable – falling back to the pure‑Python "
            "implementation.  Performance will degrade.\n"
            f"Reason: {exc}",
            RuntimeWarning,
            stacklevel=2,
        )
        from .fallback import (
            Polygon,
            Region,
            ProjectionCalculator,
            create_polygon,
            create_regular_polygon,
        )
        _USING_CPP_BACKEND = False
else:
    # Fallback forced by environment variable
    from .fallback import (
        Polygon,
        Region,
        ProjectionCalculator,
        create_polygon,
        create_regular_polygon,
    )
    _USING_CPP_BACKEND = False
    warnings.warn(
        "Using pure‑Python backend because PROJECTION_FORCE_FALLBACK is set.",
        RuntimeWarning,
        stacklevel=2,
    )

# Helper common to both backends
from .fallback import is_convex  # lightweight, pure Python

# ───────────────────────────────  Re‑exports  ────────────────────────────────
__all__ = [
    "Polygon",
    "Region",
    "ProjectionCalculator",
    "is_convex",
]

if not _USING_CPP_BACKEND:
    # expose extra helpers that only make sense with the fallback
    __all__.extend(["create_polygon", "create_regular_polygon"])

# ─────────────────────────────  Public utility  ─────────────────────────────
def using_cpp_backend() -> bool:
    """Return **True** when the high‑performance C++ backend is active."""
    return _USING_CPP_BACKEND
