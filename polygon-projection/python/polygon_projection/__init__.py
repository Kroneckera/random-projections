"""Polygon Projection Library

A library for calculating projections and distances in convex polygons.

This module provides tools for:
- Working with convex polygons
- Calculating exact and Monte Carlo projections
- Computing distances between points in polygons

Key components:
- Polygon: Represents a convex polygon
- Region: Enum for specifying interior or boundary regions
- ProjectionCalculator: Handles projection calculations
"""

import importlib.util
import sys
import warnings
import os

# Check if we should force the use of fallback implementation
_FORCE_FALLBACK = os.environ.get("POLYGON_PROJECTION_FORCE_FALLBACK", "0").lower() in ("1", "true", "yes")

# Try to import from the C++ API first (unless fallback is forced)
if not _FORCE_FALLBACK:
    try:
        from .api import Polygon, Region, ProjectionCalculator
        _USING_CPP_BACKEND = True
    except ImportError:
        # Fall back to the Python implementation if C++ backend is not available
        warnings.warn(
            "C++ backend not available, using pure Python fallback implementation. "
            "This may be significantly slower. Consider building the C++ extension "
            "for better performance.",
            UserWarning
        )
        from .fallback import Polygon, Region, ProjectionCalculator, create_polygon, create_regular_polygon
        _USING_CPP_BACKEND = False
else:
    # Deliberately use the fallback implementation
    from .fallback import Polygon, Region, ProjectionCalculator, create_polygon, create_regular_polygon
    _USING_CPP_BACKEND = False
    warnings.warn(
        "Using pure Python fallback implementation as requested via environment variable. "
        "This may be significantly slower than the C++ backend.",
        UserWarning
    )

__version__ = "0.1.0"
# Add additional helper functions for the fallback API
if not _USING_CPP_BACKEND:
    from .fallback import create_polygon, create_regular_polygon, is_convex
    __all__ = ["Polygon", "Region", "ProjectionCalculator", "create_polygon", "create_regular_polygon", "is_convex"]
else:
    # For C++ backend, define equivalent functions
    from .fallback import is_convex  # Use the fallback implementation for is_convex
    __all__ = ["Polygon", "Region", "ProjectionCalculator", "is_convex"]

def using_cpp_backend():
    """Check if the library is using the C++ backend or Python fallback.
    
    Returns:
        bool: True if using C++ backend, False if using Python fallback
    """
    return _USING_CPP_BACKEND