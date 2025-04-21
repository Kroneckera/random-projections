"""Legacy polygon projection implementations.

This package contains the legacy Python implementations of polygon projection
algorithms that are used as a fallback when the C++ backend is not available.
"""

from .polygon_projection import (
    PolygonRegion,
    PolygonProjection,
    calculate_polygon_area,
    calculate_polygon_perimeter,
    generate_regular_polygon,
    is_convex
)

from .projection_sampling import (
    projected_density,
    expected_abs_diff,
    sample_mc_pairwise_projected_distances,
    sample_point_in_convex_polygon,
    project_point_to_direction
)

from .edge_projection import (
    project_edges,
    boundary_expected_abs_diff,
    sample_point_on_boundary
)

from .piecewise_integration import (
    expected_distance_interior_exact
)

__all__ = [
    'PolygonRegion',
    'PolygonProjection',
    'calculate_polygon_area',
    'calculate_polygon_perimeter',
    'generate_regular_polygon',
    'is_convex',
    'projected_density',
    'expected_abs_diff',
    'sample_mc_pairwise_projected_distances',
    'sample_point_in_convex_polygon',
    'project_point_to_direction',
    'project_edges',
    'boundary_expected_abs_diff',
    'sample_point_on_boundary',
    'expected_distance_interior_exact'
]