#!/usr/bin/env python3
"""
Test suite for the polygon_projection module.

Tests both exact and Monte Carlo methods for interior and boundary projections.
"""

import numpy as np
import pytest
from polygon_projection import PolygonProjection, PolygonRegion, generate_regular_polygon


@pytest.fixture
def square():
    """Generate a unit square for testing."""
    return np.array([[0, 0], [1, 0], [1, 1], [0, 1]])


@pytest.fixture
def hexagon():
    """Generate a regular hexagon for testing."""
    return generate_regular_polygon(n_sides=6)


def test_interior_exact_vs_monte_carlo(square):
    """Test agreement between exact and Monte Carlo methods for interior projections."""
    directions = [
        [1, 0],     # x-axis
        [0, 1],     # y-axis 
        [1, 1],     # diagonal
    ]
    
    for direction in directions:
        proj = PolygonProjection(square, direction, region=PolygonRegion.INTERIOR, seed=42)
        comparison = proj.compare_methods(n_samples=100000)
        
        # For large sample sizes, MC should be within 1% of exact
        assert comparison['relative_error'] < 0.01
        

def test_boundary_exact_vs_monte_carlo(square):
    """Test agreement between exact and Monte Carlo methods for boundary projections."""
    # Use random directions that are not parallel to any edge
    directions = [
        [1.3, 2.7],  # Random non-axis-aligned direction
        [3.1, 1.4],  # Another random direction
        [2.2, 3.5],  # Yet another random direction
    ]
    
    for direction in directions:
        proj = PolygonProjection(square, direction, region=PolygonRegion.BOUNDARY, seed=42)
        # Use more samples for better convergence
        comparison = proj.compare_methods(n_samples=500000)
        
        # For large sample sizes, MC should be within 1% of exact
        assert comparison['relative_error'] < 0.01


def test_compare_methods_returns_correct_data(hexagon):
    """Test that compare_methods returns the expected data structure."""
    direction = [1, 1]
    n_samples = 1000
    
    proj = PolygonProjection(hexagon, direction, region=PolygonRegion.INTERIOR, seed=42)
    result = proj.compare_methods(n_samples=n_samples)
    
    # Check that all expected keys are present
    expected_keys = ['region', 'exact', 'monte_carlo', 'absolute_error', 
                    'relative_error', 'n_samples']
    
    for key in expected_keys:
        assert key in result
    
    # Check specific values
    assert result['region'] == 'interior'
    assert result['n_samples'] == n_samples
    assert isinstance(result['exact'], float)
    assert isinstance(result['monte_carlo'], float)


def test_interior_vs_boundary_different_results(hexagon):
    """Test that interior and boundary projections give different results."""
    direction = [1, 0.5]
    
    interior_proj = PolygonProjection(hexagon, direction, 
                                      region=PolygonRegion.INTERIOR, seed=42)
    boundary_proj = PolygonProjection(hexagon, direction, 
                                      region=PolygonRegion.BOUNDARY, seed=42)
    
    interior_exact = interior_proj.exact_expected_distance()
    boundary_exact = boundary_proj.exact_expected_distance()
    
    # The two methods should give different results for a non-trivial polygon
    assert abs(interior_exact - boundary_exact) > 1e-6


def test_sampling_functions(square):
    """Test that sampling functions return points with the expected properties."""
    # For interior points of a unit square
    interior_proj = PolygonProjection(square, [1, 0], region=PolygonRegion.INTERIOR, seed=42)
    for _ in range(100):
        point = interior_proj.sample_point()
        assert 0 <= point[0] <= 1
        assert 0 <= point[1] <= 1
    
    # For boundary points of a unit square
    boundary_proj = PolygonProjection(square, [1, 0], region=PolygonRegion.BOUNDARY, seed=42)
    for _ in range(100):
        point = boundary_proj.sample_point()
        # Point should be on the boundary
        on_boundary = (
            (abs(point[0]) < 1e-10 and 0 <= point[1] <= 1) or
            (abs(point[0] - 1) < 1e-10 and 0 <= point[1] <= 1) or
            (abs(point[1]) < 1e-10 and 0 <= point[0] <= 1) or
            (abs(point[1] - 1) < 1e-10 and 0 <= point[0] <= 1)
        )
        assert on_boundary


def test_regular_polygon_generation():
    """Test that generate_regular_polygon produces correct results."""
    triangle = generate_regular_polygon(n_sides=3)
    assert len(triangle) == 3
    
    square = generate_regular_polygon(n_sides=4)
    assert len(square) == 4
    
    # For a regular polygon with radius 1, all vertices should be at distance 1 from origin
    octagon = generate_regular_polygon(n_sides=8, radius=1.0)
    for vertex in octagon:
        assert abs(np.linalg.norm(vertex) - 1.0) < 1e-10


if __name__ == "__main__":
    pytest.main(["-v", "test_polygon_projection.py"])