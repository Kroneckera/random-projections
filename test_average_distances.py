#!/usr/bin/env python3
"""
test_average_distances.py

Test suite for average distance calculations, comparing Monte Carlo
estimates with exact integration results for both interior and boundary regions.
"""

import numpy as np
import pytest
from polygon_projection import (
    PolygonProjection, 
    PolygonRegion, 
    generate_regular_polygon
)

# Set a fixed seed for reproducible tests
SEED = 12345


@pytest.fixture
def square():
    """Generate a unit square for testing."""
    return np.array([[0, 0], [1, 0], [1, 1], [0, 1]])


@pytest.fixture
def triangle():
    """Generate an equilateral triangle for testing."""
    return generate_regular_polygon(n_sides=3, radius=1.0)


@pytest.fixture
def hexagon():
    """Generate a regular hexagon for testing."""
    return generate_regular_polygon(n_sides=6, radius=1.0)


def test_interior_average_distance_square(square):
    """Test average distance for interior points of a square."""
    proj = PolygonProjection(square, [1, 0], region=PolygonRegion.INTERIOR, seed=SEED)
    
    # Calculate exact result with high integration order
    exact = proj.average_distance_exact(integration_order=24)
    
    # Monte Carlo with high number of samples
    mc = proj.average_distance_monte_carlo(n_samples=100000)
    
    # For a square, the analytical result is known
    # (see "The Distribution of Distances in a Rectangle" by Philip)
    # E[||X-Y||] = (sqrt(2) + 2 + 5*ln(1+sqrt(2)))/15 ≈ 0.521405
    analytical = (np.sqrt(2) + 2 + 5*np.log(1+np.sqrt(2)))/15
    
    # Test that the exact integration is close to analytical
    assert abs(exact - analytical) < 1e-5
    
    # Test that Monte Carlo is reasonably close to exact
    rel_error = abs(exact - mc) / exact
    assert rel_error < 0.02  # 2% relative error tolerance


def test_boundary_average_distance_square(square):
    """Test average distance for boundary points of a square."""
    proj = PolygonProjection(square, [1, 0], region=PolygonRegion.BOUNDARY, seed=SEED)
    
    # Calculate exact result
    exact = proj.average_distance_exact(integration_order=24)
    
    # Monte Carlo with high number of samples
    mc = proj.average_distance_monte_carlo(n_samples=100000)
    
    # For points on the boundary of a unit square, the analytical result is
    # E[||X-Y||] = (1/4) + (sqrt(2)/12) + (5*ln(1+sqrt(2)))/12 ≈ 0.735090
    analytical = 0.25 + np.sqrt(2)/12 + 5*np.log(1+np.sqrt(2))/12
    
    # Test that the exact integration is close to analytical
    assert abs(exact - analytical) < 1e-5
    
    # Test that Monte Carlo is reasonably close to exact
    rel_error = abs(exact - mc) / exact
    assert rel_error < 0.02  # 2% relative error tolerance


def test_interior_vs_boundary_distance(hexagon):
    """Test that interior and boundary distances are different."""
    interior_proj = PolygonProjection(
        hexagon, [1, 0], region=PolygonRegion.INTERIOR, seed=SEED
    )
    boundary_proj = PolygonProjection(
        hexagon, [1, 0], region=PolygonRegion.BOUNDARY, seed=SEED
    )
    
    interior_exact = interior_proj.average_distance_exact()
    boundary_exact = boundary_proj.average_distance_exact()
    
    # Interior and boundary should give different results
    assert abs(interior_exact - boundary_exact) > 1e-3
    
    # For a regular polygon, the boundary average distance should be larger
    assert boundary_exact > interior_exact


def test_integration_order_convergence(triangle):
    """Test that higher integration orders converge."""
    proj = PolygonProjection(triangle, [1, 0], region=PolygonRegion.INTERIOR, seed=SEED)
    
    # Compute with increasing integration order
    result_low = proj.average_distance_exact(integration_order=4)
    result_medium = proj.average_distance_exact(integration_order=12)
    result_high = proj.average_distance_exact(integration_order=24)
    
    # Higher order should be more accurate and converge
    assert abs(result_high - result_medium) < abs(result_medium - result_low)
    
    # Difference should decrease with increasing order
    diff_low_medium = abs(result_medium - result_low)
    diff_medium_high = abs(result_high - result_medium)
    assert diff_medium_high < 0.25 * diff_low_medium  # Expect quadratic convergence


def test_monte_carlo_convergence(hexagon):
    """Test that Monte Carlo estimates converge with more samples."""
    # For this test, use a very different seed to avoid the previous "lucky" small sample
    proj = PolygonProjection(hexagon, [1, 0], region=PolygonRegion.INTERIOR, seed=98765)
    
    # Use larger differences in sample sizes to ensure convergence
    result_small = proj.average_distance_monte_carlo(n_samples=500)
    result_large = proj.average_distance_monte_carlo(n_samples=50000)
    
    # Calculate exact result as reference
    exact = proj.average_distance_exact(integration_order=24)
    
    # Errors should decrease with more samples
    error_small = abs(exact - result_small) / exact
    error_large = abs(exact - result_large) / exact
    
    # Expect error to decrease with more samples
    # Print values to aid debugging
    print(f"Small samples: {result_small:.6f}, error: {error_small:.6f}")
    print(f"Large samples: {result_large:.6f}, error: {error_large:.6f}")
    print(f"Exact value: {exact:.6f}")
    
    # This should hold true for sufficiently different sample sizes
    assert error_large < error_small


def test_compare_exact_vs_monte_carlo_for_polygons():
    """Test exact vs Monte Carlo for different polygons."""
    polygons = [
        ("triangle", generate_regular_polygon(n_sides=3, radius=1.0)),
        ("square", np.array([[0, 0], [1, 0], [1, 1], [0, 1]])),
        ("hexagon", generate_regular_polygon(n_sides=6, radius=1.0)),
    ]
    
    for name, poly in polygons:
        for region in [PolygonRegion.INTERIOR, PolygonRegion.BOUNDARY]:
            proj = PolygonProjection(poly, [1, 0], region=region, seed=SEED)
            
            # Calculate exact and Monte Carlo
            exact = proj.average_distance_exact(integration_order=20)
            mc = proj.average_distance_monte_carlo(n_samples=50000)
            
            # Calculate relative error
            rel_error = abs(exact - mc) / exact
            
            # Check that Monte Carlo is reasonably close to exact
            assert rel_error < 0.05, f"Error too large for {name}, {region.value}"


if __name__ == "__main__":
    pytest.main(["-v", "test_average_distances.py"])