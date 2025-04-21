#!/usr/bin/env python3
"""
Test suite for the projection module.

Tests both Python API and C++ backend functionality through the Python API.
"""

import numpy as np
import pytest
import sys
import os
from pathlib import Path

# Add the package directory to the Python path when running standalone
package_root = Path(__file__).parent.parent.parent
python_dir = package_root / "python"
sys.path.insert(0, str(python_dir))

# Now import from the package
from projection.api import Polygon, Region, ProjectionCalculator


@pytest.fixture
def square():
    """Generate a unit square for testing."""
    return np.array([[0, 0], [1, 0], [1, 1], [0, 1]])


@pytest.fixture
def hexagon():
    """Generate a regular hexagon for testing."""
    return Polygon.regular(sides=6).vertices


def test_polygon_creation():
    """Test polygon creation through different methods."""
    # Direct creation from vertices
    square_verts = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
    square = Polygon(square_verts)
    assert np.array_equal(square.vertices, square_verts)
    
    # Creation using regular polygon method
    hexagon = Polygon.regular(sides=6, radius=2.0)
    assert len(hexagon.vertices) == 6
    
    # Check radius of vertices (should all be 2.0 units from origin)
    for vertex in hexagon.vertices:
        assert abs(np.linalg.norm(vertex) - 2.0) < 1e-10


def test_polygon_properties(square, hexagon):
    """Test polygon properties like area and perimeter."""
    square_polygon = Polygon(square)
    assert abs(square_polygon.area - 1.0) < 1e-10  # Unit square has area 1.0
    assert abs(square_polygon.perimeter - 4.0) < 1e-10  # Unit square has perimeter 4.0
    
    hex_polygon = Polygon(hexagon)
    # Regular hexagon with radius 1 has area 3*sqrt(3)/2
    assert abs(hex_polygon.area - 3*np.sqrt(3)/2) < 1e-10
    # Regular hexagon with radius 1 has perimeter 6
    assert abs(hex_polygon.perimeter - 6.0) < 1e-10


def test_is_convex(square):
    """Test convexity checking."""
    convex = Polygon(square)
    assert convex.is_convex() == True
    
    # Create a non-convex shape - should raise an exception
    # Concave polygon (arrow shape)
    non_convex = np.array([[0, 0], [2, 0], [1, 1], [2, 2], [0, 2]])
    
    # The API should detect this is not convex and raise an exception
    with pytest.raises(ValueError):
        Polygon(non_convex)


def test_exact_vs_monte_carlo_interior(square):
    """Test agreement between exact and Monte Carlo methods for interior projections."""
    directions = [
        [1, 0],     # x-axis
        [0, 1],     # y-axis 
        [1, 1],     # diagonal
    ]
    
    for direction in directions:
        calc = ProjectionCalculator(
            Polygon(square), 
            direction, 
            region=Region.INTERIOR, 
            random_seed=42
        )
        
        exact = calc.projected_distance(exact=True)
        monte_carlo = calc.projected_distance(exact=False, samples=100000)
        
        # For Monte Carlo with 100k samples, allow up to 0.3% error
        relative_error = abs(exact - monte_carlo) / exact
        assert relative_error < 0.003


def test_exact_vs_monte_carlo_boundary(square):
    """Test agreement between exact and Monte Carlo methods for boundary projections."""
    # Use random directions that are not parallel to any edge
    directions = [
        [1.3, 2.7],  # Random non-axis-aligned direction
        [3.1, 1.4],  # Another random direction
        [2.2, 3.5],  # Yet another random direction
    ]
    
    for direction in directions:
        calc = ProjectionCalculator(
            Polygon(square), 
            direction, 
            region=Region.BOUNDARY, 
            random_seed=42
        )
        
        exact = calc.projected_distance(exact=True)
        monte_carlo = calc.projected_distance(exact=False, samples=500000)
        
        # For Monte Carlo with 500k samples, allow up to 0.2% error
        relative_error = abs(exact - monte_carlo) / exact
        assert relative_error < 0.002


def test_sampling_functions(square):
    """Test that sampling functions return points with the expected properties."""
    # For interior points of a unit square
    interior_calc = ProjectionCalculator(
        Polygon(square), [1, 0], 
        region=Region.INTERIOR, random_seed=42
    )
    
    # Sample multiple points at once
    interior_points = interior_calc.sample_points(count=100)
    assert interior_points.shape == (100, 2)
    
    for point in interior_points:
        assert 0 <= point[0] <= 1
        assert 0 <= point[1] <= 1
    
    # For boundary points of a unit square
    boundary_calc = ProjectionCalculator(
        Polygon(square), [1, 0], 
        region=Region.BOUNDARY, random_seed=42
    )
    
    boundary_points = boundary_calc.sample_points(count=100)
    assert boundary_points.shape == (100, 2)
    
    for point in boundary_points:
        # Point should be on the boundary (with numerical tolerance)
        on_boundary = (
            (abs(point[0]) < 1e-10 and 0 <= point[1] <= 1) or
            (abs(point[0] - 1) < 1e-10 and 0 <= point[1] <= 1) or
            (abs(point[1]) < 1e-10 and 0 <= point[0] <= 1) or
            (abs(point[1] - 1) < 1e-10 and 0 <= point[0] <= 1)
        )
        assert on_boundary


def test_monte_carlo_analysis(square):
    """Test the monte_carlo_analysis generator."""
    calc = ProjectionCalculator(
        Polygon(square), [1, 0], 
        region=Region.INTERIOR, random_seed=42
    )
    
    # Small analysis for testing
    total_samples = 1000
    batch_size = 200
    
    # Collect results from generator
    results = list(calc.monte_carlo_analysis(
        samples=total_samples,
        batch_size=batch_size,
        yield_interval=1
    ))
    
    # Check we got the expected number of results
    assert len(results) == total_samples // batch_size
    
    # Check the final result has the expected properties
    final = results[-1]
    assert final['samples_processed'] == total_samples
    assert final['total_samples'] == total_samples
    assert 'projected_distance' in final
    assert 'average_distance' in final
    assert 'progress' in final


def test_projection_functions(square):
    """Test projection functionality."""
    calc = ProjectionCalculator(
        Polygon(square), [1, 0], 
        region=Region.INTERIOR, random_seed=42
    )
    
    # Create test points
    points = np.array([
        [0, 0],  # Origin
        [1, 0],  # Right
        [0, 1],  # Top
        [1, 1],  # Top-right
        [0.5, 0.5]  # Center
    ])
    
    # Project points onto the x-axis
    projections = calc.project_points(points)
    
    # Check projections against expected values
    expected = np.array([0.0, 1.0, 0.0, 1.0, 0.5])
    assert np.allclose(projections, expected)
    
    # Test with non-axis-aligned direction
    calc2 = ProjectionCalculator(
        Polygon(square), [1, 1], 
        region=Region.INTERIOR, random_seed=42
    )
    
    # Project onto the diagonal
    diagonal_proj = calc2.project_points(points)
    
    # Expected projections (normalized by sqrt(2))
    expected2 = np.array([
        0.0,               # [0,0]
        1.0/np.sqrt(2),    # [1,0]
        1.0/np.sqrt(2),    # [0,1]
        2.0/np.sqrt(2),    # [1,1]
        1.0/np.sqrt(2)     # [0.5,0.5]
    ])
    
    assert np.allclose(diagonal_proj, expected2)


def test_compare_methods(hexagon):
    """Test that compare_methods returns the expected data structure."""
    direction = [1, 1]
    n_samples = 1000
    
    calc = ProjectionCalculator(
        Polygon(hexagon), direction, 
        region=Region.INTERIOR, random_seed=42
    )
    
    result = calc.compare_methods(samples=n_samples)
    
    # Check that all expected keys are present
    expected_keys = ['region', 'exact', 'monte_carlo', 'absolute_error', 
                    'relative_error', 'n_samples']
    
    for key in expected_keys:
        assert key in result
    
    # Check specific values
    assert result['region'] == 'INTERIOR'
    assert result['n_samples'] == n_samples
    assert isinstance(result['exact'], float)
    assert isinstance(result['monte_carlo'], float)


def test_error_handling():
    """Test error handling in the API."""
    # Invalid polygon (too few vertices)
    with pytest.raises(ValueError):
        Polygon(np.array([[0, 0], [1, 0]]))
    
    # Invalid polygon (non-convex)
    with pytest.raises(ValueError):
        Polygon(np.array([[0, 0], [2, 0], [1, 1], [2, 2], [0, 2]]))
    
    # Invalid direction (zero vector)
    square = Polygon(np.array([[0, 0], [1, 0], [1, 1], [0, 1]]))
    with pytest.raises(ValueError):
        ProjectionCalculator(square, [0, 0])
    
    # Invalid sample count
    calc = ProjectionCalculator(square, [1, 0])
    with pytest.raises(ValueError):
        calc.projected_distance(exact=False, samples=-10)


if __name__ == "__main__":
    pytest.main(["-v"])