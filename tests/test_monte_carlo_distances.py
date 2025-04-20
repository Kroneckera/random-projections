#!/usr/bin/env python3
"""
test_monte_carlo_distances.py

Test Monte Carlo distance calculations against known analytical formulas
for basic shapes.
"""

import numpy as np
import pytest
from polygon_projection import PolygonProjection, PolygonRegion

# Set a fixed seed for reproducible tests
SEED = 12345


@pytest.fixture
def unit_square():
    """Generate a unit square (side length 1) for testing."""
    return np.array([[0, 0], [1, 0], [1, 1], [0, 1]])


@pytest.fixture
def unit_disk_approximation():
    """Generate a regular polygon approximating a unit disk."""
    n_sides = 32  # High enough to approximate a circle well
    angles = np.linspace(0, 2 * np.pi, n_sides, endpoint=False)
    return np.column_stack([np.cos(angles), np.sin(angles)])


def test_square_interior_monte_carlo():
    """Test Monte Carlo for interior of unit square."""
    square = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
    
    # Create a sampler for interior points
    proj = PolygonProjection(square, [1, 0], region=PolygonRegion.INTERIOR, seed=SEED)
    
    # Run Monte Carlo with large sample size
    mc_result = proj.average_distance_monte_carlo(n_samples=500000)
    
    # Analytical value for unit square L(1,1)
    # From literature: E[||X-Y||] = (sqrt(2) + 2 + 5*ln(1+sqrt(2)))/15 ≈ 0.521405
    analytical = (np.sqrt(2) + 2 + 5*np.log(1+np.sqrt(2)))/15
    
    # Check that Monte Carlo is close to analytical
    rel_error = abs(mc_result - analytical) / analytical
    print(f"Square interior: MC = {mc_result:.6f}, Analytical = {analytical:.6f}, " 
          f"error = {rel_error:.6f}")
    
    # This should be within 1% with this many samples
    assert rel_error < 0.01


def test_square_perimeter_monte_carlo():
    """Test Monte Carlo for perimeter of unit square."""
    square = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
    perimeter = 4.0  # Perimeter of a unit square
    
    # Custom sampling function for uniform points on the perimeter
    # This bypasses our implementation to get a ground-truth comparison
    def sample_square_perimeter():
        # Choose which side (0-3)
        side = np.random.randint(0, 4)
        # Position along that side (0-1)
        t = np.random.random()
        
        if side == 0:   # Bottom
            return np.array([t, 0])
        elif side == 1: # Right
            return np.array([1, t])
        elif side == 2: # Top
            return np.array([1-t, 1])
        else:           # Left
            return np.array([0, 1-t])
    
    # Run custom Monte Carlo
    n_samples = 500000
    distances = []
    np.random.seed(SEED)
    
    for _ in range(n_samples):
        p1 = sample_square_perimeter()
        p2 = sample_square_perimeter()
        distances.append(np.linalg.norm(p1 - p2))
    
    mc_result = np.mean(distances)
    
    # For unit square perimeter, the correct formula is:
    # E[||X-Y||] = (1/4) + (sqrt(2)/12) + (5*ln(1+sqrt(2)))/12 ≈ 0.735090
    analytical = 0.25 + np.sqrt(2)/12 + 5*np.log(1+np.sqrt(2))/12
    
    # Check that Monte Carlo is close to analytical
    rel_error = abs(mc_result - analytical) / analytical
    print(f"Square perimeter: MC = {mc_result:.6f}, Analytical = {analytical:.6f}, "
          f"error = {rel_error:.6f}")
    
    # This should be within 1% with this many samples
    assert rel_error < 0.01


def test_disk_interior_monte_carlo(unit_disk_approximation):
    """Test Monte Carlo for interior of unit disk."""
    disk = unit_disk_approximation
    
    # Custom sampling function for uniform points in the unit disk
    # This bypasses our implementation to get a ground-truth comparison
    def sample_unit_disk():
        r = np.sqrt(np.random.random())  # Square root for uniform distribution
        theta = np.random.random() * 2 * np.pi
        return np.array([r * np.cos(theta), r * np.sin(theta)])
    
    # Run custom Monte Carlo
    n_samples = 500000
    distances = []
    np.random.seed(SEED)
    
    for _ in range(n_samples):
        p1 = sample_unit_disk()
        p2 = sample_unit_disk()
        distances.append(np.linalg.norm(p1 - p2))
    
    mc_result = np.mean(distances)
    
    # For unit disk, E[||X-Y||] = 128/(45*π) ≈ 0.9054
    analytical = 128/(45*np.pi)
    
    # Check that Monte Carlo is close to analytical
    rel_error = abs(mc_result - analytical) / analytical
    print(f"Disk interior: MC = {mc_result:.6f}, Analytical = {analytical:.6f}, "
          f"error = {rel_error:.6f}")
    
    # This should be within 1% with this many samples
    assert rel_error < 0.01


def test_disk_boundary_monte_carlo(unit_disk_approximation):
    """Test Monte Carlo for perimeter of unit disk (circle)."""
    circle = unit_disk_approximation
    
    # Custom sampling function for uniform points on the unit circle
    # This bypasses our implementation to get a ground-truth comparison
    def sample_unit_circle():
        theta = np.random.random() * 2 * np.pi
        return np.array([np.cos(theta), np.sin(theta)])
    
    # Run custom Monte Carlo
    n_samples = 500000
    distances = []
    np.random.seed(SEED)
    
    for _ in range(n_samples):
        p1 = sample_unit_circle()
        p2 = sample_unit_circle()
        distances.append(np.linalg.norm(p1 - p2))
    
    mc_result = np.mean(distances)
    
    # For unit circle, E[||X-Y||] = 4/π ≈ 1.273
    analytical = 4/np.pi
    
    # Check that Monte Carlo is close to analytical
    rel_error = abs(mc_result - analytical) / analytical
    print(f"Circle perimeter: MC = {mc_result:.6f}, Analytical = {analytical:.6f}, "
          f"error = {rel_error:.6f}")
    
    # This should be within 1% with this many samples
    assert rel_error < 0.01


def test_our_monte_carlo_vs_direct():
    """Test our implementation against direct Monte Carlo sampling."""
    square = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
    
    # 1. Interior tests
    # Use our implementation
    proj_interior = PolygonProjection(square, [1, 0], region=PolygonRegion.INTERIOR, seed=SEED)
    our_mc_interior = proj_interior.average_distance_monte_carlo(n_samples=100000)
    
    # Direct Monte Carlo for comparison
    np.random.seed(SEED)
    distances_interior = []
    for _ in range(100000):
        # Sample uniform points in square
        p1 = np.random.random(2)
        p2 = np.random.random(2)
        distances_interior.append(np.linalg.norm(p1 - p2))
    direct_mc_interior = np.mean(distances_interior)
    
    # 2. Boundary tests
    # Use our implementation
    proj_boundary = PolygonProjection(square, [1, 0], region=PolygonRegion.BOUNDARY, seed=SEED)
    our_mc_boundary = proj_boundary.average_distance_monte_carlo(n_samples=100000)
    
    # Direct Monte Carlo for comparison (using the square_perimeter function)
    np.random.seed(SEED)
    distances_boundary = []
    for _ in range(100000):
        # Choose sides and positions
        side1, side2 = np.random.randint(0, 4, 2)
        t1, t2 = np.random.random(2)
        
        # Map to coordinates
        if side1 == 0:   # Bottom
            p1 = np.array([t1, 0])
        elif side1 == 1: # Right
            p1 = np.array([1, t1])
        elif side1 == 2: # Top
            p1 = np.array([1-t1, 1])
        else:           # Left
            p1 = np.array([0, 1-t1])
            
        if side2 == 0:   # Bottom
            p2 = np.array([t2, 0])
        elif side2 == 1: # Right
            p2 = np.array([1, t2])
        elif side2 == 2: # Top
            p2 = np.array([1-t2, 1])
        else:           # Left
            p2 = np.array([0, 1-t2])
            
        distances_boundary.append(np.linalg.norm(p1 - p2))
    direct_mc_boundary = np.mean(distances_boundary)
    
    # Print results for comparison
    print("\nMonte Carlo implementation comparisons:")
    print(f"Interior - Our: {our_mc_interior:.6f}, Direct: {direct_mc_interior:.6f}")
    print(f"Boundary - Our: {our_mc_boundary:.6f}, Direct: {direct_mc_boundary:.6f}")
    
    # Relative differences should be small
    interior_diff = abs(our_mc_interior - direct_mc_interior) / direct_mc_interior
    boundary_diff = abs(our_mc_boundary - direct_mc_boundary) / direct_mc_boundary
    
    assert interior_diff < 0.05
    assert boundary_diff < 0.05


if __name__ == "__main__":
    pytest.main(["-v", "test_monte_carlo_distances.py"])