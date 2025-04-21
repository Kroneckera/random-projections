#!/usr/bin/env python3
"""
Compare analytical formulas with Monte Carlo calculations.

This script verifies the implementation against known analytical solutions
for simple geometric shapes like squares.
"""

import numpy as np
import sys
import os
import pytest
from pathlib import Path

# Add the package directory to the Python path when running standalone
package_root = Path(__file__).parent.parent.parent
python_dir = package_root / "python"
sys.path.insert(0, str(python_dir))

# Now import from the package
from projection.api import Polygon, Region, ProjectionCalculator


def calculate_analytical_square_interior():
    """Calculate analytical value for average distance in unit square.
    
    Returns:
        Analytical average distance value
    """
    # For a unit square, the average distance has a known formula
    # (2 + sqrt(2) + 5*ln(1+sqrt(2)))/15
    return (2 + np.sqrt(2) + 5*np.log(1+np.sqrt(2)))/15


def calculate_analytical_square_boundary():
    """Calculate analytical value for average distance on unit square boundary.
    
    Returns:
        Analytical average boundary distance value
    """
    # For a unit square perimeter, the average distance is:
    # 1/4 + sqrt(2)/12 + 5/12 * ln(1 + sqrt(2))
    return 0.25 + np.sqrt(2)/12 + 5/12 * np.log(1+np.sqrt(2))


def test_analytical_comparison():
    """Run the comparison between analytical and Monte Carlo calculations."""
    # Set up test parameters
    square = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
    seed = 12345
    n_samples = 500000
    
    print("Comparing Analytical vs Monte Carlo Results")
    print("=" * 50)
    
    # 1. Unit Square Interior - Average Euclidean distance
    interior_calc = ProjectionCalculator(
        Polygon(square), [1, 0], 
        region=Region.INTERIOR, random_seed=seed
    )
    interior_mc = interior_calc.average_distance(exact=False, samples=n_samples)
    
    # Analytical formula for unit square interior
    interior_analytical = calculate_analytical_square_interior()
    
    interior_error = abs(interior_mc - interior_analytical) / interior_analytical
    
    # Assert error is within 0.1% for large sample size
    assert interior_error < 0.001, f"Interior avg distance error too large: {interior_error:.6f}"
    
    print("\nUnit Square Interior Average Distance:")
    print(f"Monte Carlo result ({n_samples:,} samples): {interior_mc:.6f}")
    print(f"Analytical formula: {interior_analytical:.6f}")
    print(f"Relative error: {interior_error:.6f} ({interior_error*100:.4f}%)")
    
    # 2. Unit Square Boundary - Average Euclidean distance
    boundary_calc = ProjectionCalculator(
        Polygon(square), [1, 0], 
        region=Region.BOUNDARY, random_seed=seed
    )
    boundary_mc = boundary_calc.average_distance(exact=False, samples=n_samples)
    
    # Analytical formula for unit square boundary
    boundary_analytical = calculate_analytical_square_boundary()
    
    boundary_error = abs(boundary_mc - boundary_analytical) / boundary_analytical
    
    # Assert error is within 0.1% for large sample size
    assert boundary_error < 0.001, f"Boundary avg distance error too large: {boundary_error:.6f}"
    
    print("\nUnit Square Boundary Average Distance:")
    print(f"Monte Carlo result ({n_samples:,} samples): {boundary_mc:.6f}")
    print(f"Analytical formula: {boundary_analytical:.6f}")
    print(f"Relative error: {boundary_error:.6f} ({boundary_error*100:.4f}%)")
    
    # 3. Unit Square Interior - Projected distance along x-axis
    interior_proj_x = interior_calc.projected_distance(exact=True)
    interior_proj_x_mc = interior_calc.projected_distance(exact=False, samples=n_samples)
    
    # For a unit square, projected distance along x-axis is 1/3
    interior_proj_x_analytical = 1/3
    
    interior_proj_error = abs(interior_proj_x - interior_proj_x_analytical) / interior_proj_x_analytical
    interior_proj_mc_error = abs(interior_proj_x_mc - interior_proj_x_analytical) / interior_proj_x_analytical
    
    # Assert exact calculation error is negligible (numerical precision only)
    assert interior_proj_error < 1e-10, f"Exact interior projection error too large: {interior_proj_error:.10f}"
    # Assert Monte Carlo error is within 0.1%
    assert interior_proj_mc_error < 0.001, f"MC interior projection error too large: {interior_proj_mc_error:.6f}"
    
    print("\nUnit Square Interior Projected Distance (x-axis):")
    print(f"Exact calculation: {interior_proj_x:.6f}")
    print(f"Monte Carlo result ({n_samples:,} samples): {interior_proj_x_mc:.6f}")
    print(f"Analytical formula: {interior_proj_x_analytical:.6f}")
    print(f"Exact relative error: {interior_proj_error:.6f} ({interior_proj_error*100:.4f}%)")
    print(f"Monte Carlo relative error: {interior_proj_mc_error:.6f} ({interior_proj_mc_error*100:.4f}%)")
    
    # 4. Unit Square Boundary - Projected distance along x-axis
    boundary_proj_x = boundary_calc.projected_distance(exact=True)
    boundary_proj_x_mc = boundary_calc.projected_distance(exact=False, samples=n_samples)
    
    # For a unit square boundary, projected distance along x-axis is 11/24
    boundary_proj_x_analytical = 11/24
    
    boundary_proj_error = abs(boundary_proj_x - boundary_proj_x_analytical) / boundary_proj_x_analytical
    boundary_proj_mc_error = abs(boundary_proj_x_mc - boundary_proj_x_analytical) / boundary_proj_x_analytical
    
    # Assert exact calculation error is negligible (numerical precision only)
    assert boundary_proj_error < 1e-10, f"Exact boundary projection error too large: {boundary_proj_error:.10f}"
    # Assert Monte Carlo error is within 0.1%
    assert boundary_proj_mc_error < 0.001, f"MC boundary projection error too large: {boundary_proj_mc_error:.6f}"
    
    print("\nUnit Square Boundary Projected Distance (x-axis):")
    print(f"Exact calculation: {boundary_proj_x:.6f}")
    print(f"Monte Carlo result ({n_samples:,} samples): {boundary_proj_x_mc:.6f}")
    print(f"Analytical formula: {boundary_proj_x_analytical:.6f}")
    print(f"Exact relative error: {boundary_proj_error:.6f} ({boundary_proj_error*100:.4f}%)")
    print(f"Monte Carlo relative error: {boundary_proj_mc_error:.6f} ({boundary_proj_mc_error*100:.4f}%)")


if __name__ == "__main__":
    test_analytical_comparison()