#!/usr/bin/env python3
"""
Test script for the fallback implementation of polygon-projection.

This script tests the functionality of the fallback implementation by deliberately 
disabling the C++ backend and running various operations with the Python fallback.
"""

import os
import sys
import numpy as np

# Add the parent directory to sys.path
parent_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, parent_dir)

# Force the use of fallback implementation by setting an environment variable
os.environ["POLYGON_PROJECTION_FORCE_FALLBACK"] = "1"

# Import the module
import polygon_projection as pp

def main():
    """Run a simple test of the fallback implementation."""
    print("Using C++ backend:", pp.using_cpp_backend())
    
    # Create a regular hexagon
    hexagon = pp.Polygon.regular(6, 1.0)
    print("\nCreated regular hexagon with 6 sides")
    print("Vertices:", hexagon.vertices)
    print("Area:", hexagon.area)
    print("Perimeter:", hexagon.perimeter)
    print("Is convex:", hexagon.is_convex())
    
    # Test projection calculator with interior region
    direction = [1.0, 0.5]
    print("\nTesting interior projections with direction", direction)
    interior_calc = pp.ProjectionCalculator(hexagon, direction, pp.Region.INTERIOR, random_seed=42)
    
    # Calculate exact and Monte Carlo distances
    exact_distance = interior_calc.projected_distance(exact=True)
    mc_distance = interior_calc.projected_distance(exact=False, samples=10000)
    print("Exact expected projected distance:", exact_distance)
    print("Monte Carlo expected projected distance:", mc_distance)
    print("Relative error:", abs(exact_distance - mc_distance) / exact_distance)
    
    # Test average distances
    exact_avg = interior_calc.average_distance(exact=True)
    mc_avg = interior_calc.average_distance(exact=False, samples=10000)
    print("\nExact average Euclidean distance:", exact_avg)
    print("Monte Carlo average Euclidean distance:", mc_avg)
    print("Relative error:", abs(exact_avg - mc_avg) / exact_avg)
    
    # Test boundary region
    print("\nTesting boundary projections with direction", direction)
    boundary_calc = pp.ProjectionCalculator(hexagon, direction, pp.Region.BOUNDARY, random_seed=42)
    
    # Calculate exact and Monte Carlo distances for boundary
    boundary_exact = boundary_calc.projected_distance(exact=True)
    boundary_mc = boundary_calc.projected_distance(exact=False, samples=10000)
    print("Exact expected projected distance:", boundary_exact)
    print("Monte Carlo expected projected distance:", boundary_mc)
    print("Relative error:", abs(boundary_exact - boundary_mc) / boundary_exact)
    
    # Test generator-based Monte Carlo
    print("\nTesting Monte Carlo generator")
    for i, result in enumerate(interior_calc.monte_carlo_analysis(samples=50000, batch_size=1000, yield_interval=10)):
        progress = result['progress']
        proj_dist = result['projected_distance']
        avg_dist = result['average_distance']
        
        # Print progress every 20%
        if i % 10 == 0:
            print(f"Progress: {progress:.1%}, Projected distance: {proj_dist:.6f}, Average distance: {avg_dist:.6f}")
    
    print("\nTest completed successfully!")

if __name__ == "__main__":
    main()