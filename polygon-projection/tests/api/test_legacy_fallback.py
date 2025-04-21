#!/usr/bin/env python3
"""
Test script for the fallback implementation using direct imports from legacy code.

This script tests the key functionality from the legacy Python code to ensure
it works correctly before integrating with the fallback implementation.
"""

import os
import sys
import numpy as np
from enum import Enum

# Add the legacy directory to sys.path
legacy_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'legacy')
sys.path.insert(0, legacy_dir)

# Import legacy modules directly
import polygon_projection as legacy_pp
import projection_sampling as legacy_ps
import edge_projection as legacy_ep
import piecewise_integration as legacy_pi

def main():
    """Run tests directly using the legacy code."""
    print("Testing legacy code directly...")
    
    # Create a regular hexagon
    hexagon = legacy_pp.generate_regular_polygon(n_sides=6)
    print("\nCreated regular hexagon with 6 sides")
    print("Vertices shape:", hexagon.shape)
    print("Area:", legacy_pp.calculate_polygon_area(hexagon))
    print("Perimeter:", legacy_pp.calculate_polygon_perimeter(hexagon))
    
    # Test interior projections
    direction = [1.0, 0.5]
    print("\nTesting interior projections with direction", direction)
    
    # Create a legacy projection calculator
    interior_proj = legacy_pp.PolygonProjection(
        hexagon, 
        direction, 
        region=legacy_pp.PolygonRegion.INTERIOR, 
        seed=42
    )
    
    # Test exact and Monte Carlo methods
    exact = interior_proj.exact_expected_distance()
    mc = interior_proj.monte_carlo_expected_distance(n_samples=10000)
    print("Exact expected projected distance:", exact)
    print("Monte Carlo expected projected distance:", mc)
    print("Relative error:", abs(exact - mc) / exact)
    
    # Test boundary projections
    print("\nTesting boundary projections with direction", direction)
    boundary_proj = legacy_pp.PolygonProjection(
        hexagon, 
        direction, 
        region=legacy_pp.PolygonRegion.BOUNDARY, 
        seed=42
    )
    
    # Test exact and Monte Carlo methods for boundary
    boundary_exact = boundary_proj.exact_expected_distance()
    boundary_mc = boundary_proj.monte_carlo_expected_distance(n_samples=10000)
    print("Exact expected projected distance:", boundary_exact)
    print("Monte Carlo expected projected distance:", boundary_mc)
    print("Relative error:", abs(boundary_exact - boundary_mc) / boundary_exact)
    
    # Test the generator
    print("\nTesting Monte Carlo generator")
    generator = interior_proj.monte_carlo_generator(n_samples=5000, batch_size=1000, yield_interval=1)
    
    for i, result in enumerate(generator):
        if i >= 5:  # Only show first 5 results
            break
        print(f"Batch {i+1}: Projected distance: {result['projected_distance']:.6f}")
    
    print("\nLegacy code tests completed successfully!")

if __name__ == "__main__":
    main()