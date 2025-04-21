#!/usr/bin/env python3
"""
Test script for the polygon-projection fallback API.

This tests the actual API implementation using the legacy code as a backend.
"""

import os
import sys
import numpy as np

# Set environment variable to force fallback
os.environ["POLYGON_PROJECTION_FORCE_FALLBACK"] = "1"

# Add the polygon-projection module to the path
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'polygon-projection', 'python'))

# Import polygon_projection module - should use the fallback implementation
try:
    import polygon_projection as pp
    print("Using C++ backend:", pp.using_cpp_backend())
except ImportError as e:
    print(f"Failed to import polygon_projection: {e}")
    sys.exit(1)

def main():
    """Test the fallback API implementation."""
    print("\nTesting Polygon Projection Fallback API")
    
    # Create a regular hexagon
    print("\nCreating a regular hexagon...")
    try:
        hexagon = pp.Polygon.regular(6, 1.0)
        print("Created hexagon with vertices shape:", hexagon.vertices.shape)
        print("Area:", hexagon.area)
        print("Perimeter:", hexagon.perimeter)
        print("Is convex:", hexagon.is_convex())
    except Exception as e:
        print(f"Error creating polygon: {e}")
        sys.exit(1)
    
    # Test projection calculator
    print("\nTesting ProjectionCalculator...")
    try:
        direction = [1.0, 0.5]
        calculator = pp.ProjectionCalculator(hexagon, direction, pp.Region.INTERIOR, random_seed=42)
        
        # Calculate projected distances
        exact = calculator.projected_distance(exact=True)
        mc = calculator.projected_distance(exact=False, samples=5000)
        print(f"Exact projected distance: {exact}")
        print(f"Monte Carlo projected distance: {mc}")
        print(f"Relative error: {abs(exact - mc) / exact}")
        
        # Test the monte_carlo_analysis generator
        print("\nTesting monte_carlo_analysis generator...")
        for i, result in enumerate(calculator.monte_carlo_analysis(samples=2000, batch_size=500, yield_interval=1)):
            print(f"Batch {i+1}: Progress {result['progress']:.1%}, Projected distance: {result['projected_distance']:.6f}")
            if i >= 3:  # Only show 4 batches
                break
    except Exception as e:
        print(f"Error in projection calculations: {e}")
        sys.exit(1)
    
    print("\nFallback API test completed successfully!")

if __name__ == "__main__":
    main()