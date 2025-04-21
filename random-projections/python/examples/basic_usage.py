#!/usr/bin/env python3
"""
Basic usage example for the projection library.

This example creates a polygon, projects points to a direction,
and compares exact and Monte Carlo calculation methods.
"""

import numpy as np
import matplotlib.pyplot as plt
from projection import Polygon, Region, ProjectionCalculator

def main():
    # Create a regular hexagon with radius 2
    hexagon = Polygon.regular(sides=6, radius=2.0)
    
    # Define a projection direction
    direction = np.array([1.0, 0.5])
    
    # Create projection calculators for interior and boundary
    interior_calc = ProjectionCalculator(
        polygon=hexagon,
        direction=direction,
        region=Region.INTERIOR,
        random_seed=42  # For reproducibility
    )
    
    boundary_calc = ProjectionCalculator(
        polygon=hexagon,
        direction=direction,
        region=Region.BOUNDARY,
        random_seed=42  # For reproducibility
    )
    
    # Calculate exact results
    interior_exact = interior_calc.projected_distance(exact=True)
    boundary_exact = boundary_calc.projected_distance(exact=True)
    
    print(f"Exact expected distances:")
    print(f"  Interior: {interior_exact:.6f}")
    print(f"  Boundary: {boundary_exact:.6f}")
    
    # Calculate Monte Carlo results
    samples = 100000
    interior_mc = interior_calc.projected_distance(exact=False, samples=samples)
    boundary_mc = boundary_calc.projected_distance(exact=False, samples=samples)
    
    print(f"\nMonte Carlo estimates ({samples:,} samples):")
    print(f"  Interior: {interior_mc:.6f}")
    print(f"  Boundary: {boundary_mc:.6f}")
    
    # Calculate errors
    interior_error = abs(interior_exact - interior_mc) / interior_exact
    boundary_error = abs(boundary_exact - boundary_mc) / boundary_exact
    
    print(f"\nRelative errors:")
    print(f"  Interior: {interior_error:.6e}")
    print(f"  Boundary: {boundary_error:.6e}")
    
    # Calculate average distances (over all directions)
    interior_avg = interior_calc.average_distance(exact=True)
    boundary_avg = boundary_calc.average_distance(exact=True)
    
    print(f"\nAverage Euclidean distances:")
    print(f"  Interior: {interior_avg:.6f}")
    print(f"  Boundary: {boundary_avg:.6f}")
    
    # Plot the polygon and direction
    plt.figure(figsize=(8, 8))
    
    # Plot polygon
    vertices = hexagon.vertices
    plt.fill(vertices[:, 0], vertices[:, 1], alpha=0.2, edgecolor='blue', facecolor='blue')
    
    # Plot direction vector
    norm = np.linalg.norm(direction)
    if norm > 0:
        unit_dir = direction / norm
        plt.arrow(0, 0, 3*unit_dir[0], 3*unit_dir[1], 
                 head_width=0.2, head_length=0.3, fc='red', ec='red')
    
    # Plot settings
    plt.grid(True)
    plt.axis('equal')
    plt.xlim(-3, 3)
    plt.ylim(-3, 3)
    plt.title('Polygon and Projection Direction')
    
    plt.show()

if __name__ == "__main__":
    main()