#!/usr/bin/env python3
"""
Example script demonstrating how to use the Random Projections library from Python.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

# Add the Random Projections Python package to the path
sys.path.insert(0, '../../python')

# Import classes from the API
from projection import Polygon, Region, ProjectionCalculator

def main():
    # Create some polygons
    square = Polygon(np.array([[0, 0], [1, 0], [1, 1], [0, 1]]))
    hexagon = Polygon.regular(sides=6, radius=1.0)
    
    # Print information about the polygons
    print(f"Square: {len(square.vertices)} vertices, area={square.area:.2f}, perimeter={square.perimeter:.2f}")
    print(f"Hexagon: {len(hexagon.vertices)} vertices, area={hexagon.area:.2f}, perimeter={hexagon.perimeter:.2f}")
    
    # Calculate projection along x-axis for the square
    calc = ProjectionCalculator(square, [1, 0], region=Region.INTERIOR, random_seed=42)
    
    # Get the exact projected distance
    exact_distance = calc.projected_distance(exact=True)
    print(f"\nExact projected distance for square (interior): {exact_distance:.6f}")
    
    # Try Monte Carlo calculation
    mc_distance = calc.projected_distance(exact=False, samples=10000)
    print(f"Monte Carlo projected distance (10,000 samples): {mc_distance:.6f}")
    print(f"Relative error: {abs(exact_distance - mc_distance) / exact_distance:.6f}")
    
    # Sample some points
    interior_points = calc.sample_points(count=100)
    print(f"\nSampled {len(interior_points)} interior points")
    
    # Plot the polygons and sampled points
    plt.figure(figsize=(10, 5))
    
    # Plot square and interior points
    plt.subplot(1, 2, 1)
    vertices = np.vstack([square.vertices, square.vertices[0]])  # Close the polygon
    plt.plot(vertices[:, 0], vertices[:, 1], 'b-')
    plt.scatter(interior_points[:, 0], interior_points[:, 1], c='r', s=10, alpha=0.5)
    plt.title("Square with interior points")
    plt.axis('equal')
    plt.grid(True)
    
    # Plot hexagon
    plt.subplot(1, 2, 2)
    vertices = np.vstack([hexagon.vertices, hexagon.vertices[0]])  # Close the polygon
    plt.plot(vertices[:, 0], vertices[:, 1], 'g-')
    plt.title("Regular hexagon")
    plt.axis('equal')
    plt.grid(True)
    
    plt.tight_layout()
    plt.savefig("polygon_example.png")
    plt.show()

if __name__ == "__main__":
    main()