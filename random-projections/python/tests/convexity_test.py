#!/usr/bin/env python3
"""
Example demonstrating the improved convexity checking algorithm
for the polygon_projection library.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from pathlib import Path

# Add the package directory to the Python path when running standalone
package_root = Path(__file__).parent.parent.parent
python_dir = package_root / "python"
sys.path.insert(0, str(python_dir))

# Now import from the package
from projection.api import Polygon

def plot_polygon(ax, vertices, title, color='blue'):
    """Draw a polygon in the specified axis with proper formatting."""
    # Convert to numpy for easier manipulation
    vertices_np = np.array(vertices)
    
    # Close the polygon by appending the first vertex
    closed_vertices = np.vstack([vertices_np, vertices_np[0]])
    
    # Plot the polygon
    ax.plot(closed_vertices[:, 0], closed_vertices[:, 1], '-o', color=color)
    ax.set_title(title)
    ax.set_aspect('equal')
    ax.grid(True)

def main():
    """Main function to demonstrate polygon convexity testing."""
    # Define some test polygons
    square = [(0, 0), (2, 0), (2, 2), (0, 2)]
    concave_dent = [(-1, -1), (1, -1), (0, 0), (-1, 1)]
    spiral = [(0, 0), (2, 0), (2, 2), (1, 1), (0, 2)]  # self-intersecting
    hexagon = [(0, 0), (2, 0), (3, 1), (2, 2), (0, 2), (-1, 1)]  # convex

    polygons = [
        (square, "Square - Convex"),
        (concave_dent, "Concave Dent - Not Convex"),
        (spiral, "Spiral - Self-intersecting"),
        (hexagon, "Hexagon - Convex")
    ]

    # Create figure and subplots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    # Plot each polygon and test convexity
    for i, (vertices, title) in enumerate(polygons):
        try:
            # Create Polygon object to test convexity
            polygon = Polygon(vertices)
            is_convex = "Convex ✓"
            color = 'green'
        except ValueError as e:
            # If creation fails due to convexity check
            is_convex = f"Not Convex ✗ ({str(e)})"
            color = 'red'
        
        # Plot the polygon
        plot_polygon(axes[i], vertices, f"{title}\n{is_convex}", color)
    
    plt.tight_layout()
    plt.savefig("convexity_examples.png")
    plt.show()

if __name__ == "__main__":
    main()