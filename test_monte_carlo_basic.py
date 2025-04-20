#!/usr/bin/env python3
"""
test_monte_carlo_basic.py

Basic tests for the Monte Carlo distance calculations within a square.
"""

import numpy as np

# Set seed for reproducibility
np.random.seed(12345)

def sample_square_interior():
    """Sample a point uniformly from the unit square interior."""
    return np.random.random(2)

def sample_square_boundary():
    """Sample a point uniformly from the unit square boundary."""
    # Choose side (0=bottom, 1=right, 2=top, 3=left)
    side = np.random.randint(0, 4)
    # Position along the side
    t = np.random.random()
    
    if side == 0:   # Bottom
        return np.array([t, 0])
    elif side == 1: # Right
        return np.array([1, t])
    elif side == 2: # Top
        return np.array([1-t, 1])
    else:           # Left
        return np.array([0, 1-t])

def monte_carlo_average_distance(n_samples, sampling_func):
    """Calculate average Euclidean distance using Monte Carlo."""
    distances = []
    for _ in range(n_samples):
        p1 = sampling_func()
        p2 = sampling_func()
        distances.append(np.linalg.norm(p1 - p2))
    return np.mean(distances)

def main():
    # Sample size
    n_samples = 1000000
    
    # 1. Unit square interior
    interior_distance = monte_carlo_average_distance(n_samples, sample_square_interior)
    # Known analytical value: (2 + 2*sqrt(2) + 5*ln(1+sqrt(2)))/15 ≈ 0.521405
    analytical_interior = (2 + 2*np.sqrt(2) + 5*np.log(1+np.sqrt(2)))/15
    
    # 2. Unit square boundary
    boundary_distance = monte_carlo_average_distance(n_samples, sample_square_boundary)
    # Literature value: 4/3 ≈ 1.333333
    analytical_boundary = 4/3
    
    # Print results
    print("Unit Square Monte Carlo Test Results")
    print("=" * 40)
    print(f"Interior average distance:")
    print(f"  Monte Carlo result: {interior_distance:.6f}")
    print(f"  Analytical value:   {analytical_interior:.6f}")
    print(f"  Relative error:     {abs(interior_distance - analytical_interior)/analytical_interior:.6f}")
    print()
    
    print(f"Boundary average distance:")
    print(f"  Monte Carlo result: {boundary_distance:.6f}")
    print(f"  Analytical value:   {analytical_boundary:.6f}")
    print(f"  Relative error:     {abs(boundary_distance - analytical_boundary)/analytical_boundary:.6f}")

if __name__ == "__main__":
    main()