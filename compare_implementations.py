#!/usr/bin/env python3
"""
Compare our Monte Carlo implementation in PolygonProjection with the simple direct approach.
"""

import numpy as np
from polygon_projection import PolygonProjection, PolygonRegion

def sample_unit_square():
    """Sample a point uniformly from unit square: [0,1] × [0,1]"""
    return np.random.random(2)

def sample_square_boundary():
    """Sample a point uniformly from unit square boundary."""
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

def direct_monte_carlo(n_samples, sampler):
    """Direct Monte Carlo calculation."""
    distances = []
    for _ in range(n_samples):
        p1 = sampler()
        p2 = sampler()
        distances.append(np.linalg.norm(p1 - p2))
    return np.mean(distances)

def main():
    # Test parameters
    n_samples = 100000
    seed = 12345
    square = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
    
    print("Comparing Monte Carlo Implementations")
    print("=" * 50)
    
    # Set seeds for reproducibility
    np.random.seed(seed)
    
    # 1. Unit Square Interior
    print("\nUnit Square Interior:")
    
    # Our implementation
    proj_interior = PolygonProjection(square, [1, 0], region=PolygonRegion.INTERIOR, seed=seed)
    our_result = proj_interior.average_distance_monte_carlo(n_samples=n_samples)
    print(f"Our implementation: {our_result:.6f}")
    
    # Direct implementation
    np.random.seed(seed)  # Reset seed
    direct_result = direct_monte_carlo(n_samples, sample_unit_square)
    print(f"Direct computation: {direct_result:.6f}")
    
    # 2. Unit Square Boundary
    print("\nUnit Square Boundary:")
    
    # Our implementation
    proj_boundary = PolygonProjection(square, [1, 0], region=PolygonRegion.BOUNDARY, seed=seed)
    our_result = proj_boundary.average_distance_monte_carlo(n_samples=n_samples)
    print(f"Our implementation: {our_result:.6f}")
    
    # Direct implementation
    np.random.seed(seed)  # Reset seed
    direct_result = direct_monte_carlo(n_samples, sample_square_boundary)
    print(f"Direct computation: {direct_result:.6f}")

if __name__ == "__main__":
    main()