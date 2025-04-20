#!/usr/bin/env python3
"""
Simple Monte Carlo test for average distance along a unit square boundary.
No fancy math - just sampling points and computing distances.
"""

import numpy as np
import time

# Set seed for reproducibility
np.random.seed(42)

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

def main():
    # Number of sample pairs
    n_samples = 1_000_000
    
    print(f"Simple Boundary Monte Carlo Test: {n_samples:,} samples")
    print("=" * 50)
    
    # Time the computation
    start_time = time.time()
    
    # Sample point pairs and compute distances
    distances = np.zeros(n_samples)
    for i in range(n_samples):
        p1 = sample_square_boundary()
        p2 = sample_square_boundary()
        distances[i] = np.linalg.norm(p1 - p2)
    
    # Calculate average and statistics
    average_distance = np.mean(distances)
    std_dev = np.std(distances)
    max_dist = np.max(distances)
    min_dist = np.min(distances)
    
    end_time = time.time()
    elapsed = end_time - start_time
    
    # Print results
    print(f"Average distance between boundary points: {average_distance:.6f}")
    print(f"Standard deviation: {std_dev:.6f}")
    print(f"Minimum distance: {min_dist:.6f}")
    print(f"Maximum distance: {max_dist:.6f}")
    print(f"Computation time: {elapsed:.2f} seconds")
    
    # Create histogram for visualization
    hist, bins = np.histogram(distances, bins=50, density=True)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    
    print("\nHistogram data:")
    print("Distance | Frequency")
    print("-" * 25)
    for center, freq in zip(bin_centers, hist):
        print(f"{center:.4f}   | {freq:.4f}")

if __name__ == "__main__":
    main()