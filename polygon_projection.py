#!/usr/bin/env python3
"""
polygon_projection.py

A unified high-level interface for computing projections of convex polygons.
This module combines functionality from:
1. projection_sampling.py - projections of points sampled from polygon interior
2. edge_projection.py - projections of points sampled from polygon boundary

The API provides standardized methods to:
- Calculate exact expected distances between projections
- Generate Monte Carlo estimates of expected distances
- Compare exact and Monte Carlo results
"""

import numpy as np
from enum import Enum
import projection_sampling as ps
import edge_projection as ep


class PolygonRegion(Enum):
    """Enum specifying which region of the polygon to use for projections."""
    INTERIOR = 'interior'
    BOUNDARY = 'boundary'


class PolygonProjection:
    """
    Unified interface for polygon projections and expected distance calculations.
    
    This class provides a high-level API for working with projections of points
    sampled from either the interior or boundary of a convex polygon onto a
    specified direction.
    """
    
    def __init__(self, vertices, direction, region=PolygonRegion.INTERIOR, seed=None):
        """
        Initialize a PolygonProjection instance.
        
        Parameters
        ----------
        vertices : array_like, shape (n,2)
            CCW-ordered polygon vertices.
        direction : array_like, shape (2,)
            Projection direction (need not be unit; will be normalized).
        region : PolygonRegion, optional
            Which region to sample from: INTERIOR or BOUNDARY.
            Default is INTERIOR.
        seed : int, optional
            Random seed for reproducibility.
        """
        self.vertices = np.asarray(vertices, dtype=float)
        self.direction = np.asarray(direction, dtype=float)
        self.region = region
        self.rng = np.random.default_rng(seed)
        
        # Pre-compute the density representation based on the region
        if region == PolygonRegion.INTERIOR:
            self.density = ps.projected_density(vertices, direction)
        else:
            self.steps = ep.project_edges(vertices, direction)
    
    def exact_expected_distance(self):
        """
        Calculate the exact expected distance between projected points.
        
        Returns
        -------
        float
            Expected absolute difference E[|X-Y|] between projected points.
        """
        if self.region == PolygonRegion.INTERIOR:
            return ps.expected_abs_diff(self.density)
        else:
            return ep.boundary_expected_abs_diff(self.steps)
    
    def monte_carlo_expected_distance(self, n_samples=10000):
        """
        Estimate expected distance using Monte Carlo sampling.
        
        Parameters
        ----------
        n_samples : int, optional
            Number of samples to use in the estimate. Default is 10000.
            
        Returns
        -------
        float
            Monte Carlo estimate of E[|X-Y|].
        """
        if self.region == PolygonRegion.INTERIOR:
            distances = ps.sample_mc_pairwise_projected_distances(
                self.vertices, self.direction, n=n_samples, rng=self.rng)
        else:
            distances = ep.sample_mc_pairwise_projected_distances(
                self.vertices, self.direction, n=n_samples, rng=self.rng)
            
        return float(np.mean(distances))
    
    def compare_methods(self, n_samples=10000):
        """
        Compare exact and Monte Carlo methods for expected distance.
        
        Parameters
        ----------
        n_samples : int, optional
            Number of samples for Monte Carlo estimate. Default is 10000.
            
        Returns
        -------
        dict
            Dictionary containing exact value, Monte Carlo estimate, 
            absolute error, and relative error.
        """
        exact = self.exact_expected_distance()
        mc = self.monte_carlo_expected_distance(n_samples)
        abs_error = abs(exact - mc)
        rel_error = abs_error / exact if exact != 0 else float('inf')
        
        return {
            'region': self.region.value,
            'exact': exact,
            'monte_carlo': mc,
            'absolute_error': abs_error,
            'relative_error': rel_error,
            'n_samples': n_samples
        }
    
    def sample_point(self):
        """
        Sample a single point from the specified region of the polygon.
        
        Returns
        -------
        ndarray, shape (2,)
            Sampled point coordinates.
        """
        if self.region == PolygonRegion.INTERIOR:
            return ps.sample_point_in_convex_polygon(self.vertices, self.rng)
        else:
            return ep.sample_point_on_boundary(self.vertices, self.rng)
    
    def sample_point_projection(self):
        """
        Sample a point and return its projection onto the direction.
        
        Returns
        -------
        float
            Scalar projection of the sampled point onto direction.
        """
        point = self.sample_point()
        return ps.project_point_to_direction(point, self.direction)


def generate_regular_polygon(n_sides=6, radius=1.0):
    """
    Generate a regular polygon with n_sides.
    
    Parameters
    ----------
    n_sides : int, optional
        Number of sides in the polygon. Default is 6.
    radius : float, optional
        Radius of the circumscribed circle. Default is 1.0.
        
    Returns
    -------
    ndarray, shape (n_sides, 2)
        Vertices of the regular polygon.
    """
    angles = np.linspace(0, 2*np.pi, n_sides, endpoint=False)
    return radius * np.column_stack([np.cos(angles), np.sin(angles)])


if __name__ == "__main__":
    # Example usage demonstrating the unified interface
    hexagon = generate_regular_polygon(n_sides=6)
    direction = [1.0, 0.5]
    
    # Compare interior vs boundary projections
    results = []
    for region in [PolygonRegion.INTERIOR, PolygonRegion.BOUNDARY]:
        proj = PolygonProjection(hexagon, direction, region=region, seed=42)
        results.append(proj.compare_methods(n_samples=100000))
    
    # Print results 
    print("\nExpected distance comparisons:")
    print("-" * 60)
    print(f"{'Region':<10} {'Exact':<15} {'Monte Carlo':<15} {'Rel Error':<15}")
    print("-" * 60)
    for r in results:
        print(f"{r['region']:<10} {r['exact']:<15.6f} {r['monte_carlo']:<15.6f} {r['relative_error']:<15.6e}")