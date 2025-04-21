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
from . import projection_sampling as ps
from . import edge_projection as ep
from . import piecewise_integration as pi


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
        
    def monte_carlo_generator(self, n_samples=10000, batch_size=1000, yield_interval=10):
        """
        Generator that performs Monte Carlo sampling and yields intermediate results.
        
        This function samples points from the polygon (interior or boundary based on
        the region setting), computes projected distances, and periodically yields
        the current statistics, allowing for real-time UI updates.
        
        Parameters
        ----------
        n_samples : int, optional
            Total number of samples to generate. Default is 10000.
        batch_size : int, optional
            Size of batches to process. Default is 1000.
        yield_interval : int, optional
            Number of batches after which to yield results. Default is 10.
            
        Yields
        ------
        dict
            Dictionary containing:
            - 'samples_processed': Number of samples processed so far
            - 'total_samples': Total number of samples to process
            - 'projected_distance': Current mean of projected distances
            - 'average_distance': Current mean of Euclidean distances
        """
        # Setup for storing accumulated results
        projected_distances = []
        euclidean_distances = []
        samples_processed = 0
        batches_since_yield = 0
        
        # Process in batches to allow periodic updates
        for i in range(0, n_samples, batch_size):
            # Determine the actual batch size (might be smaller for the last batch)
            current_batch_size = min(batch_size, n_samples - i)
            
            # Generate sample points
            if self.region == PolygonRegion.INTERIOR:
                # Sample from interior
                points1 = [ps.sample_point_in_convex_polygon(self.vertices, self.rng) 
                           for _ in range(current_batch_size)]
                points2 = [ps.sample_point_in_convex_polygon(self.vertices, self.rng) 
                           for _ in range(current_batch_size)]
            else:
                # Sample from boundary
                points1 = [ep.sample_point_on_boundary(self.vertices, self.rng) 
                           for _ in range(current_batch_size)]
                points2 = [ep.sample_point_on_boundary(self.vertices, self.rng) 
                           for _ in range(current_batch_size)]
            
            # Calculate projections and distances for each pair
            for p1, p2 in zip(points1, points2):
                # Projected distance
                proj1 = ps.project_point_to_direction(p1, self.direction)
                proj2 = ps.project_point_to_direction(p2, self.direction)
                projected_distances.append(abs(proj1 - proj2))
                
                # Euclidean distance
                euclidean_distances.append(np.linalg.norm(p1 - p2))
            
            # Update the number of processed samples
            samples_processed += current_batch_size
            batches_since_yield += 1
            
            # Yield intermediate results periodically
            if batches_since_yield >= yield_interval or samples_processed == n_samples:
                yield {
                    'samples_processed': samples_processed,
                    'total_samples': n_samples,
                    'projected_distance': float(np.mean(projected_distances)),
                    'average_distance': float(np.mean(euclidean_distances))
                }
                batches_since_yield = 0
    
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
    
    def get_area(self):
        """
        Calculate the area of the polygon.
        
        Returns
        -------
        float
            Area of the polygon.
        """
        return calculate_polygon_area(self.vertices)
    
    def get_perimeter(self):
        """
        Calculate the perimeter of the polygon.
        
        Returns
        -------
        float
            Perimeter of the polygon.
        """
        return calculate_polygon_perimeter(self.vertices)
    
    def average_distance_exact(self, integration_order=16):
        """
        Calculate the exact average Euclidean distance E[||X-Y||] between two
        random points in the specified region (interior or boundary).
        
        This is obtained by integrating E[|X-Y|_θ] over all directions θ in [0,π]
        and dividing by 2, where E[|X-Y|_θ] is the expected absolute difference
        of projections onto direction θ.
        
        Parameters
        ----------
        integration_order : int, optional
            Order of Gauss-Legendre quadrature to use per smooth interval.
            Higher values give more accurate results. Default is 16.
            
        Returns
        -------
        float
            Expected Euclidean distance E[||X-Y||].
        """
        # Define the integrand function
        def projected_mean_abs(vertices, theta):
            direction = np.array([np.cos(theta), np.sin(theta)])
            if self.region == PolygonRegion.INTERIOR:
                pieces = ps.projected_density(vertices, direction)
                return ps.expected_abs_diff(pieces)
            else:
                steps = ep.project_edges(vertices, direction)
                return ep.boundary_expected_abs_diff(steps)
        
        # Compute the integral and divide by 2
        integral = pi.expected_distance_interior_exact(
            self.vertices, projected_mean_abs, order=integration_order
        )
        return 0.5 * integral
    
    def average_distance_monte_carlo(self, n_samples=10000):
        """
        Estimate the average Euclidean distance E[||X-Y||] between two random
        points in the specified region (interior or boundary) using Monte Carlo.
        
        Parameters
        ----------
        n_samples : int, optional
            Number of point pairs to sample. Default is 10000.
            
        Returns
        -------
        float
            Monte Carlo estimate of E[||X-Y||].
        """
        distances = []
        
        for _ in range(n_samples):
            if self.region == PolygonRegion.INTERIOR:
                p1 = ps.sample_point_in_convex_polygon(self.vertices, self.rng)
                p2 = ps.sample_point_in_convex_polygon(self.vertices, self.rng)
            else:
                p1 = ep.sample_point_on_boundary(self.vertices, self.rng)
                p2 = ep.sample_point_on_boundary(self.vertices, self.rng)
            
            distances.append(np.linalg.norm(p1 - p2))
        
        return float(np.mean(distances))


def calculate_polygon_area(vertices):
    """
    Calculate the area of a convex polygon using the Shoelace formula.
    
    Parameters
    ----------
    vertices : array_like, shape (n,2)
        CCW-ordered polygon vertices.
        
    Returns
    -------
    float
        Area of the polygon.
    """
    vertices = np.asarray(vertices, dtype=float)
    x = vertices[:, 0]
    y = vertices[:, 1]
    
    # Shoelace formula
    area = 0.5 * np.abs(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))
    return area

def calculate_polygon_perimeter(vertices):
    """
    Calculate the perimeter (sum of edge lengths) of a polygon.
    
    Parameters
    ----------
    vertices : array_like, shape (n,2)
        CCW-ordered polygon vertices.
        
    Returns
    -------
    float
        Perimeter of the polygon.
    """
    vertices = np.asarray(vertices, dtype=float)
    n = len(vertices)
    
    # Calculate the sum of distances between consecutive vertices
    perimeter = 0.0
    for i in range(n):
        p1 = vertices[i]
        p2 = vertices[(i+1) % n]
        perimeter += np.linalg.norm(p2 - p1)
        
    return perimeter

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


def is_convex(vertices):
    """
    Check if a polygon is convex by verifying that all interior angles are less than 180 degrees.
    Uses the cross product to check if all vertices are making "right turns" or all making "left turns".
    
    Parameters
    ----------
    vertices : array_like, shape (n,2)
        CCW-ordered polygon vertices.
        
    Returns
    -------
    bool
        True if the polygon is convex, False otherwise.
    """
    # Need at least 3 vertices to form a polygon
    if len(vertices) < 3:
        return True  # Degenerate case
        
    # Compute the cross product for each set of three consecutive vertices
    n = len(vertices)
    # Initialize sign of first cross product
    dx1 = vertices[1][0] - vertices[0][0]
    dy1 = vertices[1][1] - vertices[0][1]
    dx2 = vertices[2][0] - vertices[1][0]
    dy2 = vertices[2][1] - vertices[1][1]
    cross_z = dx1 * dy2 - dy1 * dx2
    sign = 1 if cross_z > 0 else -1 if cross_z < 0 else 0
    
    # Check all other cross products
    for i in range(1, n):
        dx1 = vertices[(i+1) % n][0] - vertices[i][0]
        dy1 = vertices[(i+1) % n][1] - vertices[i][1]
        dx2 = vertices[(i+2) % n][0] - vertices[(i+1) % n][0]
        dy2 = vertices[(i+2) % n][1] - vertices[(i+1) % n][1]
        cross_z = dx1 * dy2 - dy1 * dx2
        
        # If sign changes, the polygon is not convex
        curr_sign = 1 if cross_z > 0 else -1 if cross_z < 0 else 0
        if curr_sign != 0 and sign != 0 and curr_sign != sign:
            return False
        # Update sign if it was previously 0
        if sign == 0 and curr_sign != 0:
            sign = curr_sign
    
    return True


if __name__ == "__main__":
    # Example usage demonstrating the unified interface
    hexagon = generate_regular_polygon(n_sides=6)
    direction = [1.0, 0.5]
    
    # Compare interior vs boundary projections for a specific direction
    results = []
    for region in [PolygonRegion.INTERIOR, PolygonRegion.BOUNDARY]:
        proj = PolygonProjection(hexagon, direction, region=region, seed=42)
        results.append(proj.compare_methods(n_samples=100000))
    
    # Print projection results for the specific direction
    print("\nExpected projected distance comparisons (direction-specific):")
    print("-" * 60)
    print(f"{'Region':<10} {'Exact':<15} {'Monte Carlo':<15} {'Rel Error':<15}")
    print("-" * 60)
    for r in results:
        print(f"{r['region']:<10} {r['exact']:<15.6f} {r['monte_carlo']:<15.6f} {r['relative_error']:<15.6e}")
    
    # Calculate and compare average Euclidean distances (over all directions)
    print("\nAverage Euclidean distance results:")
    print("-" * 60)
    print(f"{'Region':<10} {'Exact':<15} {'Monte Carlo':<15} {'Rel Error':<15}")
    print("-" * 60)
    
    for region in [PolygonRegion.INTERIOR, PolygonRegion.BOUNDARY]:
        proj = PolygonProjection(hexagon, direction, region=region, seed=42)
        
        # Calculate exact average distance
        exact_avg = proj.average_distance_exact(integration_order=20)
        
        # Estimate with Monte Carlo
        mc_avg = proj.average_distance_monte_carlo(n_samples=50000)
        
        # Calculate relative error
        rel_error = abs(exact_avg - mc_avg) / exact_avg
        
        print(f"{region.value:<10} {exact_avg:<15.6f} {mc_avg:<15.6f} {rel_error:<15.6e}")