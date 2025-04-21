#!/usr/bin/env python3
"""
fallback.py

Fallback implementation for the polygon-projection API using the legacy Python code.
This module provides Python implementations that match the C++ API but using the
legacy Python code as a backend.
"""

import numpy as np
from enum import Enum
from typing import List, Tuple, Dict, Generator, Optional, Union, Callable, Any

# Import the legacy modules from the local package
from .legacy import (
    PolygonRegion as LegacyPolygonRegion,
    PolygonProjection,
    calculate_polygon_area,
    calculate_polygon_perimeter,
    generate_regular_polygon,
    project_point_to_direction,
    is_convex as _is_convex  # Import as private
)

# Create Region enum to match the C++ API
class Region(Enum):
    """Enum defining the region to use for projections."""
    INTERIOR = 'interior'
    BOUNDARY = 'boundary'

# Create mapping from Region enum to legacy PolygonRegion
def _to_legacy_region(region):
    """Convert Region enum to legacy PolygonRegion enum."""
    if region == Region.INTERIOR:
        return LegacyPolygonRegion.INTERIOR
    else:
        return LegacyPolygonRegion.BOUNDARY

class Polygon:
    """Represents a convex polygon with various operations.
    
    This is a fallback implementation that uses the legacy Python code.
    """
    
    def __init__(self, vertices):
        """Create a polygon from vertex array.
        
        Args:
            vertices: Array of shape (n, 2) containing CCW-ordered vertices
        """
        self._vertices = np.asarray(vertices, dtype=float)
        if len(self._vertices) < 3:
            raise ValueError("Polygon must have at least 3 vertices")
        if not self.is_convex():
            raise ValueError("Polygon must be convex")
    
    @classmethod
    def from_vertices(cls, vertices):
        """Create polygon from vertex array.
        
        Args:
            vertices: Array of shape (n, 2) containing CCW-ordered vertices
            
        Returns:
            New polygon instance
        """
        return cls(vertices)
    
    @classmethod
    def regular(cls, sides, radius=1.0):
        """Create regular polygon with specified sides and radius.
        
        Args:
            sides: Number of sides (must be >= 3)
            radius: Radius of the circumscribed circle
            
        Returns:
            Regular polygon with specified parameters
        """
        if sides < 3:
            raise ValueError("Number of sides must be at least 3")
        if radius <= 0:
            raise ValueError("Radius must be positive")
        
        # Use the legacy implementation to generate a regular polygon
        vertices = generate_regular_polygon(n_sides=sides, radius=radius)
        return cls(vertices)
    
    @property
    def vertices(self):
        """Get polygon vertices."""
        return self._vertices
    
    @property
    def area(self):
        """Get polygon area."""
        return calculate_polygon_area(self._vertices)
    
    @property
    def perimeter(self):
        """Get polygon perimeter."""
        return calculate_polygon_perimeter(self._vertices)
    
    @staticmethod
    def check_convex(vertices):
        """Static method to check if a set of vertices forms a convex polygon.
        
        Uses a simple algorithm based on cross products to check convexity.
        For a convex polygon, all cross products should have the same sign.
        
        Args:
            vertices: Array of shape (n, 2) containing CCW-ordered vertices
            
        Returns:
            bool: True if the polygon is convex, False otherwise.
        """
        # Use the imported is_convex function from the legacy module
        return _is_convex(vertices)
    
    def is_convex(self):
        """Check if polygon is convex.
        
        Uses a simple algorithm based on cross products to check convexity.
        For a convex polygon, all cross products should have the same sign.
        
        Returns:
            bool: True if the polygon is convex, False otherwise.
        """
        # Use the imported is_convex function from the legacy module
        return _is_convex(self._vertices)
    
    def contains(self, point):
        """Check if polygon contains a point.
        
        Args:
            point: Point coordinates [x, y]
            
        Returns:
            True if the point is inside the polygon
        """
        point = np.asarray(point, dtype=float)
        n = len(self._vertices)
        
        # Check if point is inside the polygon using ray casting algorithm
        inside = False
        for i in range(n):
            j = (i + 1) % n
            vi = self._vertices[i]
            vj = self._vertices[j]
            
            # Check if point is on a vertex
            if np.array_equal(point, vi) or np.array_equal(point, vj):
                return True
            
            # Check if point is inside the polygon
            if ((vi[1] > point[1]) != (vj[1] > point[1])) and \
               (point[0] < (vj[0] - vi[0]) * (point[1] - vi[1]) / (vj[1] - vi[1]) + vi[0]):
                inside = not inside
        
        return inside


class ProjectionCalculator:
    """Handles projection calculations for a polygon.
    
    This is a fallback implementation that uses the legacy Python code.
    """
    
    def __init__(self, polygon, direction, region=Region.INTERIOR, random_seed=None):
        """Initialize projection calculator.
        
        Args:
            polygon: The polygon to analyze
            direction: Direction vector for projections (2D)
            region: Region to use (INTERIOR or BOUNDARY)
            random_seed: Optional seed for random number generation
        """
        self.polygon = polygon
        self.direction = np.asarray(direction, dtype=float)
        self.region = region
        self.random_seed = 0 if random_seed is None else random_seed
        
        # Create legacy implementation
        legacy_region = _to_legacy_region(region)
        self._impl = PolygonProjection(
            polygon.vertices,
            direction,
            region=legacy_region,
            seed=random_seed
        )
    
    def projected_distance(self, exact=True, samples=10000):
        """Calculate projected distance (exact or Monte Carlo).
        
        Args:
            exact: Whether to use exact calculation (True) or Monte Carlo (False)
            samples: Number of samples for Monte Carlo (ignored if exact=True)
            
        Returns:
            Expected absolute difference between projected points
        """
        if exact:
            return self._impl.exact_expected_distance()
        else:
            return self._impl.monte_carlo_expected_distance(n_samples=samples)
    
    def average_distance(self, exact=True, samples=10000, integration_order=16):
        """Calculate average distance over all directions.
        
        Args:
            exact: Whether to use exact calculation (True) or Monte Carlo (False)
            samples: Number of samples for Monte Carlo (ignored if exact=True)
            integration_order: Order of Gauss-Legendre quadrature for exact method
            
        Returns:
            Expected Euclidean distance between points
        """
        if exact:
            return self._impl.average_distance_exact(integration_order=integration_order)
        else:
            return self._impl.average_distance_monte_carlo(n_samples=samples)
    
    def monte_carlo_analysis(self, samples=10000, batch_size=1000, yield_interval=10):
        """Generator for progressive Monte Carlo calculations.
        
        Args:
            samples: Total number of samples to use
            batch_size: Size of each batch of samples
            yield_interval: Number of batches after which to yield results
            
        Yields:
            Dict containing current results and progress information
        """
        # Use the legacy generator for Monte Carlo analysis
        generator = self._impl.monte_carlo_generator(
            n_samples=samples,
            batch_size=batch_size,
            yield_interval=yield_interval
        )
        
        for result in generator:
            # Add progress information to match the API
            result['progress'] = result['samples_processed'] / result['total_samples']
            yield result
    
    def sample_points(self, count=1):
        """Sample points from the specified region.
        
        Args:
            count: Number of points to sample
            
        Returns:
            Array of shape (count, 2) containing sampled points
        """
        if count <= 0:
            raise ValueError("Count must be positive")
        
        points = np.empty((count, 2))
        for i in range(count):
            points[i] = self._impl.sample_point()
        
        return points
    
    def project_points(self, points):
        """Project points onto the direction vector.
        
        Args:
            points: Array of shape (n, 2) containing points to project
            
        Returns:
            Array of shape (n,) containing projections
        """
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 2:
            raise ValueError("Points must be a 2D array with shape (n, 2)")
        
        # Project each point
        projections = np.empty(len(points))
        for i, point in enumerate(points):
            # Use the imported project_point_to_direction function
            projections[i] = project_point_to_direction(point, self.direction)
        
        return projections
    
    def compare_methods(self, samples=10000):
        """Compare exact and Monte Carlo methods.
        
        Args:
            samples: Number of samples for Monte Carlo estimation
            
        Returns:
            Dict containing comparison results
        """
        # Use the legacy compare_methods implementation
        result = self._impl.compare_methods(n_samples=samples)
        
        # Modify keys to match the API
        result['region'] = self.region.name
        
        return result

# Additional utility functions

def is_convex(vertices):
    """
    Check if a polygon is convex by verifying that all interior angles are less than 180 degrees.
    Uses the cross product to check if all vertices are making "right turns" or all making "left turns".
    
    Args:
        vertices: Array of shape (n, 2) containing CCW-ordered vertices
        
    Returns:
        bool: True if the polygon is convex, False otherwise.
    """
    return _is_convex(vertices)

def create_polygon(vertices):
    """Create a polygon from vertices."""
    return Polygon(vertices)

def create_regular_polygon(sides, radius=1.0):
    """Create a regular polygon with specified sides and radius."""
    return Polygon.regular(sides, radius)

def create_projection_calculator(polygon, direction, region=Region.INTERIOR, random_seed=None):
    """Create a projection calculator."""
    return ProjectionCalculator(polygon, direction, region, random_seed)