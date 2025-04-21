"""High-level Python API for polygon projection operations.

This module provides a Pythonic interface to the C++ backend.
No Python fallbacks are provided - the C++ backend is required.
"""

import os
import sys
import importlib.util
import numpy as np
from enum import Enum
from typing import List, Tuple, Dict, Generator, Optional, Union, Callable

# Import the C++ extension module
try:
    # Determine the path to the module
    module_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Find extension file by trying common suffixes
    candidates = ["_core.cpython-{}{}-{}.so", "_core.so", "_core.pyd"]
    ext_file = None
    
    for pattern in candidates:
        if sys.version_info.major == 3:
            filename = pattern.format(
                sys.version_info.major,
                sys.version_info.minor,
                sys.platform
            )
        else:
            filename = "_core.so"  # Fallback
        
        filepath = os.path.join(module_dir, filename)
        if os.path.exists(filepath):
            ext_file = filepath
            break
    
    if ext_file is None:
        raise ImportError(f"Could not find _core extension module in {module_dir}")
    
    # Import using spec to avoid package issues
    spec = importlib.util.spec_from_file_location("_core", ext_file)
    _core = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(_core)
    
except ImportError as e:
    raise ImportError(
        f"ERROR: Unable to import C++ backend '_core' module.\n"
        f"This library requires the C++ backend and does not support Python fallbacks.\n"
        f"Original error: {str(e)}"
    ) from e

# Re-export the Region enum from _core for easier use
Region = _core.Region

class Polygon:
    """Represents a convex polygon with various operations."""
    
    def __init__(self, vertices):
        """Create a polygon from vertex array.
        
        Args:
            vertices: Array of shape (n, 2) containing CCW-ordered vertices
        """
        self._vertices = np.asarray(vertices, dtype=float)
        self._impl = _core.Polygon(self._vertices)
    
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
        # Create using _core implementation
        impl = _core.Polygon.regular(sides, radius)
        
        # Create a new instance without calling __init__
        instance = cls.__new__(cls)
        instance._impl = impl
        instance._vertices = np.array(impl.vertices())
        return instance
    
    @property
    def vertices(self):
        """Get polygon vertices."""
        return self._vertices
    
    @property
    def area(self):
        """Get polygon area."""
        return self._impl.area()
    
    @property
    def perimeter(self):
        """Get polygon perimeter."""
        return self._impl.perimeter()
    
    def is_convex(self):
        """Check if polygon is convex."""
        return self._impl.is_convex()
    
    def contains(self, point):
        """Check if polygon contains a point."""
        point = np.asarray(point, dtype=float)
        return self._impl.contains(point)

class ProjectionCalculator:
    """Handles projection calculations for a polygon."""
    
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
        
        # Create C++ implementation
        self._impl = _core.ProjectionCalculator(
            polygon._impl,
            list(self.direction),
            region,
            self.random_seed
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
            return self._impl.exact_projected_distance()
        else:
            return self._impl.monte_carlo_projected_distance(samples)
    
    def average_distance(self, exact=True, samples=10000):
        """Calculate average distance over all directions.
        
        Args:
            exact: Whether to use exact calculation (True) or Monte Carlo (False)
            samples: Number of samples for Monte Carlo (ignored if exact=True)
            
        Returns:
            Expected Euclidean distance between points
        """
        if exact:
            return self._impl.exact_average_distance()
        else:
            return self._impl.monte_carlo_average_distance(samples)
    
    def monte_carlo_analysis(self, samples=10000, batch_size=1000, yield_interval=10):
        """Generator for progressive Monte Carlo calculations.
        
        Args:
            samples: Total number of samples to use
            batch_size: Size of each batch of samples
            yield_interval: Number of batches after which to yield results
            
        Yields:
            Dict containing current results and progress information
        """
        # Use C++ backend for efficient batched Monte Carlo
        samples_processed = 0
        total_samples = samples
        results_to_report = 0
        running_proj_sum = 0.0
        running_euc_sum = 0.0
        
        # Process in batches for progress updates
        while samples_processed < total_samples:
            current_batch = min(batch_size, total_samples - samples_processed)
            
            # Run the batch calculation
            result = self._impl.run_monte_carlo_batched(current_batch, min(current_batch, 100))
            
            # Update running totals
            batch_proj_dist = result["projected_distance"] * current_batch
            batch_avg_dist = result["average_distance"] * current_batch
            running_proj_sum += batch_proj_dist
            running_euc_sum += batch_avg_dist
            samples_processed += current_batch
            results_to_report += current_batch
            
            # Yield results at intervals
            if results_to_report >= yield_interval * batch_size or samples_processed >= total_samples:
                yield {
                    'samples_processed': samples_processed,
                    'total_samples': total_samples,
                    'projected_distance': running_proj_sum / samples_processed if samples_processed > 0 else 0.0,
                    'average_distance': running_euc_sum / samples_processed if samples_processed > 0 else 0.0,
                    'progress': samples_processed / total_samples
                }
                results_to_report = 0
    
    def sample_points(self, count=1):
        """Sample points from the specified region.
        
        Args:
            count: Number of points to sample
            
        Returns:
            Array of shape (count, 2) containing sampled points
        """
        if count <= 0:
            raise ValueError("Count must be positive")
        return self._impl.sample_points(count)
    
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
        return self._impl.project_points(points)
    
    def compare_methods(self, samples=10000):
        """Compare exact and Monte Carlo methods.
        
        Args:
            samples: Number of samples for Monte Carlo estimation
            
        Returns:
            Dict containing comparison results
        """
        exact = self.projected_distance(exact=True)
        mc = self.projected_distance(exact=False, samples=samples)
        
        return {
            'region': self.region.name,
            'exact': exact,
            'monte_carlo': mc,
            'absolute_error': abs(exact - mc),
            'relative_error': abs(exact - mc) / exact if exact != 0 else float('inf'),
            'n_samples': samples
        }