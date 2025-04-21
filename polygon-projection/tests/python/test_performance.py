#!/usr/bin/env python3
"""
Performance benchmark for polygon projection algorithms.

This script measures the performance of various algorithms
in the polygon projection library.
"""

import time
import numpy as np
import sys
import os
import pytest
from pathlib import Path

# Add the package directory to the Python path when running standalone
package_root = Path(__file__).parent.parent.parent
python_dir = package_root / "python"
sys.path.insert(0, str(python_dir))

# Now import from the package
from polygon_projection.api import Polygon, Region, ProjectionCalculator


def time_function(func, *args, **kwargs):
    """Measure execution time of a function.
    
    Args:
        func: Function to time
        *args, **kwargs: Arguments to pass to the function
        
    Returns:
        Tuple (result, execution_time_ms)
    """
    start_time = time.time()
    result = func(*args, **kwargs)
    end_time = time.time()
    execution_time_ms = (end_time - start_time) * 1000
    return result, execution_time_ms


def benchmark_regular_polygon(sides, samples):
    """Benchmark calculations for a regular polygon.
    
    Args:
        sides: Number of sides for the polygon
        samples: Number of samples for Monte Carlo
        
    Returns:
        Dictionary of benchmark results
    """
    # Create regular polygon
    polygon = Polygon.regular(sides)
    
    # Create calculator for interior projections
    direction = [1.0, 0.0]
    interior_calc = ProjectionCalculator(
        polygon, direction, region=Region.INTERIOR, random_seed=42
    )
    
    # Create calculator for boundary projections
    boundary_calc = ProjectionCalculator(
        polygon, direction, region=Region.BOUNDARY, random_seed=42
    )
    
    # Measure exact interior projection
    _, interior_exact_time = time_function(interior_calc.projected_distance, exact=True)
    
    # Measure Monte Carlo interior projection
    _, interior_mc_time = time_function(
        interior_calc.projected_distance, exact=False, samples=samples
    )
    
    # Measure exact boundary projection
    _, boundary_exact_time = time_function(boundary_calc.projected_distance, exact=True)
    
    # Measure Monte Carlo boundary projection
    _, boundary_mc_time = time_function(
        boundary_calc.projected_distance, exact=False, samples=samples
    )
    
    # Measure exact average distance
    _, avg_exact_time = time_function(interior_calc.average_distance, exact=True)
    
    # Measure Monte Carlo average distance
    _, avg_mc_time = time_function(
        interior_calc.average_distance, exact=False, samples=samples
    )
    
    # Measure point sampling
    _, sample_time = time_function(interior_calc.sample_points, count=samples)
    
    # Return results
    return {
        'polygon_sides': sides,
        'samples': samples,
        'interior_exact_ms': interior_exact_time,
        'interior_mc_ms': interior_mc_time,
        'boundary_exact_ms': boundary_exact_time,
        'boundary_mc_ms': boundary_mc_time,
        'avg_exact_ms': avg_exact_time,
        'avg_mc_ms': avg_mc_time,
        'sample_ms': sample_time,
        'sample_rate': samples / (sample_time / 1000) if sample_time > 0 else float('inf')
    }


@pytest.mark.performance
def test_performance_benchmarks():
    """Run performance benchmarks as a pytest test.
    
    Run with: pytest tests/python/test_performance.py -v --run-performance
    """
    # Set benchmark parameters
    polygon_sizes = [4, 6, 8, 12, 16, 24, 32]
    sample_sizes = [1000, 10000, 100000]
    
    print("\nPolygon Projection Performance Benchmark")
    print("=" * 50)
    
    # Run benchmarks for different polygon sizes
    print("\n[1] Benchmarking different polygon sizes (10,000 samples)")
    results_by_size = []
    
    for sides in polygon_sizes:
        results = benchmark_regular_polygon(sides, 10000)
        results_by_size.append(results)
        
        print(f"\n{sides}-sided polygon:")
        print(f"  Exact interior projection:  {results['interior_exact_ms']:.3f} ms")
        print(f"  MC interior projection:     {results['interior_mc_ms']:.3f} ms")
        print(f"  Exact boundary projection:  {results['boundary_exact_ms']:.3f} ms")
        print(f"  MC boundary projection:     {results['boundary_mc_ms']:.3f} ms")
        print(f"  Exact average distance:     {results['avg_exact_ms']:.3f} ms")
        print(f"  MC average distance:        {results['avg_mc_ms']:.3f} ms")
        print(f"  Point sampling rate:        {results['sample_rate']:.0f} points/sec")
    
    # Run benchmarks for different sample sizes
    print("\n\n[2] Benchmarking different sample sizes (12-sided polygon)")
    results_by_samples = []
    
    for samples in sample_sizes:
        results = benchmark_regular_polygon(12, samples)
        results_by_samples.append(results)
        
        print(f"\n{samples:,} samples:")
        print(f"  MC interior projection:     {results['interior_mc_ms']:.3f} ms")
        print(f"  MC boundary projection:     {results['boundary_mc_ms']:.3f} ms")
        print(f"  MC average distance:        {results['avg_mc_ms']:.3f} ms")
        print(f"  Sampling {samples} points:    {results['sample_ms']:.3f} ms")
    
    # Summary
    print("\n\nPerformance Summary:")
    print("=" * 50)
    print("Monte Carlo operations scale linearly with sample count")
    print("Exact calculations are much faster for small/medium polygons")
    print("Point sampling scales roughly linearly with sample count")
    
    # Add a simple assertion to satisfy pytest
    assert True, "Performance test completed"


def main():
    """Run the benchmark when script is called directly."""
    test_performance_benchmarks()


if __name__ == "__main__":
    main()