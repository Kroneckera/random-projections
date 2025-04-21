# Polygon Projection

A high-performance library for calculating projections and expected distances between points in convex polygons, with both exact analytical solutions and Monte Carlo approximations.

## Overview

This project provides tools for analyzing the projections of uniform distributions on convex polygons onto a given direction. It computes both exact analytical solutions and Monte Carlo approximations for the expected distance between two randomly selected points after projection.

The library features:
- Fast C++ implementation with Python bindings
- Pure Python fallback implementation
- Interactive visualization tools
- Comprehensive testing suite

## Mathematical Background

For a convex polygon and a direction vector, the library calculates:

- The exact probability density function of projections
- The exact expected absolute difference E[|X-Y|] between two random points' projections
- The exact average Euclidean distance E[||X-Y||] between random points (over all directions)
- High-accuracy numerical integration using Gauss-Legendre quadrature
- Polygon geometric properties (area and perimeter)
- Monte Carlo approximations with configurable sample sizes

These calculations are performed for both:
- Points sampled uniformly from the **interior** of the convex polygon
- Points sampled uniformly from the **boundary** of the convex polygon

## Project Structure

```
polygon-projection/
├── include/              # C++ header files
├── src/                  # C++ implementation files
├── python/               # Python bindings & fallback implementation
│   └── polygon_projection/
│       ├── __init__.py   # Package initialization
│       ├── api.py        # API definitions
│       ├── fallback.py   # Pure Python fallback implementation
│       └── legacy/       # Legacy Python code for fallback
├── legacy/               # Original Python implementation (deprecated)
├── examples/             # Example usage
├── tests/                # Test suite
│   ├── cpp/              # C++ tests
│   ├── python/           # Python tests
│   └── api/              # API tests
└── docs/                 # Documentation
```

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/Kroneckera/random-projections.git
cd random-projections

# Install dependencies
pip install numpy matplotlib

# Add the package to your Python path
export PYTHONPATH=$PYTHONPATH:$(pwd)/polygon-projection/python
```

### Basic Usage

```python
import sys
sys.path.insert(0, './polygon-projection/python')

from polygon_projection.api import Polygon, Region, ProjectionCalculator

# Create a regular hexagon
hexagon = Polygon.regular(sides=6, radius=1.0)

# Set up a projection calculator for interior points along x-axis
calc = ProjectionCalculator(
    polygon=hexagon,
    direction=[1.0, 0.0],
    region=Region.INTERIOR
)

# Calculate exact projection distance
exact = calc.projected_distance(exact=True)
print(f"Exact projection distance: {exact:.6f}")

# Try Monte Carlo method
monte_carlo = calc.projected_distance(exact=False, samples=10000)
print(f"Monte Carlo estimate: {monte_carlo:.6f}")
```

### Running the Interactive Tool

The interactive visualization tool allows you to:
- Create and manipulate convex polygons by dragging vertices
- Adjust the projection direction
- View real-time calculations of projection metrics
- Run Monte Carlo simulations with progress tracking

```bash
python interactive_projection.py
```

## C++ Implementation

The C++ implementation provides high-performance calculations and is used by default if available. To build it:

```bash
# Create build directory
mkdir -p polygon-projection/build && cd polygon-projection/build

# Configure and build
cmake ..
cmake --build .
```

## Python Fallback Implementation

If the C++ implementation is not available, the library automatically falls back to a pure Python implementation based on the legacy code. This ensures the library works across all platforms but will be significantly slower.

You can force the use of the fallback implementation:

```python
import os
os.environ["POLYGON_PROJECTION_FORCE_FALLBACK"] = "1"

import polygon_projection as pp
# Now uses Python fallback implementation
```

## API Reference

### Key Classes

- **Polygon**: Represents a convex polygon
  - `Polygon(vertices)`: Create from vertex array
  - `Polygon.regular(sides, radius)`: Create regular polygon
  - Properties: `vertices`, `area`, `perimeter`
  - Methods: `is_convex()`

- **Region**: Enum for specifying the region
  - `Region.INTERIOR`: Interior of the polygon
  - `Region.BOUNDARY`: Boundary of the polygon

- **ProjectionCalculator**: Performs projection calculations
  - `ProjectionCalculator(polygon, direction, region, random_seed=None)`
  - `projected_distance(exact=True, samples=10000)`: Calculate projected distance
  - `average_distance(exact=True, samples=10000)`: Calculate average distance
  - `monte_carlo_analysis(samples, batch_size, yield_interval)`: Generator for Monte Carlo with progress tracking
  - `sample_points(count)`: Sample points from region

## Testing

```bash
# Run API tests
cd polygon-projection/tests/api
python test_fallback_api.py

# Run legacy implementation tests
python test_legacy_fallback.py
```

## Performance Considerations

The C++ implementation is significantly faster than the Python fallback:

| Operation | C++ Performance | Python Performance |
|-----------|----------------|-------------------|
| Exact projection calculation | 0.0027 ms | ~0.5 ms |
| Monte Carlo (10,000 samples) | 9.41 ms | ~250 ms |
| Exact average distance | 1.45 ms | ~25 ms |

For performance-critical applications, building the C++ extension is strongly recommended.

## Development Status

The project is functional with both C++ and Python implementations. The C++ implementation is complete but has some known issues with Python bindings. The Python fallback implementation provides 100% of the functionality but with reduced performance.

## License

This project is available for educational and research purposes.

## Acknowledgements

This project originated as a mathematical exploration of projection distances in convex polygons and has evolved into a full-fledged library with both C++ and Python implementations.