# Polygon Projection API Overview

## Introduction

The Polygon Projection library provides tools for working with convex polygons and calculating projection-related metrics. It supports both exact analytic calculations and Monte Carlo approximations.

## Core Classes

### Polygon

Represents a convex polygon in 2D space.

```python
from polygon_projection import Polygon

# Create from vertices
vertices = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
square = Polygon(vertices)

# Create a regular polygon
hexagon = Polygon.regular(sides=6, radius=2.0)

# Get properties
area = hexagon.area
perimeter = hexagon.perimeter
is_convex = hexagon.is_convex()
```

### Region

Enum specifying which region of the polygon to use for calculations.

```python
from polygon_projection import Region

# Available options
interior_region = Region.INTERIOR  # For sampling from the polygon interior
boundary_region = Region.BOUNDARY  # For sampling from the polygon boundary
```

### ProjectionCalculator

Main class for projection calculations.

```python
from polygon_projection import ProjectionCalculator

# Create calculator for a polygon and direction
calc = ProjectionCalculator(
    polygon=hexagon,
    direction=np.array([1.0, 0.5]),
    region=Region.INTERIOR,
    random_seed=42  # Optional
)

# Calculate projected distances
exact_dist = calc.projected_distance(exact=True)
monte_carlo_dist = calc.projected_distance(exact=False, samples=100000)

# Calculate average distances
exact_avg = calc.average_distance(exact=True)
monte_carlo_avg = calc.average_distance(exact=False, samples=100000)

# Compare methods
comparison = calc.compare_methods(samples=100000)
```

## Detailed Method Documentation

### Polygon Methods

- `__init__(vertices)`: Create polygon from array of vertices
- `from_vertices(vertices)`: Class method to create polygon from vertices
- `regular(sides, radius)`: Class method to create regular polygon
- `vertices`: Property returning polygon vertices as numpy array
- `area`: Property returning polygon area
- `perimeter`: Property returning polygon perimeter
- `is_convex()`: Check if polygon is convex

### ProjectionCalculator Methods

- `__init__(polygon, direction, region, random_seed)`: Initialize calculator
- `projected_distance(exact, samples)`: Calculate projected distance
- `average_distance(exact, samples)`: Calculate average Euclidean distance
- `monte_carlo_analysis(samples, batch_size, yield_interval)`: Run Monte Carlo with progress reporting
- `sample_points(count)`: Sample points from specified region
- `project_points(points)`: Project points onto direction
- `compare_methods(samples)`: Compare exact and Monte Carlo methods

## Generator Interface

For progressive Monte Carlo calculations with UI updates:

```python
calculator = ProjectionCalculator(polygon, direction, region)

# Run Monte Carlo with progress updates
for result in calculator.monte_carlo_analysis(
    samples=100000,
    batch_size=1000,
    yield_interval=10
):
    # Update UI with progress information
    progress = result['samples_processed'] / result['total_samples']
    print(f"Progress: {progress:.1%}")
    
    # Get current estimates
    projected_dist = result['projected_distance']
    average_dist = result['average_distance']
```

## Performance Considerations

- Exact calculations are fast for simple polygons but may be slower for complex shapes
- Monte Carlo methods scale well with polygon complexity
- Use batched Monte Carlo for responsive UI updates during long calculations
- Setting a random seed ensures reproducible results