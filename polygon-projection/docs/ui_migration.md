# UI Migration Guide

This guide outlines how to migrate the existing UI code to use the new C++ backend with the improved Python API.

## Overview

The existing UI in `interactive_projection.py` relies on the original Python backend modules. To integrate with the new backend, we'll need to:

1. Update import statements
2. Convert the existing `PolygonRegion` enum usage to the new `Region` enum
3. Replace the `PolygonProjection` class usage with the new `ProjectionCalculator` class
4. Adapt the Monte Carlo generator pattern to work with the new API

## Step-by-Step Migration

### 1. Update Imports

Original imports:
```python
from polygon_projection import PolygonProjection, PolygonRegion, generate_regular_polygon
```

New imports:
```python
from polygon_projection import Polygon, Region, ProjectionCalculator
```

### 2. Replace Regular Polygon Generation

Original code:
```python
self.vertices = generate_regular_polygon(n_sides=6, radius=2.0)
```

New code:
```python
polygon = Polygon.regular(sides=6, radius=2.0)
self.vertices = polygon.vertices
```

### 3. Replace PolygonRegion with Region

Original code:
```python
region = PolygonRegion.INTERIOR
# or
region = PolygonRegion.BOUNDARY
```

New code:
```python
region = Region.INTERIOR
# or
region = Region.BOUNDARY
```

### 4. Replace PolygonProjection with ProjectionCalculator

Original code:
```python
interior_proj = PolygonProjection(
    self.vertices, self.direction, 
    region=PolygonRegion.INTERIOR,
    seed=42
)
```

New code:
```python
polygon = Polygon(self.vertices)
interior_calc = ProjectionCalculator(
    polygon=polygon,
    direction=self.direction,
    region=Region.INTERIOR,
    random_seed=42
)
```

### 5. Update Method Calls

#### Exact Expected Distance

Original code:
```python
exact_value = interior_proj.exact_expected_distance()
```

New code:
```python
exact_value = interior_calc.projected_distance(exact=True)
```

#### Monte Carlo Expected Distance

Original code:
```python
mc_value = interior_proj.monte_carlo_expected_distance(n_samples=10000)
```

New code:
```python
mc_value = interior_calc.projected_distance(exact=False, samples=10000)
```

#### Exact Average Distance

Original code:
```python
avg_value = interior_proj.average_distance_exact()
```

New code:
```python
avg_value = interior_calc.average_distance(exact=True)
```

#### Monte Carlo Generator

Original code:
```python
interior_generator = interior_proj.monte_carlo_generator(
    n_samples=self.n_samples,
    batch_size=batch_size,
    yield_interval=yield_interval
)

for interior_results in interior_generator:
    # Process results
    projected_distance = interior_results['projected_distance']
    average_distance = interior_results['average_distance']
    samples_processed = interior_results['samples_processed']
```

New code:
```python
for interior_results in interior_calc.monte_carlo_analysis(
    samples=self.n_samples,
    batch_size=batch_size,
    yield_interval=yield_interval
):
    # Process results - same format as before
    projected_distance = interior_results['projected_distance']
    average_distance = interior_results['average_distance']
    samples_processed = interior_results['samples_processed']
```

### 6. Additional Methods

#### Area and Perimeter

Original code:
```python
area = interior_proj.get_area()
perimeter = interior_proj.get_perimeter()
```

New code:
```python
area = polygon.area
perimeter = polygon.perimeter
```

#### Sampling Points

Original code:
```python
point = interior_proj.sample_point()
```

New code:
```python
point = interior_calc.sample_point()
```

## Complete Example: Calculator Initialization

Original code:
```python
# Initialize PolygonProjection objects
interior_proj = PolygonProjection(
    self.vertices.copy(), self.direction.copy(), 
    region=PolygonRegion.INTERIOR,
    seed=42
)
boundary_proj = PolygonProjection(
    self.vertices.copy(), self.direction.copy(), 
    region=PolygonRegion.BOUNDARY,
    seed=43
)

# Calculate exact values
interior_exact = interior_proj.exact_expected_distance()
boundary_exact = boundary_proj.exact_expected_distance()
```

New code:
```python
# Create polygon
polygon = Polygon(self.vertices.copy())

# Initialize ProjectionCalculator objects
interior_calc = ProjectionCalculator(
    polygon=polygon,
    direction=self.direction.copy(),
    region=Region.INTERIOR,
    random_seed=42
)
boundary_calc = ProjectionCalculator(
    polygon=polygon,
    direction=self.direction.copy(),
    region=Region.BOUNDARY,
    random_seed=43
)

# Calculate exact values
interior_exact = interior_calc.projected_distance(exact=True)
boundary_exact = boundary_calc.projected_distance(exact=True)
```

## Complete Example: Monte Carlo Worker

Original code:
```python
# Create generators for interior and boundary calculations
interior_generator = interior_proj.monte_carlo_generator(
    n_samples=self.n_samples, 
    batch_size=batch_size,
    yield_interval=yield_interval
)
boundary_generator = boundary_proj.monte_carlo_generator(
    n_samples=self.n_samples,
    batch_size=batch_size,
    yield_interval=yield_interval
)

# Process both generators in parallel
while not (interior_done and boundary_done):
    if not interior_done:
        try:
            interior_results = next(interior_generator)
        except StopIteration:
            interior_done = True
    
    if not boundary_done:
        try:
            boundary_results = next(boundary_generator)
        except StopIteration:
            boundary_done = True
```

New code:
```python
# Create generators for interior and boundary calculations
interior_generator = interior_calc.monte_carlo_analysis(
    samples=self.n_samples, 
    batch_size=batch_size,
    yield_interval=yield_interval
)
boundary_generator = boundary_calc.monte_carlo_analysis(
    samples=self.n_samples,
    batch_size=batch_size,
    yield_interval=yield_interval
)

# Process both generators in parallel (no change in loop structure)
while not (interior_done and boundary_done):
    if not interior_done:
        try:
            interior_results = next(interior_generator)
        except StopIteration:
            interior_done = True
    
    if not boundary_done:
        try:
            boundary_results = next(boundary_generator)
        except StopIteration:
            boundary_done = True
```

## Final Notes

1. The new API maintains the same core functionality but with cleaner, more consistent naming
2. The geometry concepts remain unchanged, only the API has been updated
3. The Monte Carlo generator returns results in the same format for easy integration
4. Performance should be significantly improved due to the C++ backend