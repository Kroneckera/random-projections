# Polygon Projection Fallback Implementation

This directory contains a fallback implementation for the Polygon Projection library. The fallback uses pure Python code derived from the legacy implementation to provide the same functionality as the C++ backend.

## Overview

The Polygon Projection library normally uses a C++ backend for high-performance calculations. However, if the C++ backend is not available (e.g., if the user hasn't built the C++ extension), the library will automatically fall back to using this pure Python implementation.

This enables users to get started with the library without having to build the C++ extension, though at the cost of significantly reduced performance.

## Implementation Details

- The fallback implementation is in `polygon_projection/fallback.py`
- It uses the legacy Python code from `/legacy/` to implement the same API as the C++ version
- It is automatically used when the C++ extension cannot be loaded
- You can force the use of the fallback implementation by setting the environment variable `POLYGON_PROJECTION_FORCE_FALLBACK=1`

## Usage

You don't need to do anything special to use the fallback implementation. The library will automatically use it if the C++ backend is not available. However, you will see a warning message when the fallback is used.

```python
import polygon_projection as pp

# Check if using C++ backend
is_using_cpp = pp.using_cpp_backend()
print(f"Using C++ backend: {is_using_cpp}")

# Create a polygon and use it normally
polygon = pp.Polygon.regular(6, 1.0)
calculator = pp.ProjectionCalculator(polygon, [1.0, 0.0], pp.Region.INTERIOR)
distance = calculator.projected_distance(exact=True)
```

## Performance Considerations

The fallback implementation is significantly slower than the C++ backend, especially for large polygons or when using Monte Carlo simulations with many samples. For production use, we strongly recommend building and using the C++ backend.

## Testing the Fallback

You can test the fallback implementation using the provided `test_fallback.py` script, which forces the use of the fallback implementation:

```bash
python test_fallback.py
```

## Building the C++ Backend

To build the C++ backend and avoid using the fallback, see the main README.md file in the project root directory for instructions.