# Python Tests for Polygon Projection Library

This directory contains Python tests for the polygon-projection library.

## Test Files

- `test_polygon_projection.py`: Core tests for the Polygon and ProjectionCalculator classes
- `test_analytical_comparison.py`: Verification against known analytical formulas
- `test_performance.py`: Performance benchmarks for different operations
- `convexity_test.py`: Visual test for convexity checking algorithm

## Running Tests

To run all tests in this directory (performance tests are skipped by default):

```bash
cd /path/to/polygon-projection
python -m pytest tests/python
```

To run a specific test file:

```bash
python -m pytest tests/python/test_polygon_projection.py -v
```

### Performance Tests

Performance tests are skipped by default to keep the test suite fast. To run performance tests:

```bash
# Run only the performance tests
python -m pytest tests/python/test_performance.py -v --run-performance

# Run all tests including performance tests
python -m pytest tests/python -v --run-performance
```

You can also run the performance test script directly:

```bash
python tests/python/test_performance.py
```

## Test Requirements

- pytest
- numpy
- matplotlib (for convexity_test.py only)

## Notes

- All tests use the Python API to test the C++ backend
- The analytical tests verify implementation against known mathematical formulas
- The performance tests benchmark different operations with varying polygon sizes and sample counts