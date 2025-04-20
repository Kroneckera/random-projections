# Random Polygon Projections

An interactive mathematical visualization tool for studying projections of points in convex polygons and computing expected distances.

## Overview

This project provides tools for analyzing the projections of uniform distributions on convex polygons onto a given direction. It computes both exact analytical solutions and Monte Carlo approximations for the expected distance between two randomly selected points after projection.

The mathematical focus is on:
1. Interior projections: points sampled from inside the polygon
2. Boundary projections: points sampled from the polygon perimeter

Both approaches yield interesting mathematical properties that can be explored through the interactive interface.

## Components

- **polygon_projection.py**: Unified high-level API for polygon projections
- **projection_sampling.py**: Analytical and sampling methods for interior points
- **edge_projection.py**: Analytical and sampling methods for boundary points
- **piecewise_integration.py**: High-accuracy integration for piecewise smooth functions with kink points
- **interactive_projection.py**: Interactive visualization tool with real-time computations
- **tests/**: Directory containing all test files
  - **test_polygon_projection.py**: Test suite for projection functionality
  - **test_average_distances.py**: Test suite comparing Monte Carlo with exact integration methods
  - Additional test files for verification and benchmarking

## Mathematical Background

For a convex polygon and a direction vector, the tool calculates:

- The exact probability density function of projections
- The exact expected absolute difference E[|X-Y|] between two random points' projections
- The exact average Euclidean distance E[||X-Y||] between random points (over all directions)
- High-accuracy numerical integration using Gauss-Legendre quadrature
- Polygon geometric properties (area and perimeter)
- Monte Carlo approximations of expected distances with configurable sample sizes

These calculations are performed for both:
- Points sampled uniformly from the **interior** of the convex polygon
- Points sampled uniformly from the **boundary** of the convex polygon

## Interactive Features

The interactive GUI allows you to:

- Create and manipulate convex polygons by dragging vertices
- Set the projection direction using two draggable endpoints
- View polygon geometric properties (area and perimeter) in real-time
- View both interior and boundary projection results simultaneously
- View exact average Euclidean distances across all directions
- Run Monte Carlo approximations with configurable sample sizes
- Compare exact analytical results with Monte Carlo approximations for both projected and average distances
- Add vertices (placed intelligently to maintain convexity)
- Remove vertices and reset the polygon
- Enforced convexity during vertex manipulation

## Installation and Usage

### Requirements
- Python 3.x
- NumPy
- Matplotlib
- Pytest (for running tests)

### Installation

```bash
git clone https://github.com/Kroneckera/random-projections.git
cd random-projections
pip install numpy matplotlib pytest
```

### Running the Interactive Tool

```bash
python interactive_projection.py
```

### Running Tests

```bash
# Run all tests
pytest tests/

# Run specific test files
pytest tests/test_polygon_projection.py
pytest tests/test_average_distances.py
```

## Demonstration

When using the interactive tool:
1. Drag the blue vertices to reshape the polygon (convexity is automatically maintained)
2. Drag both red points to change the projection direction
3. View polygon information (vertices, area, perimeter) in real-time
4. View both interior and boundary projection distances simultaneously
5. See average Euclidean distances across all directions
6. Adjust the Monte Carlo sample size using the slider
7. Click "Calculate Monte Carlo" to run simulations for both projected and average distances
8. View exact and approximate calculations in the results panel with error percentages

## License

MIT License