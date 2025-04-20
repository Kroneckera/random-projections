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
- **test_polygon_projection.py**: Comprehensive test suite validating the mathematical implementations

## Mathematical Background

For a convex polygon and a direction vector, the tool calculates:

- The exact probability density function of projections
- The exact expected absolute difference E[|X-Y|] between two random points' projections
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
- Toggle between interior and boundary projections
- Run Monte Carlo approximations with configurable sample sizes
- Compare exact analytical results with Monte Carlo approximations
- Add/remove vertices and reset the polygon

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
pytest test_polygon_projection.py
```

## Demonstration

When using the interactive tool:
1. Drag the blue vertices to reshape the polygon
2. Drag both red points to change the projection direction
3. View polygon information (vertices, area, perimeter) in real-time
4. Use the "Interior"/"Boundary" radio buttons to switch projection regions
5. Adjust the Monte Carlo sample size using the slider
6. Click "Calculate Monte Carlo" to run simulations
7. View exact and approximate calculations in the results panel

## License

MIT License