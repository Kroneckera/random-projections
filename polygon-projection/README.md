# Polygon Projection Library

A high-performance C++ library with Python bindings for performing projection calculations on convex polygons.

## Overview

This library computes expected distances between randomly chosen points within convex polygons, using both exact analytical formulas and Monte Carlo methods. These calculations are useful in computational geometry, physics simulations, and computer graphics.

### Key Features

- Create and manipulate convex polygons in 2D space
- Calculate exact analytical projections 
- Perform Monte Carlo sampling and approximations
- Work with both polygon interiors and boundaries
- Fast C++ implementation with easy-to-use Python interface

## Quick Example (Python)

```python
import sys
sys.path.insert(0, '/path/to/polygon-projection/python')

from polygon_projection.api import Polygon, Region, ProjectionCalculator

# Create a hexagon
hexagon = Polygon.regular(sides=6, radius=1.0)

# Set up a projection calculator
calc = ProjectionCalculator(
    polygon=hexagon,
    direction=[1.0, 0.0],  # Project along x-axis
    region=Region.INTERIOR
)

# Calculate exact projection distance
exact = calc.projected_distance(exact=True)
print(f"Exact projection distance: {exact:.6f}")

# Try Monte Carlo method
monte_carlo = calc.projected_distance(exact=False, samples=10000)
print(f"Monte Carlo estimate: {monte_carlo:.6f}")
```

## How to Build

<details>
<summary><b>Build Requirements (Click to expand)</b></summary>

### All Platforms
- C++ compiler with C++17 support
- CMake 3.14+
- Python 3.7+ (for Python bindings)
- NumPy (installed automatically if missing)
- pybind11 (installed automatically if missing)

### macOS
- Xcode Command Line Tools: `xcode-select --install`
- (Optional) CMake via Homebrew: `brew install cmake`

### Linux
- GCC 7+ or Clang 5+
- Python dev headers: 
  - Ubuntu/Debian: `sudo apt install python3-dev`
  - Fedora/RHEL: `sudo dnf install python3-devel`

### Windows
- Visual Studio Build Tools (recommended):
  - Download from: https://visualstudio.microsoft.com/visual-cpp-build-tools/
  - Select only "C++ build tools" during installation
- CMake from https://cmake.org/download/ (add to PATH)
- Python from https://www.python.org/
</details>

<details>
<summary><b>Building the C++ Library (Click to expand)</b></summary>

```bash
# Create build directory
mkdir -p build && cd build

# Configure
cmake ..

# Build
cmake --build .
```

### Configuration Options

- `-DBUILD_PYTHON_BINDINGS=ON/OFF` - Build Python bindings (default: ON)
- `-DBUILD_TESTING=ON/OFF` - Build tests (default: ON)
- `-DBUILD_EXAMPLES=ON/OFF` - Build examples (default: ON)
- `-DCMAKE_BUILD_TYPE=Release/Debug` - Build type (default: Release)
</details>

<details>
<summary><b>Building Python Bindings (Click to expand)</b></summary>

The simplest way to build the Python bindings is:

```bash
# From the project root
cd python

# Run the build script
python build.py
```

This will:
1. Build the C++ library
2. Build the Python extension module
3. Place the extension module in the `python/polygon_projection` directory

After building, you use the library by adding the Python directory to your path:

```python
import sys
sys.path.insert(0, '/path/to/polygon-projection/python')
from polygon_projection.api import Polygon, Region, ProjectionCalculator
```

If you make changes to the C++ code, simply run the build script again.
</details>

## API Overview

<details>
<summary><b>Python API (Click to expand)</b></summary>

### Classes

- `Polygon`: Represents a convex polygon
- `Region`: Enum for specifying regions (INTERIOR or BOUNDARY)
- `ProjectionCalculator`: Performs projection calculations

### Key Methods

- `Polygon.regular(sides, radius)`: Creates a regular polygon
- `polygon.area`: Property returning the polygon's area
- `polygon.perimeter`: Property returning the polygon's perimeter
- `polygon.is_convex()`: Checks if the polygon is convex
- `calc.projected_distance(exact=True)`: Calculates exact projected distance
- `calc.projected_distance(exact=False, samples=10000)`: Uses Monte Carlo method
- `calc.sample_points(count)`: Samples points from the specified region

### Example Usage

See the complete example in [`use_polygon.py`](/use_polygon.py) in the project root directory.
</details>

<details>
<summary><b>C++ API (Click to expand)</b></summary>

### Main Classes

- `polygon_projection::Polygon`: Represents a convex polygon
- `polygon_projection::Region`: Enum for specifying regions
- `polygon_projection::ProjectionCalculator`: Performs calculations

### Example Usage

```cpp
#include <polygon_projection/polygon.h>
#include <polygon_projection/projection.h>
#include <vector>

int main() {
    // Create a square
    std::vector<std::array<double, 2>> vertices = {
        {0, 0}, {1, 0}, {1, 1}, {0, 1}
    };
    
    // Create polygon
    polygon_projection::Polygon square(vertices);
    
    // Create a projection calculator for the interior
    std::array<double, 2> direction = {1.0, 0.0};  // Project along x-axis
    polygon_projection::ProjectionCalculator calc(
        square, direction, polygon_projection::Region::INTERIOR
    );
    
    // Calculate the projected distance
    double distance = calc.exact_projected_distance();
    
    return 0;
}
```

For more examples, see the [`examples/`](examples/) directory.
</details>

## Project Structure

```
polygon-projection/
├── include/              # C++ header files
├── src/                  # C++ implementation files
├── python/               # Python bindings
├── examples/             # Example C++ programs
├── tests/                # Test files
│   ├── cpp/              # C++ tests
│   └── python/           # Python tests
└── docs/                 # Documentation
```

## License

This project is available for educational and research purposes.