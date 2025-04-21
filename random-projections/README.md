# Random Projections

A library for calculating projections and distances in convex polygons.

## New Directory Structure

This project has been refactored to simplify the directory structure:

```
random-projections/
├── cpp/               # All C++ code 
│   ├── include/       # C++ headers (directly, without polygon_projection subdir)
│   ├── src/           # C++ implementations
│   └── tests/         # C++ tests
├── python/            # All Python code
│   ├── projection/    # Main Python package 
│   │   ├── __init__.py
│   │   ├── api.py
│   │   ├── fallback.py
│   │   └── legacy/    # Legacy code for fallback
│   ├── tests/         # Python tests 
│   └── examples/      # Python examples
├── docs/              # Documentation
└── scripts/           # Build and utility scripts
```

## Installation

### Installing from Source

```bash
# Clone the repository
git clone https://github.com/your-username/random-projections.git
cd random-projections

# Install the Python package
pip install -e python/
```

### Building C++ Only

```bash
mkdir build
cd build
cmake ..
make
```

## Usage

### Python

```python
from projection import Polygon, Region, ProjectionCalculator

# Create a polygon
polygon = Polygon.regular(5, 1.0)  # Regular pentagon with radius 1.0

# Calculate projections
calculator = ProjectionCalculator(polygon)
result = calculator.project_point([0.5, 0.5], Region.INTERIOR)
print(f"Projection distance: {result.distance}")
```

### C++

```cpp
#include <polygon.h>
#include <projection.h>

using namespace projection;

int main() {
    // Create a polygon
    auto polygon = Polygon::regular(5, 1.0);  // Regular pentagon with radius 1.0
    
    // Calculate projections
    ProjectionCalculator calculator(polygon);
    auto result = calculator.project_point({0.5, 0.5}, Region::INTERIOR);
    std::cout << "Projection distance: " << result.distance << std::endl;
    
    return 0;
}
```

## Documentation

For more detailed information, see the documentation in the [docs/](docs/) directory:

- [API Overview](docs/API_OVERVIEW.md)
- [Build System](docs/BUILD_SYSTEM_SUMMARY.md)

## License

This project is licensed under the MIT License - see the LICENSE file for details.