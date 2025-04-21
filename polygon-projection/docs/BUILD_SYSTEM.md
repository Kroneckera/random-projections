# Polygon Projection Build System

This document describes the build system for the polygon-projection library's Python bindings.

## Overview

The Python bindings for polygon-projection are built using a combination of:

1. CMake for building the C++ code
2. pybind11 for creating Python bindings 
3. setuptools for Python packaging

The build system focuses on:
- Robustness: handles errors gracefully
- Portability: works across different platforms
- Simplicity: easy to understand and maintain
- Automation: automatically manages dependencies

## Key Components

### 1. setup.py

The main build script that:
- Checks for CMake and required Python packages
- Configures and builds the C++ extension module
- Installs the Python package and its dependencies

```bash
# Install in development mode
pip install -e .

# Install normally
pip install .

# Build without installing
python setup.py build_ext
```

### 2. build_extension.py

A simplified script for directly building the extension module without installation:

```bash
# Build the extension module
python build_extension.py
```

This is useful during development when you only need to rebuild the C++ extension.

### 3. CMakeLists.txt

The CMake configuration for building the C++ code, including:
- C++17 standard requirements
- pybind11 integration
- Python detection and linking

## Dependencies

### Runtime Dependencies
- numpy: for array operations

### Build Dependencies
- CMake (≥ 3.14): for building the C++ code
- pybind11: for creating Python bindings
- C++ compiler with C++17 support

## Building for Different Platforms

The build system automatically handles different platforms:

### Linux/macOS
- Automatically uses the appropriate compiler flags
- Supports custom compiler selection via CXX environment variable

```bash
CXX=clang++ pip install .
```

### Windows
- Handles MSVC compiler specifics
- Configures for 32-bit or 64-bit builds automatically

## Best Practices

1. Always use the build system rather than manual compilation
2. Keep dependencies updated
3. Test the build on multiple platforms
4. Use virtual environments for development

## Troubleshooting

If you encounter build issues:

1. Make sure CMake is installed and in your PATH
2. Verify your C++ compiler supports C++17
3. Check for Python development headers
4. Ensure pybind11 is installed correctly