#!/usr/bin/env python3
"""
Build script for the polygon-projection Python bindings.

This script builds the C++ library and Python extension module directly,
placing the extension module in the polygon_projection directory.

Usage:
    python build.py
"""

import os
import sys
import subprocess
from pathlib import Path

def main():
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    build_dir = project_root / "build"
    python_dir = script_dir / "polygon_projection"
    
    # Make sure numpy and pybind11 are installed
    for package in ["numpy", "pybind11"]:
        try:
            __import__(package)
            print(f"Found {package}")
        except ImportError:
            print(f"Installing {package}...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])
    
    # Create directories
    build_dir.mkdir(exist_ok=True)
    python_dir.mkdir(exist_ok=True)
    
    print(f"Project root: {project_root}")
    print(f"Build directory: {build_dir}")
    print(f"Python module directory: {python_dir}")
    
    # Configure with CMake
    print("\nConfiguring with CMake...")
    
    # Configure CMake command
    cmake_cmd = [
        "cmake", "..",
        "-DBUILD_PYTHON_BINDINGS=ON",
        "-DCMAKE_BUILD_TYPE=Release",
        f"-DPYTHON_EXECUTABLE={sys.executable}",
        "-DBUILD_EXAMPLES=OFF",
        "-DBUILD_TESTING=OFF",
    ]
    
    try:
        subprocess.check_call(cmake_cmd, cwd=build_dir)
        
        # Build the library and Python module
        print("\nBuilding...")
        build_cmd = ["cmake", "--build", ".", "--config", "Release"]
        
        # Parallel build if possible
        if os.cpu_count():
            cpu_count = os.cpu_count()
            if sys.platform == "win32":
                # Use parallel build for MSVC
                build_cmd.extend(["--", f"/maxcpucount:{cpu_count}"])
            else:
                # Use parallel build for Make
                build_cmd.extend(["--", f"-j{cpu_count}"])
        
        try:
            subprocess.check_call(build_cmd, cwd=build_dir)
        except subprocess.CalledProcessError as e:
            print(f"\nBuild failed with error code {e.returncode}")
            print("Common issues:")
            print("- On Windows: Make sure you have Visual Studio Build Tools or a C++ compiler installed")
            print("- On macOS: Run 'xcode-select --install' to install the command line tools")
            print("- On Linux: Install the Python development headers (python-dev or python-devel package)")
            raise
        
        # Find the built module (with platform-specific extension)
        possible_extensions = [".pyd", ".so", ".dylib"]
        if sys.platform == "win32":
            possible_extensions = [".pyd", ".dll", ".so"]
        elif sys.platform == "darwin":
            possible_extensions = [".so", ".dylib"]
        
        module_found = False
        for ext in possible_extensions:
            module_path = python_dir / f"_core{ext}"
            if module_path.exists():
                module_found = True
                size_kb = module_path.stat().st_size / 1024
                print(f"\nSuccessfully built module: {module_path}")
                print(f"Module size: {size_kb:.1f} KB")
                break
                
        if not module_found:
            print(f"\nERROR: Failed to find the built module in {python_dir}")
            print("The build process completed, but the Python module wasn't created properly.")
            print("Check the CMake output for errors.")
            return 1
            
        # Create/update __init__.py if needed
        init_file = python_dir / "__init__.py"
        if not init_file.exists():
            with open(init_file, "w") as f:
                f.write('"""Python bindings for the polygon-projection library."""\n\n')
                f.write('from .api import Polygon, Region, ProjectionCalculator\n')
            print(f"Created {init_file}")
            
        print("\nBuild completed successfully!")
        print("You can now use the library by adding the following to your Python code:")
        print(f'import sys; sys.path.insert(0, "{script_dir}"); from polygon_projection import api')
        return 0
            
    except subprocess.CalledProcessError as e:
        print(f"\nERROR: Build failed: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())