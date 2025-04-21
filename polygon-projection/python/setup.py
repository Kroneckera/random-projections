#!/usr/bin/env python3
"""
Build script for polygon-projection Python bindings.

This script builds the C++ extension module for the polygon-projection library
using CMake and pybind11. It automatically handles dependencies and builds
the extension for the current Python version.

Requirements:
- CMake (>= 3.14)
- C++ compiler with C++17 support
- pybind11 (automatically installed if missing)
- numpy (automatically installed if missing)

Usage:
    # Install in development mode
    pip install -e .
    
    # Install normally
    pip install .
    
    # Build without installing
    python setup.py build_ext

    # Install with specific compiler
    CXX=clang++ pip install .
"""

import os
import sys
import platform
import subprocess
import setuptools
from setuptools.command.build_ext import build_ext
from setuptools import setup, Extension
import shutil

# Define version
VERSION = '0.1.0'

# Project root directory (parent of this script's directory)
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

class CMakeExtension(Extension):
    """Extension class for CMake-based builds."""
    def __init__(self, name, source_dir=""):
        Extension.__init__(self, name, sources=[])
        self.source_dir = os.path.abspath(source_dir)

class CMakeBuild(build_ext):
    """Custom build_ext command that builds using CMake."""
    
    def run(self):
        """Run the build process."""
        # Check for CMake
        try:
            cmake_version = subprocess.check_output(
                ["cmake", "--version"], 
                universal_newlines=True
            ).strip()
            cmake_version = cmake_version.split('\n')[0]
            print(f"Found {cmake_version}")
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the extension. "
                "Please install CMake from https://cmake.org/download/"
            )
            
        # Check for required packages
        self._check_dependencies()
            
        # Build extensions
        for ext in self.extensions:
            self.build_extension(ext)
    
    def _check_dependencies(self):
        """Check and install required dependencies."""
        required_packages = ['pybind11', 'numpy']
        
        for package in required_packages:
            try:
                __import__(package)
                version = subprocess.check_output(
                    [sys.executable, '-c', f'import {package}; print({package}.__version__)'],
                    universal_newlines=True
                ).strip()
                print(f"Found {package} {version}")
            except (ImportError, subprocess.CalledProcessError):
                print(f"Installing {package}...")
                subprocess.check_call([sys.executable, '-m', 'pip', 'install', package])
    
    def build_extension(self, ext):
        """Build a single extension."""
        # Determine output directory
        ext_dir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        
        # Print information about the build
        print(f"Building extension {ext.name} for Python {platform.python_version()} on {platform.system()}")
        
        # CMake configuration
        source_dir = ext.source_dir if ext.source_dir else ROOT_DIR
        build_temp = os.path.join(self.build_temp, ext.name)
        os.makedirs(build_temp, exist_ok=True)
        
        # Determine the CMake configuration based on the build type
        config = 'Debug' if self.debug else 'Release'
        
        # CMake args
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={ext_dir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DBUILD_PYTHON_BINDINGS=ON",
            f"-DCMAKE_BUILD_TYPE={config}",
            # Disable components not needed for Python bindings
            "-DBUILD_EXAMPLES=OFF",
        ]
        
        # Build args
        build_args = ['--config', config]
        
        # Handle platform-specific settings
        if platform.system() == "Windows":
            if sys.maxsize > 2**32:  # 64-bit Windows
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']  # Parallel build on Windows
        else:  # Unix-like (Linux, Mac)
            build_args += ['--', '-j', str(os.cpu_count() or 2)]
            
            # Add C++17 flag directly if on macOS or Linux
            if platform.system() in ["Darwin", "Linux"]:
                env = os.environ.copy()
                env['CXXFLAGS'] = f"{env.get('CXXFLAGS', '')} -std=c++17"
                
                # Set compiler if specified by the user
                if 'CXX' in env:
                    print(f"Using compiler: {env['CXX']}")
                    cmake_args.append(f"-DCMAKE_CXX_COMPILER={env['CXX']}")
        
        # Run CMake configure step
        print("Configuring CMake...")
        subprocess.check_call(['cmake', source_dir] + cmake_args, cwd=build_temp)
        
        # Run the build
        print("Building extension...")
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=build_temp)
        
        # Verify the extension module was built correctly
        extension_path = next(
            (f for f in os.listdir(ext_dir) if f.startswith(ext.name.split('.')[-1])),
            None
        )
        
        if extension_path:
            full_path = os.path.join(ext_dir, extension_path)
            print(f"Extension built successfully: {full_path} ({os.path.getsize(full_path)} bytes)")
        else:
            raise RuntimeError(f"Failed to build extension {ext.name}. No extension file found in {ext_dir}.")
        
        print(f"Build successful for {ext.name}")

setup(
    name="polygon-projection",
    version=VERSION,
    description="Polygon projection calculations library",
    packages=setuptools.find_packages(),
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.17.0",
    ],
    ext_modules=[
        CMakeExtension("polygon_projection._core"),
    ],
    cmdclass={
        "build_ext": CMakeBuild,
    },
    zip_safe=False,
)