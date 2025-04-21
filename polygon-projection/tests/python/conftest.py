"""
Configuration file for pytest.

This file sets up the Python path so that the tests can find the polygon_projection package.
"""

import sys
import os
from pathlib import Path
import pytest

# Add the package root directory to the Python path
package_root = Path(__file__).parent.parent.parent
python_dir = package_root / "python"
sys.path.insert(0, str(python_dir))

# Define command line options
def pytest_addoption(parser):
    parser.addoption(
        "--run-performance",
        action="store_true",
        default=False,
        help="Run performance tests (normally skipped)"
    )

# Register the custom markers
def pytest_configure(config):
    config.addinivalue_line("markers", "performance: mark test as a performance benchmark")

# Skip performance tests unless --run-performance is provided    
def pytest_collection_modifyitems(config, items):
    if not config.getoption("--run-performance"):
        skip_performance = pytest.mark.skip(reason="Use --run-performance to run")
        for item in items:
            if "performance" in item.keywords:
                item.add_marker(skip_performance)