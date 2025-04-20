"""
Configuration file for pytest.
This file is automatically recognized by pytest and allows configuring its behavior.
"""
import os
import sys

# Add root directory to Python path
# This ensures tests can import modules from the root directory
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))