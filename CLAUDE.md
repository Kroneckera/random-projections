# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build/Test Commands
- Run a Python script: `python script_name.py`
- Run tests (if added): `pytest` or `pytest test_script.py`
- Run with specific seed: `python script_name.py --seed=42`

## Code Style Guidelines
- Use NumPy style docstrings with parameter descriptions and return types
- Line length: Keep under 88 characters 
- Imports: Group standard library, third-party (numpy), and local imports
- Variable naming: Use snake_case for variables and functions
- Error handling: Use explicit error messages with ValueError
- Type hints: Not currently used, but can be added if desired
- Arrays: Use np.asarray() to convert inputs with explicit dtype=float
- Constants: Use ALL_CAPS for constants
- No Unicode: All code uses ASCII only (no unicode symbols)
- Vectorized operations: Prefer NumPy vectorized operations over loops