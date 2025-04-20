#!/usr/bin/env python3
"""
Verify the analytical formula for average distance along square boundary
"""

import numpy as np

# The correct formula for square boundary is:
# (1/4) + (sqrt(2)/12) + (5*ln(1+sqrt(2)))/12

# Calculate step by step
sqrt2 = np.sqrt(2)
ln_term = np.log(1 + sqrt2)

term1 = 1/4
term2 = sqrt2/12
term3 = (5*ln_term)/12

result = term1 + term2 + term3

print("Calculation of analytical formula for square boundary average distance:")
print(f"sqrt(2) = {sqrt2}")
print(f"ln(1+sqrt(2)) = {ln_term}")
print(f"Term 1: 1/4 = {term1}")
print(f"Term 2: sqrt(2)/12 = {term2}")
print(f"Term 3: 5*ln(1+sqrt(2))/12 = {term3}")
print(f"Sum: {result}")

# Compare with Monte Carlo and integration results
monte_carlo = 0.735073  # From our simple_square_boundary_test.py
integration = 0.735090  # From our test_updated.py

print(f"\nMonte Carlo result: {monte_carlo}")
print(f"Integration result: {integration}")
print(f"Analytical formula: {result}")

print(f"\nError formula vs Monte Carlo: {abs(result - monte_carlo)/monte_carlo*100:.6f}%")
print(f"Error formula vs Integration: {abs(result - integration)/integration*100:.6f}%")