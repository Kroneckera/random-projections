#!/usr/bin/env python3
"""
Calculate the correct analytical formula for unit square average distance
"""

import numpy as np

# The correct analytical formula for a unit square is:
# (sqrt(2) + 2 + 5*ln(1+sqrt(2)))/15

# Let's calculate it step by step
sqrt2 = np.sqrt(2)
ln_term = np.log(1 + sqrt2)

numerator = sqrt2 + 2 + 5*ln_term
denominator = 15

result = numerator / denominator

print("Calculation of CORRECT analytical formula for unit square average distance:")
print(f"sqrt(2) = {sqrt2}")
print(f"ln(1+sqrt(2)) = {ln_term}")
print(f"Numerator: sqrt(2) + 2 + 5*ln(1+sqrt(2)) = {numerator}")
print(f"Denominator: 15")
print(f"Result: {result}")

# Compare with our Monte Carlo and integration results
monte_carlo = 0.522023
integration = 0.521405

print(f"\nMonte Carlo result: {monte_carlo}")
print(f"Integration result: {integration}")
print(f"Analytical formula: {result}")

print(f"\nError formula vs Monte Carlo: {abs(result - monte_carlo)/monte_carlo*100:.6f}%")
print(f"Error formula vs Integration: {abs(result - integration)/integration*100:.6f}%")