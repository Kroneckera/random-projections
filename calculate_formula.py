#!/usr/bin/env python3
"""
Calculate the analytical formula for unit square average distance
"""

import numpy as np

# The analytical formula for a unit square is:
# (2 + 2*sqrt(2) + 5*ln(1+sqrt(2)))/15

# Let's break it down step by step
sqrt2 = np.sqrt(2)
ln_term = np.log(1 + sqrt2)

numerator = 2 + 2*sqrt2 + 5*ln_term
denominator = 15

result = numerator / denominator

print("Calculation of analytical formula for unit square average distance:")
print(f"sqrt(2) = {sqrt2}")
print(f"ln(1+sqrt(2)) = {ln_term}")
print(f"Numerator: 2 + 2*sqrt(2) + 5*ln(1+sqrt(2)) = {numerator}")
print(f"Denominator: 15")
print(f"Result: {result}")

# Let's try another formula version
result2 = (2 + sqrt2 + 5*ln_term) / 6
print(f"\nAlternative formula (1/6)*(2 + sqrt(2) + 5*ln(1+sqrt(2))): {result2}")

# Let's try a numerical approximation of the formula
approx = 0.521405
print(f"\nNumerical approximation from literature: 0.521405")
print(f"Error from first formula: {abs(result - approx)/approx*100:.6f}%")
print(f"Error from second formula: {abs(result2 - approx)/approx*100:.6f}%")