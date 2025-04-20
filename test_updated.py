#!/usr/bin/env python3
"""
Test Monte Carlo vs Analytical calculations with correct formulas
"""

import numpy as np
from polygon_projection import PolygonProjection, PolygonRegion

def main():
    # Set up test parameters
    square = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
    seed = 12345
    n_samples = 500000
    
    print("Comparing Monte Carlo vs Analytical Results (Corrected)")
    print("=" * 50)
    
    # 1. Unit Square Interior
    interior_proj = PolygonProjection(square, [1, 0], region=PolygonRegion.INTERIOR, seed=seed)
    interior_mc = interior_proj.average_distance_monte_carlo(n_samples=n_samples)
    
    # CORRECT analytical formula for unit square interior
    interior_analytical = (2 + 2*np.sqrt(2) + 5*np.log(1+np.sqrt(2)))/15
    
    interior_error = abs(interior_mc - interior_analytical) / interior_analytical
    
    print("\nUnit Square Interior:")
    print(f"Monte Carlo result ({n_samples:,} samples): {interior_mc:.6f}")
    print(f"Analytical formula (correct): {interior_analytical:.6f}")
    print(f"Relative error: {interior_error:.6f} ({interior_error*100:.4f}%)")
    
    # 2. Unit Square Boundary
    boundary_proj = PolygonProjection(square, [1, 0], region=PolygonRegion.BOUNDARY, seed=seed)
    boundary_mc = boundary_proj.average_distance_monte_carlo(n_samples=n_samples)
    
    # Monte Carlo gives consistently ~0.735, which matches literature
    print("\nUnit Square Boundary:")
    print(f"Monte Carlo result ({n_samples:,} samples): {boundary_mc:.6f}")
    
    # 3. Test the integration method too
    print("\nTesting integration method:")
    interior_exact = interior_proj.average_distance_exact(integration_order=24)
    boundary_exact = boundary_proj.average_distance_exact(integration_order=24)
    
    print(f"Interior - Exact Integration: {interior_exact:.6f}, MC: {interior_mc:.6f}")
    print(f"Boundary - Exact Integration: {boundary_exact:.6f}, MC: {boundary_mc:.6f}")
    
    # Calculate errors for integration vs exact formula and Monte Carlo
    int_analytical_error = abs(interior_exact - interior_analytical) / interior_analytical
    int_mc_error = abs(interior_exact - interior_mc) / interior_mc
    
    print(f"\nInterior Integration vs Analytical: {int_analytical_error:.6f} ({int_analytical_error*100:.4f}%)")
    print(f"Interior Integration vs Monte Carlo: {int_mc_error:.6f} ({int_mc_error*100:.4f}%)")

if __name__ == "__main__":
    main()