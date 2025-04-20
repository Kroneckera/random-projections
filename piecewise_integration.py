#!/usr/bin/env python3
"""
piecewise_integration.py

High‑accuracy integration of a \pi‑periodic integrand that is
piecewise smooth but may have finitely many kink points (first–derivative
discontinuities) determined by a convex polygon.

We expose a single helper

    expected_distance_interior_exact(vertices, f, order=16)

where
    * vertices  – (n,2) ndarray of CCW polygon vertices
    * f(v,theta) – callable returning the value of the integrand at angle theta
                   for those vertices (example: the exact projected
                   mean absolute distance implemented in
                   projection_sampling.expected_abs_diff).
    * order     – Gauss–Legendre order used inside every smooth interval.

The routine automatically:
  1. enumerates all break angles where projections of two vertices coincide;
  2. subdivides [0,pi] into smooth pieces;
  3. applies an order‑`order` Gauss–Legendre rule on each sub‑interval;
  4. returns the integral as a float.

For the particular application of E||X-Y|| = 0.5 * int_{0}^{pi} E_proj(theta) dtheta
simply multiply the returned integral by 0.5.
"""

import numpy as np
from numpy.polynomial.legendre import leggauss

__all__ = [
    "break_angles",
    "expected_distance_interior_exact",
]


def break_angles(vertices, tol=1e-12):
    """Return sorted unique angles theta in (0,pi) where projection order changes."""
    V = np.asarray(vertices, dtype=float)
    m = len(V)
    S = set()
    for i in range(m):
        for j in range(i + 1, m):
            dx, dy = V[i] - V[j]
            if abs(dx) < tol and abs(dy) < tol:
                continue  # duplicate vertex
            phi = np.arctan2(dy, dx) + np.pi / 2  # perpendicular direction
            phi = (phi + np.pi) % np.pi  # wrap into [0,pi)
            # round to avoid numerically duplicated entries near 0 or pi
            S.add(round(phi, 12))
    return np.sort(list(S))


def expected_distance_interior_exact(vertices, f, order=16):
    """
    Integrate f(vertices, theta) over theta in [0,pi] accounting for kinks.

    Parameters
    ----------
    vertices : ndarray (n,2)
    f        : callable (vertices, theta) -> float
    order    : int, Gauss–Legendre nodes per smooth piece (>=4 recommended)

    Returns
    -------
    integral : float  (\int_0^{pi} f(theta) dtheta)
    """
    # 1. break points
    bks = break_angles(vertices)
    if bks.size == 0:
        bks = np.array([0.0, np.pi])
    else:
        bks = np.concatenate(([0.0], bks, [np.pi]))

    # 2. Gauss nodes & weights on [-1,1]
    xi, wi = leggauss(order)

    total = 0.0
    for a, b in zip(bks[:-1], bks[1:]):
        length = b - a
        if length < 1e-14:
            continue
        # affine map xi -> theta = c0 + c1*xi
        c1 = length / 2
        c0 = (a + b) / 2
        theta = c0 + c1 * xi
        vals = np.array([f(vertices, th) for th in theta])
        total += c1 * np.dot(wi, vals)
    return total


# -----------------------------------------------------
# DEMONSTRATION (requires projection_sampling.py nearby)
# -----------------------------------------------------
if __name__ == "__main__":
    from projection_sampling import projected_density, expected_abs_diff

    def projected_mean_abs(vertices, theta):
        direction = np.array([np.cos(theta), np.sin(theta)])
        pieces = projected_density(vertices, direction)
        return expected_abs_diff(pieces)

    # regular pentagon
    ang = np.linspace(0, 2 * np.pi, 6)[:-1]
    pent = np.stack([np.cos(ang), np.sin(ang)], axis=1)

    integral = expected_distance_interior_exact(pent, projected_mean_abs, order=20)
    E_dist = 0.5 * integral
    print("E[||X-Y||] for pentagon ≈", E_dist)