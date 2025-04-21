#!/usr/bin/env python3
"""
edge_projection.py

Compute the projection of each edge of a convex polygon onto a given direction.
Each edge contributes a piecewise-constant "step" density on its projected interval.
The total (normalized) boundary projection density is the sum of these steps.
Compute the exact expected distance E[|X-Y|] under this boundary projection.

Also supports sampling uniformly random points from the polygon boundary.

All code uses ASCII only.
"""

import numpy as np


def project_edge(v0, v1, direction):
    """
    Project a single edge [v0, v1] onto 'direction'.

    v0, v1: array_like, shape (2,)
        Endpoints of the edge.
    direction: array_like, shape (2,)
        The projection direction (need not be unit).

    Returns:
        lo: float, lower scalar bound of projection interval
        hi: float, upper scalar bound of projection interval
        density: float, constant density on [lo, hi] from this edge
    """
    p0 = np.asarray(v0, dtype=float)
    p1 = np.asarray(v1, dtype=float)
    u = np.asarray(direction, dtype=float)
    u = u / np.linalg.norm(u)

    x0 = float(np.dot(p0, u))
    x1 = float(np.dot(p1, u))
    if x0 == x1:
        return None

    lo, hi = (x0, x1) if x0 < x1 else (x1, x0)
    edge_length = np.linalg.norm(p1 - p0)
    density = edge_length / (hi - lo)
    return lo, hi, density


def project_edges(vertices, direction):
    """
    Project all edges of a polygon, returning normalized step densities.

    vertices: array_like, shape (n,2)
        CCW-ordered polygon vertices.
    direction: array_like, shape (2,)
        Projection direction.

    Returns:
        steps: list of (lo, hi, c_norm)
            Step density c_norm on [lo, hi]; sum integrates to 1.
    """
    V = np.asarray(vertices, dtype=float)
    n = len(V)
    raw = []
    for i in range(n):
        proj = project_edge(V[i], V[(i+1) % n], direction)
        if proj is not None:
            raw.append(proj)
    total_length = sum((hi - lo) * d for lo, hi, d in raw)
    if total_length <= 0:
        raise ValueError("Polygon perimeter must be positive.")
    steps = [(lo, hi, d / total_length) for lo, hi, d in raw]
    return steps


def combine_steps(steps):
    """
    Combine overlapping step contributions into disjoint constant segments.

    steps: list of (lo, hi, c)
    Returns:
        combined: list of (x0, x1, f_const)
    """
    pts = []
    for lo, hi, _ in steps:
        pts.append(lo)
        pts.append(hi)
    xs = np.sort(np.unique(pts))
    combined = []
    for x0, x1 in zip(xs[:-1], xs[1:]):
        f = sum(c for lo, hi, c in steps if lo <= x0 < hi)
        combined.append((x0, x1, f))
    return combined


def boundary_expected_abs_diff(steps):
    """
    Compute E[|X - Y|] exactly for X,Y iid with step density from boundary.

    steps: list of (lo, hi, c_norm) from project_edges.
    Returns:
        E: float, expected absolute difference.
    """
    pieces = combine_steps(steps)
    total = 0.0
    F0 = 0.0
    for x0, x1, c in pieces:
        L = x1 - x0
        A = F0*(1 - F0)
        B = c*(1 - 2*F0)
        D = -c*c
        total += A*L + B*(L*L)/2 + D*(L*L*L)/3
        F0 += c*L
    return 2 * total


def sample_point_on_boundary(vertices, rng=None):
    """
    Sample a point uniformly at random from the boundary of a convex polygon.

    vertices: array_like, shape (n,2)
        CCW-ordered polygon vertices.
    rng: np.random.Generator or None
    Returns:
        point: ndarray of shape (2,), sampled boundary point.
    """
    if rng is None:
        rng = np.random.default_rng()
    V = np.asarray(vertices, dtype=float)
    # compute edge lengths and cumulative
    lengths = []
    for i in range(len(V)):
        p0 = V[i]
        p1 = V[(i+1) % len(V)]
        lengths.append(np.linalg.norm(p1 - p0))
    lengths = np.array(lengths)
    cum = np.cumsum(lengths)
    total = cum[-1]
    # select random position along perimeter
    u = rng.random() * total
    edge_idx = np.searchsorted(cum, u)
    # find distance along that edge
    prev = cum[edge_idx-1] if edge_idx > 0 else 0.0
    t = (u - prev) / lengths[edge_idx]
    # interpolate
    p0 = V[edge_idx]
    p1 = V[(edge_idx+1) % len(V)]
    return p0 + t * (p1 - p0)

def project_point_to_direction(point, direction):
    """
    Project a 2D point onto a direction, returning the scalar coordinate.
    """
    p = np.asarray(point, dtype=float)
    u = np.asarray(direction, dtype=float)
    u = u / np.linalg.norm(u)
    return float(np.dot(p, u))


def projected_distance(pt1, pt2, direction):
    """
    Compute the absolute difference of projections of pt1 and pt2 onto direction.
    """
    x1 = project_point_to_direction(pt1, direction)
    x2 = project_point_to_direction(pt2, direction)
    return abs(x1 - x2)


def sample_mc_pairwise_projected_distances(vertices, direction, n=1, rng=None):
    """
    Monte Carlo: sample n distances |proj(p1)-proj(p2)| from random points in polygon.
    """
    if rng is None:
        rng = np.random.default_rng()
    distances = np.empty(n)
    u = np.asarray(direction, dtype=float)
    u = u / np.linalg.norm(u)
    for i in range(n):
        p1 = sample_point_on_boundary(vertices, rng)
        p2 = sample_point_on_boundary(vertices, rng)
        distances[i] = abs(np.dot(p1, u) - np.dot(p2, u))
    return distances if n > 1 else float(distances[0])

def convex_hull(points):
    """
    Compute the convex hull of a set of 2D points using the monotone chain algorithm.
    Returns hull vertices in CCW order (no duplicate start/end).
    """
    P = sorted({(x, y) for x, y in points})
    if len(P) <= 1:
        return np.array(P)
    # build lower hull
    lower = []
    for x, y in P:
        while len(lower) >= 2:
            x0, y0 = lower[-2]
            x1, y1 = lower[-1]
            if (x1 - x0)*(y - y0) - (y1 - y0)*(x - x0) <= 0:
                lower.pop()
            else:
                break
        lower.append((x, y))
    # build upper hull
    upper = []
    for x, y in reversed(P):
        while len(upper) >= 2:
            x0, y0 = upper[-2]
            x1, y1 = upper[-1]
            if (x1 - x0)*(y - y0) - (y1 - y0)*(x - x0) <= 0:
                upper.pop()
            else:
                break
        upper.append((x, y))
    # concatenate lower+upper, removing duplicate endpoints
    hull = lower[:-1] + upper[:-1]
    return np.array(hull, dtype=float)


def generate_random_convex_polygon(n, rmin=0.5, rmax=2.0, rng=None):
    """
    Generate a random convex polygon with n vertices in CCW order.
    Does so by sampling more points in an annulus, computing their convex hull,
    and then selecting n evenly spaced hull vertices.
    """
    if rng is None:
        rng = np.random.default_rng()
    # sample extra random points in annulus
    m = max(n*5, 50)
    thetas = rng.random(m) * 2 * np.pi
    rs = rng.uniform(rmin, rmax, size=m)
    pts = np.column_stack([rs * np.cos(thetas), rs * np.sin(thetas)])
    hull = convex_hull(pts)
    # if hull too small, resample recursively
    if len(hull) < n:
        return generate_random_convex_polygon(n, rmin, rmax, rng)
    # select n vertices evenly around hull
    idx = np.linspace(0, len(hull), n, endpoint=False, dtype=int)
    poly = hull[idx]
    return poly

if __name__ == "__main__":
    # Example: square boundary projection and sampling
    square = np.array([[0.0,0.0],[1.0,0.0],[1.0,1.0],[0.0,1.0]])
    direction = [1.0, 3.0]
    steps = project_edges(square, direction)
    print("Steps (lo, hi, c_norm):")
    for lo, hi, c in steps:
        print(f"[{lo:.3f}, {hi:.3f}]: {c:.3f}")
    E = boundary_expected_abs_diff(steps)
    print("Exact boundary E[|X - Y|]:", E)
    # sample some boundary points
    rng = np.random.default_rng(42)
    dists = sample_mc_pairwise_projected_distances(square, direction, n=1000000, rng=rng)
    print("Monte Carlo E[|X - Y|]:", dists.mean())
