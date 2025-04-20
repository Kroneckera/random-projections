#!/usr/bin/env python3
"""
projection_sampling.py

Compute the projection of a uniform distribution on a convex polygon
onto a given direction, sample from that projected density, and compute
exact and Monte Carlo estimates of pairwise distances.

All code uses ASCII only (no unicode symbols).
"""

import numpy as np


def projected_density(vertices, direction, tol=1e-12):
    """
    Compute piecewise-linear density f(x) of the projection of the uniform
    distribution on a convex polygon onto a given direction.

    vertices: array_like, shape (n,2)
        CCW-ordered polygon vertices.
    direction: array_like, shape (2,)
        Projection direction (need not be unit; will be normalized).
    tol: float
        Threshold below which dx is treated as zero (skip vertical edges).

    Returns:
        pieces: list of tuples (x0, x1, m, c)
            On each [x0, x1], f(x) = m*x + c.
    """
    V = np.asarray(vertices, dtype=float)
    u = np.asarray(direction, dtype=float)
    u = u / np.linalg.norm(u)

    # 1) rotate so u aligns with the x-axis
    phi = np.arctan2(u[1], u[0])
    R = np.array([[np.cos(phi), np.sin(phi)],
                  [-np.sin(phi), np.cos(phi)]])
    V2 = (R @ V.T).T
    xs = V2[:, 0]
    ys = V2[:, 1]

    # 2) compute polygon area via shoelace formula
    A = 0.5 * np.sum(xs * np.roll(ys, -1) - np.roll(xs, -1) * ys)
    if A <= 0:
        raise ValueError("Polygon area must be positive; check vertex order.")

    # 3) collect each edge's contribution
    events = []
    n = len(V2)
    for i in range(n):
        x0, y0 = xs[i], ys[i]
        x1, y1 = xs[(i+1) % n], ys[(i+1) % n]
        dx, dy = x1 - x0, y1 - y0
        if abs(dx) < tol:
            continue  # vertical edge, no continuous contribution
        k = dy / dx
        b = y0 - k * x0
        sign = -np.sign(dx)  # +1 if dx<0, -1 if dx>0
        lo = x0 if x0 < x1 else x1
        hi = x1 if x1 > x0 else x0
        events.append((lo, hi, sign*k, sign*b))

    # 4) build piecewise sum over sorted breakpoints
    breaks = np.sort(np.unique(xs))
    pieces = []
    for lo, hi in zip(breaks[:-1], breaks[1:]):
        a_sum = 0.0
        b_sum = 0.0
        for (l, h, ak, bk) in events:
            if l <= lo < h:
                a_sum += ak
                b_sum += bk
        pieces.append((lo, hi, a_sum/A, b_sum/A))

    return pieces


def eval_density(x, pieces):
    """
    Evaluate the piecewise-linear density at x (scalar or array).
    """
    x_arr = np.asarray(x, dtype=float)
    f = np.zeros_like(x_arr)
    for (lo, hi, m, c) in pieces:
        mask = (x_arr >= lo) & (x_arr <= hi)
        f[mask] = m * x_arr[mask] + c
    if f.shape == ():
        return float(f)
    return f


def sample_point_in_convex_polygon(vertices, rng=None):
    """
    Sample a single point uniformly from a convex polygon.

    vertices: array_like, shape (n,2), CCW-ordered
    rng:       instance of np.random.Generator (optional)
    """
    if rng is None:
        rng = np.random.default_rng()

    V = np.asarray(vertices, dtype=float)
    n = len(V)
    if n < 3:
        raise ValueError("Need at least 3 vertices")

    # Build triangles (fan from V[0]) of shape (n-2, 3, 2)
    tris = np.stack([np.repeat(V[0][None,:], n-2, axis=0),
                     V[1:-1],
                     V[2:]], axis=1)

    # Compute areas and select one triangle
    v0 = tris[:, 0, :]
    v1 = tris[:, 1, :]
    v2 = tris[:, 2, :]
    cross = (v1 - v0)[:, 0] * (v2 - v0)[:, 1] - (v1 - v0)[:, 1] * (v2 - v0)[:, 0]
    areas = 0.5 * np.abs(cross)
    cum_areas = np.cumsum(areas)
    total = cum_areas[-1]
    u = rng.random() * total
    tri_idx = np.searchsorted(cum_areas, u)
    p0, p1, p2 = tris[tri_idx]

    # Sample uniformly in that triangle via barycentric coords
    r1, r2 = rng.random(), rng.random()
    if r1 + r2 > 1:
        r1, r2 = 1 - r1, 1 - r2
    return p0 + r1 * (p1 - p0) + r2 * (p2 - p0)


def sample_two_points_in_polygon(vertices, rng=None):
    """
    Sample two independent points uniformly from the polygon.
    Returns: (pt1, pt2), each is a length-2 array.
    """
    p1 = sample_point_in_convex_polygon(vertices, rng)
    p2 = sample_point_in_convex_polygon(vertices, rng)
    return p1, p2


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
        p1, p2 = sample_two_points_in_polygon(vertices, rng)
        distances[i] = abs(np.dot(p1, u) - np.dot(p2, u))
    return distances if n > 1 else float(distances[0])


def expected_abs_diff(pieces):
    """
    Compute E[|X - Y|] exactly for X,Y iid with piecewise-linear density
    defined by `pieces = [(x0,x1,m,c), …]` on each [x0,x1], f(x)=m*x+c.
    """
    # 1) accumulate F0 at each segment start
    segs = []
    F0 = 0.0
    for (x0, x1, m, c) in pieces:
        segs.append({'x0': x0, 'x1': x1, 'm': m, 'c': c, 'F0': F0})
        mass = 0.5*m*(x1**2 - x0**2) + c*(x1 - x0)
        F0 += mass
    if abs(F0 - 1.0) > 1e-8:
        raise ValueError(f"Total mass = {F0}, but should be 1")

    # 2) integrate 2*∫ F(x)*(1-F(x)) dx piecewise
    total = 0.0
    for s in segs:
        x0, x1 = s['x0'], s['x1']
        m, c, Fstart = s['m'], s['c'], s['F0']

        # F(x) = alpha*x^2 + beta*x + gamma
        alpha = 0.5 * m
        beta  = c
        gamma = Fstart - (alpha*x0**2 + beta*x0)

        # Expand g(x) = F(x)*(1 - F(x)) = -F(x)^2 + F(x)
        #   = p4*x^4 + p3*x^3 + p2*x^2 + p1*x + p0
        p4 = -alpha * alpha
        p3 = -2*alpha * beta
        p2 = -(2*alpha*gamma + beta*beta) + alpha
        p1 = -2*beta*gamma + beta
        p0 = -gamma*gamma + gamma

        # ∫ x^n dx = (x1^n - x0^n)/n
        def I(n): 
            return (x1**n - x0**n) / n

        total += (
            p4 * I(5) +
            p3 * I(4) +
            p2 * I(3) +
            p1 * I(2) +
            p0 * I(1)
        )

    return 2 * total

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
    # Example usage: a more nontrivial random convex polygon
    rng = np.random.default_rng(1234)
    poly = generate_random_convex_polygon(10, rmin=0.5, rmax=2.0, rng=rng)
    direction = np.array([2.1, 2.0])

    pieces = projected_density(poly, direction)
    print("Pieces:", pieces)

    # exact expected distance
    E_exact = expected_abs_diff(pieces)
    print("Exact E[|X - Y|]:", E_exact)

    # Monte Carlo estimate
    dists = sample_mc_pairwise_projected_distances(poly, direction, n=1000000, rng=rng)
    print("Monte Carlo E[|X - Y|]:", dists.mean())
