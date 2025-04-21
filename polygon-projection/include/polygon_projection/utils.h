#ifndef POLYGON_PROJECTION_UTILS_H
#define POLYGON_PROJECTION_UTILS_H

#include "polygon.h"
#include "projection.h"

namespace polygon_projection {

// ===== Polygon Generation and Geometric Utilities =====

/**
 * Generate a regular polygon with given number of sides and radius.
 * 
 * @param n_sides Number of sides (must be >= 3)
 * @param radius Radius of the circumscribed circle
 * @return List of polygon vertices
 */
Polygon::VertexList generate_regular_polygon(int n_sides, double radius = 1.0);

/**
 * Check if a polygon is convex.
 * 
 * Determines if vertices form a simple (non-self-intersecting),
 * weakly convex polygon. Collinear edges are allowed; reflex angles are not.
 * 
 * @param vertices List of polygon vertices in order (clockwise or counter-clockwise)
 * @param check_ccw If true, requires vertices to be in counter-clockwise order
 * @param eps Tolerance for floating-point comparisons (optional, default=1e-12)
 * @return True if the polygon is convex (and CCW if check_ccw is true)
 */
bool is_convex(const Polygon::VertexList& vertices, bool check_ccw = false, double eps = 1e-12);

/**
 * Calculate the area of a polygon.
 * 
 * @param vertices List of polygon vertices in CCW order
 * @return Area of the polygon
 */
double calculate_polygon_area(const Polygon::VertexList& vertices);

/**
 * Calculate the perimeter of a polygon.
 * 
 * @param vertices List of polygon vertices in CCW order
 * @return Perimeter of the polygon
 */
double calculate_polygon_perimeter(const Polygon::VertexList& vertices);

/**
 * Normalize a direction vector.
 * 
 * @param direction Direction vector to normalize
 * @return Normalized direction vector
 */
ProjectionCalculator::Direction normalize_direction(const ProjectionCalculator::Direction& direction);

/**
 * Generate a random convex polygon.
 * 
 * @param n Number of vertices
 * @param rmin Minimum radius for sampling points
 * @param rmax Maximum radius for sampling points
 * @param rng Random number generator
 * @return List of polygon vertices
 */
Polygon::VertexList generate_random_convex_polygon(
    int n, double rmin, double rmax, std::mt19937& rng);

// ===== Projection and Monte Carlo Calculation Functions =====

/**
 * Calculate the exact expected projected distance for interior points.
 * 
 * @param calc The projection calculator instance
 * @return The exact expected distance
 */
double interior_exact_projected_distance(const ProjectionCalculator& calc);

/**
 * Calculate the exact expected projected distance for boundary points.
 * 
 * @param calc The projection calculator instance
 * @return The exact expected distance
 */
double boundary_exact_projected_distance(const ProjectionCalculator& calc);

/**
 * Calculate Monte Carlo estimate of projected distance for interior points.
 * 
 * @param calc The projection calculator instance
 * @param n_samples Number of point pairs to sample
 * @return The estimated expected distance
 */
double interior_monte_carlo_projected_distance(const ProjectionCalculator& calc, int n_samples);

/**
 * Calculate Monte Carlo estimate of projected distance for boundary points.
 * 
 * @param calc The projection calculator instance
 * @param n_samples Number of point pairs to sample
 * @return The estimated expected distance
 */
double boundary_monte_carlo_projected_distance(const ProjectionCalculator& calc, int n_samples);

// ===== Sampling Functions =====

/**
 * Sample a point from the interior of a polygon.
 * 
 * @param polygon The polygon to sample from
 * @param rng Random number generator
 * @return A randomly sampled point
 */
Polygon::Point sample_interior_point(const Polygon& polygon, std::mt19937& rng);

/**
 * Sample a point from the boundary of a polygon.
 * 
 * @param polygon The polygon to sample from
 * @param rng Random number generator
 * @return A randomly sampled point
 */
Polygon::Point sample_boundary_point(const Polygon& polygon, std::mt19937& rng);

// ===== Validation Functions =====

/**
 * Validate inputs for Monte Carlo calculations.
 * 
 * @param n_samples Number of samples (must be positive)
 * @param batch_size Batch size (must be positive and not larger than n_samples)
 * @throws std::invalid_argument if the inputs are invalid
 */
void validate_monte_carlo_inputs(int n_samples, int batch_size);

/**
 * Validate a polygon for projection calculations.
 * 
 * @param polygon The polygon to validate
 * @throws std::invalid_argument if the polygon is invalid
 */
void validate_polygon(const Polygon& polygon);

} // namespace polygon_projection

#endif // POLYGON_PROJECTION_UTILS_H