#ifndef PROJECTION_PROJECTION_H
#define PROJECTION_PROJECTION_H

#include <vector>
#include <array>
#include <random>
#include <functional>
#include "polygon.h"

namespace projection {

/**
 * Enum defining the region to use for projections.
 */
enum class Region {
    INTERIOR,  ///< Interior of the polygon
    BOUNDARY   ///< Boundary of the polygon
};

/**
 * Class for handling projections of a polygon onto a direction.
 */
class ProjectionCalculator {
public:
    using Direction = std::array<double, 2>;
    using ProgressCallback = std::function<void(double)>;
    
    /**
     * Create a projection calculator for a polygon and direction.
     * 
     * @param polygon The polygon to project
     * @param direction The direction to project onto (2D vector, does not need to be normalized)
     * @param region The region to use for calculations (INTERIOR or BOUNDARY)
     * @param seed Optional random seed for reproducibility
     */
    ProjectionCalculator(const Polygon& polygon, 
                         const Direction& direction,
                         Region region = Region::INTERIOR,
                         unsigned int seed = 0);
    
    /**
     * Calculate the exact expected distance between projected points.
     * 
     * @return The exact expected absolute difference
     */
    double exact_projected_distance() const;
    
    /**
     * Calculate the exact average Euclidean distance between pairs of points.
     * 
     * @param integration_order Order of Gauss-Legendre quadrature (higher = more accurate)
     * @return The exact average Euclidean distance
     */
    double exact_average_distance(int integration_order = 16) const;
    
    /**
     * Calculate Monte Carlo estimate of projected distance.
     * 
     * @param n_samples Number of point pairs to sample
     * @return Estimated expected absolute difference
     */
    double monte_carlo_projected_distance(int n_samples = 10000) const;
    
    /**
     * Calculate Monte Carlo estimate of average Euclidean distance.
     * 
     * @param n_samples Number of point pairs to sample
     * @return Estimated average Euclidean distance
     */
    double monte_carlo_average_distance(int n_samples = 10000) const;
    
    /**
     * Run Monte Carlo simulation with progress updates.
     * 
     * @param n_samples Total number of samples to use
     * @param batch_size Size of each batch of samples
     * @param progress_callback Function to call with progress (0.0 to 1.0)
     * @return Results of the Monte Carlo simulation
     */
    struct MonteCarloResults {
        double projected_distance;
        double average_distance;
        int samples_processed;
    };
    
    MonteCarloResults run_monte_carlo_batched(
        int n_samples,
        int batch_size,
        const ProgressCallback& progress_callback = nullptr) const;
    
    /**
     * Sample a point from the specified region of the polygon.
     * 
     * @return A randomly sampled point
     */
    Polygon::Point sample_point() const;
    
    /**
     * Project a point onto the direction vector.
     * 
     * @param point The point to project
     * @return The scalar projection value
     */
    double project_point(const Polygon::Point& point) const;
    
    /**
     * Get the polygon associated with this calculator.
     * 
     * @return Reference to the polygon
     */
    const Polygon& polygon() const { return polygon_; }
    
    /**
     * Get the projection direction.
     * 
     * @return The projection direction vector
     */
    const Direction& direction() const { return direction_; }
    
 private:
    const Polygon& polygon_;
    Direction direction_;
    Direction normalized_direction_;
    Region region_;
    mutable std::mt19937 rng_;
    
    // Internal structures and functions for implementation

public:
    // Represents piecewise-linear density for interior projections
    struct ProjectedDensity {
        // A piece of the density function: on [x0, x1], f(x) = m*x + c
        struct Piece {
            double x0;  // Lower bound
            double x1;  // Upper bound
            double m;   // Slope
            double c;   // Intercept
        };
        
        std::vector<Piece> pieces;
    };
    
    // Represents piecewise-constant density for boundary projections
    using BoundarySteps = std::vector<std::tuple<double, double, double>>;
    
    // Density calculation methods
    ProjectedDensity compute_interior_density() const;
    std::vector<std::tuple<double, double, double>> compute_boundary_steps() const;
    
    // Sampling methods
    Polygon::Point sample_interior_point() const;
    Polygon::Point sample_boundary_point() const;
    
private:
};

} // namespace projection

#endif // PROJECTION_PROJECTION_H