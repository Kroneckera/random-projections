#include "projection.h"
#include "utils.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace projection {

// Constructor implementation
ProjectionCalculator::ProjectionCalculator(
    const Polygon& polygon,
    const Direction& direction,
    Region region,
    unsigned int seed)
    : polygon_(polygon), 
      direction_(direction),
      region_(region),
      rng_(seed)
{
    // Validate polygon
    validate_polygon(polygon);
    
    // Normalize direction vector
    double norm = std::sqrt(direction_[0] * direction_[0] + direction_[1] * direction_[1]);
    if (norm < 1e-10) {
        throw std::invalid_argument("Direction vector cannot be zero");
    }
    normalized_direction_ = {direction_[0] / norm, direction_[1] / norm};
}

// Project a point onto the direction vector
double ProjectionCalculator::project_point(const Polygon::Point& point) const {
    // Project the point onto the normalized direction
    return point[0] * normalized_direction_[0] + 
           point[1] * normalized_direction_[1];
}

// We'll use the sample_interior_point and sample_boundary_point functions
// from utils.h instead of implementing our own here

// Sample a point from the specified region of the polygon
Polygon::Point ProjectionCalculator::sample_point() const {
    if (region_ == Region::INTERIOR) {
        return projection::sample_interior_point(polygon_, rng_);
    } else {
        return projection::sample_boundary_point(polygon_, rng_);
    }
}

// Calculate exact expected projected distance
double ProjectionCalculator::exact_projected_distance() const {
    if (region_ == Region::INTERIOR) {
        return interior_exact_projected_distance(*this);
    } else {
        return boundary_exact_projected_distance(*this);
    }
}

// Calculate Monte Carlo estimate of projected distance
double ProjectionCalculator::monte_carlo_projected_distance(int n_samples) const {
    // Validate inputs
    if (n_samples <= 0) {
        throw std::invalid_argument("Number of samples must be positive");
    }
    
    if (region_ == Region::INTERIOR) {
        return interior_monte_carlo_projected_distance(*this, n_samples);
    } else {
        return boundary_monte_carlo_projected_distance(*this, n_samples);
    }
}

// Validation functions are implemented in sampling.cpp

} // namespace projection