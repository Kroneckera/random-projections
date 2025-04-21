#include "polygon_projection/projection.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <numeric>
#include <map>
#include <stdexcept>
#include <iostream>  // For debugging output

namespace polygon_projection {

// Project a single edge onto the direction
// Returns (lo, hi, mass) instead of density for better numerical behavior
std::tuple<double, double, double> project_edge(
    const Polygon::Point& v0, 
    const Polygon::Point& v1, 
    const ProjectionCalculator::Direction& direction,
    double tol = 1e-9)
{
    // Normalize direction
    double norm = std::sqrt(direction[0] * direction[0] + direction[1] * direction[1]);
    double ux = direction[0] / norm;
    double uy = direction[1] / norm;
    
    // Project endpoints
    double x0 = v0[0] * ux + v0[1] * uy;
    double x1 = v1[0] * ux + v1[1] * uy;
    
    // If they project to the same point, we need special handling
    // Return the point with its mass (edge length) directly
    if (std::abs(x0 - x1) < tol) {
        // Calculate edge length
        double dx = v1[0] - v0[0];
        double dy = v1[1] - v0[1];
        double edge_length = std::sqrt(dx * dx + dy * dy);
        
        // Return point projection with its mass
        return {x0, x0, edge_length}; // When lo==hi, we interpret the third value as mass
    }
    
    // Calculate bounds
    double lo = std::min(x0, x1);
    double hi = std::max(x0, x1);
    
    // Calculate edge length (the mass)
    double dx = v1[0] - v0[0];
    double dy = v1[1] - v0[1];
    double edge_length = std::sqrt(dx * dx + dy * dy);
    
    // Return the interval bounds and the mass directly
    return {lo, hi, edge_length};
}

// Implementation of compute_boundary_steps from the ProjectionCalculator class
// Returns a vector of tuples (lo, hi, normalized_mass)
std::vector<std::tuple<double, double, double>> 
ProjectionCalculator::compute_boundary_steps() const {
    const auto& vertices = polygon_.vertices();
    const auto& direction = normalized_direction_;
    
    std::vector<std::tuple<double, double, double>> raw;
    const size_t n = vertices.size();
    
    for (size_t i = 0; i < n; ++i) {
        auto proj = project_edge(vertices[i], vertices[(i+1) % n], direction, 1e-9);
        
        // The third value is now the mass, not density
        double lo = std::get<0>(proj);
        double hi = std::get<1>(proj);
        double mass = std::get<2>(proj);
        
        // Check if it's a valid projection with non-zero mass
        if (mass > 0.0) {
            raw.emplace_back(lo, hi, mass);
        }
    }
    
    // Calculate total mass (sum of edge lengths)
    double total_mass = 0.0;
    for (const auto& [lo, hi, mass] : raw) {
        total_mass += mass;
    }
    
    if (total_mass <= 0.0) {
        throw std::runtime_error("Polygon perimeter must be positive.");
    }
    
    // Normalize masses (divide by total perimeter)
    std::vector<std::tuple<double, double, double>> steps;
    for (const auto& [lo, hi, mass] : raw) {
        steps.emplace_back(lo, hi, mass / total_mass);
    }
    
    return steps;
}

// Combine steps by decomposing into elementary intervals and point masses
std::vector<std::tuple<double, double, double>> combine_steps(
    const std::vector<std::tuple<double, double, double>>& steps)
{
    // Collect all unique boundary points
    std::vector<double> points;
    for (const auto& [lo, hi, _] : steps) {
        points.push_back(lo);
        points.push_back(hi);
    }
    
    // Sort and remove duplicates with a small tolerance
    std::sort(points.begin(), points.end());
    
    // Use tolerance-based deduplication
    const double tol = 1e-9;
    std::vector<double> unique_points;
    for (double p : points) {
        if (unique_points.empty() || p - unique_points.back() > tol) {
            unique_points.push_back(p);
        }
    }
    
    const size_t N = unique_points.size();
    
    // Initialize arrays for point masses and interval masses
    std::vector<double> point_masses(N, 0.0);       // Mass at each point
    std::vector<double> interval_masses(N-1, 0.0);  // Mass for each interval between adjacent points
    
    // Process each step to update the elementary masses
    for (const auto& [lo, hi, mass] : steps) {
        // Find indices of lo and hi in unique_points
        int lo_idx = -1;
        int hi_idx = -1;
        
        for (size_t i = 0; i < unique_points.size(); ++i) {
            if (std::abs(unique_points[i] - lo) < tol) {
                lo_idx = static_cast<int>(i);
            }
            if (std::abs(unique_points[i] - hi) < tol) {
                hi_idx = static_cast<int>(i);
            }
        }
        
        if (lo_idx == -1 || hi_idx == -1) {
            throw std::runtime_error("Could not find step boundaries in unique points");
        }
        
        // If this is a point mass (lo == hi), add to the point mass array
        if (lo_idx == hi_idx) {
            point_masses[lo_idx] += mass;
            continue;
        }
        
        // For an interval, distribute mass proportionally across elementary intervals
        double total_length = unique_points[hi_idx] - unique_points[lo_idx];
        for (int i = lo_idx; i < hi_idx; ++i) {
            double interval_length = unique_points[i+1] - unique_points[i];
            double fraction = interval_length / total_length;
            interval_masses[i] += mass * fraction;
        }
    }
    
    // Create the combined result with all non-zero masses
    std::vector<std::tuple<double, double, double>> combined;
    
    // Add all non-zero masses (both point and interval) in sorted order
    for (size_t i = 0; i < N; ++i) {
        // Add point mass at this point if non-zero
        if (point_masses[i] > tol) {
            combined.emplace_back(unique_points[i], unique_points[i], point_masses[i]);
        }
        
        // Add interval mass starting at this point if non-zero and not the last point
        if (i < N-1 && interval_masses[i] > tol) {
            combined.emplace_back(unique_points[i], unique_points[i+1], interval_masses[i]);
        }
    }
    
    return combined;
}

// Calculate boundary expected absolute difference
// Uses the formula 2 * integral of F(x) * (1 - F(x)) dx
// where F(x) is the cumulative distribution function
double boundary_expected_abs_diff(
    const std::vector<std::tuple<double, double, double>>& steps)
{
    const double tol = 1e-9;
    
    // First decompose steps into elementary intervals
    std::vector<std::tuple<double, double, double>> pieces = combine_steps(steps);
    
    // Verify the total mass is normalized to 1.0
    double total_mass = 0.0;
    for (const auto& [x0, x1, mass] : pieces) {
        total_mass += mass;
    }
    
    if (std::abs(total_mass - 1.0) > 1e-6) {
        throw std::runtime_error("Total mass should be 1.0");
    }
    
    // Calculate the CDF points by accumulating normalized masses
    // For each interval, store (x0, x1, F_start, F_end)
    std::vector<std::tuple<double, double, double, double>> cdf_pieces;
    double F0 = 0.0;  // Cumulative probability at current point
    
    for (const auto& [x0, x1, mass] : pieces) {
        // For point masses, update the CDF but don't create an interval piece
        if (std::abs(x1 - x0) < tol) {
            // Just update the CDF value
            F0 += mass;
            continue;
        }
        
        // For intervals, create a CDF piece with start and end values
        // When we integrate from x0 to x1, the CDF increases from F0 to F0+mass
        double F1 = F0 + mass;  // CDF value at the end of this interval
        cdf_pieces.emplace_back(x0, x1, F0, F1);
        
        // Update the cumulative probability
        F0 = F1;
    }
    
    // Integrate F(x) * (1 - F(x)) over each interval
    double integral = 0.0;
    
    for (const auto& [x0, x1, F_start, F_end] : cdf_pieces) {
        // Skip zero-width intervals
        if (std::abs(x1 - x0) < tol) {
            continue;
        }
        
        // Calculate parameters for this interval
        double interval_length = x1 - x0;
        
        // For this interval, the CDF is F(x) = F_start + (x - x0) * slope
        // where slope = (F_end - F_start) / interval_length
        double slope = (F_end - F_start) / interval_length;
        
        // Integrate F(x) * (1 - F(x)) over [x0, x1]
        double A = F_start * (1.0 - F_start);          // Term: F_start - F_start²
        double B = (1.0 - 2.0 * F_start) * slope;      // Term: (1 - 2*F_start) * slope
        double C = -slope * slope;                     // Term: -slope²
        
        // ∫ A dx = A * (x1 - x0)
        double term1 = A * interval_length;
        
        // ∫ B * (x-x0) dx = B * (x1-x0)²/2
        double term2 = B * std::pow(interval_length, 2) / 2.0;
        
        // ∫ C * (x-x0)² dx = C * (x1-x0)³/3
        double term3 = C * std::pow(interval_length, 3) / 3.0;
        
        double interval_integral = term1 + term2 + term3;
        integral += interval_integral;
    }
    
    // The final result is 2 * integral
    return 2.0 * integral;
}

// Implement sample_boundary_point from ProjectionCalculator
Polygon::Point ProjectionCalculator::sample_boundary_point() const {
    const auto& vertices = polygon_.vertices();
    const size_t n = vertices.size();
    
    // Compute edge lengths
    std::vector<double> lengths;
    lengths.reserve(n);
    
    for (size_t i = 0; i < n; ++i) {
        const auto& p0 = vertices[i];
        const auto& p1 = vertices[(i + 1) % n];
        
        double dx = p1[0] - p0[0];
        double dy = p1[1] - p0[1];
        lengths.push_back(std::sqrt(dx * dx + dy * dy));
    }
    
    // Compute cumulative lengths
    std::vector<double> cum_lengths(n);
    std::partial_sum(lengths.begin(), lengths.end(), cum_lengths.begin());
    double total_length = cum_lengths.back();
    
    // Select random position along perimeter
    std::uniform_real_distribution<double> dist(0.0, total_length);
    double u = dist(rng_);
    
    // Find which edge this position falls on
    size_t edge_idx = 0;
    while (edge_idx < cum_lengths.size() && u > cum_lengths[edge_idx]) {
        edge_idx++;
    }
    
    // Calculate distance along the edge
    double prev = (edge_idx > 0) ? cum_lengths[edge_idx - 1] : 0.0;
    double t = (u - prev) / lengths[edge_idx];
    
    // Get edge endpoints
    const auto& p0 = vertices[edge_idx];
    const auto& p1 = vertices[(edge_idx + 1) % n];
    
    // Interpolate
    Polygon::Point result = {
        p0[0] + t * (p1[0] - p0[0]),
        p0[1] + t * (p1[1] - p0[1])
    };
    
    return result;
}

// Helper function for boundary Monte Carlo
double boundary_monte_carlo_projected_distance(const ProjectionCalculator& calc, int n_samples) {
    double sum = 0.0;
    
    for (int i = 0; i < n_samples; ++i) {
        Polygon::Point p1 = calc.sample_boundary_point();
        Polygon::Point p2 = calc.sample_boundary_point();
        
        double proj1 = calc.project_point(p1);
        double proj2 = calc.project_point(p2);
        
        sum += std::abs(proj1 - proj2);
    }
    
    return sum / n_samples;
}

// Forward declare the interior expected_abs_diff function
double expected_abs_diff(const ProjectionCalculator::ProjectedDensity& density);

// Helper function for boundary exact calculation
double boundary_exact_projected_distance(const ProjectionCalculator& calc) {
    // Compute boundary steps with mass-based representation
    auto steps = calc.compute_boundary_steps();
    
    // Calculate expected absolute difference using steps
    return boundary_expected_abs_diff(steps);
}

} // namespace polygon_projection