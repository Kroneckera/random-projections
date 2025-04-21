#include "projection.h"
#include "integration.h"
#include "utils.h"
#include <cmath>
#include <algorithm>
#include <set>
#include <vector>
#include <stdexcept>

namespace projection {

std::vector<double> calculate_break_angles(const Polygon::VertexList& vertices, double tol) {
    const size_t m = vertices.size();
    std::set<double> breakpoints;
    
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = i + 1; j < m; ++j) {
            double dx = vertices[i][0] - vertices[j][0];
            double dy = vertices[i][1] - vertices[j][1];
            
            // Skip duplicate vertices
            if (std::abs(dx) < tol && std::abs(dy) < tol) {
                continue;
            }
            
            // Calculate perpendicular direction
            double phi = std::atan2(dy, dx) + M_PI / 2.0;
            
            // Wrap into [0, π)
            phi = std::fmod(phi + M_PI, M_PI);
            
            // Round to avoid numerical duplicates
            double rounded_phi = std::round(phi * 1e12) / 1e12;
            
            breakpoints.insert(rounded_phi);
        }
    }
    
    // Convert to sorted vector
    std::vector<double> result(breakpoints.begin(), breakpoints.end());
    std::sort(result.begin(), result.end());
    
    return result;
}

// Helper function for Gauss-Legendre quadrature
void gauss_legendre_coefficients(int n, std::vector<double>& x, std::vector<double>& w) {
    // This is a simple implementation for n <= 5, using precomputed values
    
    if (n <= 0) {
        throw std::invalid_argument("Order must be positive");
    }
    
    if (n > 16) {
        throw std::invalid_argument("Quadrature order cannot exceed 16");
    }
    
    if (n == 1) {
        x = {0.0};
        w = {2.0};
    }
    else if (n == 2) {
        x = {-0.5773502691896257, 0.5773502691896257};
        w = {1.0, 1.0};
    }
    else if (n == 3) {
        x = {-0.7745966692414834, 0.0, 0.7745966692414834};
        w = {0.5555555555555556, 0.8888888888888888, 0.5555555555555556};
    }
    else if (n == 4) {
        x = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
        w = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538};
    }
    else if (n == 5) {
        x = {-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640};
        w = {0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891};
    }
    else if (n == 16) {
        x = {
            -0.9894009349916499, -0.9445750230732326, -0.8656312023878318, -0.7554044083550030,
            -0.6178762444026438, -0.4580167776572274, -0.2816035507792589, -0.0950125098376374,
             0.0950125098376374,  0.2816035507792589,  0.4580167776572274,  0.6178762444026438,
             0.7554044083550030,  0.8656312023878318,  0.9445750230732326,  0.9894009349916499
        };
        w = {
            0.0271524594117541, 0.0622535239386479, 0.0951585116824928, 0.1246289712555339,
            0.1495959888165767, 0.1691565193950025, 0.1826034150449236, 0.1894506104550685,
            0.1894506104550685, 0.1826034150449236, 0.1691565193950025, 0.1495959888165767,
            0.1246289712555339, 0.0951585116824928, 0.0622535239386479, 0.0271524594117541
        };
    }
    else {
        throw std::invalid_argument("Unsupported quadrature order; use 1-5 or 16");
    }
}

double integrate_over_directions(
    const Polygon::VertexList& vertices,
    const std::function<double(const Polygon::VertexList&, double)>& f,
    int order)
{
    // Step 1: Calculate breakpoints
    auto bks = calculate_break_angles(vertices);
    
    // Ensure we have endpoints at 0 and π
    if (bks.empty()) {
        bks = {0.0, M_PI};
    } else {
        if (bks.front() > 0.0) {
            bks.insert(bks.begin(), 0.0);
        }
        if (bks.back() < M_PI) {
            bks.push_back(M_PI);
        }
    }
    
    // Step 2: Get Gauss-Legendre nodes and weights
    std::vector<double> xi, wi;
    gauss_legendre_coefficients(order, xi, wi);
    
    // Step 3: Integrate over each subinterval
    double total = 0.0;
    
    for (size_t i = 0; i < bks.size() - 1; ++i) {
        double a = bks[i];
        double b = bks[i + 1];
        double length = b - a;
        
        if (length < 1e-14) {
            continue;  // Skip tiny intervals
        }
        
        // Affine map from [-1, 1] to [a, b]
        double c0 = (a + b) / 2.0;  // Midpoint
        double c1 = length / 2.0;   // Half-length
        
        // Evaluate at Gauss points
        double interval_sum = 0.0;
        for (size_t j = 0; j < xi.size(); ++j) {
            double theta = c0 + c1 * xi[j];
            double value = f(vertices, theta);
            interval_sum += wi[j] * value;
        }
        
        // Scale by the half-length
        total += c1 * interval_sum;
    }
    
    return total;
}

// Implementation of exact_average_distance method in ProjectionCalculator
double ProjectionCalculator::exact_average_distance(int integration_order) const {
    // Define the integrand function for the expected projected distance
    auto projected_mean_abs = [this](const Polygon::VertexList& /* verts */, double theta) -> double {
        Direction dir = {std::cos(theta), std::sin(theta)};
        
        // Create a temporary ProjectionCalculator with the same settings but for this direction
        ProjectionCalculator temp_calc(polygon_, dir, region_, 0);
        return temp_calc.exact_projected_distance();
    };
    
    // Validate integration_order parameter
    if (integration_order <= 0) {
        throw std::invalid_argument("Integration order must be positive");
    }
    
    // Use the specified integration_order for more accurate results
    // Higher values provide more accuracy but take longer to compute
    int order = std::max(4, integration_order); // Ensure minimum of 4 for reasonable accuracy
    
    // Integrate over all directions [0, π]
    double integral = integrate_over_directions(
        polygon_.vertices(), projected_mean_abs, order
    );
    
    // The average distance is 0.5 * integral
    return 0.5 * integral;
}

// Implementation of monte_carlo_average_distance
double ProjectionCalculator::monte_carlo_average_distance(int n_samples) const {
    double sum = 0.0;
    
    for (int i = 0; i < n_samples; ++i) {
        Polygon::Point p1, p2;
        
        // Sample points from the appropriate region
        if (region_ == Region::INTERIOR) {
            p1 = sample_interior_point();
            p2 = sample_interior_point();
        } else {
            p1 = sample_boundary_point();
            p2 = sample_boundary_point();
        }
        
        // Calculate Euclidean distance
        double dx = p2[0] - p1[0];
        double dy = p2[1] - p1[1];
        sum += std::sqrt(dx * dx + dy * dy);
    }
    
    return sum / n_samples;
}

// Implementation of batched Monte Carlo analysis
ProjectionCalculator::MonteCarloResults 
ProjectionCalculator::run_monte_carlo_batched(
    int n_samples,
    int batch_size,
    const ProgressCallback& progress_callback) const 
{
    // Validate inputs
    validate_monte_carlo_inputs(n_samples, batch_size);
    
    // Initialize accumulators for projected and Euclidean distances
    double proj_dist_sum = 0.0;
    double euclidean_dist_sum = 0.0;
    int total_processed = 0;
    
    while (total_processed < n_samples) {
        // Determine current batch size
        int current_batch = std::min(batch_size, n_samples - total_processed);
        
        // Process batch
        for (int i = 0; i < current_batch; ++i) {
            Polygon::Point p1 = sample_point();
            Polygon::Point p2 = sample_point();
            
            // Calculate projected distance
            double proj1 = project_point(p1);
            double proj2 = project_point(p2);
            proj_dist_sum += std::abs(proj1 - proj2);
            
            // Calculate Euclidean distance
            double dx = p2[0] - p1[0];
            double dy = p2[1] - p1[1];
            euclidean_dist_sum += std::sqrt(dx * dx + dy * dy);
        }
        
        total_processed += current_batch;
        
        // Report progress if callback provided
        if (progress_callback) {
            double progress = static_cast<double>(total_processed) / n_samples;
            progress_callback(progress);
        }
    }
    
    // Return final results
    return {
        proj_dist_sum / n_samples,       // projected_distance
        euclidean_dist_sum / n_samples,  // average_distance
        total_processed                  // samples_processed
    };
}

} // namespace projection