#include "polygon_projection/projection.h"
#include "polygon_projection/utils.h"
#include <cmath>
#include <algorithm>
#include <random>
#include <set>
#include <vector>
#include <stdexcept>

namespace polygon_projection {

// This file contains common sampling utilities and helpers that aren't
// specific to either interior or boundary regions.

// Implementation functions for sampling points
namespace {
    // Implementation of interior point sampling
    Polygon::Point impl_sample_interior_point(const Polygon& polygon, std::mt19937& rng) {
        // Use rejection sampling for interior points
        const auto& vertices = polygon.vertices();
        
        // Find bounding box
        double min_x = vertices[0][0];
        double min_y = vertices[0][1];
        double max_x = vertices[0][0];
        double max_y = vertices[0][1];
        
        for (const auto& v : vertices) {
            min_x = std::min(min_x, v[0]);
            min_y = std::min(min_y, v[1]);
            max_x = std::max(max_x, v[0]);
            max_y = std::max(max_y, v[1]);
        }
        
        // Uniform distributions for sampling within bounding box
        std::uniform_real_distribution<double> dist_x(min_x, max_x);
        std::uniform_real_distribution<double> dist_y(min_y, max_y);
        
        // Rejection sampling loop
        Polygon::Point point;
        do {
            point = {dist_x(rng), dist_y(rng)};
        } while (!polygon.contains(point));
        
        return point;
    }

    // Implementation of boundary point sampling
    Polygon::Point impl_sample_boundary_point(const Polygon& polygon, std::mt19937& rng) {
        const auto& vertices = polygon.vertices();
        const size_t n = vertices.size();
        
        // Calculate edge lengths
        std::vector<double> edge_lengths;
        double total_length = 0.0;
        
        for (size_t i = 0; i < n; ++i) {
            const auto& v1 = vertices[i];
            const auto& v2 = vertices[(i + 1) % n];
            
            double dx = v2[0] - v1[0];
            double dy = v2[1] - v1[1];
            double length = std::sqrt(dx * dx + dy * dy);
            
            edge_lengths.push_back(length);
            total_length += length;
        }
        
        // Sample edge proportional to length
        std::uniform_real_distribution<double> dist(0.0, total_length);
        double target = dist(rng);
        
        double cumulative = 0.0;
        size_t edge_idx = 0;
        
        for (size_t i = 0; i < n; ++i) {
            cumulative += edge_lengths[i];
            if (target <= cumulative) {
                edge_idx = i;
                break;
            }
        }
        
        // Sample point along the selected edge
        const auto& v1 = vertices[edge_idx];
        const auto& v2 = vertices[(edge_idx + 1) % n];
        
        // Calculate parameter t along edge
        double t_max = edge_lengths[edge_idx];
        double t = target - (cumulative - t_max);
        double normalized_t = t / t_max;
        
        // Linearly interpolate between vertices
        return {
            v1[0] + normalized_t * (v2[0] - v1[0]),
            v1[1] + normalized_t * (v2[1] - v1[1])
        };
    }
} // anonymous namespace

// Public sample point function for interior and boundary points
Polygon::Point sample_interior_point(const Polygon& polygon, std::mt19937& rng) {
    return impl_sample_interior_point(polygon, rng);
}

Polygon::Point sample_boundary_point(const Polygon& polygon, std::mt19937& rng) {
    return impl_sample_boundary_point(polygon, rng);
}

// Validation functions for Monte Carlo calculations
void validate_monte_carlo_inputs(int n_samples, int batch_size) {
    if (n_samples <= 0) {
        throw std::invalid_argument("Number of samples must be positive");
    }
    if (batch_size <= 0) {
        throw std::invalid_argument("Batch size must be positive");
    }
    if (batch_size > n_samples) {
        throw std::invalid_argument("Batch size cannot be larger than number of samples");
    }
}

// Validation function for polygons
void validate_polygon(const Polygon& polygon) {
    if (polygon.vertices().size() < 3) {
        throw std::invalid_argument("Polygon must have at least 3 vertices");
    }
    if (!polygon.is_convex()) {
        throw std::invalid_argument("Polygon must be convex");
    }
    if (polygon.area() < 1e-10) {
        throw std::invalid_argument("Polygon area must be positive");
    }
}

// Note: monte_carlo_projected_distance is implemented in projection.cpp

// Helper function to generate a random convex polygon (useful for testing)
Polygon::VertexList generate_random_convex_polygon(
    int n, double rmin, double rmax, std::mt19937& rng)
{
    // Sample more points than needed in an annulus
    int m = std::max(n * 5, 50);
    std::vector<Polygon::Point> pts;
    
    std::uniform_real_distribution<double> angle_dist(0.0, 2.0 * M_PI);
    std::uniform_real_distribution<double> radius_dist(rmin, rmax);
    
    for (int i = 0; i < m; ++i) {
        double theta = angle_dist(rng);
        double r = radius_dist(rng);
        pts.push_back({r * std::cos(theta), r * std::sin(theta)});
    }
    
    // Compute convex hull
    // Note: This is a simplified implementation, not a full convex hull algorithm
    // For a complete implementation, you would use Graham scan or similar
    
    // Sort by x-coordinate
    std::sort(pts.begin(), pts.end(), [](const Polygon::Point& a, const Polygon::Point& b) {
        return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]);
    });
    
    // Remove duplicates
    pts.erase(std::unique(pts.begin(), pts.end(), [](const Polygon::Point& a, const Polygon::Point& b) {
        return a[0] == b[0] && a[1] == b[1];
    }), pts.end());
    
    // Extract convex hull
    std::vector<Polygon::Point> hull;
    
    // Build lower hull
    for (const auto& p : pts) {
        while (hull.size() >= 2) {
            auto& p1 = hull[hull.size() - 2];
            auto& p2 = hull[hull.size() - 1];
            
            double cross = (p2[0] - p1[0]) * (p[1] - p1[1]) - (p2[1] - p1[1]) * (p[0] - p1[0]);
            if (cross <= 0) {
                hull.pop_back();
            } else {
                break;
            }
        }
        hull.push_back(p);
    }
    
    // Save size of lower hull
    size_t lower_size = hull.size();
    
    // Build upper hull
    for (auto it = pts.rbegin(); it != pts.rend(); ++it) {
        const auto& p = *it;
        while (hull.size() > lower_size) {
            auto& p1 = hull[hull.size() - 2];
            auto& p2 = hull[hull.size() - 1];
            
            double cross = (p2[0] - p1[0]) * (p[1] - p1[1]) - (p2[1] - p1[1]) * (p[0] - p1[0]);
            if (cross <= 0) {
                hull.pop_back();
            } else {
                break;
            }
        }
        hull.push_back(p);
    }
    
    // Remove duplicate endpoints
    hull.pop_back();
    
    // If hull has too few points, try again recursively with more points
    if (hull.size() < static_cast<size_t>(n)) {
        return generate_random_convex_polygon(n, rmin, rmax, rng);
    }
    
    // Extract n evenly spaced points from the hull
    Polygon::VertexList result;
    for (int i = 0; i < n; ++i) {
        size_t idx = static_cast<size_t>(i * hull.size() / n);
        result.push_back(hull[idx]);
    }
    
    return result;
}

} // namespace polygon_projection