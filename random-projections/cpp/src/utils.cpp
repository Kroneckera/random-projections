#include "polygon.h"
#include "projection.h"
#include "utils.h"
#include <cmath>
#include <vector>

namespace projection {

// Helper function to generate a regular polygon with n sides
Polygon::VertexList generate_regular_polygon(int n_sides, double radius) {
    if (n_sides < 3) {
        throw std::invalid_argument("Regular polygon must have at least 3 sides");
    }
    
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    
    Polygon::VertexList vertices;
    vertices.reserve(n_sides);
    
    for (int i = 0; i < n_sides; ++i) {
        double angle = 2.0 * M_PI * i / n_sides;
        vertices.push_back({radius * std::cos(angle), radius * std::sin(angle)});
    }
    
    return vertices;
}

// Helper function to calculate squared length between two points
inline double _lensq(const Polygon::Point& p, const Polygon::Point& q) {
    return (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]);
}

// Helper function for cross product calculation
inline double _cross(const Polygon::Point& a, const Polygon::Point& b, const Polygon::Point& c) {
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
}

/**
 * Check if a polygon is convex according to the following definition:
 * 
 * A polygon is convex if:
 * 1. It has at least 3 vertices
 * 2. All interior angles are less than or equal to 180 degrees
 * 3. All vertices "turn" in the same direction (all clockwise or all counter-clockwise)
 * 4. The polygon is simple (non-self-intersecting) 
 * 5. Collinear edges are allowed as long as they don't form reflex angles
 * 
 * @param vertices List of polygon vertices in order
 * @param check_ccw If true, requires vertices to be in counter-clockwise order
 * @param eps Tolerance for floating-point comparisons
 * @return True if the polygon is convex (and CCW if check_ccw is true)
 */
bool is_convex(const Polygon::VertexList& vertices, bool check_ccw, double eps) {
    size_t n = vertices.size();
    if (n < 3) {
        return false;
    }

    int orient = 0;            // +1 = CCW, −1 = CW, 0 = not set yet
    double signed_turn_sum = 0.0;

    for (size_t i = 0; i < n; ++i) {
        const auto& prev = vertices[(i + n - 1) % n];
        const auto& curr = vertices[i];
        const auto& next = vertices[(i + 1) % n];

        // Reject zero-length edges (duplicate vertices)
        if (_lensq(prev, curr) < eps || _lensq(curr, next) < eps) {
            return false;
        }

        // Cross product to determine orientation (negative = CW, positive = CCW)
        double z = _cross(curr, next, prev);

        if (std::abs(z) >= eps) {  // Non-collinear points (genuine corner)
            int s = (z > 0.0) ? 1 : -1;
            if (orient == 0) {
                orient = s;  // Remember first valid turn direction
            } else if (s != orient) {
                return false;  // Sign flip => concave/reflex angle
            }
        }
        // else: collinear points - allowed in a convex polygon

        // Accumulate signed exterior angle for total turn calculation
        double e_i_x = curr[0] - prev[0];
        double e_i_y = curr[1] - prev[1];
        double e_ip1_x = next[0] - curr[0];
        double e_ip1_y = next[1] - curr[1];
        
        double cross_e = e_i_x * e_ip1_y - e_i_y * e_ip1_x;
        double dot_e = e_i_x * e_ip1_x + e_i_y * e_ip1_y;
        signed_turn_sum += std::atan2(cross_e, dot_e);
    }

    // If all vertices are collinear, not a proper polygon
    if (orient == 0) {
        return false;
    }

    // Simple polygon must have total turn of ±2π
    // This test ensures the polygon is not self-intersecting
    if (std::abs(std::abs(signed_turn_sum) - 2 * M_PI) > 1e-8) {
        return false;
    }

    // If check_ccw is true, reject clockwise polygons
    if (check_ccw && orient == -1) {
        return false;  // Polygon is convex but clockwise
    }

    return true;
}

// Helper function to compute a polygon's area
double calculate_polygon_area(const Polygon::VertexList& vertices) {
    // Shoelace formula
    double sum = 0.0;
    size_t n = vertices.size();
    
    for (size_t i = 0; i < n; ++i) {
        const auto& v1 = vertices[i];
        const auto& v2 = vertices[(i + 1) % n];
        sum += v1[0] * v2[1] - v2[0] * v1[1];
    }
    
    return 0.5 * std::abs(sum);
}

// Helper function to compute a polygon's perimeter
double calculate_polygon_perimeter(const Polygon::VertexList& vertices) {
    double sum = 0.0;
    size_t n = vertices.size();
    
    for (size_t i = 0; i < n; ++i) {
        const auto& v1 = vertices[i];
        const auto& v2 = vertices[(i + 1) % n];
        
        double dx = v2[0] - v1[0];
        double dy = v2[1] - v1[1];
        sum += std::sqrt(dx * dx + dy * dy);
    }
    
    return sum;
}

// Helper function to normalize a direction vector
ProjectionCalculator::Direction normalize_direction(const ProjectionCalculator::Direction& dir) {
    double len = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1]);
    if (len < 1e-10) {
        throw std::invalid_argument("Direction vector cannot be zero");
    }
    return {dir[0] / len, dir[1] / len};
}

} // namespace projection