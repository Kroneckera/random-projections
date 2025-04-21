#include "polygon_projection/polygon.h"
#include "polygon_projection/utils.h"
#include <cmath>
#include <numeric>
#include <algorithm>

namespace polygon_projection {

Polygon::Polygon(const VertexList& vertices) : vertices_(vertices) {
    if (vertices_.size() < 3) {
        throw std::invalid_argument("Polygon must have at least 3 vertices");
    }
    
    if (!is_convex()) {
        throw std::invalid_argument("Polygon must be convex");
    }
}

Polygon Polygon::regular(int sides, double radius) {
    if (sides < 3) {
        throw std::invalid_argument("Regular polygon must have at least 3 sides");
    }
    
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    
    VertexList vertices;
    vertices.reserve(sides);
    
    for (int i = 0; i < sides; ++i) {
        double angle = 2.0 * M_PI * i / sides;
        vertices.push_back({radius * cos(angle), radius * sin(angle)});
    }
    
    return Polygon(vertices);
}

double Polygon::area() const {
    // Shoelace formula
    double sum = 0.0;
    for (size_t i = 0; i < vertices_.size(); ++i) {
        const auto& v1 = vertices_[i];
        const auto& v2 = vertices_[(i + 1) % vertices_.size()];
        sum += (v1[0] * v2[1]) - (v2[0] * v1[1]);
    }
    return 0.5 * std::abs(sum);
}

double Polygon::perimeter() const {
    double sum = 0.0;
    for (size_t i = 0; i < vertices_.size(); ++i) {
        const auto& v1 = vertices_[i];
        const auto& v2 = vertices_[(i + 1) % vertices_.size()];
        double dx = v2[0] - v1[0];
        double dy = v2[1] - v1[1];
        sum += std::sqrt(dx*dx + dy*dy);
    }
    return sum;
}

bool Polygon::is_convex() const {
    // Delegate to the standalone implementation
    // In Polygon context, we need the vertices to be counter-clockwise
    return polygon_projection::is_convex(vertices_, true);
}

bool Polygon::contains(const Point& point) const {
    // First, check if point is exactly on any edge - consider these inside
    size_t n = vertices_.size();
    for (size_t i = 0; i < n; i++) {
        const auto& v1 = vertices_[i];
        const auto& v2 = vertices_[(i + 1) % n];
        
        // Vector from v1 to v2
        double edge_x = v2[0] - v1[0];
        double edge_y = v2[1] - v1[1];
        
        // Vector from v1 to point
        double point_x = point[0] - v1[0];
        double point_y = point[1] - v1[1];
        
        // Edge length squared
        double edge_length_squared = edge_x * edge_x + edge_y * edge_y;
        
        // If edge is zero length, check if point equals vertex
        if (edge_length_squared < 1e-10) {
            if (std::abs(point_x) < 1e-10 && std::abs(point_y) < 1e-10) {
                return true; // Point is at vertex
            }
            continue;
        }
        
        // Projection of point onto edge
        double dot = point_x * edge_x + point_y * edge_y;
        double proj = dot / edge_length_squared;
        
        // Check if projection is within edge
        if (proj >= 0 && proj <= 1) {
            // Calculate distance from point to projected point
            double proj_x = v1[0] + proj * edge_x;
            double proj_y = v1[1] + proj * edge_y;
            double dist_x = point[0] - proj_x;
            double dist_y = point[1] - proj_y;
            double dist_squared = dist_x * dist_x + dist_y * dist_y;
            
            // If distance is small enough, point is on edge
            if (dist_squared < 1e-10) {
                return true;
            }
        }
    }
    
    // If not on edge, use ray casting algorithm
    bool inside = false;
    
    for (size_t i = 0, j = n - 1; i < n; j = i++) {
        const auto& vi = vertices_[i];
        const auto& vj = vertices_[j];
        
        if (((vi[1] > point[1]) != (vj[1] > point[1])) &&
            (point[0] < (vj[0] - vi[0]) * (point[1] - vi[1]) / (vj[1] - vi[1]) + vi[0])) {
            inside = !inside;
        }
    }
    
    return inside;
}

} // namespace polygon_projection