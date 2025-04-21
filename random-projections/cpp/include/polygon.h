#ifndef PROJECTION_POLYGON_H
#define PROJECTION_POLYGON_H

#include <vector>
#include <array>
#include <stdexcept>

namespace projection {

/**
 * Represents a convex polygon in 2D space.
 * All operations assume and maintain the convexity of the polygon.
 */
class Polygon {
public:
    // Type definitions
    using Point = std::array<double, 2>;
    using VertexList = std::vector<Point>;

    /**
     * Create a polygon from a list of vertices.
     * Vertices must be in counter-clockwise order.
     * 
     * @param vertices List of vertices in counter-clockwise order
     * @throws std::invalid_argument if fewer than 3 vertices are provided or polygon is not convex
     */
    explicit Polygon(const VertexList& vertices);

    /**
     * Create a regular polygon with n sides and given radius.
     * 
     * @param sides Number of sides (must be >= 3)
     * @param radius Radius of the circumscribed circle
     * @return A regular polygon
     * @throws std::invalid_argument if sides < 3 or radius <= 0
     */
    static Polygon regular(int sides, double radius = 1.0);

    /**
     * Get the vertices of the polygon.
     * 
     * @return List of vertices in counter-clockwise order
     */
    const VertexList& vertices() const { return vertices_; }

    /**
     * Get the area of the polygon.
     * 
     * @return Area value
     */
    double area() const;

    /**
     * Get the perimeter of the polygon.
     * 
     * @return Perimeter value
     */
    double perimeter() const;

    /**
     * Check if the polygon is convex.
     * 
     * @return true if the polygon is convex
     */
    bool is_convex() const;

    /**
     * Check if a point is inside the polygon.
     * 
     * @param point The point to check
     * @return true if the point is inside the polygon
     */
    bool contains(const Point& point) const;

private:
    VertexList vertices_;
};

} // namespace projection

#endif // PROJECTION_POLYGON_H