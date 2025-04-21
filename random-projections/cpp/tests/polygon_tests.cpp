#include <gtest/gtest.h>
#include "polygon.h"

using namespace projection;

TEST(PolygonTest, RegularPolygonCreation) {
    // Test regular triangle creation
    Polygon triangle = Polygon::regular(3, 1.0);
    EXPECT_EQ(triangle.vertices().size(), static_cast<size_t>(3));
    
    // Test regular square creation
    Polygon square = Polygon::regular(4, 1.0);
    EXPECT_EQ(square.vertices().size(), static_cast<size_t>(4));
    
    // Test with invalid parameters
    EXPECT_THROW(Polygon::regular(2, 1.0), std::invalid_argument); // Too few sides
    EXPECT_THROW(Polygon::regular(3, -1.0), std::invalid_argument); // Negative radius
}

TEST(PolygonTest, AreaCalculation) {
    // Test square with side length 2
    Polygon::VertexList square_vertices = {{{-1, -1}}, {{1, -1}}, {{1, 1}}, {{-1, 1}}};
    Polygon square(square_vertices);
    EXPECT_NEAR(square.area(), 4.0, 1e-10);
    
    // Test equilateral triangle with side length 1
    double h = std::sqrt(3.0) / 2.0; // Height of equilateral triangle
    Polygon::VertexList triangle_vertices = {{{0, 0}}, {{1, 0}}, {{0.5, h}}};
    Polygon triangle(triangle_vertices);
    EXPECT_NEAR(triangle.area(), 0.5 * 1.0 * h, 1e-10);
}

TEST(PolygonTest, PerimeterCalculation) {
    // Test square with side length 2
    Polygon::VertexList square_vertices = {{{-1, -1}}, {{1, -1}}, {{1, 1}}, {{-1, 1}}};
    Polygon square(square_vertices);
    EXPECT_NEAR(square.perimeter(), 8.0, 1e-10); // 4 sides of length 2
    
    // Test equilateral triangle with side length 1
    double h = std::sqrt(3.0) / 2.0; // Height of equilateral triangle
    Polygon::VertexList triangle_vertices = {{{0, 0}}, {{1, 0}}, {{0.5, h}}};
    Polygon triangle(triangle_vertices);
    EXPECT_NEAR(triangle.perimeter(), 3.0, 1e-10); // 3 sides of length 1
}

TEST(PolygonTest, ConvexityChecking) {
    // Convex square
    Polygon::VertexList square_vertices = {{{-1, -1}}, {{1, -1}}, {{1, 1}}, {{-1, 1}}};
    Polygon square(square_vertices);
    EXPECT_TRUE(square.is_convex());
    
    // Actually a convex shape with collinear points
    Polygon::VertexList old_vertices = {{{-1, -1}}, {{1, -1}}, {{0, 0}}, {{-1, 1}}};
    
    // Create a truly non-convex shape (concave pentagon)
    Polygon::VertexList non_convex_vertices = {
        {{0.0, 0.0}}, {{2.0, 0.0}}, {{1.0, 1.0}}, {{2.0, 2.0}}, {{0.0, 2.0}}
    };
    
    // Custom check - manually calculate cross products
    bool convex = true;
    int last_sign = 0;
    
    std::cout << "Debug - Testing truly non-convex pentagon:" << std::endl;
    for (size_t i = 0; i < non_convex_vertices.size(); ++i) {
        const auto& p1 = non_convex_vertices[i];
        const auto& p2 = non_convex_vertices[(i + 1) % non_convex_vertices.size()];
        const auto& p3 = non_convex_vertices[(i + 2) % non_convex_vertices.size()];
        
        // Vectors from p1->p2 and p2->p3
        double v1x = p2[0] - p1[0];
        double v1y = p2[1] - p1[1];
        double v2x = p3[0] - p2[0];
        double v2y = p3[1] - p2[1];
        
        // z-component of the cross product
        double cross_z = v1x * v2y - v1y * v2x;
        
        std::cout << "  Point " << i << ": (" << p1[0] << "," << p1[1] << ")" << std::endl;
        std::cout << "  Vector " << i << "->" << (i+1)%non_convex_vertices.size() << ": (" << v1x << "," << v1y << ")" << std::endl;
        std::cout << "  Vector " << (i+1)%non_convex_vertices.size() << "->" << (i+2)%non_convex_vertices.size() << ": (" << v2x << "," << v2y << ")" << std::endl;
        std::cout << "  Cross product: " << cross_z << std::endl;
        
        // Skip near-zero cross products (collinear points)
        if (std::abs(cross_z) < 1e-10) {
            std::cout << "  Skipping (collinear)" << std::endl;
            continue;
        }
        
        int sign = (cross_z > 0) ? 1 : -1;
        std::cout << "  Sign: " << sign << std::endl;
        
        if (last_sign == 0) {
            last_sign = sign;
        } else if (sign != last_sign) {
            // If the sign of the cross product changes, the polygon is not convex
            std::cout << "  Sign changed from " << last_sign << " to " << sign << " -> NOT CONVEX" << std::endl;
            convex = false;
            break;
        }
    }
    
    std::cout << "Manual check result: " << (convex ? "CONVEX" : "NOT CONVEX") << std::endl;
    
    // Now try to construct the polygon, which should throw
    EXPECT_THROW({
        Polygon p(non_convex_vertices);
    }, std::invalid_argument);
}

TEST(PolygonTest, PointContainment) {
    // Square with vertices at (±1, ±1)
    Polygon::VertexList square_vertices = {{{-1, -1}}, {{1, -1}}, {{1, 1}}, {{-1, 1}}};
    Polygon square(square_vertices);
    
    // Test points
    EXPECT_TRUE(square.contains({{0, 0}}));      // Center
    EXPECT_TRUE(square.contains({{0.5, 0.5}}));  // Inside
    EXPECT_FALSE(square.contains({{2, 0}}));     // Outside
    EXPECT_FALSE(square.contains({{-2, -2}}));   // Outside
    
    // Edge cases are considered inside by this algorithm
    EXPECT_TRUE(square.contains({{1, 0}}));      // On edge
    EXPECT_TRUE(square.contains({{0, 1}}));      // On edge
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}