#include <iostream>
#include "polygon_projection/polygon.h"
#include "polygon_projection/utils.h"

using namespace polygon_projection;

int main() {
    try {
        // Create a regular hexagon with radius 2
        Polygon::VertexList vertices = generate_regular_polygon(6, 2.0);
        
        // Create polygon object
        Polygon hexagon(vertices);
        
        // Print polygon properties
        std::cout << "Hexagon area: " << hexagon.area() << std::endl;
        std::cout << "Hexagon perimeter: " << hexagon.perimeter() << std::endl;
        std::cout << "Is convex: " << (is_convex(vertices) ? "yes" : "no") << std::endl;
        
        // Print vertices
        std::cout << "\nVertices:" << std::endl;
        const auto& verts = hexagon.vertices();
        for (size_t i = 0; i < verts.size(); ++i) {
            std::cout << "  " << i << ": (" << verts[i][0] << ", " << verts[i][1] << ")" << std::endl;
        }
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}