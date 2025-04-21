#include <iostream>
#include "polygon.h"
#include "utils.h"

using namespace projection;

int main() {
    try {
        // Create a regular hexagon
        Polygon::VertexList vertices = generate_regular_polygon(6, 2.0);
        Polygon hexagon(vertices);
        
        // Print basic polygon info
        std::cout << "Polygon with " << vertices.size() << " vertices" << std::endl;
        std::cout << "Area: " << hexagon.area() << std::endl;
        std::cout << "Perimeter: " << hexagon.perimeter() << std::endl;
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}