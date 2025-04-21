#include <iostream>
#include <iomanip>
#include "polygon_projection/polygon.h"
#include "polygon_projection/projection.h"
#include "polygon_projection/utils.h"

using namespace polygon_projection;

int main() {
    try {
        // Create a regular hexagon with radius 2
        Polygon::VertexList vertices = generate_regular_polygon(6, 2.0);
        Polygon hexagon(vertices);
        
        // Define projection direction
        ProjectionCalculator::Direction direction = {1.0, 0.5};
        
        std::cout << "Testing projections for a regular hexagon (radius 2.0)\n";
        std::cout << "Direction: [" << direction[0] << ", " << direction[1] << "]\n\n";
        
        // Test interior projections
        std::cout << "Interior projections:\n";
        std::cout << "--------------------\n";
        
        ProjectionCalculator interior_calc(hexagon, direction, Region::INTERIOR, 42);
        
        // Exact calculations
        double interior_exact = interior_calc.exact_projected_distance();
        std::cout << "Exact expected distance: " << std::fixed << std::setprecision(6) 
                  << interior_exact << "\n";
        
        // Monte Carlo calculations
        std::cout << "\nMonte Carlo approximations:\n";
        for (int samples : {1000, 10000, 100000}) {
            double monte_carlo = interior_calc.monte_carlo_projected_distance(samples);
            double error = std::abs(monte_carlo - interior_exact) / interior_exact;
            
            std::cout << std::setw(7) << samples << " samples: " 
                      << std::fixed << std::setprecision(6) << monte_carlo
                      << " (relative error: " << std::scientific << std::setprecision(4) 
                      << error << ")\n";
        }
        
        // Test boundary projections
        std::cout << "\nBoundary projections:\n";
        std::cout << "--------------------\n";
        
        ProjectionCalculator boundary_calc(hexagon, direction, Region::BOUNDARY, 42);
        
        // Exact calculations
        double boundary_exact = boundary_calc.exact_projected_distance();
        std::cout << "Exact expected distance: " << std::fixed << std::setprecision(6) 
                  << boundary_exact << "\n";
        
        // Monte Carlo calculations
        std::cout << "\nMonte Carlo approximations:\n";
        for (int samples : {1000, 10000, 100000}) {
            double monte_carlo = boundary_calc.monte_carlo_projected_distance(samples);
            double error = std::abs(monte_carlo - boundary_exact) / boundary_exact;
            
            std::cout << std::setw(7) << samples << " samples: " 
                      << std::fixed << std::setprecision(6) << monte_carlo
                      << " (relative error: " << std::scientific << std::setprecision(4) 
                      << error << ")\n";
        }
        
        // Test average distance calculation (exact)
        std::cout << "\nAverage Euclidean distances:\n";
        std::cout << "---------------------------\n";
        double interior_avg = interior_calc.exact_average_distance();
        double boundary_avg = boundary_calc.exact_average_distance();
        
        std::cout << "Interior average distance: " << std::fixed << std::setprecision(6) 
                  << interior_avg << "\n";
        std::cout << "Boundary average distance: " << std::fixed << std::setprecision(6) 
                  << boundary_avg << "\n";
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}