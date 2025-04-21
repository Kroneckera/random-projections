#include <iostream>
#include <vector>
#include <cmath>
#include "polygon.h"
#include "projection.h"
#include "utils.h"

using namespace projection;

int main() {
    // Create a regular hexagon
    Polygon hexagon(generate_regular_polygon(6, 2.0));
    
    // Print polygon info
    std::cout << "Hexagon area: " << hexagon.area() << std::endl;
    std::cout << "Hexagon perimeter: " << hexagon.perimeter() << std::endl;
    
    // Define projection direction
    ProjectionCalculator::Direction direction = {1.0, 0.5};
    
    // Create calculators for interior and boundary regions
    ProjectionCalculator interior_calc(hexagon, direction, Region::INTERIOR, 42);
    ProjectionCalculator boundary_calc(hexagon, direction, Region::BOUNDARY, 43);
    
    // Calculate exact projections
    double interior_exact = interior_calc.exact_projected_distance();
    double boundary_exact = boundary_calc.exact_projected_distance();
    
    std::cout << "Exact expected projected distances:" << std::endl;
    std::cout << "  Interior: " << interior_exact << std::endl;
    std::cout << "  Boundary: " << boundary_exact << std::endl;
    
    // Calculate Monte Carlo estimates
    int n_samples = 100000;
    double interior_mc = interior_calc.monte_carlo_projected_distance(n_samples);
    double boundary_mc = boundary_calc.monte_carlo_projected_distance(n_samples);
    
    std::cout << "\nMonte Carlo estimates (" << n_samples << " samples):" << std::endl;
    std::cout << "  Interior: " << interior_mc << std::endl;
    std::cout << "  Boundary: " << boundary_mc << std::endl;
    
    // Calculate errors
    double interior_error = std::abs(interior_exact - interior_mc) / interior_exact;
    double boundary_error = std::abs(boundary_exact - boundary_mc) / boundary_exact;
    
    std::cout << "\nRelative errors:" << std::endl;
    std::cout << "  Interior: " << interior_error << std::endl;
    std::cout << "  Boundary: " << boundary_error << std::endl;
    
    // Calculate average distances
    double interior_avg = interior_calc.exact_average_distance();
    double boundary_avg = boundary_calc.exact_average_distance();
    
    std::cout << "\nAverage Euclidean distances:" << std::endl;
    std::cout << "  Interior: " << interior_avg << std::endl;
    std::cout << "  Boundary: " << boundary_avg << std::endl;
    
    // Test batched Monte Carlo calculation
    std::cout << "\nRunning batched Monte Carlo calculation..." << std::endl;
    
    ProjectionCalculator::MonteCarloResults results = interior_calc.run_monte_carlo_batched(
        50000, 1000, 
        [](double progress) {
            std::cout << "Progress: " << int(progress * 100) << "%" << std::endl;
        }
    );
    
    std::cout << "Batched MC results:" << std::endl;
    std::cout << "  Projected distance: " << results.projected_distance << std::endl;
    std::cout << "  Average distance: " << results.average_distance << std::endl;
    std::cout << "  Samples processed: " << results.samples_processed << std::endl;
    
    return 0;
}