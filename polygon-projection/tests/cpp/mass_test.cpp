#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <cmath>
#include "polygon_projection/polygon.h"
#include "polygon_projection/projection.h"
#include "polygon_projection/utils.h"

using namespace polygon_projection;

// Simple utility to create polygon objects
Polygon create_square() {
    Polygon::VertexList vertices = {{{0, 0}}, {{1, 0}}, {{1, 1}}, {{0, 1}}};
    return Polygon(vertices);
}

int main() {
    std::cout << "Detailed debugging for C++ mass calculation issues\n" << std::endl;
    
    // Create unit square and projection directions
    Polygon square = create_square();
    ProjectionCalculator::Direction x_direction = {1.0, 0.0};
    ProjectionCalculator::Direction y_direction = {0.0, 1.0};
    ProjectionCalculator::Direction diag_direction = {1.0, 1.0};

    // 1. Test interior projection along x-axis
    std::cout << "=== Interior Projection (X-axis) ===" << std::endl;
    ProjectionCalculator calc_x(square, x_direction, Region::INTERIOR);
    auto interior_density_x = calc_x.compute_interior_density();
    
    std::cout << "Number of pieces: " << interior_density_x.pieces.size() << std::endl;
    for (size_t i = 0; i < interior_density_x.pieces.size(); ++i) {
        const auto& p = interior_density_x.pieces[i];
        std::cout << "Piece " << i << ": [" << std::fixed << std::setprecision(10) << p.x0 
                  << ", " << p.x1 << "], m=" << p.m << ", c=" << p.c << std::endl;
    }
    
    // Verify total mass integration
    double total_mass = 0.0;
    for (const auto& p : interior_density_x.pieces) {
        double mass = 0.5 * p.m * (p.x1 * p.x1 - p.x0 * p.x0) + p.c * (p.x1 - p.x0);
        total_mass += mass;
        std::cout << "Piece mass: " << std::fixed << std::setprecision(10) << mass 
                  << ", Running total: " << total_mass << std::endl;
    }
    std::cout << "Total mass should be 1.0, actual: " << std::fixed << std::setprecision(10) 
              << total_mass << "\n" << std::endl;
    
    // 2. Test interior projection along y-axis
    std::cout << "=== Interior Projection (Y-axis) ===" << std::endl;
    ProjectionCalculator calc_y(square, y_direction, Region::INTERIOR);
    auto interior_density_y = calc_y.compute_interior_density();
    
    std::cout << "Number of pieces: " << interior_density_y.pieces.size() << std::endl;
    for (size_t i = 0; i < interior_density_y.pieces.size(); ++i) {
        const auto& p = interior_density_y.pieces[i];
        std::cout << "Piece " << i << ": [" << std::fixed << std::setprecision(10) << p.x0 
                  << ", " << p.x1 << "], m=" << p.m << ", c=" << p.c << std::endl;
    }
    
    // Verify total mass integration
    total_mass = 0.0;
    for (const auto& p : interior_density_y.pieces) {
        double mass = 0.5 * p.m * (p.x1 * p.x1 - p.x0 * p.x0) + p.c * (p.x1 - p.x0);
        total_mass += mass;
        std::cout << "Piece mass: " << std::fixed << std::setprecision(10) << mass 
                  << ", Running total: " << total_mass << std::endl;
    }
    std::cout << "Total mass should be 1.0, actual: " << std::fixed << std::setprecision(10) 
              << total_mass << "\n" << std::endl;
    
    // 3. Test interior projection along diagonal
    std::cout << "=== Interior Projection (Diagonal) ===" << std::endl;
    ProjectionCalculator calc_diag(square, diag_direction, Region::INTERIOR);
    auto interior_density_diag = calc_diag.compute_interior_density();
    
    std::cout << "Number of pieces: " << interior_density_diag.pieces.size() << std::endl;
    for (size_t i = 0; i < interior_density_diag.pieces.size(); ++i) {
        const auto& p = interior_density_diag.pieces[i];
        std::cout << "Piece " << i << ": [" << std::fixed << std::setprecision(10) << p.x0 
                  << ", " << p.x1 << "], m=" << p.m << ", c=" << p.c << std::endl;
    }
    
    // Verify total mass integration
    total_mass = 0.0;
    for (const auto& p : interior_density_diag.pieces) {
        double mass = 0.5 * p.m * (p.x1 * p.x1 - p.x0 * p.x0) + p.c * (p.x1 - p.x0);
        total_mass += mass;
        std::cout << "Piece mass: " << std::fixed << std::setprecision(10) << mass 
                  << ", Running total: " << total_mass << std::endl;
    }
    std::cout << "Total mass should be 1.0, actual: " << std::fixed << std::setprecision(10) 
              << total_mass << "\n" << std::endl;
    
    // 4. Test boundary projection along x-axis
    std::cout << "=== Boundary Projection (X-axis) ===" << std::endl;
    ProjectionCalculator calc_x_boundary(square, x_direction, Region::BOUNDARY);
    auto boundary_steps_x = calc_x_boundary.compute_boundary_steps();
    
    std::cout << "Number of steps: " << boundary_steps_x.size() << std::endl;
    for (size_t i = 0; i < boundary_steps_x.size(); ++i) {
        const auto& step = boundary_steps_x[i];
        double lo = std::get<0>(step);
        double hi = std::get<1>(step);
        double c = std::get<2>(step);
        std::cout << "Step " << i << ": [" << std::fixed << std::setprecision(10) << lo 
                  << ", " << hi << "], c=" << c << std::endl;
    }
    
    // Verify total mass integration
    total_mass = 0.0;
    for (const auto& step : boundary_steps_x) {
        double lo = std::get<0>(step);
        double hi = std::get<1>(step);
        double c = std::get<2>(step);
        double mass = c * (hi - lo);
        total_mass += mass;
        std::cout << "Step mass: " << std::fixed << std::setprecision(10) << mass 
                  << ", Running total: " << total_mass << std::endl;
    }
    std::cout << "Total mass should be 1.0, actual: " << std::fixed << std::setprecision(10) 
              << total_mass << "\n" << std::endl;
    
    // 5. Create regular hexagon and test
    Polygon hexagon = Polygon::regular(6, 1.0);
    ProjectionCalculator::Direction hex_direction = {0.7, 0.4};
    
    // Test interior projection of hexagon
    std::cout << "=== Hexagon Interior Projection ===" << std::endl;
    ProjectionCalculator calc_hex(hexagon, hex_direction, Region::INTERIOR);
    auto interior_density_hex = calc_hex.compute_interior_density();
    
    std::cout << "Number of pieces: " << interior_density_hex.pieces.size() << std::endl;
    for (size_t i = 0; i < interior_density_hex.pieces.size(); ++i) {
        const auto& p = interior_density_hex.pieces[i];
        std::cout << "Piece " << i << ": [" << std::fixed << std::setprecision(10) << p.x0 
                  << ", " << p.x1 << "], m=" << p.m << ", c=" << p.c << std::endl;
    }
    
    // Verify total mass integration
    total_mass = 0.0;
    for (const auto& p : interior_density_hex.pieces) {
        double mass = 0.5 * p.m * (p.x1 * p.x1 - p.x0 * p.x0) + p.c * (p.x1 - p.x0);
        total_mass += mass;
        std::cout << "Piece mass: " << std::fixed << std::setprecision(10) << mass 
                  << ", Running total: " << total_mass << std::endl;
    }
    std::cout << "Total mass should be 1.0, actual: " << std::fixed << std::setprecision(10) 
              << total_mass << "\n" << std::endl;
    
    // 6. Test boundary projection of hexagon
    std::cout << "=== Hexagon Boundary Projection ===" << std::endl;
    ProjectionCalculator calc_hex_boundary(hexagon, hex_direction, Region::BOUNDARY);
    auto boundary_steps_hex = calc_hex_boundary.compute_boundary_steps();
    
    std::cout << "Number of steps: " << boundary_steps_hex.size() << std::endl;
    for (size_t i = 0; i < boundary_steps_hex.size(); ++i) {
        const auto& step = boundary_steps_hex[i];
        double lo = std::get<0>(step);
        double hi = std::get<1>(step);
        double c = std::get<2>(step);
        std::cout << "Step " << i << ": [" << std::fixed << std::setprecision(10) << lo 
                  << ", " << hi << "], c=" << c << std::endl;
    }
    
    // Verify total mass integration
    total_mass = 0.0;
    for (const auto& step : boundary_steps_hex) {
        double lo = std::get<0>(step);
        double hi = std::get<1>(step);
        double c = std::get<2>(step);
        double mass = c * (hi - lo);
        total_mass += mass;
        std::cout << "Step mass: " << std::fixed << std::setprecision(10) << mass 
                  << ", Running total: " << total_mass << std::endl;
    }
    std::cout << "Total mass should be 1.0, actual: " << std::fixed << std::setprecision(10) 
              << total_mass << "\n" << std::endl;
    
    // 7. Test expected_abs_diff calculations for each case
    try {
        std::cout << "=== Expected Absolute Difference Calculations ===" << std::endl;
        
        // Square interior projections
        double square_x = calc_x.exact_projected_distance();
        std::cout << "Square Interior (X-axis): " << square_x << std::endl;
        
        double square_y = calc_y.exact_projected_distance();
        std::cout << "Square Interior (Y-axis): " << square_y << std::endl;
        
        double square_diag = calc_diag.exact_projected_distance();
        std::cout << "Square Interior (Diagonal): " << square_diag << std::endl;
        
        // Square boundary projection
        double square_boundary = calc_x_boundary.exact_projected_distance();
        std::cout << "Square Boundary (X-axis): " << square_boundary << std::endl;
        
        // Hexagon projections
        double hex_interior = calc_hex.exact_projected_distance();
        std::cout << "Hexagon Interior: " << hex_interior << std::endl;
        
        double hex_boundary = calc_hex_boundary.exact_projected_distance();
        std::cout << "Hexagon Boundary: " << hex_boundary << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error in exact distance calculations: " << e.what() << std::endl;
    }
    
    return 0;
}