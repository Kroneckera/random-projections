#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <cmath>
#include "polygon_projection/polygon.h"
#include "polygon_projection/projection.h"
#include "polygon_projection/utils.h"

using namespace polygon_projection;

int main() {
    std::cout << "===== C++ IMPLEMENTATION DETAILS =====\n" << std::endl;
    
    // Create a unit square and diagonal projection direction
    Polygon::VertexList square_vertices = {{{0, 0}}, {{1, 0}}, {{1, 1}}, {{0, 1}}};
    Polygon square(square_vertices);
    ProjectionCalculator::Direction diag_direction = {1.0, 1.0};
    
    std::cout << "Step 1: Create a projection calculator with a diagonal direction\n";
    // Create a calculator without actually calling the density computation yet
    ProjectionCalculator calc(square, diag_direction, Region::INTERIOR);
    
    std::cout << "\nStep 2: Extract the implementation details from compute_interior_density\n";
    // Get the interior density
    auto density = calc.compute_interior_density();
    
    // Display the computed pieces
    std::cout << "\nResults from compute_interior_density:" << std::endl;
    std::cout << "Number of pieces: " << density.pieces.size() << std::endl;
    
    double total_mass = 0.0;
    for (size_t i = 0; i < density.pieces.size(); ++i) {
        const auto& piece = density.pieces[i];
        
        // Output piece information
        std::cout << "Piece " << i << ": [" << std::fixed << std::setprecision(10) 
                  << piece.x0 << ", " << piece.x1 << "], m=" << piece.m 
                  << ", c=" << piece.c << std::endl;
        
        // Calculate mass contribution
        double mass = 0.5 * piece.m * (piece.x1 * piece.x1 - piece.x0 * piece.x0) + 
                      piece.c * (piece.x1 - piece.x0);
        total_mass += mass;
        
        std::cout << "  Mass: " << mass << ", Running total: " << total_mass << std::endl;
    }
    
    std::cout << "\nTotal mass: " << total_mass << std::endl;
    
    // Manually verify by recomputing with our own formula
    std::cout << "\nManual verification:" << std::endl;
    
    // Calculate projection for the diagonal direction
    std::cout << "Direction: (" << diag_direction[0] << ", " << diag_direction[1] << ")" << std::endl;
    
    // Normalize direction
    double norm = std::sqrt(diag_direction[0] * diag_direction[0] + diag_direction[1] * diag_direction[1]);
    double direction_x = diag_direction[0] / norm;
    double direction_y = diag_direction[1] / norm;
    
    std::cout << "Normalized direction: (" << direction_x << ", " << direction_y << ")" << std::endl;
    
    // Rotation angle
    double phi = std::atan2(direction_y, direction_x);
    std::cout << "Rotation angle: " << phi << " rad" << std::endl;
    
    // Rotation matrix
    double cos_phi = std::cos(phi);
    double sin_phi = std::sin(phi);
    std::cout << "cos(phi): " << cos_phi << ", sin(phi): " << sin_phi << std::endl;
    
    // Rotate the vertices
    std::vector<std::array<double, 2>> rotated_vertices;
    for (const auto& vertex : square_vertices) {
        double x = cos_phi * vertex[0] + sin_phi * vertex[1];
        double y = -sin_phi * vertex[0] + cos_phi * vertex[1];
        rotated_vertices.push_back({x, y});
        std::cout << "  Rotated: (" << x << ", " << y << ")" << std::endl;
    }
    
    // Calculate area
    double area = 0.0;
    for (size_t i = 0; i < rotated_vertices.size(); ++i) {
        size_t j = (i + 1) % rotated_vertices.size();
        double term = rotated_vertices[i][0] * rotated_vertices[j][1] - 
                     rotated_vertices[j][0] * rotated_vertices[i][1];
        area += term;
        std::cout << "Area term " << i << "->" << j << ": " << term << std::endl;
    }
    area = 0.5 * std::abs(area);
    std::cout << "Area: " << area << std::endl;
    
    // Collect edge contributions
    std::vector<std::tuple<double, double, double, double>> events;
    const double tol = 1e-12;
    
    std::cout << "\nEdge contributions:" << std::endl;
    for (size_t i = 0; i < rotated_vertices.size(); ++i) {
        double x0 = rotated_vertices[i][0];
        double y0 = rotated_vertices[i][1];
        double x1 = rotated_vertices[(i+1) % rotated_vertices.size()][0];
        double y1 = rotated_vertices[(i+1) % rotated_vertices.size()][1];
        
        double dx = x1 - x0;
        double dy = y1 - y0;
        
        std::cout << "Edge " << i << ": (" << x0 << "," << y0 << ") -> (" 
                  << x1 << "," << y1 << "), dx=" << dx << ", dy=" << dy << std::endl;
        
        if (std::abs(dx) < tol) {
            std::cout << "  Skipping vertical edge" << std::endl;
            continue;
        }
        
        double k = dy / dx;
        double b = y0 - k * x0;
        int sign = (dx < 0) ? 1 : -1;
        double lo = std::min(x0, x1);
        double hi = std::max(x0, x1);
        
        std::cout << "  Adding event: lo=" << lo << ", hi=" << hi 
                  << ", k=" << k << ", b=" << b << ", sign=" << sign << std::endl;
        events.emplace_back(lo, hi, sign * k, sign * b);
    }
    
    // Display collected events
    std::cout << "\nCollected events:" << std::endl;
    for (size_t i = 0; i < events.size(); ++i) {
        const auto& event = events[i];
        std::cout << "Event " << i << ": [" << std::get<0>(event) << ", " << std::get<1>(event) 
                  << "], ak=" << std::get<2>(event) << ", bk=" << std::get<3>(event) << std::endl;
    }
    
    // Get breakpoints
    std::vector<double> xs;
    for (const auto& vertex : rotated_vertices) {
        xs.push_back(vertex[0]);
    }
    
    std::cout << "\nBreakpoints before sorting: ";
    for (double x : xs) {
        std::cout << x << " ";
    }
    std::cout << std::endl;
    
    std::sort(xs.begin(), xs.end());
    xs.erase(std::unique(xs.begin(), xs.end(), 
            [&tol](double a, double b) { return std::abs(a - b) < tol; }), 
            xs.end());
    
    std::cout << "Breakpoints after sorting and deduplication: ";
    for (double x : xs) {
        std::cout << x << " ";
    }
    std::cout << std::endl;
    
    // Build piecewise density
    using Piece = ProjectionCalculator::ProjectedDensity::Piece;
    std::vector<Piece> manual_pieces;
    
    std::cout << "\nPiecewise density calculation:" << std::endl;
    for (size_t i = 0; i < xs.size() - 1; ++i) {
        double lo = xs[i];
        double hi = xs[i + 1];
        
        std::cout << "Piece " << i << ": [" << lo << ", " << hi << "]" << std::endl;
        
        double a_sum = 0.0;
        double b_sum = 0.0;
        
        std::cout << "  Checking events:" << std::endl;
        for (size_t j = 0; j < events.size(); ++j) {
            const auto& event = events[j];
            double l = std::get<0>(event);
            double h = std::get<1>(event);
            double ak = std::get<2>(event);
            double bk = std::get<3>(event);
            
            bool contributes = (l <= lo && lo < h);
            std::cout << "    Event " << j << ": [" << l << ", " << h 
                      << "], ak=" << ak << ", bk=" << bk 
                      << " - " << (contributes ? "ACTIVE" : "inactive") << std::endl;
            
            if (contributes) {
                a_sum += ak;
                b_sum += bk;
            }
        }
        
        double m = a_sum / area;
        double c = b_sum / area;
        manual_pieces.push_back({lo, hi, m, c});
        
        std::cout << "  Final: a_sum=" << a_sum << ", b_sum=" << b_sum 
                  << ", m=" << m << ", c=" << c << std::endl;
    }
    
    // Verify total mass from manual calculation
    double manual_total_mass = 0.0;
    std::cout << "\nManual pieces:" << std::endl;
    for (size_t i = 0; i < manual_pieces.size(); ++i) {
        const auto& piece = manual_pieces[i];
        
        // Output piece information
        std::cout << "Piece " << i << ": [" << piece.x0 << ", " << piece.x1 
                  << "], m=" << piece.m << ", c=" << piece.c << std::endl;
        
        // Calculate mass contribution
        double mass = 0.5 * piece.m * (piece.x1 * piece.x1 - piece.x0 * piece.x0) + 
                      piece.c * (piece.x1 - piece.x0);
        manual_total_mass += mass;
        
        std::cout << "  Mass: " << mass << ", Running total: " << manual_total_mass << std::endl;
    }
    
    std::cout << "\nTotal mass from manual calculation: " << manual_total_mass << std::endl;
    
    // Check differences between automated and manual calculations
    std::cout << "\nDifferences between automated and manual calculations:" << std::endl;
    std::cout << "Number of pieces - Automated: " << density.pieces.size() 
              << ", Manual: " << manual_pieces.size() << std::endl;
    std::cout << "Total mass - Automated: " << total_mass 
              << ", Manual: " << manual_total_mass << std::endl;
              
    return 0;
}