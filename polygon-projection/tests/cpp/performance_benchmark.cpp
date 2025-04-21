#include <iostream>
#include <iomanip>
#include <chrono>
#include "polygon_projection/polygon.h"
#include "polygon_projection/projection.h"
#include "polygon_projection/utils.h"

using namespace polygon_projection;
using Clock = std::chrono::high_resolution_clock;
using Duration = std::chrono::duration<double>;

// Simple benchmark function
template<typename Func>
double benchmark(Func&& func, int repeats = 1) {
    auto start = Clock::now();
    
    for (int i = 0; i < repeats; ++i) {
        func();
    }
    
    auto end = Clock::now();
    Duration elapsed = end - start;
    return elapsed.count() / repeats;
}

int main() {
    try {
        // Create a complex polygon with many vertices
        Polygon::VertexList vertices = generate_regular_polygon(24, 2.0);
        Polygon polygon(vertices);
        
        // Define projection direction
        ProjectionCalculator::Direction direction = {1.0, 0.5};
        
        std::cout << "Performance benchmark for a " << vertices.size() << "-sided polygon\n";
        std::cout << "============================================\n\n";
        
        // Interior calculator
        ProjectionCalculator interior_calc(polygon, direction, Region::INTERIOR, 42);
        
        // Exact projection benchmark
        double time_exact = benchmark([&]() {
            interior_calc.exact_projected_distance();
        }, 100);
        
        std::cout << "Exact projection calculation (100 repeats):\n";
        std::cout << "  Time per calculation: " << std::fixed << std::setprecision(6) 
                  << time_exact * 1000 << " ms\n";
        
        // Monte Carlo benchmark with different sample sizes
        std::cout << "\nMonte Carlo calculation times:\n";
        
        for (int samples : {1000, 10000, 100000}) {
            double time_mc = benchmark([&]() {
                interior_calc.monte_carlo_projected_distance(samples);
            }, 5);
            
            std::cout << "  " << std::setw(7) << samples << " samples: " 
                      << std::fixed << std::setprecision(6) 
                      << time_mc * 1000 << " ms\n";
        }
        
        // Average distance benchmark
        double time_avg = benchmark([&]() {
            interior_calc.exact_average_distance();
        }, 5);
        
        std::cout << "\nExact average distance calculation (5 repeats):\n";
        std::cout << "  Time per calculation: " << std::fixed << std::setprecision(6) 
                  << time_avg * 1000 << " ms\n";
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}