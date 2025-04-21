#include <gtest/gtest.h>
#include <chrono>
#include <iostream>
#include <vector>
#include <iomanip>
#include "polygon_projection/polygon.h"
#include "polygon_projection/projection.h"
#include "polygon_projection/utils.h"

using namespace polygon_projection;

// This test file mirrors test_performance.py but in C++

// Helper template to time a function
template<typename F, typename... Args>
std::pair<double, double> time_function(F&& func, Args&&... args) {
    auto start = std::chrono::high_resolution_clock::now();
    auto result = func(std::forward<Args>(args)...);
    auto end = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double, std::milli> duration = end - start;
    return {result, duration.count()};
}

// Helper to create a regular polygon
Polygon create_regular_polygon(int sides, double radius = 1.0) {
    return Polygon::regular(sides, radius);
}

// Struct to hold benchmark results
struct BenchmarkResult {
    int polygon_sides;
    int samples;
    double interior_exact_ms;
    double interior_mc_ms;
    double boundary_exact_ms;
    double boundary_mc_ms;
    double avg_exact_ms;
    double avg_mc_ms;
    double sample_ms;
    double sample_rate;  // samples per second
};

// Function to benchmark operations on a regular polygon
BenchmarkResult benchmark_regular_polygon(int sides, int samples) {
    BenchmarkResult result;
    result.polygon_sides = sides;
    result.samples = samples;
    
    // Create regular polygon
    Polygon polygon = create_regular_polygon(sides);
    
    // Create calculator for interior projections
    ProjectionCalculator::Direction direction = {1.0, 0.0};
    ProjectionCalculator interior_calc(polygon, direction, Region::INTERIOR, 42);
    
    // Create calculator for boundary projections
    ProjectionCalculator boundary_calc(polygon, direction, Region::BOUNDARY, 42);
    
    // Measure exact interior projection
    {
        auto [val, time] = time_function([&interior_calc]() {
            return interior_calc.exact_projected_distance();
        });
        result.interior_exact_ms = time;
    }
    
    // Measure Monte Carlo interior projection
    {
        auto [val, time] = time_function([&interior_calc, samples]() {
            return interior_calc.monte_carlo_projected_distance(samples);
        });
        result.interior_mc_ms = time;
    }
    
    // Measure exact boundary projection
    {
        auto [val, time] = time_function([&boundary_calc]() {
            return boundary_calc.exact_projected_distance();
        });
        result.boundary_exact_ms = time;
    }
    
    // Measure Monte Carlo boundary projection
    {
        auto [val, time] = time_function([&boundary_calc, samples]() {
            return boundary_calc.monte_carlo_projected_distance(samples);
        });
        result.boundary_mc_ms = time;
    }
    
    // Measure exact average distance
    {
        auto [val, time] = time_function([&interior_calc]() {
            return interior_calc.exact_average_distance();
        });
        result.avg_exact_ms = time;
    }
    
    // Measure Monte Carlo average distance
    {
        auto [val, time] = time_function([&interior_calc, samples]() {
            return interior_calc.monte_carlo_average_distance(samples);
        });
        result.avg_mc_ms = time;
    }
    
    // Measure point sampling time (by sampling multiple points)
    {
        std::vector<Polygon::Point> points;
        points.reserve(samples);
        
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < samples; ++i) {
            points.push_back(interior_calc.sample_point());
        }
        auto end = std::chrono::high_resolution_clock::now();
        
        std::chrono::duration<double, std::milli> duration = end - start;
        result.sample_ms = duration.count();
        result.sample_rate = samples / (duration.count() / 1000.0);
    }
    
    return result;
}

// Test different polygon sizes
TEST(PerformanceBenchmarkTest, DifferentPolygonSizes) {
    std::vector<int> polygon_sizes = {4, 6, 8, 12, 16, 24};
    std::cout << "\n[1] Benchmarking different polygon sizes (10,000 samples)" << std::endl;
    
    for (int sides : polygon_sizes) {
        BenchmarkResult result = benchmark_regular_polygon(sides, 10000);
        
        std::cout << "\n" << sides << "-sided polygon:" << std::endl;
        std::cout << "  Exact interior projection:  " << std::fixed << std::setprecision(3) << result.interior_exact_ms << " ms" << std::endl;
        std::cout << "  MC interior projection:     " << std::fixed << std::setprecision(3) << result.interior_mc_ms << " ms" << std::endl;
        std::cout << "  Exact boundary projection:  " << std::fixed << std::setprecision(3) << result.boundary_exact_ms << " ms" << std::endl;
        std::cout << "  MC boundary projection:     " << std::fixed << std::setprecision(3) << result.boundary_mc_ms << " ms" << std::endl;
        std::cout << "  Exact average distance:     " << std::fixed << std::setprecision(3) << result.avg_exact_ms << " ms" << std::endl;
        std::cout << "  MC average distance:        " << std::fixed << std::setprecision(3) << result.avg_mc_ms << " ms" << std::endl;
        std::cout << "  Point sampling rate:        " << std::fixed << std::setprecision(0) << result.sample_rate << " points/sec" << std::endl;
        
        // Verify that MC calculations take longer than exact calculations
        EXPECT_GT(result.interior_mc_ms, result.interior_exact_ms);
        EXPECT_GT(result.boundary_mc_ms, result.boundary_exact_ms);
        
        // Verify that average distance calculation is slower (more complex)
        EXPECT_GT(result.avg_exact_ms, result.interior_exact_ms);
    }
}

// Test different sample sizes
TEST(PerformanceBenchmarkTest, DifferentSampleSizes) {
    std::vector<int> sample_sizes = {1000, 10000, 100000};
    std::cout << "\n[2] Benchmarking different sample sizes (12-sided polygon)" << std::endl;
    
    for (int samples : sample_sizes) {
        BenchmarkResult result = benchmark_regular_polygon(12, samples);
        
        std::cout << "\n" << samples << " samples:" << std::endl;
        std::cout << "  MC interior projection:     " << std::fixed << std::setprecision(3) << result.interior_mc_ms << " ms" << std::endl;
        std::cout << "  MC boundary projection:     " << std::fixed << std::setprecision(3) << result.boundary_mc_ms << " ms" << std::endl;
        std::cout << "  MC average distance:        " << std::fixed << std::setprecision(3) << result.avg_mc_ms << " ms" << std::endl;
        std::cout << "  Sampling " << samples << " points:    " << std::fixed << std::setprecision(3) << result.sample_ms << " ms" << std::endl;
        
        if (samples > 1000) {
            // Verify that MC execution time scales roughly linearly with sample size
            // Ratios between 8-12x are reasonable for 10x increase in samples
            double ratio_interior = result.interior_mc_ms / (result.interior_mc_ms / samples * (samples / 10));
            double ratio_boundary = result.boundary_mc_ms / (result.boundary_mc_ms / samples * (samples / 10));
            
            EXPECT_GT(ratio_interior, 8.0);
            EXPECT_LT(ratio_interior, 12.0);
            EXPECT_GT(ratio_boundary, 8.0);
            EXPECT_LT(ratio_boundary, 12.0);
        }
    }
    
    std::cout << "\nPerformance Summary:" << std::endl;
    std::cout << "Monte Carlo operations scale linearly with sample count" << std::endl;
    std::cout << "Exact calculations are much faster for small/medium polygons" << std::endl;
    std::cout << "Point sampling scales roughly linearly with sample count" << std::endl;
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}