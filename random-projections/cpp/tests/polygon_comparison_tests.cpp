#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <iostream>
#include "polygon.h"
#include "projection.h"
#include "utils.h"

using namespace projection;

// This test file mirrors test_polygon_projection.py but in C++

// Helper function to create a unit square
Polygon create_unit_square() {
    Polygon::VertexList vertices = {{{0, 0}}, {{1, 0}}, {{1, 1}}, {{0, 1}}};
    return Polygon(vertices);
}

// Helper function to create a regular hexagon
Polygon create_regular_hexagon() {
    return Polygon::regular(6, 1.0);
}

// Test polygon creation
TEST(PolygonComparisonTest, PolygonCreation) {
    // Direct creation from vertices
    Polygon::VertexList square_vertices = {{{0, 0}}, {{1, 0}}, {{1, 1}}, {{0, 1}}};
    Polygon square(square_vertices);
    EXPECT_EQ(square.vertices().size(), static_cast<size_t>(4));
    
    // Creation using regular polygon method
    Polygon hexagon = Polygon::regular(6, 2.0);
    EXPECT_EQ(hexagon.vertices().size(), static_cast<size_t>(6));
    
    // Check radius of vertices (should all be 2.0 units from origin)
    for (const auto& vertex : hexagon.vertices()) {
        double dist = std::sqrt(vertex[0] * vertex[0] + vertex[1] * vertex[1]);
        EXPECT_NEAR(dist, 2.0, 1e-10);
    }
}

// Test polygon properties
TEST(PolygonComparisonTest, PolygonProperties) {
    Polygon square = create_unit_square();
    EXPECT_NEAR(square.area(), 1.0, 1e-10);  // Unit square has area 1.0
    EXPECT_NEAR(square.perimeter(), 4.0, 1e-10);  // Unit square has perimeter 4.0
    
    Polygon hexagon = create_regular_hexagon();
    // Regular hexagon with radius 1 has area 3*sqrt(3)/2
    EXPECT_NEAR(hexagon.area(), 3.0 * std::sqrt(3.0) / 2.0, 1e-10);
    // Regular hexagon with radius 1 has perimeter 6
    EXPECT_NEAR(hexagon.perimeter(), 6.0, 1e-10);
}

// Test exact vs Monte Carlo interior projections
TEST(PolygonComparisonTest, ExactVsMonteCarloInterior) {
    Polygon square = create_unit_square();
    
    std::vector<ProjectionCalculator::Direction> directions = {
        {1.0, 0.0},    // x-axis
        {0.0, 1.0},    // y-axis
        {1.0, 1.0}     // diagonal
    };
    
    for (const auto& direction : directions) {
        ProjectionCalculator calc(square, direction, Region::INTERIOR, 42);
        
        double exact = calc.exact_projected_distance();
        double monte_carlo = calc.monte_carlo_projected_distance(100000);
        
        // For large sample sizes, MC should be within 1% of exact
        double relative_error = std::abs(exact - monte_carlo) / exact;
        EXPECT_LT(relative_error, 0.01);
    }
}

// Test exact vs Monte Carlo boundary projections
TEST(PolygonComparisonTest, ExactVsMonteCarloBoundary) {
    Polygon square = create_unit_square();
    
    std::vector<ProjectionCalculator::Direction> directions = {
        {1.3, 2.7},    // Random direction
        {3.1, 1.4},    // Another random direction
        {2.2, 3.5}     // Yet another random direction
    };
    
    for (const auto& direction : directions) {
        ProjectionCalculator calc(square, direction, Region::BOUNDARY, 42);
        
        double exact = calc.exact_projected_distance();
        double monte_carlo = calc.monte_carlo_projected_distance(500000);
        
        // For large sample sizes, MC should be within 1% of exact
        double relative_error = std::abs(exact - monte_carlo) / exact;
        EXPECT_LT(relative_error, 0.01);
    }
}

// Test sampling functions
TEST(PolygonComparisonTest, SamplingFunctions) {
    Polygon square = create_unit_square();
    
    // Test interior sampling
    ProjectionCalculator interior_calc(square, {1.0, 0.0}, Region::INTERIOR, 42);
    
    for (int i = 0; i < 100; ++i) {
        Polygon::Point point = interior_calc.sample_point();
        
        // Point should be inside the unit square
        EXPECT_GE(point[0], 0.0);
        EXPECT_LE(point[0], 1.0);
        EXPECT_GE(point[1], 0.0);
        EXPECT_LE(point[1], 1.0);
    }
    
    // Test boundary sampling
    ProjectionCalculator boundary_calc(square, {1.0, 0.0}, Region::BOUNDARY, 42);
    
    for (int i = 0; i < 100; ++i) {
        Polygon::Point point = boundary_calc.sample_point();
        
        // Point should be on the boundary
        bool on_boundary = 
            (std::abs(point[0]) < 1e-10 && 0.0 <= point[1] && point[1] <= 1.0) ||
            (std::abs(point[0] - 1.0) < 1e-10 && 0.0 <= point[1] && point[1] <= 1.0) ||
            (std::abs(point[1]) < 1e-10 && 0.0 <= point[0] && point[0] <= 1.0) ||
            (std::abs(point[1] - 1.0) < 1e-10 && 0.0 <= point[0] && point[0] <= 1.0);
            
        EXPECT_TRUE(on_boundary);
    }
}

// Test Monte Carlo batched calculation
TEST(PolygonComparisonTest, MonteCarloAnalysis) {
    Polygon square = create_unit_square();
    
    ProjectionCalculator calc(square, {1.0, 0.0}, Region::INTERIOR, 42);
    
    // Set up parameters
    int total_samples = 1000;
    int batch_size = 200;
    
    // Track progress
    std::vector<double> progress_values;
    auto progress_callback = [&progress_values](double progress) {
        progress_values.push_back(progress);
    };
    
    // Run batched Monte Carlo
    auto results = calc.run_monte_carlo_batched(total_samples, batch_size, progress_callback);
    
    // Check results
    EXPECT_EQ(results.samples_processed, total_samples);
    EXPECT_GT(results.projected_distance, 0.0);
    EXPECT_GT(results.average_distance, 0.0);
    
    // Check progress reporting
    EXPECT_EQ(progress_values.size(), static_cast<size_t>(total_samples / batch_size));
    EXPECT_NEAR(progress_values.back(), 1.0, 1e-10);
}

// Test projection functions
TEST(PolygonComparisonTest, ProjectionFunctions) {
    Polygon square = create_unit_square();
    
    // Test projection onto x-axis
    ProjectionCalculator x_proj(square, {1.0, 0.0}, Region::INTERIOR, 42);
    
    std::vector<Polygon::Point> points = {
        {0.0, 0.0},    // Origin
        {1.0, 0.0},    // Right
        {0.0, 1.0},    // Top
        {1.0, 1.0},    // Top-right
        {0.5, 0.5}     // Center
    };
    
    std::vector<double> expected_x = {0.0, 1.0, 0.0, 1.0, 0.5};
    
    for (size_t i = 0; i < points.size(); ++i) {
        double proj = x_proj.project_point(points[i]);
        EXPECT_NEAR(proj, expected_x[i], 1e-10);
    }
    
    // Test projection onto diagonal
    ProjectionCalculator diag_proj(square, {1.0, 1.0}, Region::INTERIOR, 42);
    
    std::vector<double> expected_diag = {
        0.0,                    // Origin
        1.0 / std::sqrt(2.0),   // Right
        1.0 / std::sqrt(2.0),   // Top
        2.0 / std::sqrt(2.0),   // Top-right
        1.0 / std::sqrt(2.0)    // Center
    };
    
    for (size_t i = 0; i < points.size(); ++i) {
        double proj = diag_proj.project_point(points[i]);
        EXPECT_NEAR(proj, expected_diag[i], 1e-10);
    }
}

// Test error handling
TEST(PolygonComparisonTest, ErrorHandling) {
    // Test invalid polygon creation (too few vertices)
    Polygon::VertexList too_few_vertices = {{{0, 0}}, {{1, 0}}};
    EXPECT_THROW({
        Polygon p(too_few_vertices);
    }, std::invalid_argument);
    
    // Test non-convex polygon
    Polygon::VertexList non_convex = {{{0, 0}}, {{2, 0}}, {{1, 1}}, {{2, 2}}, {{0, 2}}};
    EXPECT_THROW({
        Polygon p(non_convex);
    }, std::invalid_argument);
    
    // Test invalid direction (zero vector)
    Polygon square = create_unit_square();
    EXPECT_THROW({
        ProjectionCalculator calc(square, {0.0, 0.0}, Region::INTERIOR, 42);
    }, std::invalid_argument);
    
    // Test invalid sample count
    ProjectionCalculator calc(square, {1.0, 0.0}, Region::INTERIOR, 42);
    EXPECT_THROW({
        calc.monte_carlo_projected_distance(-10);
    }, std::invalid_argument);
    
    // Test invalid batch size
    EXPECT_THROW({
        calc.run_monte_carlo_batched(100, 0, nullptr);
    }, std::invalid_argument);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}