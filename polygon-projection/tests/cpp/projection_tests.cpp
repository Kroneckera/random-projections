#include <gtest/gtest.h>
#include <cmath>
#include "polygon_projection/polygon.h"
#include "polygon_projection/projection.h"
#include "polygon_projection/utils.h"

using namespace polygon_projection;

// Helper function to create a square with vertices at (0,0), (1,0), (1,1), (0,1)
Polygon create_unit_square() {
    Polygon::VertexList square_vertices = {{{0, 0}}, {{1, 0}}, {{1, 1}}, {{0, 1}}};
    return Polygon(square_vertices);
}

// Helper function to create a regular hexagon with specified radius
Polygon create_regular_hexagon(double radius = 1.0) {
    return Polygon::regular(6, radius);
}

TEST(ProjectionTest, RegularPolygonProjections) {
    // Create a unit square
    Polygon square = create_unit_square();
    
    // Create projection calculator for interior
    ProjectionCalculator x_proj(square, {1.0, 0.0}, Region::INTERIOR, 42);
    ProjectionCalculator y_proj(square, {0.0, 1.0}, Region::INTERIOR, 42);
    ProjectionCalculator diag_proj(square, {1.0, 1.0}, Region::INTERIOR, 42);
    
    // Calculate exact projected distances
    double x_dist = x_proj.exact_projected_distance();
    double y_dist = y_proj.exact_projected_distance();
    double diag_dist = diag_proj.exact_projected_distance();
    
    // For a unit square, both x and y projections should be the same
    EXPECT_NEAR(x_dist, y_dist, 1e-10);
    
    // For a unit square, the diagonal projection has a known analytical value: 7*sqrt(2)/30
    double expected_diag = 7.0 * std::sqrt(2.0) / 30.0;
    EXPECT_NEAR(diag_dist, expected_diag, 1e-10);
}

TEST(ProjectionTest, ExactVsMonteCarloInterior) {
    // Create a unit square
    Polygon square = create_unit_square();
    
    // Test directions
    std::vector<ProjectionCalculator::Direction> directions = {
        {1.0, 0.0},    // x-axis
        {0.0, 1.0},    // y-axis
        {1.0, 1.0}     // diagonal
    };
    
    for (const auto& dir : directions) {
        ProjectionCalculator calc(square, dir, Region::INTERIOR, 42);
        
        // Calculate exact and Monte Carlo distances
        double exact = calc.exact_projected_distance();
        double monte_carlo = calc.monte_carlo_projected_distance(100000);
        
        // For large sample sizes, Monte Carlo should be within 1% of exact
        double relative_error = std::abs(exact - monte_carlo) / exact;
        EXPECT_LT(relative_error, 0.01);
    }
}

TEST(ProjectionTest, ExactVsMonteCarloBoundary) {
    // Create a unit square
    Polygon square = create_unit_square();
    
    // Test random directions
    std::vector<ProjectionCalculator::Direction> directions = {
        {1.3, 2.7},    // Random direction
        {3.1, 1.4},    // Another random direction
        {2.2, 3.5}     // Yet another random direction
    };
    
    for (const auto& dir : directions) {
        ProjectionCalculator calc(square, dir, Region::BOUNDARY, 42);
        
        // Calculate exact and Monte Carlo distances
        double exact = calc.exact_projected_distance();
        double monte_carlo = calc.monte_carlo_projected_distance(500000);
        
        // For large sample sizes, Monte Carlo should be within 1% of exact
        double relative_error = std::abs(exact - monte_carlo) / exact;
        EXPECT_LT(relative_error, 0.01);
    }
}

TEST(ProjectionTest, SamplingFunctions) {
    // Create a unit square
    Polygon square = create_unit_square();
    
    // Test interior sampling
    ProjectionCalculator interior_calc(square, {1.0, 0.0}, Region::INTERIOR, 42);
    
    // Sample 100 interior points
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
    
    // Sample 100 boundary points
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

TEST(ProjectionTest, ProjectionFunctions) {
    // Create a unit square
    Polygon square = create_unit_square();
    
    // Test point projection onto x-axis
    ProjectionCalculator x_proj(square, {1.0, 0.0}, Region::INTERIOR, 42);
    
    // Test points
    std::vector<Polygon::Point> points = {
        {0.0, 0.0},    // Origin
        {1.0, 0.0},    // Right
        {0.0, 1.0},    // Top
        {1.0, 1.0},    // Top-right
        {0.5, 0.5}     // Center
    };
    
    // Expected x-projections
    std::vector<double> expected_x = {0.0, 1.0, 0.0, 1.0, 0.5};
    
    for (size_t i = 0; i < points.size(); ++i) {
        double proj = x_proj.project_point(points[i]);
        EXPECT_NEAR(proj, expected_x[i], 1e-10);
    }
    
    // Test point projection onto diagonal
    ProjectionCalculator diag_proj(square, {1.0, 1.0}, Region::INTERIOR, 42);
    
    // Expected diagonal projections (normalized by sqrt(2))
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

TEST(ProjectionTest, BatchedMonteCarloCalculation) {
    // Create a unit square
    Polygon square = create_unit_square();
    
    // Create projection calculator
    ProjectionCalculator calc(square, {1.0, 0.0}, Region::INTERIOR, 42);
    
    // Test parameters
    int n_samples = 1000;
    int batch_size = 200;
    
    // Track progress
    std::vector<double> progress_values;
    auto progress_callback = [&progress_values](double progress) {
        progress_values.push_back(progress);
    };
    
    // Run batched Monte Carlo
    auto results = calc.run_monte_carlo_batched(n_samples, batch_size, progress_callback);
    
    // Check results
    EXPECT_EQ(results.samples_processed, n_samples);
    EXPECT_GE(results.projected_distance, 0.0);
    EXPECT_GE(results.average_distance, 0.0);
    
    // Check progress reporting
    EXPECT_EQ(static_cast<int>(progress_values.size()), n_samples / batch_size);
    if (!progress_values.empty()) {
        EXPECT_NEAR(progress_values.back(), 1.0, 1e-10);
    }
}

TEST(ProjectionTest, ExactAverageDistance) {
    // Create a unit square
    Polygon square = create_unit_square();
    
    // Create projection calculator
    ProjectionCalculator calc(square, {1.0, 0.0}, Region::INTERIOR, 42);
    
    // Calculate exact average distance
    double avg_dist = calc.exact_average_distance();
    
    // For a unit square, the average distance has a known analytical value
    // (2 + sqrt(2) + 5*ln(1+sqrt(2)))/15 ≈ 0.5214...
    double analytical = (2.0 + std::sqrt(2.0) + 5.0*std::log(1.0+std::sqrt(2.0)))/15.0;
    
    EXPECT_NEAR(avg_dist, analytical, 1e-3);  // Allow modest error due to numerical integration
}

TEST(ProjectionTest, ErrorHandling) {
    // Create a unit square
    Polygon square = create_unit_square();
    
    // Test invalid direction (zero vector)
    EXPECT_THROW(
        ProjectionCalculator(square, {0.0, 0.0}, Region::INTERIOR, 42),
        std::invalid_argument
    );
    
    // Test invalid sample count
    ProjectionCalculator calc(square, {1.0, 0.0}, Region::INTERIOR, 42);
    EXPECT_THROW(
        calc.monte_carlo_projected_distance(-10),
        std::invalid_argument
    );
    
    // Test invalid batch size
    EXPECT_THROW(
        calc.run_monte_carlo_batched(100, 0, nullptr),
        std::invalid_argument
    );
    
    // Test invalid integration order
    EXPECT_THROW(
        calc.exact_average_distance(-5),
        std::invalid_argument
    );
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}