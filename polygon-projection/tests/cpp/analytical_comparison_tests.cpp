#include <gtest/gtest.h>
#include <cmath>
#include <iostream>
#include "polygon_projection/polygon.h"
#include "polygon_projection/projection.h"
#include "polygon_projection/utils.h"

using namespace polygon_projection;

// This test file mirrors test_analytical_comparison.py but in C++

// Helper function to create a unit square
Polygon create_unit_square() {
    Polygon::VertexList vertices = {{{0, 0}}, {{1, 0}}, {{1, 1}}, {{0, 1}}};
    return Polygon(vertices);
}

// Calculate analytical value for average distance in unit square
double calculate_analytical_square_interior() {
    // For a unit square, the average distance has a known formula
    // (2 + sqrt(2) + 5*ln(1+sqrt(2)))/15
    return (2.0 + std::sqrt(2.0) + 5.0 * std::log(1.0 + std::sqrt(2.0))) / 15.0;
}

// Calculate analytical value for average distance on unit square boundary
double calculate_analytical_square_boundary() {
    // For a unit square perimeter, the average distance is correctly:
    // 1/4 + sqrt(2)/12 + 5/12 * ln(1 + sqrt(2))
    return 0.25 + std::sqrt(2.0)/12.0 + 5.0/12.0 * std::log(1.0 + std::sqrt(2.0));
}

// Unit Square Interior - Average Euclidean distance
TEST(AnalyticalComparisonTest, UnitSquareInteriorAverageDistance) {
    Polygon square = create_unit_square();
    
    ProjectionCalculator calc(square, {1.0, 0.0}, Region::INTERIOR, 12345);
    
    // Monte Carlo with high sample count
    double monte_carlo = calc.monte_carlo_average_distance(500000);
    
    // Analytical formula for unit square interior
    double analytical = calculate_analytical_square_interior();
    
    // Check that our Monte Carlo estimate is close to the analytical value
    // Allow for 1% error with large sample size
    double relative_error = std::abs(monte_carlo - analytical) / analytical;
    EXPECT_LT(relative_error, 0.01);
    
    std::cout << "Unit Square Interior Average Distance:" << std::endl;
    std::cout << "Monte Carlo result: " << monte_carlo << std::endl;
    std::cout << "Analytical formula: " << analytical << std::endl;
    std::cout << "Relative error: " << relative_error << std::endl;
}

// Unit Square Boundary - Average Euclidean distance
TEST(AnalyticalComparisonTest, UnitSquareBoundaryAverageDistance) {
    Polygon square = create_unit_square();
    
    ProjectionCalculator calc(square, {1.0, 0.0}, Region::BOUNDARY, 12345);
    
    // Monte Carlo with high sample count
    double monte_carlo = calc.monte_carlo_average_distance(500000);
    
    // Analytical formula for unit square boundary
    double analytical = calculate_analytical_square_boundary();
    
    // Check that our Monte Carlo estimate is close to the analytical value
    // Allow for 1% error with large sample size
    double relative_error = std::abs(monte_carlo - analytical) / analytical;
    EXPECT_LT(relative_error, 0.01);
    
    std::cout << "Unit Square Boundary Average Distance:" << std::endl;
    std::cout << "Monte Carlo result: " << monte_carlo << std::endl;
    std::cout << "Analytical formula: " << analytical << std::endl;
    std::cout << "Relative error: " << relative_error << std::endl;
}

// Unit Square Interior - Projected distance along x-axis
TEST(AnalyticalComparisonTest, UnitSquareInteriorProjectedDistance) {
    Polygon square = create_unit_square();
    
    ProjectionCalculator calc(square, {1.0, 0.0}, Region::INTERIOR, 12345);
    
    // Exact calculation
    double exact = calc.exact_projected_distance();
    
    // Monte Carlo with high sample count
    double monte_carlo = calc.monte_carlo_projected_distance(500000);
    
    // For a unit square, projected distance along x-axis is 1/3
    double analytical = 1.0 / 3.0;
    
    // Check that our exact calculation matches the analytical value
    double exact_error = std::abs(exact - analytical) / analytical;
    EXPECT_LT(exact_error, 1e-6);  // Should be very close
    
    // Check that our Monte Carlo estimate is close to the analytical value
    double monte_carlo_error = std::abs(monte_carlo - analytical) / analytical;
    EXPECT_LT(monte_carlo_error, 0.01);  // Allow for 1% error
    
    std::cout << "Unit Square Interior Projected Distance (x-axis):" << std::endl;
    std::cout << "Exact calculation: " << exact << std::endl;
    std::cout << "Monte Carlo result: " << monte_carlo << std::endl;
    std::cout << "Analytical formula: " << analytical << std::endl;
    std::cout << "Exact relative error: " << exact_error << std::endl;
    std::cout << "Monte Carlo relative error: " << monte_carlo_error << std::endl;
}

// Unit Square Boundary - Projected distance along x-axis
TEST(AnalyticalComparisonTest, UnitSquareBoundaryProjectedDistance) {
    Polygon square = create_unit_square();
    
    ProjectionCalculator calc(square, {1.0, 0.0}, Region::BOUNDARY, 12345);
    
    // Exact calculation
    double exact = calc.exact_projected_distance();
    
    // Monte Carlo with high sample count
    double monte_carlo = calc.monte_carlo_projected_distance(500000);
    
    // For a unit square boundary, projected distance along x-axis is 11/24
    double analytical = 11.0 / 24.0;
    
    // Check that our exact calculation matches the analytical value
    double exact_error = std::abs(exact - analytical) / analytical;
    EXPECT_LT(exact_error, 1e-6);  // Should be very close
    
    // Check that our Monte Carlo estimate is close to the analytical value
    double monte_carlo_error = std::abs(monte_carlo - analytical) / analytical;
    EXPECT_LT(monte_carlo_error, 0.01);  // Allow for 1% error
    
    std::cout << "Unit Square Boundary Projected Distance (x-axis):" << std::endl;
    std::cout << "Exact calculation: " << exact << std::endl;
    std::cout << "Monte Carlo result: " << monte_carlo << std::endl;
    std::cout << "Analytical formula: " << analytical << std::endl;
    std::cout << "Exact relative error: " << exact_error << std::endl;
    std::cout << "Monte Carlo relative error: " << monte_carlo_error << std::endl;
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}