#include "polygon_projection/projection.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <set>
#include <stdexcept>

namespace polygon_projection {

// Implementation of compute_interior_density method from ProjectionCalculator class
ProjectionCalculator::ProjectedDensity 
ProjectionCalculator::compute_interior_density() const {
    // Get polygon vertices
    const auto& vertices = polygon_.vertices();
    const auto& direction = normalized_direction_;
    
    // Step 1: Rotate so direction aligns with x-axis
    double phi = std::atan2(direction[1], direction[0]);
    double cos_phi = std::cos(phi);
    double sin_phi = std::sin(phi);
    
    // Rotation matrix
    std::array<std::array<double, 2>, 2> R = {{
        {cos_phi, sin_phi},
        {-sin_phi, cos_phi}
    }};
    
    // Rotate vertices
    std::vector<std::array<double, 2>> rotated_vertices;
    rotated_vertices.reserve(vertices.size());
    
    std::vector<double> xs;
    xs.reserve(vertices.size());
    
    for (const auto& vertex : vertices) {
        double x = R[0][0] * vertex[0] + R[0][1] * vertex[1];
        double y = R[1][0] * vertex[0] + R[1][1] * vertex[1];
        rotated_vertices.push_back({x, y});
        xs.push_back(x);
    }
    
    // Step 2: Compute polygon area
    double area = 0.0;
    for (size_t i = 0; i < rotated_vertices.size(); ++i) {
        size_t j = (i + 1) % rotated_vertices.size();
        double term = rotated_vertices[i][0] * rotated_vertices[j][1] - 
                      rotated_vertices[j][0] * rotated_vertices[i][1];
        area += term;
    }
    area = 0.5 * std::abs(area);
    
    if (area <= 0) {
        throw std::runtime_error("Polygon area must be positive; check vertex order.");
    }
    
    // Step 3: Collect each edge's contribution
    std::vector<std::tuple<double, double, double, double>> events;
    
    const double tol = 1e-12;  // Tolerance for vertical edges
    const size_t n = rotated_vertices.size();
    
    for (size_t i = 0; i < n; ++i) {
        double x0 = rotated_vertices[i][0];
        double y0 = rotated_vertices[i][1];
        double x1 = rotated_vertices[(i+1) % n][0];
        double y1 = rotated_vertices[(i+1) % n][1];
        
        double dx = x1 - x0;
        double dy = y1 - y0;
        
        if (std::abs(dx) < tol) {
            continue;  // Vertical edge, no continuous contribution
        }
        
        double k = dy / dx;
        double b = y0 - k * x0;
        
        // Key fix: This is +1 if dx<0, -1 if dx>0 (exact match to Python)
        int sign = (dx < 0) ? 1 : -1;
        
        double lo = std::min(x0, x1);
        double hi = std::max(x0, x1);
        
        events.emplace_back(lo, hi, sign * k, sign * b);
    }
    
    // Step 4: Build piecewise sum over sorted breakpoints
    // First, sort xs and remove duplicates
    std::sort(xs.begin(), xs.end());
    
    // Match Python's np.sort(np.unique(xs)) behavior exactly
    std::set<double> unique_set(xs.begin(), xs.end());
    std::vector<double> unique_xs(unique_set.begin(), unique_set.end());
    
    // Replace xs with unique_xs
    xs = std::move(unique_xs);
    
    ProjectedDensity result;
    
    // Now build the piecewise linear density
    for (size_t i = 0; i < xs.size() - 1; ++i) {
        double lo = xs[i];
        double hi = xs[i + 1];
        
        double a_sum = 0.0;
        double b_sum = 0.0;
        
        for (const auto& e : events) {
            double l = std::get<0>(e);
            double h = std::get<1>(e);
            double ak = std::get<2>(e);
            double bk = std::get<3>(e);
            
            // Python's simple rule is: if the interval start (lo) is contained in [l, h)
            bool is_active = (l <= lo) && (lo < h);
                
            // For floating point endpoints, use tolerance
            if ((std::abs(l - lo) < tol) || (std::abs(h - lo) < tol)) {
                if (std::abs(h - lo) < tol) { 
                    // If lo == h, then NOT active (interval start is at edge end)
                    is_active = false;
                } else {
                    // If lo == l, still active (interval start is at edge start)
                    is_active = true;
                }
            }
            
            if (is_active) {
                a_sum += ak;
                b_sum += bk;
            }
        }
        
        // Scale by area
        result.pieces.push_back({lo, hi, a_sum / area, b_sum / area});
    }
    
    return result;
}

// Implementation of evaluate_density function
double evaluate_density(double x, const ProjectionCalculator::ProjectedDensity& density) {
    double result = 0.0;
    
    for (const auto& piece : density.pieces) {
        if (x >= piece.x0 && x <= piece.x1) {
            result = piece.m * x + piece.c;
            break;
        }
    }
    
    return result;
}

// Implementation of expected_abs_diff from projection_sampling.py
double expected_abs_diff(const ProjectionCalculator::ProjectedDensity& density) {
    // 1) Accumulate F0 at each segment start
    struct SegmentInfo {
        double x0;
        double x1;
        double m;
        double c;
        double F0;
    };
    
    std::vector<SegmentInfo> segs;
    double F0 = 0.0;
    
    // Strict match to Python - need to ensure mass calculation is accurate
    for (const auto& piece : density.pieces) {
        segs.push_back({
            piece.x0,
            piece.x1,
            piece.m,
            piece.c,
            F0
        });
        
        // Calculate mass exactly as in Python: 0.5*m*(x1**2 - x0**2) + c*(x1 - x0)
        double mass = 0.5 * piece.m * (piece.x1 * piece.x1 - piece.x0 * piece.x0) + 
                      piece.c * (piece.x1 - piece.x0);
        F0 += mass;
    }
    
    // Check if total mass is reasonably close to 1.0
    const double tol = 1e-8;  // Use same tolerance as Python (1e-8)
    if (std::abs(F0 - 1.0) > tol) {
        throw std::runtime_error("Total mass should be 1.0");
    }
    
    // 2) Integrate 2*∫ F(x)*(1-F(x)) dx piecewise - match Python code exactly
    double total = 0.0;
    
    for (const auto& s : segs) {
        double x0 = s.x0;
        double x1 = s.x1;
        double m = s.m;
        double c = s.c;
        double Fstart = s.F0;
        
        // F(x) = alpha*x^2 + beta*x + gamma
        double alpha = 0.5 * m;
        double beta = c;
        double gamma = Fstart - (alpha * x0 * x0 + beta * x0);
        
        // Expand g(x) = F(x)*(1 - F(x)) = -F(x)^2 + F(x)
        //   = p4*x^4 + p3*x^3 + p2*x^2 + p1*x + p0
        double p4 = -alpha * alpha;
        double p3 = -2 * alpha * beta;
        double p2 = -(2 * alpha * gamma + beta * beta) + alpha;
        double p1 = -2 * beta * gamma + beta;
        double p0 = -gamma * gamma + gamma;
        
        // Helper for computing integral of power function - important for correctness
        auto I = [x0, x1](int n) -> double {
            return (std::pow(x1, n) - std::pow(x0, n)) / n;
        };
        
        total += p4 * I(5) + p3 * I(4) + p2 * I(3) + p1 * I(2) + p0 * I(1);
    }
    
    return 2.0 * total;
}

// Helper function for interior projection
double interior_exact_projected_distance(const ProjectionCalculator& calc) {
    auto density = calc.compute_interior_density();
    return expected_abs_diff(density);
}

// Implementation of sample_interior_point method
Polygon::Point ProjectionCalculator::sample_interior_point() const {
    const auto& vertices = polygon_.vertices();
    if (vertices.size() < 3) {
        throw std::invalid_argument("Need at least 3 vertices for sampling");
    }
    
    // Build triangles (fan from vertices[0])
    std::vector<std::array<Polygon::Point, 3>> triangles;
    for (size_t i = 1; i < vertices.size() - 1; ++i) {
        triangles.push_back({vertices[0], vertices[i], vertices[i+1]});
    }
    
    // Compute areas and build cumulative distribution
    std::vector<double> areas;
    areas.reserve(triangles.size());
    
    for (const auto& tri : triangles) {
        const auto& v0 = tri[0];
        const auto& v1 = tri[1];
        const auto& v2 = tri[2];
        
        double cross = (v1[0] - v0[0]) * (v2[1] - v0[1]) - 
                      (v1[1] - v0[1]) * (v2[0] - v0[0]);
        areas.push_back(0.5 * std::abs(cross));
    }
    
    std::vector<double> cum_areas(areas.size());
    std::partial_sum(areas.begin(), areas.end(), cum_areas.begin());
    double total_area = cum_areas.back();
    
    // Select random triangle
    std::uniform_real_distribution<double> dist(0.0, total_area);
    double u = dist(rng_);
    
    size_t tri_idx = 0;
    while (tri_idx < cum_areas.size() && u > cum_areas[tri_idx]) {
        tri_idx++;
    }
    
    // Sample uniformly in triangle using barycentric coordinates
    std::uniform_real_distribution<double> unit_dist(0.0, 1.0);
    double r1 = unit_dist(rng_);
    double r2 = unit_dist(rng_);
    
    // Ensure point is in triangle
    if (r1 + r2 > 1.0) {
        r1 = 1.0 - r1;
        r2 = 1.0 - r2;
    }
    
    const auto& tri = triangles[tri_idx];
    const auto& p0 = tri[0];
    const auto& p1 = tri[1];
    const auto& p2 = tri[2];
    
    // Compute point using barycentric coordinates
    Polygon::Point result = {
        p0[0] + r1 * (p1[0] - p0[0]) + r2 * (p2[0] - p0[0]),
        p0[1] + r1 * (p1[1] - p0[1]) + r2 * (p2[1] - p0[1])
    };
    
    return result;
}

// Helper function for interior Monte Carlo
double interior_monte_carlo_projected_distance(const ProjectionCalculator& calc, int n_samples) {
    double sum = 0.0;
    
    for (int i = 0; i < n_samples; ++i) {
        Polygon::Point p1 = calc.sample_interior_point();
        Polygon::Point p2 = calc.sample_interior_point();
        
        double proj1 = calc.project_point(p1);
        double proj2 = calc.project_point(p2);
        
        sum += std::abs(proj1 - proj2);
    }
    
    return sum / n_samples;
}

} // namespace polygon_projection