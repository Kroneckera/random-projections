#ifndef POLYGON_PROJECTION_INTEGRATION_H
#define POLYGON_PROJECTION_INTEGRATION_H

#include <vector>
#include <functional>
#include "polygon.h"

namespace polygon_projection {

/**
 * Calculate break angles where projection order changes.
 * 
 * @param vertices Polygon vertices
 * @param tol Tolerance for detecting coincident points
 * @return Vector of angles in [0,π] where projections coincide
 */
std::vector<double> calculate_break_angles(const Polygon::VertexList& vertices, double tol = 1e-12);

/**
 * Integrate a function over angles [0,π] using Gauss-Legendre quadrature,
 * accounting for breakpoints where the function may have kinks.
 * 
 * @param vertices Polygon vertices
 * @param f Function to integrate, takes (vertices, angle) as arguments
 * @param order Gauss-Legendre quadrature order per smooth interval
 * @return The integral value
 */
double integrate_over_directions(
    const Polygon::VertexList& vertices,
    const std::function<double(const Polygon::VertexList&, double)>& f,
    int order = 16
);

} // namespace polygon_projection

#endif // POLYGON_PROJECTION_INTEGRATION_H