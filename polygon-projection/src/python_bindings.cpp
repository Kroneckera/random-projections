#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include "polygon_projection/polygon.h"
#include "polygon_projection/projection.h"

namespace py = pybind11;
using namespace polygon_projection;

// Convert NumPy array to vertex list
Polygon::VertexList vertices_from_numpy(py::array_t<double> array) {
    if (array.ndim() != 2 || array.shape(1) != 2) {
        throw std::invalid_argument("Vertices must be a 2D array with shape (n, 2)");
    }
    
    const size_t n_points = array.shape(0);
    Polygon::VertexList vertices;
    vertices.reserve(n_points);
    
    auto r = array.unchecked<2>();
    for (size_t i = 0; i < n_points; ++i) {
        vertices.push_back({r(i, 0), r(i, 1)});
    }
    
    return vertices;
}

// Convert vertex list to NumPy array
py::array_t<double> vertices_to_numpy(const Polygon::VertexList& vertices) {
    py::ssize_t n = static_cast<py::ssize_t>(vertices.size());
    py::array_t<double> result(std::vector<py::ssize_t>{n, 2});
    auto r = result.mutable_unchecked<2>();
    
    for (size_t i = 0; i < vertices.size(); ++i) {
        r(i, 0) = vertices[i][0];
        r(i, 1) = vertices[i][1];
    }
    
    return result;
}

// Module definition
PYBIND11_MODULE(_core, m) {
    m.doc() = "C++ backend for polygon projection calculations";
    m.attr("__version__") = "0.1.0";
    
    // Enum: Region
    py::enum_<Region>(m, "Region")
        .value("INTERIOR", Region::INTERIOR)
        .value("BOUNDARY", Region::BOUNDARY);
    
    // Class: Polygon
    py::class_<Polygon>(m, "Polygon")
        // Constructors
        .def(py::init([](py::array_t<double> vertices) {
            return Polygon(vertices_from_numpy(vertices));
        }))
        
        // Static methods
        .def_static("regular", &Polygon::regular,
                   py::arg("sides"), py::arg("radius") = 1.0)
        
        // Properties
        .def("vertices", &Polygon::vertices)
        .def("area", &Polygon::area)
        .def("perimeter", &Polygon::perimeter)
        .def("is_convex", &Polygon::is_convex)
        .def("contains", &Polygon::contains);
    
    // Class: ProjectionCalculator
    py::class_<ProjectionCalculator>(m, "ProjectionCalculator")
        // Constructor
        .def(py::init([](const Polygon& polygon, py::list direction_list, Region region, unsigned int seed) {
            // Convert py::list to Direction
            if (py::len(direction_list) != 2)
                throw std::invalid_argument("Direction must be a 2D vector");
                
            ProjectionCalculator::Direction direction = {
                direction_list[0].cast<double>(),
                direction_list[1].cast<double>()
            };
            
            return ProjectionCalculator(polygon, direction, region, seed);
        }),
        py::arg("polygon"),
        py::arg("direction"),
        py::arg("region") = Region::INTERIOR,
        py::arg("seed") = 0)
        
        // Methods for distance calculation
        .def("exact_projected_distance", &ProjectionCalculator::exact_projected_distance)
        .def("monte_carlo_projected_distance", &ProjectionCalculator::monte_carlo_projected_distance,
             py::arg("n_samples") = 10000)
        .def("exact_average_distance", &ProjectionCalculator::exact_average_distance,
             py::arg("integration_order") = 16)
        .def("monte_carlo_average_distance", &ProjectionCalculator::monte_carlo_average_distance,
             py::arg("n_samples") = 10000)
        
        // Sample points
        .def("sample_point", &ProjectionCalculator::sample_point)
        .def("project_point", &ProjectionCalculator::project_point)
        
        // Sample multiple points
        .def("sample_points", [](const ProjectionCalculator& self, int count) {
            py::array_t<double> result(std::vector<py::ssize_t>{count, 2});
            auto r = result.mutable_unchecked<2>();
            
            for (int i = 0; i < count; ++i) {
                auto point = self.sample_point();
                r(i, 0) = point[0];
                r(i, 1) = point[1];
            }
            
            return result;
        }, py::arg("count") = 1)
        
        // Project multiple points
        .def("project_points", [](const ProjectionCalculator& self, py::array_t<double> points) {
            if (points.ndim() != 2 || points.shape(1) != 2) {
                throw std::invalid_argument("Points must be a 2D array with shape (n, 2)");
            }
            
            const size_t n_points = points.shape(0);
            py::array_t<double> result(std::vector<py::ssize_t>{static_cast<py::ssize_t>(n_points)});
            
            auto r_points = points.unchecked<2>();
            auto r_result = result.mutable_unchecked<1>();
            
            for (size_t i = 0; i < n_points; ++i) {
                Polygon::Point point = {r_points(i, 0), r_points(i, 1)};
                r_result(i) = self.project_point(point);
            }
            
            return result;
        })
        
        // Monte Carlo batched with results dict
        .def("run_monte_carlo_batched", [](const ProjectionCalculator& self, 
                                          int n_samples, 
                                          int batch_size) {
            auto results = self.run_monte_carlo_batched(n_samples, batch_size);
            
            py::dict result_dict;
            result_dict["projected_distance"] = results.projected_distance;
            result_dict["average_distance"] = results.average_distance;
            result_dict["samples_processed"] = results.samples_processed;
            
            return result_dict;
        }, py::arg("n_samples"), py::arg("batch_size") = 100);
}