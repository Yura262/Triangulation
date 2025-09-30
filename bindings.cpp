// bindings.cpp
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "3.cpp"  // point this to your real implementation file

namespace py = pybind11;

PYBIND11_MODULE(triangulation, m) {
    m.doc() = "Delaunay triangulation module";

    // Bind Point class
    py::class_<Point>(m, "Point")
        .def(py::init<>())
        .def(py::init<double, double>())
        .def_readwrite("x", &Point::x)
        .def_readwrite("y", &Point::y)
        .def("__repr__", [](const Point &p) {
            return "<Point(" + std::to_string(p.x) + ", " + std::to_string(p.y) + ")>";
        });

    // Bind Triangle class
    py::class_<Triangle>(m, "Triangle")
        .def(py::init<>())
        .def(py::init<Point, Point, Point>())
        .def_readwrite("a", &Triangle::a)
        .def_readwrite("b", &Triangle::b)
        .def_readwrite("c", &Triangle::c)
        .def("centroid", &Triangle::centroid)
        .def("__repr__", [](const Triangle &t) {
            return "<Triangle((" + std::to_string(t.a.x) + "," + std::to_string(t.a.y) + "), (" +
                   std::to_string(t.b.x) + "," + std::to_string(t.b.y) + "), (" +
                   std::to_string(t.c.x) + "," + std::to_string(t.c.y) + "))>";
        });

    // Bind Triangulation class
    py::class_<Triangulation>(m, "Triangulation")
        .def(py::init<>())
        .def("addTriangle", &Triangulation::addTriangle)
        .def("InsertPoint", &Triangulation::InsertPoint)
        .def("getTriangles", &Triangulation::getTriangles)
        .def("DeleteTriangle", &Triangulation::DeleteTriangle)
        .def("RemoveSuperTriangleTriangles", &Triangulation::RemoveSuperTriangleTriangles)
        // expose setters so UI can set quality parameters
        .def("setMinAngle", &Triangulation::setMinAngle)
        .def("setMaxArea", &Triangulation::setMaxArea)
        // direct binding if present in your C++ class
        .def("enforceQuality", &Triangulation::enforceQuality)
        // convenience wrapper: set params then call enforceQuality
        .def("enforceQualityWith", [](Triangulation &self, double angle, double area) {
            self.setMinAngle(angle);
            self.setMaxArea(area);
            self.enforceQuality();
        })
        // boundary extraction: returns std::vector<Point> -> becomes list[Point]
        .def("getBoundaryVertices", &Triangulation::getBoundaryVertices)
        ;

    // Bind the createSuperTriangle function
    m.def("createSuperTriangle", &createSuperTriangle, "Create a super triangle from points");
}
