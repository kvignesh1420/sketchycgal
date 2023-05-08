#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "sketchycgal/cc/cgal.h"

namespace py = pybind11;

PYBIND11_MODULE(_libsketchycgal, m){
    m.doc() = "api wrapper for the sketchycgal c++ library.";

    // the py::dynamic_attr() tag has been set to enable dynamic attributes on the
    // ApproximateCholesky class at the python layer. However, this can lead to a small
    // runtime cost as a __dict__ is now added to this python class, and the garbage
    // collection becomes a bit expensive.
    py::class_<SketchyCGAL>(m, "SketchyCGAL", py::dynamic_attr())
        .def(py::init<>())
        .def("setup", &SketchyCGAL::setup, py::arg("filepath"), py::arg("max_iters"), 
            py::arg("sketch_rank"), py::arg("tolerance"), 
            "setup the edge_info matrix and precondition the laplacian")
        .def("run", &SketchyCGAL::run,
            "run the randomized algorithm to find the max-cut")
        .def("get_adjacency", &SketchyCGAL::getAdjacencyMatrix,
            "retrieve the sparse adjacency matrix after setup (if necessary)")
        .def("get_laplacian", &SketchyCGAL::getLaplacian,
            "retrieve the sparse laplacian matrix after setup (if necessary)")
        .def("__repr__", 
            [](const SketchyCGAL& a){ return "SketchyCGAL()";}
        );
}
