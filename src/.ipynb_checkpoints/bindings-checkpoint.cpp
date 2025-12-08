#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "nonlinfunc.hpp"
#include "Newton.hpp"
#include <vector.hpp>
#include <matrix.hpp>

namespace py = pybind11;
using namespace ASC_ode;
using namespace nanoblas;

PYBIND11_MODULE(asc_ode, m) {
    m.doc() = "ASC ODE solver";

    // Vector
    py::class_<Vector<double>>(m, "Vector")
        .def(py::init<size_t>())
        .def(py::init<std::vector<double>>())
        .def("__len__", [](const Vector<double>& v) { return v.size(); })
        .def("__getitem__", [](Vector<double>& v, size_t i) -> double& { return v(i); })
        .def("__setitem__", [](Vector<double>& v, size_t i, double val) { v(i) = val; })
        .def("__repr__", [](const Vector<double>& v) {
            std::stringstream ss; ss << "[";
            for(size_t i=0; i<v.size(); i++) { ss << v(i); if(i+1<v.size()) ss << ", "; }
            ss << "]"; return ss.str();
        });

    // NonlinearFunction base + Python subclass
    py::class_<NonlinearFunction, PyNonlinearFunction, std::shared_ptr<NonlinearFunction>>(m, "NonlinearFunction", py::dynamic_attr(), py::module_local())
        .def(py::init<>())
        .def("dimX", &NonlinearFunction::dimX)
        .def("dimF", &NonlinearFunction::dimF)
        .def("evaluate", [](NonlinearFunction& self, Vector<double>& x, Vector<double>& f) { self.evaluate(x,f); })
        .def("evaluateDeriv", [](NonlinearFunction& self, Vector<double>& x, Matrix<double>& df) { self.evaluateDeriv(x,df); });

    // NewtonSolver
    m.def("NewtonSolver", [](std::shared_ptr<NonlinearFunction> func, Vector<double>& x,
                             double tol = 1e-12, int maxit = 50, bool linesearch = true) {
        return NewtonSolver(func, x, tol, maxit, linesearch);
    });
}
