#include <pybind11/pybind11.h>
namespace py = pybind11;

PYBIND11_MODULE(asc_ode, m) {
    m.doc() = "tiny asc_ode stub bindings to satisfy CMake configure";
    // no bindings required for the C++ demos; this is only a placeholder
}
