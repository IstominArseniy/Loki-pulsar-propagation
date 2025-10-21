#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../lib/CppInterfaceClass.h"

namespace py = pybind11;

// PYBIND11_MODULE(cpp_interface, m, py::mod_gil_not_used()) {
//     py::class_<CppInterface>(m, "Interface")
//         .def(py::init<const std::map<std::string, double>, const std::map<std::string, double> >());
// }

int add(int i, int j) {
    return i+j;
}

PYBIND11_MODULE(cpp_interface, m, py::mod_gil_not_used()) {
    m.doc() = "pybind11 example plugin"; // optional module docstring
    py::class_<CppInterface>(m, "ProfileCalculator")
        .def(py::init<const std::map<std::string, double>, const std::map<std::string, double> >())
        .def("calculate_profile", &CppInterface::calculate_profile);
}