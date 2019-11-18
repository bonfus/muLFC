#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

void init_lattice(py::module &);
void init_dipolesum(py::module &);

PYBIND11_MODULE(lfclib, m) {
    m.doc() = "lfclib python bindings";
    
    init_lattice(m);
    init_dipolesum(m);
}
