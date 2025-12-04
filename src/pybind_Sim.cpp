#include <pybind11/pybind11.h>
#include <string>

#include "Eigen/Dense"
#include "random_mars.h"
#include "nlohmann/json.hpp"
#include "prof_timer.cpp"

#include "Bead.h"
#include "Bond.h"
#include "DSS_Bond.h"
#include "Cell.cpp"
#include "Grid.h"
#include "Sim.h"

#include "Sim.cpp"
#include "Grid.cpp"
#include "random_mars.cpp"

/*
 * crude implementation for pybind
 * usage: create python module in command line:
 *		>make pybind
 *
 * ...will generate pyticg.cpythonxxx.so file that can be imported in python as:
 *		import pyticg
 *		sim = pytig.Sim()
 *		sim.run()
 *
 * todo: move to Sim.cpp source code and integrate with build system
 */

namespace py = pybind11;

PYBIND11_MODULE(pyticg, m) {
    py::class_<Sim>(m, "Sim")
        .def(py::init<>())
        .def(py::init<const std::string &>())
        .def("run", &Sim::run);
}
