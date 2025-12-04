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
 * Pybind11 wrapper for SCRIBE simulation engine
 * 
 * Build from src/ directory:
 *     make pybind
 *
 * This generates scribe_engine.cpython-*.so which can be imported as:
 *     from scribe.scribe_engine import Sim
 *     sim = Sim()
 *     sim.run()
 */

namespace py = pybind11;

PYBIND11_MODULE(scribe_engine, m) {
    py::class_<Sim>(m, "Sim")
        .def(py::init<>())
        .def(py::init<const std::string &>())
        .def("run", &Sim::run);
}
