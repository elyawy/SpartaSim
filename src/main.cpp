#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../SpartaSRC/simulator.h"


namespace py = pybind11;

PYBIND11_MODULE(Sparta, m) {
    m.doc() = R"pbdoc(
        Sparta indel simulator
        -----------------------

        .. currentmodule:: Sparta

        .. autosummary::
           :toctree: _generate

           init_sim
           simulate_based_on_tree
    )pbdoc";

    py::class_<Simulator>(m, "Sim")
        .def(py::init<const std::string&>())
        .def("init_sim", &Simulator::InitSimulator)
        .def("simulate_based_on_tree", &Simulator::simulateBasedOnTree);

    py::class_<MSA>(m, "Msa")
        .def(py::init<std::string>())
        .def(py::init<const std::vector<std::string> &>())
        .def("get_sum_stats", &MSA::getStatVec)
        .def("print_msa", &MSA::printMSA)
        .def("calc_stats", &MSA::recomputeStats)
        .def("get_msa_length", &MSA::getMSAlength)
        .def("get_flag", &MSA::getStatFlag);

}
