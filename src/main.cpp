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

           Simulator
           Msa
    )pbdoc";

    py::class_<Simulator>(m, "Sim")
        .def(py::init<const std::string&>())
        .def("init_sim", &Simulator::InitSimulator)
        .def("set_seed", &Simulator::setSeed)
        .def("set_max_indel_size", &Simulator::setMaxIndelSize)
        .def("run_sim", &Simulator::simulateBasedOnTree);

    py::class_<MSA>(m, "Msa")
        .def(py::init<std::string>())
        .def(py::init<const std::vector<std::string> &>())
        .def("get_sum_stats", &MSA::getStatVec, "Return the summary statistics of the MSA, initialized to 0 except MSA length.")
        .def("print_msa", &MSA::printMSA)
        .def("calc_stats", &MSA::recomputeStats, "Compute the summary statistics of the MSA.")
        .def("get_longest_seq_length", &MSA::getMSALongestSeqLength)
        .def("get_shortest_seq_length", &MSA::getMSAshortestSeqLength)
        .def("get_seq_with_indels", &MSA::getAlignedSeqs)
        .def("get_seq_without_indels", &MSA::getUnalignedSeqs)
        .def("get_num_seq", &MSA::getNumberOfSequences);

}
