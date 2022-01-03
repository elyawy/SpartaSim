#include <iostream>
#include "simulator.h"


int main() {

    Vdouble deletions = {0.3};
    Vdouble insertions = {0.4};

    size_t seq_len = 200;

    double insertionRate = 0.05;
    double deletionRate = 0.02;
    // SpartaABC_options::setRLPrior(300, 300);
    // SpartaABC_options::_lengthDistribution = "zipf";

    Simulator sim = Simulator("/home/elyawy/Desktop/alan-2.1.1/RAxML_tree.tree");

    sim.InitSimulator(seq_len, "geometric", deletions, insertions, deletionRate, insertionRate);
    MSA msa = sim.simulateBasedOnTree();

    msa.printMSA();

    for(const auto& value: msa.getStatVec()) {
        std::cout << value << "\t";
    }
    
    return 0;
}