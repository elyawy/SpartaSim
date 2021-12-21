#ifndef _SUMM_STATS_WRAPPER
#define _SUMM_STATS_WRAPPER
#include <vector>
#include "MSA.h"
#include "MSA_decomposer.h"
#include "SpartaABC_options.h"
#include "read_seqs.h"
#include "definitions.h"

using namespace std;


vector<double> getStatVec(MSA &currMSA);// This function gets an MSA as input and return a vector of summary statistics.
vector<double> getStatVec(string inputFile); // this function return the summary statistics vector from an input MSA read from a file
string getSpecificHeader();
string NiceHeader(); // prints a nice header in the output file
Vdouble generateLengthDistributionParameters();
vector<double> getWeightsVector(); // gets a weight vector either from an input file or using simulations.

#endif