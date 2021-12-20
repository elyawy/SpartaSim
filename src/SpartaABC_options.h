#ifndef _SPARTAABC_OPTIONS
#define _SPARTAABC_OPTIONS

#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cfloat>
#include "errorMsg.h"
#include "Parameters.h"

using namespace std;

class SpartaABC_options {
		
public:
    static void initOptions(const string& paramFileName);
	virtual ~SpartaABC_options();	
	static void initDefault();
	static void getParamsFromFile(const string& paramFileName);
	static void setRLPrior(int msa_min_len, int msa_max_len);
	static void decreaseAlignmentParam();

public:
	//specific output config
	static bool _only_real_stats;
	static bool _predefinedMode;
	static string _predefinedSimsFile;


	static int _alignments_output;
	static string _outputAlignmnetsFile;

	static int _numSimulations;
	static int _numBurnIn;

	//summary statistics weigths (0 = not used, 1 for all = unweighted, any positive double = given weight, -1 = estimate [not implmented yet])
	static double _wAvgGapSize; 
	static double _wAvgUniqueGapSize;
	static double _wMSALen;
	static double _wMSAMax;
	static double _wMSAMin; 
	static double _wTotNumGaps; // total number of gaps
	static double _wNumGapsLenOne; // number of gaps of length one
	static double _wNumGapsLenTwo; // number of gaps of length two
	static double _wNumGapsLenThree; // number of gaps of length three
	static double _wNumGapsLenAtLeastFour;
	static double _wTotNumUniqueGaps; // total number of unique gaps

	static double _wNumGapsLenOneIn1Pos;
	static double _wNumGapsLenOneIn2Pos;
	static double _wNumGapsLenOneInNMinus1Pos;
	static double _wNumGapsLenTwoIn1Pos;
	static double _wNumGapsLenTwoIn2Pos;
	static double _wNumGapsLenTwoInNMinus1Pos;
	static double _wNumGapsLenThreeIn1Pos;
	static double _wNumGapsLenThreeIn2Pos;
	static double _wNumGapsLenThreeInNMinus1Pos;
	static double _wNumGapsLenAtLeastFourIn1Pos;
	static double _wNumGapsLenAtLeastFourIn2Pos;
	static double _wNumGapsLenAtLeastFourInNMinus1Pos;

	static double _wNumberOfMSA_position_with_0_gaps;
	static double _wNumberOfMSA_position_with_1_gaps;
	static double _wNumberOfMSA_position_with_2_gaps;
	static double _wNumberOfMSA_position_with_n_minus_1_gaps;

	//template control
	static string _indelibleTemplateControlFile;
	static string _dawgTemplateControlFile;

	//input real MSA
	static string _inputRealMSAFile;

	//input tree file name
	static string _inputTreeFileName;

	//output collected samples
	static string _outputGoodParamsFile;

	//parameters prior options
	static string _priorDistTypeRL;
	static string _priorDistTypeA;
	static string _priorDistTypeIR;

	static int _minRLVal;
	static int _maxRLVal;

	static string _modelType;

	static int _numDistributionParams;

	static double _minAIVal;
	static double _maxAIVal;

	static double _minADVal;
	static double _maxADVal;


	static double _minALavalette;
	static double _maxALavalette;

	static double _minCLavalette;
	static double _maxCLavalette;

	static double _minLambdaPoisson;
	static double _maxLambdaPoisson;

	static double _minPGeometric;
	static double _maxPGeometric;


	static double _minIRVal;
	static double _maxIRVal;
	static double _minDRVal;
	static double _maxDRVal;



	//ABC Reject params
	static int _numberOfSamplesToKeep;
	static double _distanceCutOff;

	//pairwise or multiple
	static int _alignmentMode;

	//similarity - simple / nucleotide / aa
	static int _similarityMode;

	// should dawg simulator be used (defaults to 'false')
	static bool _dawgSimulator;

	// temporary- can be deleted.
	static int  boomed_times;
	static double AVG_time1;
	static double AVG_time2;

	//indel length distribution
	static string _lengthDistribution;

};


#endif