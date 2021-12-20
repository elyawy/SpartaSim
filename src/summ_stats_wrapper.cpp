
#include "summ_stats_wrapper.h"
#include "param_priors.h"


vector<double> simulateSequencesAndReturnSummaryStatistics(size_t randRootLength,
	Vdouble randInsertonParams,
	Vdouble randDeletionParams,
	double randInsertRatio,
	double randDeletionRatio,
	bool isBurnIn);

vector<double> getStatVec(MSA &currMSA) {
	vector<double> statVals;

	statVals.push_back(currMSA.getStatValByType(AVG_GAP_SIZE));
	statVals.push_back(currMSA.getStatValByType(MSA_LEN));
	statVals.push_back(currMSA.getStatValByType(MSA_MAX_LEN));
	statVals.push_back(currMSA.getStatValByType(MSA_MIN_LEN));
	statVals.push_back(currMSA.getStatValByType(TOT_NUM_GAPS));
	statVals.push_back(currMSA.getStatValByType(NUM_GAPS_LEN_ONE));
	statVals.push_back(currMSA.getStatValByType(NUM_GAPS_LEN_TWO));
	statVals.push_back(currMSA.getStatValByType(NUM_GAPS_LEN_THREE));
	statVals.push_back(currMSA.getStatValByType(NUM_GAPS_LEN_AT_LEAST_FOUR));
	statVals.push_back(currMSA.getStatValByType(AVG_UNIQUE_GAP_SIZE));
	statVals.push_back(currMSA.getStatValByType(TOT_NUM_UNIQUE_GAPS));
	statVals.push_back(currMSA.getStatValByType(NUM_GAPS_LEN_ONE_POS_1_GAPS));
	statVals.push_back(currMSA.getStatValByType(NUM_GAPS_LEN_ONE_POS_2_GAPS));
	statVals.push_back(currMSA.getStatValByType(NUM_GAPS_LEN_ONE_POS_N_MINUS_1_GAPS));
	statVals.push_back(currMSA.getStatValByType(NUM_GAPS_LEN_TWO_POS_1_GAPS));
	statVals.push_back(currMSA.getStatValByType(NUM_GAPS_LEN_TWO_POS_2_GAPS));
	statVals.push_back(currMSA.getStatValByType(NUM_GAPS_LEN_TWO_POS_N_MINUS_1_GAPS));
	statVals.push_back(currMSA.getStatValByType(NUM_GAPS_LEN_THREE_POS_1_GAPS));
	statVals.push_back(currMSA.getStatValByType(NUM_GAPS_LEN_THREE_POS_2_GAPS));
	statVals.push_back(currMSA.getStatValByType(NUM_GAPS_LEN_THREE_POS_N_MINUS_1_GAPS));
	statVals.push_back(currMSA.getStatValByType(NUM_GAPS_LEN_AT_LEAST_FOUR_POS_1_GAPS));
	statVals.push_back(currMSA.getStatValByType(NUM_GAPS_LEN_AT_LEAST_FOUR_POS_2_GAPS));
	statVals.push_back(currMSA.getStatValByType(NUM_GAPS_LEN_AT_LEAST_FOUR_POS_N_MINUS_1_GAPS));
	statVals.push_back(currMSA.getStatValByType(MSA_POSITION_WITH_0_GAPS));
	statVals.push_back(currMSA.getStatValByType(MSA_POSITION_WITH_1_GAPS));
	statVals.push_back(currMSA.getStatValByType(MSA_POSITION_WITH_2_GAPS));
	statVals.push_back(currMSA.getStatValByType(MSA_POSITION_WITH_N_MINUS_1_GAPS));

	return statVals;
}

vector<double> getStatVec(string inputFile) { // this function return the summary statistics vector from an input MSA read from a file
	MSA curr_MSA(inputFile);
	return getStatVec(curr_MSA);
}

string getSpecificHeader() {
	string dist = SpartaABC_options::_lengthDistribution;
	if(dist == "zipf") {
		return "AIR\tADR\t";
	}
	if(dist == "lavalette"){
		return "aLavaletteI\tCLavaletteI\taLavaletteD\tCLavaletteD\t";
	}
	if(dist == "poisson"){
		return "LambdaPoissonI\tLambdaPoissonD\t";
	}
	if(dist == "geometric"){
		return "pGeometricI\tpGeometricD\t";
	}
	return "AIR\tADR\t";
}


string NiceHeader() {
	string niceHeader = "DISTANCE\tRL\t";
	niceHeader += getSpecificHeader();
	niceHeader += "IR\tDR\tAVG_GAP_SIZE\tMSA_LEN\tMSA_MAX_LEN\tMSA_MIN_LEN\tTOT_NUM_GAPS\tNUM_GAPS_LEN_ONE\tNUM_GAPS_LEN_TWO\tNUM_GAPS_LEN_THREE\tNUM_GAPS_LEN_AT_LEAST_FOUR\tAVG_UNIQUE_GAP_SIZE\tTOT_NUM_UNIQE_GAPS\tNUM_GAPS_LEN_ONE_POS_1_GAPS\tNUM_GAPS_LEN_ONE_POS_2_GAPS\tNUM_GAPS_LEN_ONE_POS_N_MINUS_1_GAPS\tNUM_GAPS_LEN_TWO_POS_1_GAPS\tNUM_GAPS_LEN_TWO_POS_2_GAPS\tNUM_GAPS_LEN_TWO_POS_N_MINUS_1_GAPS\tNUM_GAPS_LEN_THREE_POS_1_GAPS\tNUM_GAPS_LEN_THREE_POS_2_GAPS\tNUM_GAPS_LEN_THREE_POS_N_MINUS_1_GAPS\tNUM_GAPS_LEN_AT_LEAST_FOUR_POS_1_GAPS\tNUM_GAPS_LEN_AT_LEAST_FOUR_POS_2_GAPS\tNUM_GAPS_LEN_AT_LEAST_FOUR_POS_N_MINUS_1_GAPS\tMSA_POSITION_WITH_0_GAPS\tMSA_POSITION_WITH_1_GAPS\tMSA_POSITION_WITH_2_GAPS\tMSA_POSITION_WITH_N_MINUS_1_GAPS\n";
	
	return niceHeader;
}

Vdouble getRandLengthDistributionParameters(){
	Vdouble distParams;
	if (SpartaABC_options::_lengthDistribution == "zipf"){
		double aParam = 1.0;
		aParam = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeA, SpartaABC_options::_minAIVal, SpartaABC_options::_maxAIVal);
		distParams.push_back(aParam);
	}
	if (SpartaABC_options::_lengthDistribution == "lavalette"){
		double aParam = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeA, SpartaABC_options::_minALavalette, SpartaABC_options::_maxALavalette);
		double cParam = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeA, SpartaABC_options::_minCLavalette, SpartaABC_options::_maxCLavalette);
		distParams.push_back(aParam);
		distParams.push_back(cParam);
	}
	if (SpartaABC_options::_lengthDistribution == "poisson"){
		double lambdaParam = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeA, SpartaABC_options::_minLambdaPoisson, SpartaABC_options::_maxLambdaPoisson);
		distParams.push_back(lambdaParam);
	}
	if (SpartaABC_options::_lengthDistribution == "geometric"){
		double pParam = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeA, SpartaABC_options::_minPGeometric, SpartaABC_options::_maxPGeometric);
		distParams.push_back(pParam);
	}

	return distParams;
}


vector<double> getWeightsVectorUsingSimulations() {
	// simulate datasets from the prior
	const size_t numberOfSimulations = SpartaABC_options::_numBurnIn;//temp, suppose to be 10,000
	size_t real_numberOfSimulations = numberOfSimulations;
	vector<double> summaryStatisticsSum;
	vector<double> summaryStatisticsSquareSum;
	vector<double> summStatistics;

	vector<double> summaryStatisticsSum_upgraded;
	vector<double> summaryStatisticsSquareSum_upgraded;
	vector<double> summStatistics_upgraded;
	for (size_t i = 0; i < numberOfSimulations; ++i) {
		
		if ((i % 1000) == 0) {
			cout << "burn in simulation number " << i << endl;
		}
		int randRootLength = getRandIntParamVal(SpartaABC_options::_priorDistTypeRL, SpartaABC_options::_minRLVal, SpartaABC_options::_maxRLVal);
	
		Vdouble	lengthDistParamsInsertion = getRandLengthDistributionParameters();
		Vdouble	lengthDistParamsDeletion = getRandLengthDistributionParameters();
	
		double randInsertRatio = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeIR, SpartaABC_options::_minIRVal, SpartaABC_options::_maxIRVal);
		double randDeletionRatio = randInsertRatio;//getRandDoubleParamVal(SpartaABC_options::_priorDistTypeIR, SpartaABC_options::_minIRVal, SpartaABC_options::_maxIRVal);
		while (randInsertRatio == 0/* && randDeletionRatio == 0*/) {
			randInsertRatio = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeIR, SpartaABC_options::_minIRVal, SpartaABC_options::_maxIRVal);
			randDeletionRatio = randInsertRatio;// getRandDoubleParamVal(SpartaABC_options::_priorDistTypeIR, SpartaABC_options::_minIRVal, SpartaABC_options::_maxIRVal);
		}
			summStatistics = simulateSequencesAndReturnSummaryStatistics(randRootLength,
				lengthDistParamsInsertion,
				lengthDistParamsDeletion,
				randInsertRatio,
				randDeletionRatio,
				true);
			if (summStatistics[0]==-1)
			{
				real_numberOfSimulations--;
			}
			else
			{
				if (i == 0) {
					summaryStatisticsSum.resize(summStatistics.size());
					summaryStatisticsSquareSum.resize(summStatistics.size());


				}
				for (size_t j = 0; j < summStatistics.size(); ++j) {
					summaryStatisticsSum[j] += summStatistics[j];
					summaryStatisticsSquareSum[j] += (summStatistics[j] * summStatistics[j]);

				}
			}
		
	}// end of simulations
	for (size_t j = 0; j < summStatistics.size(); ++j) {
		summaryStatisticsSum[j] /= real_numberOfSimulations;
		summaryStatisticsSquareSum[j] /= real_numberOfSimulations;

		
	}
	
	vector<double> variance(summaryStatisticsSum.size());
	for (size_t j = 0; j < variance.size(); ++j) {
		variance[j] = summaryStatisticsSquareSum[j] - summaryStatisticsSum[j] * summaryStatisticsSum[j];
		variance[j] = sqrt(variance[j]); // now it is standard deviation
		if (variance[j] != 0) {
			variance[j] = 1.0 / variance[j]; // now it is weights.
		}
	}
	//here we pushing our variance into the variance vector will be deleted //29.3 OI

	return variance; // now the weights are updated.
}


vector<double> getWeightsVector() {
	double wAvgGapSize = SpartaABC_options::_wAvgGapSize;
	double wMSALen = SpartaABC_options::_wMSALen;
	double wMSAMax = SpartaABC_options::_wMSAMax;
	double wMSAMin = SpartaABC_options::_wMSAMin;
	double wTotNumGaps = SpartaABC_options::_wTotNumGaps;
	double wNumGapsLenOne = SpartaABC_options::_wNumGapsLenOne;
	double wNumGapsLenTwo = SpartaABC_options::_wNumGapsLenTwo;
	double wNumGapsLenThree = SpartaABC_options::_wNumGapsLenThree;
	double wNumGapsLenAtLeastFour = SpartaABC_options::_wNumGapsLenAtLeastFour;
	double wAvgUniqueGapSize = SpartaABC_options::_wAvgUniqueGapSize;
	double wTotNumUniqueGaps = SpartaABC_options::_wTotNumUniqueGaps;

	double wNumGapsLenOneIn1Pos = SpartaABC_options::_wNumGapsLenOneIn1Pos;
	double wNumGapsLenOneIn2Pos = SpartaABC_options::_wNumGapsLenOneIn2Pos;
	double wNumGapsLenOneInNMinus1Pos = SpartaABC_options::_wNumGapsLenOneInNMinus1Pos;
	double wNumGapsLenTwoIn1Pos = SpartaABC_options::_wNumGapsLenTwoIn1Pos;
	double wNumGapsLenTwoIn2Pos = SpartaABC_options::_wNumGapsLenTwoIn2Pos;
	double wNumGapsLenTwoInNMinus1Pos = SpartaABC_options::_wNumGapsLenTwoInNMinus1Pos;
	double wNumGapsLenThreeIn1Pos = SpartaABC_options::_wNumGapsLenThreeIn1Pos;
	double wNumGapsLenThreeIn2Pos = SpartaABC_options::_wNumGapsLenThreeIn2Pos;
	double wNumGapsLenThreeInNMinus1Pos = SpartaABC_options::_wNumGapsLenThreeInNMinus1Pos;
	double wNumGapsLenAtLeastFourIn1Pos = SpartaABC_options::_wNumGapsLenAtLeastFourIn1Pos;
	double wNumGapsLenAtLeastFourIn2Pos = SpartaABC_options::_wNumGapsLenAtLeastFourIn2Pos;
	double wNumGapsLenAtLeastFourInNMinus1Pos = SpartaABC_options::_wNumGapsLenAtLeastFourInNMinus1Pos;

	double wNumberOfMSA_position_with_0_gaps = SpartaABC_options::_wNumberOfMSA_position_with_0_gaps;
	double wNumberOfMSA_position_with_1_gaps = SpartaABC_options::_wNumberOfMSA_position_with_1_gaps;
	double wNumberOfMSA_position_with_2_gaps = SpartaABC_options::_wNumberOfMSA_position_with_2_gaps;
	double wNumberOfMSA_position_with_n_minus_1_gaps = SpartaABC_options::_wNumberOfMSA_position_with_n_minus_1_gaps;

	if ((wAvgGapSize == -1.0) && (wMSALen == -1.0) && (wMSAMax == -1.0) && (wMSAMin == -1.0) && (wTotNumGaps == -1.0)
		&& (wNumGapsLenOne == -1.0) && (wNumGapsLenTwo == -1.0) && (wNumGapsLenThree == -1.0) && (wNumGapsLenAtLeastFour == -1.0)
		&& (wAvgUniqueGapSize == -1.0) && (wTotNumUniqueGaps == -1.0) && (wNumGapsLenOneIn1Pos == -1.0) && (wNumGapsLenOneIn2Pos == -1.0)
		&& (wNumGapsLenOneInNMinus1Pos == -1.0) && (wNumGapsLenTwoIn1Pos == -1.0)&& (wNumGapsLenTwoIn2Pos == -1.0) 
		&& (wNumGapsLenTwoInNMinus1Pos == -1.0) && (wNumGapsLenThreeIn1Pos == -1.0) && (wNumGapsLenThreeIn2Pos == -1.0) 
		&& (wNumGapsLenThreeInNMinus1Pos == -1.0) && (wNumGapsLenAtLeastFourIn1Pos == -1.0) && (wNumGapsLenAtLeastFourIn2Pos == -1.0) 
		&& (wNumGapsLenAtLeastFourInNMinus1Pos == -1.0) &&(wNumberOfMSA_position_with_0_gaps == -1.0)
		&& (wNumberOfMSA_position_with_1_gaps == -1.0) && (wNumberOfMSA_position_with_2_gaps == -1.0) && (wNumberOfMSA_position_with_n_minus_1_gaps == -1.0))
	{
		return getWeightsVectorUsingSimulations();
	}
	else if ((wAvgGapSize == -1.0) || (wMSALen == -1.0) || (wMSAMax == -1.0) || (wMSAMin == -1.0) || (wTotNumGaps == -1.0)
		|| (wNumGapsLenOne == -1.0) || (wNumGapsLenTwo == -1.0) || (wNumGapsLenThree == -1.0) || (wNumGapsLenAtLeastFour == -1.0)
		|| (wAvgUniqueGapSize == -1.0) || (wTotNumUniqueGaps == -1.0) || (wNumGapsLenOneIn1Pos == -1.0) || (wNumGapsLenOneIn2Pos == -1.0)
		|| (wNumGapsLenOneInNMinus1Pos == -1.0) || (wNumGapsLenTwoIn1Pos == -1.0) || (wNumGapsLenTwoIn2Pos == -1.0)
		|| (wNumGapsLenTwoInNMinus1Pos == -1.0) || (wNumGapsLenThreeIn1Pos == -1.0) || (wNumGapsLenThreeIn2Pos == -1.0)
		|| (wNumGapsLenThreeInNMinus1Pos == -1.0) || (wNumGapsLenAtLeastFourIn1Pos == -1.0) || (wNumGapsLenAtLeastFourIn2Pos == -1.0)
		|| (wNumGapsLenAtLeastFourInNMinus1Pos == -1.0)||(wNumberOfMSA_position_with_0_gaps == -1.0)
		|| (wNumberOfMSA_position_with_1_gaps == -1.0) || (wNumberOfMSA_position_with_2_gaps == -1.0) || (wNumberOfMSA_position_with_n_minus_1_gaps == -1.0))
	{
		cout << " error. Either all weights are given or all of them not. You cannot give only some of the weights";
		exit(4);
	}
	else { // read from the input file
		vector<double> _summStatWeights;
		_summStatWeights.push_back(wAvgGapSize);
		_summStatWeights.push_back(wMSALen);
		_summStatWeights.push_back(wMSAMax);
		_summStatWeights.push_back(wMSAMin);
		_summStatWeights.push_back(wTotNumGaps);
		_summStatWeights.push_back(wNumGapsLenOne);
		_summStatWeights.push_back(wNumGapsLenTwo);
		_summStatWeights.push_back(wNumGapsLenThree);
		_summStatWeights.push_back(wNumGapsLenAtLeastFour);
		_summStatWeights.push_back(wAvgUniqueGapSize);
		_summStatWeights.push_back(wTotNumUniqueGaps);

		_summStatWeights.push_back(wNumGapsLenOneIn1Pos);
		_summStatWeights.push_back(wNumGapsLenOneIn2Pos);
		_summStatWeights.push_back(wNumGapsLenOneInNMinus1Pos);
		_summStatWeights.push_back(wNumGapsLenTwoIn1Pos);
		_summStatWeights.push_back(wNumGapsLenTwoIn2Pos);
		_summStatWeights.push_back(wNumGapsLenTwoInNMinus1Pos);
		_summStatWeights.push_back(wNumGapsLenThreeIn1Pos);
		_summStatWeights.push_back(wNumGapsLenThreeIn2Pos);
		_summStatWeights.push_back(wNumGapsLenThreeInNMinus1Pos);
		_summStatWeights.push_back(wNumGapsLenAtLeastFourIn1Pos);
		_summStatWeights.push_back(wNumGapsLenAtLeastFourIn2Pos);
		_summStatWeights.push_back(wNumGapsLenAtLeastFourInNMinus1Pos);

		_summStatWeights.push_back(wNumberOfMSA_position_with_0_gaps);
		_summStatWeights.push_back(wNumberOfMSA_position_with_1_gaps);
		_summStatWeights.push_back(wNumberOfMSA_position_with_2_gaps);
		_summStatWeights.push_back(wNumberOfMSA_position_with_n_minus_1_gaps);

		return _summStatWeights;
	}

}
