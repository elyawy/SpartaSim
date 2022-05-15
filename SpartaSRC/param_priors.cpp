#include "param_priors.h"
#include <chrono>


int seed = static_cast<int>(chrono::system_clock::now().time_since_epoch().count());
default_random_engine generator(seed);
mt19937 mt_rand(seed);


double getRandDoubleParamVal(string distributionName, double minVal, double maxVal)
{
	if(distributionName == "uniform")
	{
		double rand_unif_between_0_and_1 = ((double)rand() / (double)RAND_MAX);
		double rand_unif_in_range = minVal + rand_unif_between_0_and_1 * (maxVal - minVal);
		double rand_unif_in_range_rounded = roundParamVal(rand_unif_in_range,9);
		return rand_unif_in_range_rounded;
	}
	else
	{
		//at the moment - no other implementation
		return -1.0;
	}
}

int getRandIntParamVal(string distributionName, int minVal, int maxVal)
{
	if(distributionName == "uniform")
	{
		int rand_int_0_to_max = rand();
		int rand_int_unif_in_range = minVal + (rand_int_0_to_max % (int)(maxVal - minVal + 1));
		return rand_int_unif_in_range;
	}
	else
	{
		//at the moment - no other implementation
		return -1;
	}
}

double roundParamVal(double paramToRound, int numDigAfterDecPoint)
{
	double factor = pow(10.0,numDigAfterDecPoint);
	// double roundedVal = (int)(paramToRound * (int)factor) / factor;
	double roundedVal = (floor(paramToRound * factor) / factor);
	return roundedVal;
}


void resetSeed(size_t new_seed) {
	srand(new_seed);
	generator.seed(new_seed);
	mt_rand.seed(new_seed);
}


double drawExp(double lambda) {
	
	exponential_distribution<double> distribution(lambda);
	double number = distribution(generator);
	return number;
}
double uniform() { // uniform between 0 and 1
	uniform_real_distribution<double> distribution(0.0, 1.0);
	double number = distribution(generator);
	return number;
}

int uniform(int a, int b) { // a random number between a and b, including a and b.
	uniform_int_distribution<int> distribution(a, b);
	int number = distribution(generator);
	return number;
}
