#include "RandomGenerators.h"
// was danaRandomGenerators.h

using namespace std;

int RandomGenerators::seed;
mt19937 RandomGenerators::mt_rand;

void RandomGenerators::initRandomGenerator(){
	int chrono_seed = static_cast<int>(chrono::system_clock::now().time_since_epoch().count());
	srand((unsigned)(chrono_seed));
	mt_rand.seed(chrono_seed);
}

void RandomGenerators::setSeed(int seed){
	srand((unsigned)(seed));
	mt_rand.seed(seed);
}



double RandomGenerators::getRandDoubleParamVal(string distributionName, double minVal, double maxVal)
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

int RandomGenerators::getRandIntParamVal(string distributionName, int minVal, int maxVal)
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

double RandomGenerators::roundParamVal(double paramToRound, int numDigAfterDecPoint)
{
	double factor = pow(10.0,numDigAfterDecPoint);
	// double roundedVal = (int)(paramToRound * (int)factor) / factor;
	double roundedVal = (floor(paramToRound * factor) / factor);
	return roundedVal;
}


double RandomGenerators::drawExp(double lambda) {
	
	exponential_distribution<double> distribution(lambda);
	double number = distribution(mt_rand);
	return number;
}

double RandomGenerators::uniform() { // uniform between 0 and 1
	uniform_real_distribution<double> distribution(0.0, 1.0);
	double number = distribution(mt_rand);
	return number;
}

int RandomGenerators::uniform(int a, int b) { // a random number between a and b, including a and b.
	uniform_int_distribution<int> distribution(a, b);
	int number = distribution(mt_rand);
	return number;
}

mt19937 RandomGenerators::getEngine(){
	return mt_rand;
}


RandomGenerators::~RandomGenerators()
{
}
