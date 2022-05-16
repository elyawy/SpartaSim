#pragma once

#include <algorithm>
#include <cmath>
#include <chrono>
#include <random>
#include <cassert>
using namespace std;


class RandomGenerators
{
private:
	static int seed;
	static mt19937 mt_rand;
public:
	static void initRandomGenerator();
	static void setSeed(int seed);

	static double getRandDoubleParamVal(string distributionName, double minVal, double maxVal);
	static int getRandIntParamVal(string distributionName, int minVal, int maxVal);
	static double roundParamVal(double paramToRound, int numDigAfterDecPoint);


	static double drawExp(double lambda);
	static double uniform();
	static int uniform(int a, int b);

	static mt19937 getEngine();

	virtual ~RandomGenerators();
};



