#include "lengthDistribution.h"
#include "param_priors.h"


class poissonDist: public lengthDistribution 
{
public:
	poissonDist(Vdouble params ,int max): lengthDistribution(params, max) {};
	void generateLengthDistribution();
	int drawLength();
	~poissonDist();
private:
    std::mt19937 generator;
	std::poisson_distribution<int> distribution;
};


