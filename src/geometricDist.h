#include "lengthDistribution.h"
#include "param_priors.h"


class geometricDist: public lengthDistribution 
{
public:
	geometricDist(Vdouble params ,int max): lengthDistribution(params, max) {};
	void generateLengthDistribution();
	int drawLength();
	~geometricDist();
private:
    std::mt19937 generator;
	std::geometric_distribution<int> distribution;
};


