#include "lengthDistribution.h"
#include "RandomGenerators.h"


class poissonDist: public lengthDistribution 
{
public:
	poissonDist(Vdouble params ,int max): lengthDistribution(params, max) {};
	void generateLengthDistribution();
	int drawLength();
	~poissonDist();
private:
	std::poisson_distribution<int> distribution;
};


