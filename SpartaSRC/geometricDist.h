#include "lengthDistribution.h"
#include "RandomGenerators.h"


class geometricDist: public lengthDistribution 
{
public:
	geometricDist(Vdouble params ,int max): lengthDistribution(params, max) {};
	void generateLengthDistribution();
	int drawLength();
	~geometricDist();
private:
	std::geometric_distribution<int> distribution;
};


