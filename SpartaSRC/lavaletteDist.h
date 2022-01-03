#include "lengthDistribution.h"
#include "param_priors.h"


class lavaletteDist: public lengthDistribution 
{
public:
	lavaletteDist(Vdouble params ,int max): lengthDistribution(params, max) {};
	void generateLengthDistribution();
	int drawLength();
	~lavaletteDist();
private:
	// for truncated distribution:
	MDOUBLE left_limit;
	MDOUBLE right_limit;

};


