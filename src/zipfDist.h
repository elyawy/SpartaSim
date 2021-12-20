#include "lengthDistribution.h"
#include "param_priors.h"
#include "FastZip.h"

class zipfDist: public lengthDistribution 
{
public:
	zipfDist(Vdouble params ,int max): lengthDistribution(params, max), _fastZip(params[0], max){};
	void generateLengthDistribution();
	int drawLength();
	~zipfDist();
private:
    FastZip _fastZip;
};


