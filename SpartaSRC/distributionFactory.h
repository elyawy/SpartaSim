#include "geometricDist.h"
#include "poissonDist.h"
#include "zipfDist.h"
// #include "lavaletteDist.h"

using namespace std;


static void getDistribution(lengthDistribution* &l, string distributionName ,Vdouble parameters, int maxLength){
	if(l != nullptr){
		// cout << "l is nullptr" << endl;
		delete l;
		l = nullptr;
	}
	if(distributionName == "geometric") l =  new geometricDist(parameters, maxLength);
	if(distributionName == "poisson") l = new poissonDist(parameters, maxLength);
	if(distributionName == "zipf") l = new zipfDist(parameters, maxLength);
	// if(distributionName == "lavalette") l = new lavaletteDist(parameters, maxLength);

	l->generateLengthDistribution();

}