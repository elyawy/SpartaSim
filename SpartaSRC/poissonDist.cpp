#include "poissonDist.h"

using namespace std;

void poissonDist::generateLengthDistribution() {
    poisson_distribution<int> d(this->parameters[0]);
	distribution = d;
}


int poissonDist::drawLength() {
	mt19937 generator = RandomGenerators::getEngine();
	int length  = distribution(generator)+1;
	while (length > this->maxLength){
		length  = distribution(generator)+1;
	}
	return length;
}

poissonDist::~poissonDist() {
	
}