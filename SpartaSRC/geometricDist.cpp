#include "geometricDist.h"

using namespace std;

void geometricDist::generateLengthDistribution() {
    geometric_distribution<int> d(this->parameters[0]);
	distribution = d;
}


int geometricDist::drawLength() {
	mt19937 generator = RandomGenerators::getEngine();
	int length  = distribution(generator)+1;
	while (length > this->maxLength){
		length  = distribution(generator)+1;
	}
	return length;
}

geometricDist::~geometricDist() {
	
}