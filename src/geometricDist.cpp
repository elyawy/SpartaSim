#include "geometricDist.h"

using namespace std;

void geometricDist::generateLengthDistribution() {
	mt19937 gen(std::rand());
	generator = gen;
    geometric_distribution<int> d(this->parameters[0]);
	distribution = d;
}


int geometricDist::drawLength() {
	int length  = distribution(generator)+1;
	while (length > this->maxLength){
		length  = distribution(generator)+1;
	}
	return length;
}

geometricDist::~geometricDist() {
	
}