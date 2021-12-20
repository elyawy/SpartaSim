#include "poissonDist.h"

using namespace std;

void poissonDist::generateLengthDistribution() {
	mt19937 gen(std::rand());
	generator = gen;
    poisson_distribution<int> d(this->parameters[0]);
	distribution = d;
}


int poissonDist::drawLength() {
	int length  = distribution(generator)+1;
	while (length > this->maxLength){
		length  = distribution(generator)+1;
	}
	return length;
}

poissonDist::~poissonDist() {
	
}