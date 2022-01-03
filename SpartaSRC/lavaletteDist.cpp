#include "lavaletteDist.h"

using namespace std;

void lavaletteDist::generateLengthDistribution() {
	this->left_limit = 1.0 - 1.0/(1.0+pow((1-0.5)/this->parameters[1],1.0/this->parameters[0]));
	this->right_limit = 1.0 - 1.0/(1.0+pow((maxLength+0.5)/this->parameters[1],1.0/this->parameters[0]));
}


int lavaletteDist::drawLength() {
	MDOUBLE y = getRandDoubleParamVal("uniform", this->left_limit, this->right_limit);
	int x = floor(this->parameters[1]*pow(y/(1-y),this->parameters[0]) + 0.5);
	return x;
}

lavaletteDist::~lavaletteDist() {

}