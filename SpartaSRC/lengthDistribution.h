#pragma once

#include <vector>
#include "definitions.h"
#include "SpartaABC_options.h"


class lengthDistribution
{
public:
	lengthDistribution(Vdouble params,int max): parameters(params), maxLength(max) {};
	virtual void generateLengthDistribution() = 0;
	virtual int drawLength() = 0;
	virtual ~lengthDistribution() {} ;
protected:
	int maxLength;
	Vdouble parameters;
};
