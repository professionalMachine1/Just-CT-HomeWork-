#pragma once

#ifndef _OPWITHFUNCTIONS_H_
#define _OPWITHFUNCTIONS_H_

#include <iostream>
#include <math.h>

using namespace std;

#define fVAL1(s) s
#define fVAL2 (1+x+sin(x))

class OpWithFucntions
{
private:
	string s;
public:
	OpWithFucntions() {}
	double fValue(double x);
	~OpWithFucntions() {}
};

#endif