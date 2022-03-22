#pragma once

#ifndef DIFFERENCIALEQUATIONS_H
#define DIFFERENCIALEQUATIONS_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <tuple>
#include <list>
#include "NextGenMatrix.h"

struct ContinuosSystemParametrs
{
	Matrix<double> A, BTheta, LC, ABTh;
	Vector<double> theta, L, B, C;
};

struct SystemState 
{
	Vector<double> x;
	double t = 0;
};

class DifferencialEquations
{
private:
	int Number;
	double eps, tprev, h;
	SystemState state;
	Vector<double> X, xprev, Xprev, p;
	ContinuosSystemParametrs SP;
public:
	DifferencialEquations();
	DifferencialEquations(double t_, Vector<double> x_, double h_, double eps_, int number);

	SystemState& ComplitWithoutControl();
	SystemState& ComplitWithControl(bool boolf = true);

	Vector<double> RunKut4(double t_, Vector<double>& v_, double h_);
	Vector<double> systValue(double t_, Vector<double>& v_);

	void SetVecParam(Vector<double>& par) { p = par; }
	void SetStructParam(ContinuosSystemParametrs& par) { SP = par; }

	void SetCoordinate(double t_) { state.t = t_;; }
	void SetStep(double h_) { h = h_; }
	void SetSystemCondition(Vector<double>& x_) { state.x = x_; }

	SystemState& GetData() { return state; }
	void BackStep(bool f = true);

	void FillInTheFile(double min_t, double max_t, std::list<SystemState>& v, int MaxIt = 100000);
};

#endif