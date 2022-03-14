#pragma once

#ifndef DIFFERENCIALEQUATIONS_H
#define DIFFERENCIALEQUATIONS_H

#include <iostream>
#include <math.h>
#include <tuple>
#include "NextGenMatrix.h"

struct SystParametrs
{
	Matrix<double> A, BTheta, LC, ABTh;
	Vector<double> theta;
};

class DifferencialEquations
{
private:
	double eps, tprev;
	double t, h;
	Vector<double> x, X, xprev, Xprev, p;
	SystParametrs SP;
	int Number;
public:
	DifferencialEquations(double t_, Vector<double> x_, double h_, double eps_, int number);
	std::tuple<Vector<double>, double> ComplitWithoutControl();
	std::tuple<Vector<double>, double> ComplitWithControl(bool boolf = true);
	Vector<double> RunKut4(double t_, Vector<double>& v_, double h_);
	Vector<double> Function(double t_, Vector<double>& v_);
	void TakeVecParam(Vector<double>& par) { p = par; }
	void TakeStructParam(SystParametrs& par) { SP = par; }
	void ChangeX(double t_) { t = t_; }
	void ChangeH(double h_) { h = h_; }
	void ChangeV(Vector<double>& x_) { x = x_; }
	std::tuple<Vector<double>, double> GetData();
	void BackStep(bool f = true);
};

#endif