#pragma once

#ifndef DIFFERENCIALEQUATIONS_H
#define DIFFERENCIALEQUATIONS_H

#include <iostream>
#include <math.h>
#include <tuple>
#include "TMatrix.h"

class DifferencialEquations
{
private:
	double eps, xprev;
	double x, h;
	TVector<double> v, V, prev, preV, p;
	TMatrix<double> A;
	TVector<int> count;
	int Number;
public:
	DifferencialEquations(double X, TVector<double> U, double H, double e, int number);
	tuple<TVector<double>, double> ComplitWithoutControl();
	tuple<TVector<double>, double> ComplitWithControl(bool boolf = true);
	TVector<double> RunKut4(double X, TVector<double>& V_, double H);
	TVector<double> Function(double X, TVector<double>& V_);
	void TakeVecParam(TVector<double> par) { p = par; }
	void TakeMatrParam(TMatrix<double> matr) { A = matr; }
	double max(TVector<double>& Vector);
	TVector<int> GetCounts() { return count; }
	void ChangeX(double X) { x = X; }
	void ChangeH(double H) { h = H; }
	void ChangeCounts(TVector<int> Count) { count = Count; }
	void ChangeV(TVector<double>& V) { v = V; }
	tuple<TVector<double>, double> GetData();
	void BackStep(bool f = true);
	~DifferencialEquations();
};

#endif