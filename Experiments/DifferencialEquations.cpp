#include "DifferencialEquations.h"

DifferencialEquations::DifferencialEquations(double X, TVector<double> U, double H, double e, int number)
{
	x = X;
	v = U;
	V = TVector<double>(v.Size());
	prev = v;
	preV = V;
	h = H;
	Number = number;
	eps = e;
	count = TVector<int>(2);
	p = TVector<double>(0);
}

TVector<double> DifferencialEquations::RunKut4(double X, TVector<double>& V_, double H)
{
	int k = v.Size();
	TVector<double> k1(k), k2(k), k3(k), k4(k), res(k);
	k1 = Function(X, V_);
	k2 = V_ + k1 * H / 2, k2 = Function(X + H / 2, k2);
	k3 = V_ + k2 * H / 2, k3 = Function(X + H / 2, k3);
	k4 = V_ + k3 * H, k4 = Function(X + H, k4);
	res = V_ + (k1 + k2 * 2 + k3 * 2 + k4) * H / 6;
	return res;
}

TVector<double> DifferencialEquations::Function(double X, TVector<double>& V_)
{
	TVector<double> res(V_.Size());
	switch (Number)
	{
	case 1:
		res[0] = V_[2];
		res[1] = V_[3];
		res[2] = -131.432907 * V_[0] + 14.4427874 * V_[1]
			- 26.3830839 * V_[2] + 25.3923040 * V_[3];
		res[3] = -40.1635050 * V_[0] + 3.33546013 * V_[1]
			- 5.98670166 * V_[2] + 5.86417393 * V_[3];
		break;
	case 2:
		double u, denum;
		u = (p[12] * V_[0] + p[13] * V_[1] + p[14] * V_[2] + p[15] * V_[3]);
		denum = p[10] - p[11] * pow(cos(V_[0]), 2);
		res[0] = V_[2];
		res[1] = V_[3];
		res[2] = (p[0] * cos(V_[0]) * u - p[1] * cos(V_[0]) * sin(V_[0]) * pow(V_[2], 2)
			- p[2] * cos(V_[0]) * V_[3] + p[3] * sin(V_[0]) - p[4] * V_[2]) / denum;
		res[3] = (p[5] * u - p[6] * sin(V_[0]) * pow(V_[2], 2) - p[7] * V_[3]
			+ p[8] * cos(V_[0]) * sin(V_[0]) - p[9] * cos(V_[0]) * V_[2]) / denum;
		break;
	case 3:
		res[0] = V_[2];
		res[1] = V_[3];
		res[2] = p[8] * V_[0] + p[9] * V_[1] + p[10] * V_[2] + p[11] * V_[3];
		res[3] = p[12] * V_[0] + p[13] * V_[1] + p[14] * V_[2] + p[15] * V_[3];
		break;
	default:
		return 0;
	}
	return res;
}

double DifferencialEquations::max(TVector<double>& Vector)
{
	double res = Vector[0];
	for (int i = 1; i < v.Size(); i++)
	{
		if (Vector[i] > res)
			res = Vector[i];
	}
	return res;
}

tuple<TVector<double>, double> DifferencialEquations::ComplitWithoutControl()
{
	v = RunKut4(x, v, h);
	x = x + h;
	return make_tuple(v, x);
}

tuple<TVector<double>, double> DifferencialEquations::ComplitWithControl(bool boolf)
{
	//count[0] - число делений, count[1] - число умножений
	TVector<double> S(v.Size());
	prev = v;
	preV = V;
	// Vn+1/2
	V = RunKut4(x, v, h / 2.0);
	//V - Vn+1 с двумя крышечками (полученное с шагом h/2)
	V = RunKut4(x + h / 2.0, V, h / 2.0);
	//v - Vn+1 вычисленное с шагом h
	v = RunKut4(x, v, h);
	S = (V - v) / 15.0;
	//Берем модуль погрешности
	for (int i = 0; i < v.Size(); i++)
	{
		S[i] = fabs(S[i]);
	}
	//Находим max S
	S[0] = max(S);

	xprev = x;
	x = x + h;
	if ((boolf) && (S[0] < eps / 32.0))
	{
		h = h * 2.0;
		count[1]++;
	}
	else
		if (S[0] > eps)
			BackStep();

	return make_tuple(v, x);
}

tuple<TVector<double>, double> DifferencialEquations::GetData()
{
	return make_tuple(v, x);
}

void DifferencialEquations::BackStep(bool f)
{
	v = prev;
	V = preV;
	x = xprev;
	if (f)
	{
		h = h / 2.0;
		count[0]++;
	}
}

DifferencialEquations::~DifferencialEquations()
{

}