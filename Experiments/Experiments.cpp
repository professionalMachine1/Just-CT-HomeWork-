#include <iostream>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <list>
#include "TMatrix.h"
#include "DifferencialEquations.h"

using namespace std;

void SLAE_Solve(TMatrix<double>& Matr, TVector<double>& Vec);
double Determinant(TMatrix<double>& mt);
TMatrix<double> Addit(TMatrix<double>& mt, int row, int cell);
TMatrix<double> MatrixMult(TVector<double>& a, TVector<double>& b);
TVector<double> LinearSystem(TVector<double>& p, bool LinorNonlin);
void NonLinearSystem(TVector<double>& p, TVector<double> theta);
void GetImage(TVector<double>& pL, TVector<double>& pNon);

int main()
{   
	system("chcp 1251");
	system("cls");

	TVector<double> pL, pNon;
	LinearSystem(pL, true);
	NonLinearSystem(pNon, LinearSystem(pNon, true));
	GetImage(pL, pNon);

	return 0;
}

void SLAE_Solve(TMatrix<double>& Matr, TVector<double>& Vec) 
{

	int k = Vec.Size();
	double temp = 0, num, denum;

	ofstream file("dat.txt");
	file << "Изначальная система:\n";
	for (int i1 = 0; i1 < k; i1++)
	{
		for (int j1 = 0; j1 < k; j1++)
			file << round(Matr[i1][j1] * 1000) / 1000 << setw(6);
		file << "| " << Vec[i1] << endl;
	}
	file << "\n\n";
	file << "temp      num      denum\n";

	for (int j = 0; j < k - 1; j++)
	{
		num = Matr[j + 1][j], denum = Matr[j][j];
		for (int i = j + 1; i < k; i++)
		{
			temp = Matr[i][j] / Matr[j][j];
			Matr[i] = Matr[i] - Matr[j] * temp;
			Vec[i] = Vec[i] - Vec[j] * temp;
		}
		temp = num / denum;
		file << temp << setw(10) << num << setw(10) << denum << endl;
		for (int i1 = 0; i1 < k; i1++)
		{
			for (int j1 = 0; j1 < k; j1++)
				file << round(Matr[i1][j1] * 1000) / 1000 << setw(6);
			file << "| " << Vec[i1] << endl;
		}
		file << "\n\n";
	}

	for (int j = k - 1; j > 0; j--)
	{
		num = Matr[j - 1][j], denum = Matr[j][j];
		for (int i = j - 1; i >= 0; i--)
		{
			temp = Matr[i][j] / Matr[j][j];
			Matr[i] = Matr[i] - Matr[j] * temp;
			Vec[i] = Vec[i] - Vec[j] * temp;
		}
		temp = num / denum;
		file << temp << setw(10) << num << setw(10) << denum << endl;
		for (int i1 = 0; i1 < k; i1++)
		{
			for (int j1 = 0; j1 < k; j1++)
				file << round(Matr[i1][j1] * 1000) / 1000 << setw(6);
			file << "| " << Vec[i1] << endl;
		}
		file << "\n\n";
	}

}

double Determinant(TMatrix<double>& mt)
{
	double det = 0;
	int n = mt.Size();

	if (n == 2)
		return mt[0][0] * mt[1][1] - mt[1][0] * mt[0][1];

	TMatrix<double> sup(n - 1);
	for (int i = 0; i < n; i++)
	{
		sup = Addit(mt, 0, i);
		det += pow(-1, i) * Determinant(sup) * mt[0][i];
	}

	return det;
}

TMatrix<double> Addit(TMatrix<double>& mt, int row, int cell)
{
	int n = mt.Size();
	TMatrix<double> res(n - 1);

	for (int i = 0, i1 = 0; i < n; i++)
	{
		if (i != row)
		{
			for (int j = 0, j1 = 0; j < n; j++)
			{
				if (j != cell)
				{
					res[i1][j1] = mt[i][j];
					j1++;
				}
			}
			i1++;
		}
	}

	return res;

}

TMatrix<double> MatrixMult(TVector<double>& a, TVector<double>& b)
{
	int n = a.Size();
	TMatrix<double> Q(n);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			Q[i][j] = a[i] * b[j];

	return Q;
}\

TVector<double> LinearSystem(TVector<double>& p, bool LinOrNonlin)
{
	double m, M, I, l, Bp, Beq, g, Kf, Ks;
	m = 0.127, M = 1.206, I = 1.2 * pow(10, -3), l = 0.1778, Kf = 1.726, Ks = 4.487;
	Beq = 5.4, Bp = 2.4 * pow(10, -3), g = 9.81;

	double temp = (m + M) * (I + m * pow(l, 2)) - pow(m, 2) * pow(l, 2);
	TVector<double> a(5), b(5);
	a[0] = m * g * l * (m + M) / temp;
	a[1] = 0;
	a[2] = -Bp * (m + M) / temp;
	a[3] = -m * l * (Kf * Ks + Beq) / temp;
	a[4] = m * l * Kf / temp;

	b[0] = pow(m, 2) * pow(l, 2) * g / temp;
	b[1] = 0;
	b[2] = -m * l * Bp / temp;
	b[3] = -(I + m * pow(l, 2)) * (Kf * Ks + Beq) / temp;
	b[4] = Kf * (I + m * pow(l, 2)) / temp;

	//Составили матрицы A и B
	TMatrix<double> A(4), C(4);
	TVector<double> B(4);
	A[0][2] = 1, A[1][3] = 1;
	for (int i = 0; i < 4; i++)
		A[2][i] = a[i], A[3][i] = b[i];
	B[2] = a[4], B[3] = b[4];

	//Матрица управляемости
	for (int i = 0; i < 4; i++)
		C[i] = A.POW(i) * B;
	C.Transponir();

	//Определитель матрицы управляемости
	//cout << "Det C = " << Determinant(C) << endl;

	//Находим матрицу замкнутой системы, управление ищем с помощью формулы Аккермана
	TMatrix<double> E(4);
	TVector<double> e_n(4), theta(4);
	E.SingleMatrix();
	e_n[3] = 1;

	E = C.ReverseMatrix() * (A + E * 11.1856) * (A + E * 6.33331) * (A + E) * (A + E * 2);
	E.Transponir();

	theta = E * e_n * (-1);

	A = A + MatrixMult(B, theta);

	//Перевели все коэффициенты в вектор
	if (LinOrNonlin)
	{
		p = TVector<double>(16);
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				p[4 * i + j] = A[i][j];
	}

	return theta;

}

void NonLinearSystem(TVector<double>& p, TVector<double> theta)
{
	double m, M, I, l, Bp, Beq, g, Kf, Ks;
	m = 0.127, M = 1.206, I = 1.2 * pow(10, -3), l = 0.1778, Kf = 1.726, Ks = 4.487;
	Beq = 5.4, Bp = 2.4 * pow(10, -3), g = 9.81;

	p = TVector<double>(16);
	p[0] = m * l * Kf;
	p[1] = pow(m, 2) * pow(l, 2);
	p[2] = m * l * (Kf * Ks + Beq);
	p[3] = m * g * l * (m + M);
	p[4] = Bp * (m + M);
	p[5] = (I + m * pow(l, 2)) * Kf;
	p[6] = m * l * (I + m * pow(l, 2));
	p[7] = (Kf * Ks + Beq) * (I + m * pow(l, 2));
	p[8] = pow(m, 2) * g * pow(l, 2);
	p[9] = m * l * Bp;
	p[10] = (m + M) * (I + m * pow(l, 2));
	p[11] = pow(m, 2) * pow(l, 2);
	for (int i = 0; i < 4; i++)
		p[12 + i] = theta[i];

}

void GetImage(TVector<double>& pL, TVector<double>& pNon)
{
	TVector<double> x_start(4), v(4);
	double t = 0, eps = pow(10, -10), h = pow(10, -3);
	int MaxIt = 10000;
	double tmax = 10;
	//Начальные условия
	x_start[0] = 0.1, x_start[1] = 0.5, x_start[2] = 0.01, x_start[3] = 0.01;
	
	//Заполняем файлы линеаризованной системы
	ofstream fi("Lfi.txt"), x("Lx.txt"), dfi("Ldfi.txt"), dx("Ldx.txt");

	DifferencialEquations solve(t, x_start, h, eps, 3);
	solve.TakeVecParam(pL);

	fi << t << " " << x_start[0] << endl;
	x << t << " " << x_start[1] << endl;
	dfi << t << " " << x_start[2] << endl;
	dx << t << " " << x_start[3] << endl;

	for (int i = 0; (i < MaxIt) && (t < tmax); i++)
	{
		tie(v, t) = solve.ComplitWithControl();
		fi << t << " " << v[0] << endl;
		x << t << " " << v[1] << endl;
		dfi << t << " " << v[2] << endl;
		dx << t << " " << v[3] << endl;
	}

	fi.close(), x.close(), dfi.close(), dx.close();

	//Заполняем файлы нелинейной системы 
	fi.open("Nonfi.txt"), x.open("Nonx.txt"), dfi.open("Nondfi.txt"), dx.open("Nondx.txt");

	t = 0;
	solve = DifferencialEquations(t, x_start, h, eps, 2);
	solve.TakeVecParam(pNon);

	fi << t << " " << x_start[0] << endl;
	x << t << " " << x_start[1] << endl;
	dfi << t << " " << x_start[2] << endl;
	dx << t << " " << x_start[3] << endl;

	for (int i = 0; (i < MaxIt) && (t < tmax); i++)
	{
		tie(v, t) = solve.ComplitWithControl();
		fi << t << " " << v[0] << endl;
		x << t << " " << v[1] << endl;
		dfi << t << " " << v[2] << endl;
		dx << t << " " << v[3] << endl;
	}

	fi.close(), x.close(), dfi.close(), dx.close();

	FILE* gpipe = _popen("C:\\Users\\gnuplot\\bin\\gnuplot.exe -persist", "w");
	string s, way;
	way = "C:\\Users\\HOME\\source\\repos\\Experiments\\Experiments\\";
	s = "set terminal png size 1920,1080 enhanced font 'Times New Roman,12'\n";

	s += "set output '1fi.png'\n";
	s += "plot '" + way + "Lfi.txt' with lines lc 'red' title 'Lfi', ";
	s += "'" + way + "Nonfi.txt' with lines lc 'blue' title 'Nonfi'\n";

	s += "set output '1x.png'\n";
	s += "plot '" + way + "Lx.txt' with lines lc 'red' title 'Lx', ";
	s += "'" + way + "Nonx.txt' with lines lc 'blue' title 'Nonx'\n";

	s += "set output '1dfi.png'\n";
	s += "plot '" + way + "Ldfi.txt' with lines lc 'red' title 'Ldfi', ";
	s += "'" + way + "Nondfi.txt' with lines lc 'blue' title 'Nondfi'\n";

	s += "set output '1dx.png'\n";
	s += "plot '" + way + "Ldx.txt' with lines lc 'red' title 'Ldx', ";
	s += "'" + way + "Nondx.txt' with lines lc 'blue' title 'Nondx'\n";
	fprintf(gpipe, s.c_str());
	_pclose(gpipe);
	
}
