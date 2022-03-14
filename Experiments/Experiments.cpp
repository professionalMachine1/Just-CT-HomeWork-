#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include "NextGenMatrix.h"
#include "DifferencialEquations.h"

using namespace std;

struct TaskParametrs
{
	double m, M, I, l, Bp, Beq, g, Kf, Ks;
};

struct DiscreteStruct
{
	Matrix<double> Ad, ABTd, BTd, LCd;
	Vector<double> Bd, Ld, thetad;
	double h = 0;
};

TaskParametrs P;
Vector<double> EigV(4);
bool REAL1, REAL2;

const string setting = "set terminal pngcairo size 1080,720 enhanced font 'Times New Roman,12'\n";
const string way = "C:\\Users\\HOME\\source\\repos\\Experiments\\Experiments\\";
void PLOT(const string& commands)
{
	FILE* gpipe = _popen("C:\\Users\\gnuplot\\bin\\gnuplot.exe -persist", "w");
	fprintf(gpipe, commands.c_str());
	_pclose(gpipe);
}

void ContinuosSystem();
void DiscreteSystem();
Matrix<double> MatrixExp(Matrix<double> A);
Matrix<double> IntegrateMatrixExp(Matrix<double>& A, double h);
Vector<double>& GetControl(Matrix<double>& A, const Vector<double>& B, Vector<double>& theta);
void InicialSyst(Matrix<double>& A, Vector<double>& B);
Vector<double>& NonLinSystCoef(Vector<double>& p);
void GetImageDiscreteSystem(DiscreteStruct& DS, Vector<double>& x_start);
void GetImageDiscrSystemObs(DiscreteStruct& DS, Vector<double>& x_start, Vector<double>& ksi_start);
void GetImageOriginalSystem(SystParametrs& SP, Vector<double>& pNon, Vector<double>& x_start);
void GetImageLinObserver(SystParametrs& SP, Vector<double>& x_start);
void GetImageNonLinObserver(SystParametrs& SP, Vector<double>& pNon, Vector<double>& x_start);

#define PRINT(name) cout << #name << " = " << name << endl
#define GETNAME(name) cout << #name << endl

int main()
{
	system("chcp 1251");
	system("cls");

	P.m = 0.127, P.M = 1.206, P.I = 1.2 * pow(10, -3), P.l = 0.1778, P.Kf = 1.726, P.Ks = 4.487;
	P.Beq = 5.4, P.Bp = 2.4 * pow(10, -3), P.g = 9.81;

	//Решение для непрерывной системы
	GETNAME(Непрерывная система);
	ContinuosSystem();

	//Решение для дискретной системы
	GETNAME(Дискретная система);
	DiscreteSystem();
	
	return 0;
}

void ContinuosSystem()
{
	//Инициализируем необходимые структуры
	SystParametrs SP;

	int n = 4;
	Matrix<double> A(n), AT(n);
	Vector<double> B(n), C(n), L(n), pNon, x_start(n);
	C[0] = 1, C[1] = 1;

	/*

		Без наблюдателя

	*/

	//Инициализируем A, B
	InicialSyst(A, B);

	//Получаем theta с заданными собственными числами и выводим его на экран
	REAL1 = true, REAL2 = true;
	EigV[0] = -11.1856, EigV[1] = -6.33331, EigV[2] = -1, EigV[3] = -2;
	GetControl(A, B, SP.theta);
	PRINT(SP.theta);

	//Заполняем структуру, которая потом будет применятся в решении ДУ
	SP.A = A, SP.BTheta = CellRowMultiply(B, SP.theta), SP.ABTh = SP.A + SP.BTheta;
	//Заполняем коэффициенты для нелинейной системы
	NonLinSystCoef(pNon);

	//Начальные условия
	x_start[0] = 0.1, x_start[1] = 0.5, x_start[2] = 0.01, x_start[3] = -0.01;

	//Строим графики исходной системы и сохраняем их в формате картинок
	GetImageOriginalSystem(SP, pNon, x_start);

	/*

		Наблюдатель

	*/

	//Ищем L = -theta^T для наблюдателя и выводим его
	REAL1 = true, REAL2 = true;
	EigV[0] = -11.1856, EigV[1] = -6.33331, EigV[2] = -10, EigV[3] = -10;
	AT = TransposedMatrix(A);
	GetControl(AT, C, L);
	L = -1.0 * L;
	PRINT(L);

	//Заполняем ещё один элемент структуры
	SP.LC = CellRowMultiply(L, C);

	x_start = Vector<double>(2 * n);
	//Начальные условия для наблюдателя
	x_start[0] = 0.1, x_start[1] = 0.5, x_start[2] = 0.1, x_start[3] = -0.1;
	//Начальные условия для исходной системы
	x_start[4] = 0.1, x_start[5] = 0.5, x_start[6] = 0.01, x_start[7] = -0.01;

	//Строим графики линейной и нелинейной систем с наблюдателем и сохраняем их в формате картинок
	GetImageLinObserver(SP, x_start);
	GetImageNonLinObserver(SP, pNon, x_start);
}

void DiscreteSystem()
{
	DiscreteStruct DS;
	int n = 4;
	Matrix<double> A(n), AT(n), I = Single<double>(n);
	Vector<double> B(n), C(n), x_start(n), ksi_start(n);
	C[0] = 1, C[1] = 1;
	DS.h = pow(10, -1);

	/*
		
		Без наблюдателя

	*/

	//Инициализируем A, B
	InicialSyst(A, B);
	DS.Ad = MatrixExp(A * DS.h);
	DS.Bd = IntegrateMatrixExp(A, DS.h) * B;

	//Получаем theta с заданными собственными числами
	REAL1 = true, REAL2 = true;
	EigV[0] = exp(-11.1856 * DS.h), EigV[1] = exp(-6.33331 * DS.h);
	EigV[2] = exp(-1 * DS.h), EigV[3] = exp(-2 * DS.h);
	GetControl(DS.Ad, DS.Bd, DS.thetad);
	PRINT(DS.thetad);

	DS.BTd = CellRowMultiply(DS.Bd, DS.thetad);
	DS.ABTd = DS.Ad + DS.BTd;

	//Начальные условия
	x_start[0] = 0.1, x_start[1] = 0.5, x_start[2] = 0.01, x_start[3] = -0.01;
	GetImageDiscreteSystem(DS, x_start);

	/*
	
		Наблюдатель
	
	*/

	//Получаем L = -theta^T с заданными собственными числами
	REAL1 = true, REAL2 = true;
	EigV[0] = exp(-11.1856 * DS.h), EigV[1] = exp(-6.33331 * DS.h);
	EigV[2] = exp(-10 * DS.h), EigV[3] = exp(-10 * DS.h);
	AT = TransposedMatrix(DS.Ad);
	GetControl(AT, C, DS.Ld);
	DS.Ld = -1.0 * DS.Ld;
	PRINT(DS.Ld);

	DS.LCd = CellRowMultiply(DS.Ld, C);

	//Начальные условия
	x_start[0] = 0.1, x_start[1] = 0.5, x_start[2] = 0.01, x_start[3] = -0.01;
	ksi_start[0] = 0.1, ksi_start[1] = 0.5, ksi_start[2] = 0.1, ksi_start[3] = -0.1;
	GetImageDiscrSystemObs(DS, x_start, ksi_start);

}

Vector<double>& GetControl(Matrix<double>& A, const Vector<double>& B, Vector<double>& theta)
{
	int n = A.size();
	//Матрица управляемости
	Matrix<double> ControlMatrix(n);
	//Единичная матрица и вектор e_n = (0, 0, ..., 1)^T
	Matrix<double> I = Single<double>(n);
	Vector<double> en(n);
	double temp;
	en[3] = 1;

	//Заполняем матрицу управляемости
	/*ControlMatrix.SetCell(0, B);
	for (int i = 1; i < n; i++)
		ControlMatrix.SetCell(i, A * ControlMatrix.GetCell(i - 1));*/
	ControlMatrix[0] = B;
	for (int i = 1; i < n; i++)
		ControlMatrix[i] = A * ControlMatrix[i - 1];
	ControlMatrix = TransposedMatrix(ControlMatrix);

	//Определитель матрицы управляемости
	PRINT(Determinant(ControlMatrix));
	//PRINT(InverseMatrix(ControlMatrix));

	theta = -1.0 * en * InverseMatrix(ControlMatrix);
	//Находим матрицу замкнутой системы, управление ищем с помощью формулы Аккермана
	if (REAL1 && REAL2)
	{
		theta = theta * (A - EigV[0] * I) * (A - EigV[1] * I)
			* (A - EigV[2] * I) * (A - EigV[3] * I);
	}
	if (REAL1 && !REAL2)
	{
		theta = theta * (A - EigV[0] * I) * (A - EigV[1] * I)
			* (A * A - 2 * EigV[2] * A + (pow(EigV[2], 2) + pow(EigV[3], 2)) * I);
	}
	if (!REAL1 && REAL2)
	{
		theta = theta * (A - EigV[2] * I) * (A - EigV[3] * I)
			* (A * A - 2 * EigV[0] * A + (pow(EigV[0], 2) + pow(EigV[1], 2)) * I);
	}
	if (!REAL1 && !REAL2)
	{
		theta = theta * (A * A - 2 * EigV[2] * A + (pow(EigV[2], 2) + pow(EigV[3], 2)) * I)
			* (A * A - 2 * EigV[0] * A + (pow(EigV[0], 2) + pow(EigV[1], 2)) * I);
	}

	return theta;
}

Matrix<double> MatrixExp(Matrix<double> A)
{
	Matrix<double> GapMt = Single<double>(A.size()), res(A.size());
	double GapCoef = 1;
	for (int i = 0; i < 14; i++)
	{
		res += GapMt / GapCoef;
		GapMt *= A;
		GapCoef *= (i + 1.0);
	}
	return res;
}

Matrix<double> IntegrateMatrixExp(Matrix<double>& A, double h)
{
	Matrix<double> GapMt = Single<double>(A.size()), res(A.size());
	double GapCoef = h;
	for (int i = 1; i < 14; i++)
	{
		res += GapMt * GapCoef;
		GapMt *= A;
		GapCoef *= h / (i + 1.0);
	}
	return res;
}

void InicialSyst(Matrix<double>& A, Vector<double>& B)
{
	//Заполнили матрицу A и матрицу-столбец B
	double temp = (P.m + P.M) * (P.I + P.m * pow(P.l, 2)) - pow(P.m, 2) * pow(P.l, 2);
	A[2][0] = P.m * P.g * P.l * (P.m + P.M) / temp;
	A[2][1] = 0;
	A[2][2] = -P.Bp * (P.m + P.M) / temp;
	A[2][3] = -P.m * P.l * (P.Kf * P.Ks + P.Beq) / temp;
	B[2] = P.m * P.l * P.Kf / temp;

	A[3][0] = pow(P.m, 2) * pow(P.l, 2) * P.g / temp;
	A[3][1] = 0;
	A[3][2] = -P.m * P.l * P.Bp / temp;
	A[3][3] = -(P.I + P.m * pow(P.l, 2)) * (P.Kf * P.Ks + P.Beq) / temp;
	B[3] = P.Kf * (P.I + P.m * pow(P.l, 2)) / temp;

	A[0][2] = 1, A[1][3] = 1;

}

Vector<double>& NonLinSystCoef(Vector<double>& p)
{
	p = Vector<double>(12);
	p[0] = P.m * P.l * P.Kf;
	p[1] = pow(P.m, 2) * pow(P.l, 2);
	p[2] = P.m * P.l * (P.Kf * P.Ks + P.Beq);
	p[3] = P.m * P.g * P.l * (P.m + P.M);
	p[4] = P.Bp * (P.m + P.M);
	p[5] = (P.I + P.m * pow(P.l, 2)) * P.Kf;
	p[6] = P.m * P.l * (P.I + P.m * pow(P.l, 2));
	p[7] = (P.Kf * P.Ks + P.Beq) * (P.I + P.m * pow(P.l, 2));
	p[8] = pow(P.m, 2) * P.g * pow(P.l, 2);
	p[9] = P.m * P.l * P.Bp;
	p[10] = (P.m + P.M) * (P.I + P.m * pow(P.l, 2));
	p[11] = pow(P.m, 2) * pow(P.l, 2);
	return p;
}

void GetImageDiscreteSystem(DiscreteStruct& DS, Vector<double>& x_start)
{
	Vector<double> x(4);
	double t = 0, tmax = 10;
	int MaxIt = 20000;

	//Заполняем файлы линеаризованной системы
	ofstream fi("Lfi.txt"), xg("Lx.txt"), dfi("Ldfi.txt"), dx("Ldx.txt"), Lu("Lu.txt");

	x = x_start;

	for (int i = 0; (i < MaxIt) && (t < tmax); i++)
	{
		fi << t << " " << x[0] << endl;
		xg << t << " " << x[1] << endl;
		dfi << t << " " << x[2] << endl;
		dx << t << " " << x[3] << endl;
		Lu << t << " " << x * DS.thetad << endl;
		x = DS.ABTd * x, t += DS.h;
	}

	fi.close(), xg.close(), dfi.close(), dx.close(), Lu.close();

	string s;
	s = setting;
	s += "set grid\n";

	s += "set output '1DRfi.png'\n";
	s += "plot '" + way + "Lfi.txt' with lines lc 'red' lw 2 title 'Lfi\n";

	s += "set output '1DRx.png'\n";
	s += "plot '" + way + "Lx.txt' with lines lc 'red' lw 2 title 'Lx'\n";

	s += "set output '1DRdfi.png'\n";
	s += "plot '" + way + "Ldfi.txt' with lines lc 'red' lw 2 title 'Ldfi'\n";

	s += "set output '1DRdx.png'\n";
	s += "plot '" + way + "Ldx.txt' with lines lc 'red' lw 2 title 'Ldx'\n";

	s += "set output '1DRu.png'\n";
	s += "plot '" + way + "Lu.txt' with lines lc 'green' lw 2 title 'Lu'\n";
	PLOT(s);
}

void GetImageDiscrSystemObs(DiscreteStruct& DS, Vector<double>& x_start, Vector<double>& ksi_start)
{
	Vector<double> x(4), ksi(4), ksiprev(4), xprev(4);
	double t = 0, tmax = 10;
	int MaxIt = 20000;

	//Заполняем файлы линеаризованной системы
	ofstream fi("Lfi.txt"), xg("Lx.txt"), dfi("Ldfi.txt"), dx("Ldx.txt"), Lu("Lu.txt");
	ofstream obsfi("LObsfi.txt"), obsxg("LObsx.txt"), obsdfi("LObsdfi.txt"), obsdx("LObsdx.txt");

	x = x_start, ksi = ksi_start;

	for (int i = 0; (i < MaxIt) && (t < tmax); i++)
	{
		xprev = x, ksiprev = ksi;
		fi << t << " " << x[0] << endl;
		xg << t << " " << x[1] << endl;
		dfi << t << " " << x[2] << endl;
		dx << t << " " << x[3] << endl;
		obsfi << t << " " << ksi[0] << endl;
		obsxg << t << " " << ksi[1] << endl;
		obsdfi << t << " " << ksi[2] << endl;
		obsdx << t << " " << ksi[3] << endl;
		Lu << t << " " << ksi * DS.thetad << endl;
		ksi = DS.ABTd * ksiprev + DS.LCd * (xprev - ksiprev);
		x = DS.Ad * xprev + DS.BTd * ksiprev;
		t += DS.h;
	}

	fi.close(), xg.close(), dfi.close(), dx.close(), Lu.close();
	obsfi.close(), obsxg.close(), obsdfi.close(), obsdx.close();

	string s;
	s = setting;
	s += "set grid\n";

	s += "set output '1DOBSfi.png'\n";
	s += "plot '" + way + "Lfi.txt' with lines lc 'red' lw 2 title 'Lfi', ";
	s += "'" + way + "LObsfi.txt' with lines lc 'blue' lw 2 title 'LObsfi\n";

	s += "set output '1DOBSx.png'\n";
	s += "plot '" + way + "Lx.txt' with lines lc 'red' lw 2 title 'Lx', ";
	s += "'" + way + "LObsx.txt' with lines lc 'blue' lw 2 title 'LObsx\n";

	s += "set output '1DOBSdfi.png'\n";
	s += "plot '" + way + "Ldfi.txt' with lines lc 'red' lw 2 title 'Ldfi', ";
	s += "'" + way + "LObsdfi.txt' with lines lc 'blue' lw 2 title 'LObsdfi\n";

	s += "set output '1DOBSdx.png'\n";
	s += "plot '" + way + "Ldx.txt' with lines lc 'red' lw 2 title 'Ldx', ";
	s += "'" + way + "LObsdx.txt' with lines lc 'blue' lw 2 title 'LObsdx\n";

	s += "set output '1DOBSu.png'\n";
	s += "plot '" + way + "Lu.txt' with lines lc 'green' lw 2 title 'Lu'\n";
	PLOT(s);
}

void GetImageOriginalSystem(SystParametrs& SP, Vector<double>& pNon, Vector<double>& x_start)
{
	Vector<double> x(4);
	double t = 0, eps = pow(10, -10), h = pow(10, -3), tmax = 10;
	int MaxIt = 20000;

	//Заполняем файлы линеаризованной системы
	ofstream fi("Lfi.txt"), xg("Lx.txt"), dfi("Ldfi.txt"), dx("Ldx.txt"), Lu("Lu.txt");
	ofstream Nfi("NLfi.txt"), Nxg("NLx.txt"), Ndfi("NLdfi.txt"), Ndx("NLdx.txt"), Nu("NLu.txt");

	x = x_start;
	DifferencialEquations solve(t, x_start, h, eps, 1);
	solve.TakeStructParam(SP);

	for (int i = 0; (i < MaxIt) && (t < tmax); i++)
	{
		fi << t << " " << x[0] << endl;
		xg << t << " " << x[1] << endl;
		dfi << t << " " << x[2] << endl;
		dx << t << " " << x[3] << endl;
		Lu << t << " " << x * SP.theta << endl;
		tie(x, t) = solve.ComplitWithControl();
	}

	fi.close(), xg.close(), dfi.close(), dx.close(), Lu.close();

	//Заполняем файлы нелинейной системы 
	t = 0;
	x = x_start;
	solve = DifferencialEquations(t, x_start, h, eps, 2);
	solve.TakeStructParam(SP);
	solve.TakeVecParam(pNon);

	for (int i = 0; (i < MaxIt) && (t < tmax); i++)
	{
		Nfi << t << " " << x[0] << endl;
		Nxg << t << " " << x[1] << endl;
		Ndfi << t << " " << x[2] << endl;
		Ndx << t << " " << x[3] << endl;
		Nu << t << " " << x * SP.theta << endl;
		tie(x, t) = solve.ComplitWithControl();
	}

	Nfi.close(), Nxg.close(), Ndfi.close(), Ndx.close(), Nu.close();

	string s;
	s = setting;
	s += "set grid\n";

	s += "set output '1Rfi.png'\n";
	s += "plot '" + way + "Lfi.txt' with lines lc 'red' lw 2 title 'Lfi', ";
	s += "'" + way + "NLfi.txt' with lines lc 'blue' lw 2 title 'NLfi'\n";

	s += "set output '1Rx.png'\n";
	s += "plot '" + way + "Lx.txt' with lines lc 'red' lw 2 title 'Lx', ";
	s += "'" + way + "NLx.txt' with lines lc 'blue' lw 2 title 'NLx'\n";

	s += "set output '1Rdfi.png'\n";
	s += "plot '" + way + "Ldfi.txt' with lines lc 'red' lw 2 title 'Ldfi', ";
	s += "'" + way + "NLdfi.txt' with lines lc 'blue' lw 2 title 'NLdfi'\n";

	s += "set output '1Rdx.png'\n";
	s += "plot '" + way + "Ldx.txt' with lines lc 'red' lw 2 title 'Ldx', ";
	s += "'" + way + "NLdx.txt' with lines lc 'blue' lw 2 title 'NLdx'\n";

	s += "set output '1Ru.png'\n";
	s += "plot '" + way + "Lu.txt' with lines lc 'red' lw 2 title 'Lu', ";
	s += "'" + way + "NLu.txt' with lines lc 'blue' lw 2 title 'NLu'\n";
	PLOT(s);

}

void GetImageLinObserver(SystParametrs& SP, Vector<double>& x_start)
{
	Vector<double> x(8);
	double t = 0, eps = pow(10, -10), h = pow(10, -3), tmax = 10, uVal;
	int MaxIt = 20000;

	//Заполняем файлы линеаризованной системы
	ofstream fi("Lfi.txt"), xg("Lx.txt"), dfi("Ldfi.txt"), dx("Ldx.txt"), u("Lu.txt");
	ofstream obsfi("LObsfi.txt"), obsxg("LObsx.txt"), obsdfi("LObsdfi.txt"), obsdx("LObsdx.txt");

	x = x_start;
	DifferencialEquations solve(t, x_start, h, eps, 3);
	solve.TakeStructParam(SP);

	for (int i = 0; (i < MaxIt) && (t < tmax); i++)
	{
		obsfi << t << " " << x[0] << endl;
		obsxg << t << " " << x[1] << endl;
		obsdfi << t << " " << x[2] << endl;
		obsdx << t << " " << x[3] << endl;
		fi << t << " " << x[4] << endl;
		xg << t << " " << x[5] << endl;
		dfi << t << " " << x[6] << endl;
		dx << t << " " << x[7] << endl;
		uVal = 0;
		for (int i = 0; i < 4; i++)
			uVal += SP.theta[i] * x[i];
		u << t << " " << uVal << endl;
		tie(x, t) = solve.ComplitWithControl();
	}

	fi.close(), xg.close(), dfi.close(), dx.close(), u.close();
	obsfi.close(), obsxg.close(), obsdfi.close(), obsdx.close();
	
	string s;
	s = setting;
	s += "set grid\n";

	s += "set output '1LOBSfi.png'\n";
	s += "plot '" + way + "Lfi.txt' with lines lc 'red' lw 2 title 'Lfi', ";
	s += "'" + way + "LObsfi.txt' with lines lc 'blue' lw 2 title 'LObsfi'\n";

	s += "set output '1LOBSx.png'\n";
	s += "plot '" + way + "Lx.txt' with lines lc 'red' lw 2 title 'Lx', ";
	s += "'" + way + "LObsx.txt' with lines lc 'blue' lw 2 title 'LObsx'\n";

	s += "set output '1LOBSdfi.png'\n";
	s += "plot '" + way + "Ldfi.txt' with lines lc 'red' lw 2 title 'Ldfi', ";
	s += "'" + way + "LObsdfi.txt' with lines lc 'blue' lw 2 title 'LObsdfi'\n";

	s += "set output '1LOBSdx.png'\n";
	s += "plot '" + way + "Ldx.txt' with lines lc 'red' lw 2 title 'Ldx', ";
	s += "'" + way + "LObsdx.txt' with lines lc 'blue' lw 2 title 'LObsdx'\n";

	s += "set output '1LOBSu.png'\n";
	s += "plot '" + way + "Lu.txt' with lines lc 'green' lw 2 title 'Lu'\n";
	PLOT(s);

}

void GetImageNonLinObserver(SystParametrs& SP, Vector<double>& pNon, Vector<double>& x_start)
{
	Vector<double> x(8);
	double t = 0, eps = pow(10, -10), h = pow(10, -3), tmax = 10, uVal;
	int MaxIt = 20000;

	//Заполняем файлы линеаризованной системы
	ofstream fi("NLfi.txt"), xg("NLx.txt"), dfi("NLdfi.txt"), dx("NLdx.txt"), u("NLu.txt");
	ofstream obsfi("NLObsfi.txt"), obsxg("NLObsx.txt"), obsdfi("NLObsdfi.txt"), obsdx("NLObsdx.txt");

	x = x_start;
	DifferencialEquations solve(t, x_start, h, eps, 4);
	solve.TakeStructParam(SP);
	solve.TakeVecParam(pNon);

	for (int i = 0; (i < MaxIt) && (t < tmax); i++)
	{
		obsfi << t << " " << x[0] << endl;
		obsxg << t << " " << x[1] << endl;
		obsdfi << t << " " << x[2] << endl;
		obsdx << t << " " << x[3] << endl;
		fi << t << " " << x[4] << endl;
		xg << t << " " << x[5] << endl;
		dfi << t << " " << x[6] << endl;
		dx << t << " " << x[7] << endl;
		uVal = 0;
		for (int i = 0; i < 4; i++)
			uVal += SP.theta[i] * x[i];
		u << t << " " << uVal << endl;
		tie(x, t) = solve.ComplitWithControl();
	}

	fi.close(), xg.close(), dfi.close(), dx.close(), u.close();
	obsfi.close(), obsxg.close(), obsdfi.close(), obsdx.close();
	
	string s;
	s = setting;
	s += "set grid\n";

	s += "set output '1NLOBSfi.png'\n";
	s += "plot '" + way + "NLfi.txt' with lines lc 'red' lw 2 title 'NLfi', ";
	s += "'" + way + "NLObsfi.txt' with lines lc 'blue' lw 2 title 'NLObsfi'\n";

	s += "set output '1NLOBSx.png'\n";
	s += "plot '" + way + "NLx.txt' with lines lc 'red' lw 2 title 'NLx', ";
	s += "'" + way + "NLObsx.txt' with lines lc 'blue' lw 2 title 'NLObsx'\n";

	s += "set output '1NLOBSdfi.png'\n";
	s += "plot '" + way + "NLdfi.txt' with lines lc 'red' lw 2 title 'NLdfi', ";
	s += "'" + way + "NLObsdfi.txt' with lines lc 'blue' lw 2 title 'NLObsdfi'\n";

	s += "set output '1NLOBSdx.png'\n";
	s += "plot '" + way + "NLdx.txt' with lines lc 'red' lw 2 title 'NLdx', ";
	s += "'" + way + "NLObsdx.txt' with lines lc 'blue' lw 2 title 'NLObsdx'\n";

	s += "set output '1NLOBSu.png'\n";
	s += "plot '" + way + "NLu.txt' with lines lc 'green' lw 2 title 'NLu'\n";
	PLOT(s);

}