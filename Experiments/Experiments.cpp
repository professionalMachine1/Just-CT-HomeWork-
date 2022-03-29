#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include "NextGenMatrix.h"
#include "DifferencialEquations.h"

using namespace std;

struct TaskParametrs
{
	double m, M, I, l, Bp, Beq, g, Kf, Ks, tmax;
};

struct DiscreteStruct
{
	Matrix<double> Ad, ABTdc, BTdc, ABTdp, BTdp, LCd;
	Vector<double> Bdc, Bdp, Ld, thetadc, thetadp;
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
void GetImageOriginalSystem(ContinuosSystemParametrs& SP, Vector<double>& pNon, Vector<double>& x_start);
void GetImageLinObserver(ContinuosSystemParametrs& SP, Vector<double>& x_start);
void GetImageNonLinObserver(ContinuosSystemParametrs& SP, Vector<double>& pNon, Vector<double>& x_start);

#define PRINT(name) cout << #name << " = " << name << endl
#define GETNAME(name) cout << #name << endl

int main()
{
	system("chcp 1251");
	system("cls");

	P.tmax = 16;
	P.m = 0.127, P.M = 1.206, P.I = 1.2 * pow(10, -3), P.l = 0.1778, P.Kf = 1.726, P.Ks = 4.487;
	P.Beq = 5.4, P.Bp = 2.4 * pow(10, -3), P.g = 9.81;

	//Решение для непрерывной системы
	//GETNAME(Непрерывная система);
	//ContinuosSystem();

	//Решение для дискретной системы
	//GETNAME(Дискретная система);
	//DiscreteSystem();
	
	return 0;
}

void ContinuosSystem()
{
	//Инициализируем необходимые структуры
	ContinuosSystemParametrs CSP;

	int n = 4;

	Matrix<double> AT(n);
	CSP.A = Matrix<double>(n);
	CSP.B = Vector<double>(n), CSP.C = Vector<double>(n), CSP.L = Vector<double>(n);
	Vector<double> pNon, x_start(n);
	CSP.C[0] = 1, CSP.C[1] = 1;

	/*

		Без наблюдателя

	*/

	//Инициализируем A, B
	InicialSyst(CSP.A, CSP.B);

	//Получаем theta с заданными собственными числами и выводим его на экран
	REAL1 = true, REAL2 = true;
	EigV[0] = -11.1856, EigV[1] = -6.33331, EigV[2] = -1, EigV[3] = -1;
	GetControl(CSP.A, CSP.B, CSP.theta);
	PRINT(CSP.theta);

	//Заполняем структуру, которая потом будет применятся в решении ДУ
	CSP.BTheta = CellRowMultiply(CSP.B, CSP.theta), CSP.ABTh = CSP.A + CSP.BTheta;
	//Заполняем коэффициенты для нелинейной системы
	NonLinSystCoef(pNon);

	//Начальные условия
	x_start[0] = 1, x_start[1] = 0.5, x_start[2] = 0.3, x_start[3] = 0.5;

	//Строим графики исходной системы и сохраняем их в формате картинок
	GetImageOriginalSystem(CSP, pNon, x_start);

	/*

		Наблюдатель

	*/

	//Ищем L = -theta^T для наблюдателя и выводим его
	REAL1 = true, REAL2 = true;
	EigV[0] = -11.1856, EigV[1] = -6.33331, EigV[2] = -10, EigV[3] = -10;
	AT = TransposedMatrix(CSP.A);
	GetControl(AT, CSP.C, CSP.L);
	CSP.L = -1.0 * CSP.L;
	PRINT(CSP.L);

	//Заполняем ещё один элемент структуры
	CSP.LC = CellRowMultiply(CSP.L, CSP.C);

	x_start = Vector<double>(2 * n);
	//Начальные условия для наблюдателя
	x_start[0] = 1, x_start[1] = 0.5, x_start[2] = 0.3, x_start[3] = 0.5;
	//Начальные условия для исходной системы
	x_start[4] = 1, x_start[5] = 0.5, x_start[6] = 0.1, x_start[7] = 0.2;

	//Строим графики линейной и нелинейной систем с наблюдателем и сохраняем их в формате картинок
	GetImageLinObserver(CSP, x_start);
	GetImageNonLinObserver(CSP, pNon, x_start);
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

	//Инициализируем Ad, Bd
	InicialSyst(A, B);
	DS.Ad = MatrixExp(A * DS.h);

	//Считаем Bd в зависимости от того, каким является управление - импульсным или кусочно-постоянным
	DS.Bdp = DS.Ad * B;
	DS.Bdc = IntegrateMatrixExp(A, DS.h) * B;

	//Получаем theta с заданными собственными числами
	REAL1 = true, REAL2 = true;
	EigV[0] = exp(-11.1856 * DS.h), EigV[1] = exp(-6.33331 * DS.h);
	EigV[2] = exp(-1 * DS.h), EigV[3] = exp(-2 * DS.h);
	GetControl(DS.Ad, DS.Bdc, DS.thetadc);
	PRINT(DS.thetadc);
	GetControl(DS.Ad, DS.Bdp, DS.thetadp);
	PRINT(DS.thetadp);

	DS.BTdc = CellRowMultiply(DS.Bdc, DS.thetadc);
	DS.ABTdc = DS.Ad + DS.BTdc;

	DS.BTdp = CellRowMultiply(DS.Bdp, DS.thetadp);
	DS.ABTdp = DS.Ad + DS.BTdp;

	//Начальные условия
	x_start[0] = 1, x_start[1] = 0.5, x_start[2] = 0.4, x_start[3] = -0.2;
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
	ksi_start[0] = 1, ksi_start[1] = 0.5, ksi_start[2] = 0.4, ksi_start[3] = -0.2;
	x_start[0] = 1, x_start[1] = 0.5, x_start[2] = -0.4, x_start[3] = 0.5;
	
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
	Vector<double> x(x_start.size());
	int MaxIt = 20000;
	double t;

	//Заполняем файлы линеаризованной системы
	ofstream cfi("cfi.txt"), cxg("cx.txt"), cdfi("cdfi.txt"), cdx("cdx.txt"), cu("cu.txt");
	ofstream pfi("pfi.txt"), pxg("px.txt"), pdfi("pdfi.txt"), pdx("pdx.txt"), pu("pu.txt");

	/*
	
		С кусочно-постоянным управлением
	
	*/

	t = 0, x = x_start;

	for (int i = 0; (i < MaxIt) && (t < P.tmax); i++)
	{
		cfi << t << " " << x[0] << endl;
		cxg << t << " " << x[1] << endl;
		cdfi << t << " " << x[2] << endl;
		cdx << t << " " << x[3] << endl;
		cu << t << " " << x * DS.thetadc << endl;
		x = DS.ABTdc * x, t += DS.h;
	}

	cfi.close(), cxg.close(), cdfi.close(), cdx.close(), cu.close();

	/*
	
		С импульсным управлением
	
	*/

	t = 0, x = x_start;

	for (int i = 0; (i < MaxIt) && (t < P.tmax); i++)
	{
		pfi << t << " " << x[0] << endl;
		pxg << t << " " << x[1] << endl;
		pdfi << t << " " << x[2] << endl;
		pdx << t << " " << x[3] << endl;
		pu << t << " " << x * DS.thetadp << endl;
		x = DS.ABTdp * x, t += DS.h;
	}

	pfi.close(), pxg.close(), pdfi.close(), pdx.close(), pu.close();


	string s;
	s = setting;
	s += "set grid\n";

	s += "set output '1D_R_fi.png'\n";
	s += "plot '" + way + "cfi.txt' with lines lc 'red' lw 2 title 'CONSTfi', ";
	s += "'" + way + "pfi.txt' with lines lc 'blue' lw 2 title 'IMPULSEfi'\n";

	s += "set output '1D_R_x.png'\n";
	s += "plot '" + way + "cx.txt' with lines lc 'red' lw 2 title 'CONSTx', ";
	s += "'" + way + "px.txt' with lines lc 'blue' lw 2 title 'IMPUSLEx'\n";

	s += "set output '1D_R_dfi.png'\n";
	s += "plot '" + way + "cdfi.txt' with lines lc 'red' lw 2 title 'CONSTdfi', ";
	s += "'" + way + "pdfi.txt' with lines lc 'blue' lw 2 title 'IMPUSLEdfi'\n";

	s += "set output '1D_R_dx.png'\n";
	s += "plot '" + way + "cdx.txt' with lines lc 'red' lw 2 title 'CONSTdx', ";
	s += "'" + way + "pdx.txt' with lines lc 'blue' lw 2 title 'IMPUSLEdx'\n";

	s += "set output '1D_R_u.png'\n";
	s += "plot '" + way + "cu.txt' with lines lc 'green' lw 2 title 'CONSTu', ";
	s += "'" + way + "pu.txt' with lines lc 'blue' lw 2 title 'IMPUSLEu'\n";
	PLOT(s);
}

void GetImageDiscrSystemObs(DiscreteStruct& DS, Vector<double>& x_start, Vector<double>& ksi_start)
{
	int n = 4;
	Vector<double> x(n), ksi(n), ksiprev(n), xprev(n);
	int MaxIt = 20000;
	double t;

	//Заполняем файлы линеаризованной системы
	Vector<ofstream> xt(n), ksit(n);
	ofstream u("u.txt", ios_base::out);
	string name;

	/*
	
		Для кусочно-постоянного управления
	
	*/
	
	for (int i = 0; i < n; i++)
	{
		name = to_string(i) + "x.txt";
		xt[i].open(name, ios_base::out);
		name = to_string(i) + "ksi.txt";
		ksit[i].open(name, ios_base::out);
	}

	t = 0, x = x_start, ksi = ksi_start;

	for (int i = 0; (i < MaxIt) && (t < P.tmax); i++)
	{
		xprev = x, ksiprev = ksi;
		for (int j = 0; j < n; j++)
			xt[j] << t << " " << x[j] << endl;
		for (int j = 0; j < n; j++)
			ksit[j] << t << " " << ksi[j] << endl;
		u << t << " " << ksi * DS.thetadc << endl;
		ksi = DS.ABTdc * ksiprev + DS.LCd * (xprev - ksiprev);
		x = DS.Ad * xprev + DS.BTdc * ksiprev;
		t += DS.h;
	}

	string s;
	s = setting;
	s += "set grid\n";

	s += "set output '1D_OBS_CONST_fi.png'\n";
	s += "plot '" + way + "0x.txt' with lines lc 'red' lw 2 title 'CONSTfi', ";
	s += "'" + way + "0ksi.txt' with lines lc 'blue' lw 2 title 'CONSTObsfi\n";

	s += "set output '1D_OBS_CONST_x.png'\n";
	s += "plot '" + way + "1x.txt' with lines lc 'red' lw 2 title 'CONSTx', ";
	s += "'" + way + "1ksi.txt' with lines lc 'blue' lw 2 title 'CONSTObsx\n";

	s += "set output '1D_OBS_CONST_dfi.png'\n";
	s += "plot '" + way + "2x.txt' with lines lc 'red' lw 2 title 'CONSTdfi', ";
	s += "'" + way + "2ksi.txt' with lines lc 'blue' lw 2 title 'CONSTObsdfi\n";

	s += "set output '1D_OBS_CONST_dx.png'\n";
	s += "plot '" + way + "3x.txt' with lines lc 'red' lw 2 title 'CONSTdx', ";
	s += "'" + way + "3ksi.txt' with lines lc 'blue' lw 2 title 'CONSTObsdx\n";

	s += "set output '1D_OBS_CONST_u.png'\n";
	s += "plot '" + way + "u.txt' with lines lc 'green' lw 2 title 'CONSTu'\n";
	PLOT(s);

	/*
	
		Для импульсного управления
	
	*/

	for (int i = 0; i < n; i++)
		xt[i].close(), ksit[i].close();
	u.close();

	for (int i = 0; i < n; i++)
	{
		name = to_string(i) + "x.txt";
		xt[i].open(name, ios_base::out);
		name = to_string(i) + "ksi.txt";
		ksit[i].open(name, ios_base::out);
	}
	u.open("u.txt", ios_base::out);

	t = 0, x = x_start, ksi = ksi_start;

	for (int i = 0; (i < MaxIt) && (t < P.tmax); i++)
	{
		xprev = x, ksiprev = ksi;
		for (int j = 0; j < n; j++)
			xt[j] << t << " " << x[j] << endl;
		for (int j = 0; j < n; j++)
			ksit[j] << t << " " << ksi[j] << endl;
		u << t << " " << ksi * DS.thetadp << endl;
		ksi = DS.ABTdp * ksiprev + DS.LCd * (xprev - ksiprev);
		x = DS.Ad * xprev + DS.BTdp * ksiprev;
		t += DS.h;
	}
	
	s = setting;
	s += "set grid\n";

	s += "set output '1D_OBS_IMPUSLE_fi.png'\n";
	s += "plot '" + way + "0x.txt' with lines lc 'red' lw 2 title 'IMPUSLEfi', ";
	s += "'" + way + "0ksi.txt' with lines lc 'blue' lw 2 title 'IMPUSLEObsfi\n";

	s += "set output '1D_OBS_IMPUSLE_x.png'\n";
	s += "plot '" + way + "1x.txt' with lines lc 'red' lw 2 title 'IMPUSLEx', ";
	s += "'" + way + "1ksi.txt' with lines lc 'blue' lw 2 title 'IMPUSLEObsx\n";

	s += "set output '1D_OBS_IMPUSLE_dfi.png'\n";
	s += "plot '" + way + "2x.txt' with lines lc 'red' lw 2 title 'IMPUSLEdfi', ";
	s += "'" + way + "2ksi.txt' with lines lc 'blue' lw 2 title 'IMPUSLEObsdfi\n";

	s += "set output '1D_OBS_IMPUSLE_dx.png'\n";
	s += "plot '" + way + "3x.txt' with lines lc 'red' lw 2 title 'IMPUSLEdx', ";
	s += "'" + way + "3ksi.txt' with lines lc 'blue' lw 2 title 'IMPUSLEObsdx\n";

	s += "set output '1D_OBS_IMPUSLE_u.png'\n";
	s += "plot '" + way + "u.txt' with lines lc 'green' lw 2 title 'IMPUSLEu'\n";
	PLOT(s);

}

void GetImageOriginalSystem(ContinuosSystemParametrs& SP, Vector<double>& pNon, Vector<double>& x_start)
{
	SystemState state;
	double eps = pow(10, -10), h = pow(10, -3);
	int MaxIt = 20000;

	//Заполняем файлы линеаризованной системы
	ofstream fi("Lfi.txt"), xg("Lx.txt"), dfi("Ldfi.txt"), dx("Ldx.txt"), Lu("Lu.txt");
	ofstream Nfi("NLfi.txt"), Nxg("NLx.txt"), Ndfi("NLdfi.txt"), Ndx("NLdx.txt"), Nu("NLu.txt");

	state.t = 0, state.x = x_start;
	DifferencialEquations solve(state.t, x_start, h, eps, 1);
	solve.SetStructParam(SP);

	for (int i = 0; (i < MaxIt) && (state.t < P.tmax); i++)
	{
		fi << state.t << " " << state.x[0] << endl;
		xg << state.t << " " << state.x[1] << endl;
		dfi << state.t << " " << state.x[2] << endl;
		dx << state.t << " " << state.x[3] << endl;
		Lu << state.t << " " << state.x * SP.theta << endl;
		state = solve.ComplitWithControl();
	}

	fi.close(), xg.close(), dfi.close(), dx.close(), Lu.close();

	//Заполняем файлы нелинейной системы 
	state.t = 0, state.x = x_start;
	solve = DifferencialEquations(state.t, x_start, h, eps, 2);
	solve.SetStructParam(SP);
	solve.SetVecParam(pNon);

	for (int i = 0; (i < MaxIt) && (state.t < P.tmax); i++)
	{
		Nfi << state.t << " " << state.x[0] << endl;
		Nxg << state.t << " " << state.x[1] << endl;
		Ndfi << state.t << " " << state.x[2] << endl;
		Ndx << state.t << " " << state.x[3] << endl;
		Nu << state.t << " " << state.x * SP.theta << endl;
		state = solve.ComplitWithControl();
	}

	Nfi.close(), Nxg.close(), Ndfi.close(), Ndx.close(), Nu.close();

	string s;
	s = setting;
	s += "set grid\n";

	s += "set output '1R_fi.png'\n";
	s += "plot '" + way + "Lfi.txt' with lines lc 'red' lw 2 title 'Lfi', ";
	s += "'" + way + "NLfi.txt' with lines lc 'blue' lw 2 title 'NLfi'\n";

	s += "set output '1R_x.png'\n";
	s += "plot '" + way + "Lx.txt' with lines lc 'red' lw 2 title 'Lx', ";
	s += "'" + way + "NLx.txt' with lines lc 'blue' lw 2 title 'NLx'\n";

	s += "set output '1R_dfi.png'\n";
	s += "plot '" + way + "Ldfi.txt' with lines lc 'red' lw 2 title 'Ldfi', ";
	s += "'" + way + "NLdfi.txt' with lines lc 'blue' lw 2 title 'NLdfi'\n";

	s += "set output '1R_dx.png'\n";
	s += "plot '" + way + "Ldx.txt' with lines lc 'red' lw 2 title 'Ldx', ";
	s += "'" + way + "NLdx.txt' with lines lc 'blue' lw 2 title 'NLdx'\n";

	s += "set output '1R_u.png'\n";
	s += "plot '" + way + "Lu.txt' with lines lc 'red' lw 2 title 'Lu', ";
	s += "'" + way + "NLu.txt' with lines lc 'blue' lw 2 title 'NLu'\n";
	PLOT(s);

}

void GetImageLinObserver(ContinuosSystemParametrs& SP, Vector<double>& x_start)
{
	int n = 8;
	SystemState state;
	double eps = pow(10, -10), h = pow(10, -3), uVal;
	int MaxIt = 20000;

	//Заполняем файлы линеаризованной системы
	Vector<ofstream> data(n);
	ofstream u("u.txt");
	string name;

	for (int i = 0; i < n; i++)
	{
		name = to_string(i) + ".txt";
		data[i].open(name);
	}

	state.t = 0, state.x = x_start;
	DifferencialEquations solve(state.t, x_start, h, eps, 3);
	solve.SetStructParam(SP);

	for (int i = 0; (i < MaxIt) && (state.t < P.tmax); i++)
	{
		for (int j = 0; j < n; j++)
			data[j] << state.t << " " << state.x[j] << endl;
		uVal = 0;
		for (int i = 0; i < 4; i++)
			uVal += SP.theta[i] * state.x[i];
		u << state.t << " " << uVal << endl;
		state = solve.ComplitWithControl();
	}

	for (int i = 0; i < n; i++)
		data[i].close();
	u.close();

	string s;
	s = setting;
	s += "set grid\n";

	s += "set output '1L_OBS_fi.png'\n";
	s += "plot '" + way + "4.txt' with lines lc 'red' lw 2 title 'Lfi', ";
	s += "'" + way + "0.txt' with lines lc 'blue' lw 2 title 'LObsfi'\n";

	s += "set output '1L_OBS_x.png'\n";
	s += "plot '" + way + "5.txt' with lines lc 'red' lw 2 title 'Lx', ";
	s += "'" + way + "1.txt' with lines lc 'blue' lw 2 title 'LObsx'\n";

	s += "set output '1L_OBS_dfi.png'\n";
	s += "plot '" + way + "6.txt' with lines lc 'red' lw 2 title 'Ldfi', ";
	s += "'" + way + "2.txt' with lines lc 'blue' lw 2 title 'LObsdfi'\n";

	s += "set output '1L_OBS_dx.png'\n";
	s += "plot '" + way + "7.txt' with lines lc 'red' lw 2 title 'Ldx', ";
	s += "'" + way + "3.txt' with lines lc 'blue' lw 2 title 'LObsdx'\n";

	s += "set output '1L_OBS_u.png'\n";
	s += "plot '" + way + "u.txt' with lines lc 'green' lw 2 title 'Lu'\n";
	PLOT(s);

}

void GetImageNonLinObserver(ContinuosSystemParametrs& SP, Vector<double>& pNon, Vector<double>& x_start)
{
	int n = 8;
	SystemState state;
	double eps = pow(10, -10), h = pow(10, -3), uVal;
	int MaxIt = 20000;

	//Заполняем файлы линеаризованной системы
	Vector<ofstream> data(n);
	ofstream u("u.txt");
	string name;

	for (int i = 0; i < n; i++)
	{
		name = to_string(i) + ".txt";
		data[i].open(name);
	}

	state.t = 0, state.x = x_start;
	DifferencialEquations solve(state.t, x_start, h, eps, 4);
	solve.SetStructParam(SP);
	solve.SetVecParam(pNon);

	for (int i = 0; (i < MaxIt) && (state.t < P.tmax); i++)
	{
		for (int j = 0; j < n; j++)
			data[j] << state.t << " " << state.x[j] << endl;
		uVal = 0;
		for (int i = 0; i < 4; i++)
			uVal += SP.theta[i] * state.x[i];
		u << state.t << " " << uVal << endl;
		state = solve.ComplitWithControl();
	}

	for (int i = 0; i < n; i++)
		data[i].close();
	u.close();
	
	string s;
	s = setting;
	s += "set grid\n";

	s += "set output '1NL_OBS_fi.png'\n";
	s += "plot '" + way + "4.txt' with lines lc 'red' lw 2 title 'NLfi', ";
	s += "'" + way + "0.txt' with lines lc 'blue' lw 2 title 'NLObsfi'\n";

	s += "set output '1NL_OBS_x.png'\n";
	s += "plot '" + way + "5.txt' with lines lc 'red' lw 2 title 'NLx', ";
	s += "'" + way + "1.txt' with lines lc 'blue' lw 2 title 'NLObsx'\n";

	s += "set output '1NL_OBS_dfi.png'\n";
	s += "plot '" + way + "6.txt' with lines lc 'red' lw 2 title 'NLdfi', ";
	s += "'" + way + "2.txt' with lines lc 'blue' lw 2 title 'NLObsdfi'\n";

	s += "set output '1NL_OBS_dx.png'\n";
	s += "plot '" + way + "7.txt' with lines lc 'red' lw 2 title 'NLdx', ";
	s += "'" + way + "3.txt' with lines lc 'blue' lw 2 title 'NLObsdx'\n";

	s += "set output '1NL_OBS_u.png'\n";
	s += "plot '" + way + "u.txt' with lines lc 'green' lw 2 title 'NLu'\n";
	PLOT(s);

}
