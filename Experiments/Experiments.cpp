#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include "ODE.h"
#include "ORE.h"

using namespace std;

struct TaskParametrs
{
	double m, M, I, l, Bp, Beq, g, Kf, Ks, tmax;
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
void GetImageDiscrSystem(DiscreteSystemParametrs& DS, Vector<double>& x_start);
void GetImageDiscrSystemObs(DiscreteSystemParametrs& DS, Vector<double>& x_start, Vector<double>& ksi_start);
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
	DiscreteSystemParametrs DS;
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
	DS.Bd_impulse = DS.Ad * B;
	DS.Bd_const = IntegrateMatrixExp(A, DS.h) * B;

	//Получаем theta с заданными собственными числами
	REAL1 = true, REAL2 = true;
	EigV[0] = exp(-11.1856 * DS.h), EigV[1] = exp(-6.33331 * DS.h);
	EigV[2] = exp(-1 * DS.h), EigV[3] = exp(-2 * DS.h);
	GetControl(DS.Ad, DS.Bd_const, DS.thetad_const);
	PRINT(DS.thetad_const);
	GetControl(DS.Ad, DS.Bd_impulse, DS.thetad_impulse);
	PRINT(DS.thetad_impulse);

	DS.BTd_const = CellRowMultiply(DS.Bd_const, DS.thetad_const);
	DS.ABTd_const = DS.Ad + DS.BTd_const;

	DS.BTd_impulse = CellRowMultiply(DS.Bd_impulse, DS.thetad_impulse);
	DS.ABTd_impulse = DS.Ad + DS.BTd_impulse;

	//Начальные условия
	x_start[0] = 1, x_start[1] = 0.5, x_start[2] = 0.4, x_start[3] = -0.2;
	GetImageDiscrSystem(DS, x_start);

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

void GetImageDiscrSystem(DiscreteSystemParametrs& DS, Vector<double>& x_start)
{
	int MaxIt = 20000;
	double t;

	/*
	
		С кусочно-постоянным управлением
	
	*/

	t = 0;
	ORE solve(t, x_start, DS.h, 2);
	solve.SetSystParametrs(DS);
	solve.WriteToFileWithControl(P.tmax, "const", MaxIt);

	/*
	
		С импульсным управлением
	
	*/

	t = 0;
	solve = ORE(t, x_start, DS.h, 1);
	solve.SetSystParametrs(DS);
	solve.WriteToFileWithControl(P.tmax, "impulse", MaxIt);

	string s;
	s = setting;
	s += "set grid\n";

	s += "set output '1D_R_fi.png'\n";
	s += "plot '" + way + "const0.txt' with lines lc 'red' lw 2 title 'CONSTfi', ";
	s += "'" + way + "impulse0.txt' with lines lc 'blue' lw 2 title 'IMPULSEfi'\n";

	s += "set output '1D_R_x.png'\n";
	s += "plot '" + way + "const1.txt' with lines lc 'red' lw 2 title 'CONSTx', ";
	s += "'" + way + "impulse1.txt' with lines lc 'blue' lw 2 title 'IMPUSLEx'\n";

	s += "set output '1D_R_dfi.png'\n";
	s += "plot '" + way + "const2.txt' with lines lc 'red' lw 2 title 'CONSTdfi', ";
	s += "'" + way + "impulse2.txt' with lines lc 'blue' lw 2 title 'IMPUSLEdfi'\n";

	s += "set output '1D_R_dx.png'\n";
	s += "plot '" + way + "const3.txt' with lines lc 'red' lw 2 title 'CONSTdx', ";
	s += "'" + way + "impulse3.txt' with lines lc 'blue' lw 2 title 'IMPUSLEdx'\n";

	s += "set output '1D_R_u.png'\n";
	s += "plot '" + way + "constcontrol.txt' with lines lc 'red' lw 2 title 'CONSTu', ";
	s += "'" + way + "impulsecontrol.txt' with lines lc 'blue' lw 2 title 'IMPUSLEu'\n";
	PLOT(s);
}

void GetImageDiscrSystemObs(DiscreteSystemParametrs& DS, Vector<double>& x_start, Vector<double>& ksi_start)
{
	Vector<double> x(2 * x_start.size());
	int MaxIt = 20000;
	double t;

	/*
	
		Для кусочно-постоянного управления
	
	*/
	

	t = 0;
	for (int i = 0; i < x_start.size(); i++)
	{
		x[i] = ksi_start[i];
		x[i + 4] = x_start[i];
	}
	ORE solve(t, x, DS.h, 4);
	solve.SetSystParametrs(DS);
	solve.WriteToFileWithControl(P.tmax, "constOBS", MaxIt);

	string s;
	s = setting;
	s += "set grid\n";

	s += "set output '1D_OBS_CONST_fi.png'\n";
	s += "plot '" + way + "constOBS4.txt' with lines lc 'red' lw 2 title 'CONSTfi', ";
	s += "'" + way + "constOBS0.txt' with lines lc 'blue' lw 2 title 'CONSTObsfi\n";

	s += "set output '1D_OBS_CONST_x.png'\n";
	s += "plot '" + way + "constOBS5.txt' with lines lc 'red' lw 2 title 'CONSTx', ";
	s += "'" + way + "constOBS1.txt' with lines lc 'blue' lw 2 title 'CONSTObsx\n";

	s += "set output '1D_OBS_CONST_dfi.png'\n";
	s += "plot '" + way + "constOBS6.txt' with lines lc 'red' lw 2 title 'CONSTdfi', ";
	s += "'" + way + "constOBS2.txt' with lines lc 'blue' lw 2 title 'CONSTObsdfi\n";

	s += "set output '1D_OBS_CONST_dx.png'\n";
	s += "plot '" + way + "constOBS7.txt' with lines lc 'red' lw 2 title 'CONSTdx', ";
	s += "'" + way + "constOBS3.txt' with lines lc 'blue' lw 2 title 'CONSTObsdx\n";

	s += "set output '1D_OBS_CONST_u.png'\n";
	s += "plot '" + way + "constOBScontrol.txt' with lines lc 'green' lw 2 title 'CONSTu'\n";
	PLOT(s);

	/*
	
		Для импульсного управления
	
	*/

	t = 0;
	for (int i = 0; i < x_start.size(); i++)
	{
		x[i] = ksi_start[i];
		x[i + 4] = x_start[i];
	}
	solve = ORE(t, x, DS.h, 3);
	solve.SetSystParametrs(DS);
	solve.WriteToFileWithControl(P.tmax, "impulseOBS", MaxIt);
	
	s = setting;
	s += "set grid\n";

	s += "set output '1D_OBS_IMPUSLE_fi.png'\n";
	s += "plot '" + way + "impulseOBS4.txt' with lines lc 'red' lw 2 title 'IMPUSLEfi', ";
	s += "'" + way + "impulseOBS0.txt' with lines lc 'blue' lw 2 title 'IMPUSLEObsfi\n";

	s += "set output '1D_OBS_IMPUSLE_x.png'\n";
	s += "plot '" + way + "impulseOBS5.txt' with lines lc 'red' lw 2 title 'IMPUSLEx', ";
	s += "'" + way + "impulseOBS1.txt' with lines lc 'blue' lw 2 title 'IMPUSLEObsx\n";

	s += "set output '1D_OBS_IMPUSLE_dfi.png'\n";
	s += "plot '" + way + "impulseOBS6.txt' with lines lc 'red' lw 2 title 'IMPUSLEdfi', ";
	s += "'" + way + "impulseOBS2.txt' with lines lc 'blue' lw 2 title 'IMPUSLEObsdfi\n";

	s += "set output '1D_OBS_IMPUSLE_dx.png'\n";
	s += "plot '" + way + "impulseOBS7.txt' with lines lc 'red' lw 2 title 'IMPUSLEdx', ";
	s += "'" + way + "impulseOBS3.txt' with lines lc 'blue' lw 2 title 'IMPUSLEObsdx\n";

	s += "set output '1D_OBS_IMPUSLE_u.png'\n";
	s += "plot '" + way + "impulseOBScontrol.txt' with lines lc 'green' lw 2 title 'IMPUSLEu'\n";
	PLOT(s);

}

void GetImageOriginalSystem(ContinuosSystemParametrs& SP, Vector<double>& pNon, Vector<double>& x_start)
{
	SystemState state;
	double eps = pow(10, -10), h = pow(10, -3);
	int MaxIt = 20000;

	state.t = 0, state.x = x_start;
	ODE4 solve(state.t, x_start, h, eps, 1);
	solve.SetStructParam(SP);
	solve.WriteToFileWithControl(P.tmax, "L", MaxIt);

	//Заполняем файлы нелинейной системы 
	state.t = 0, state.x = x_start;
	solve = ODE4(state.t, x_start, h, eps, 2);
	solve.SetStructParam(SP);
	solve.SetVecParam(pNon);
	solve.WriteToFileWithControl(P.tmax, "NL", MaxIt);

	string s;
	s = setting;
	s += "set grid\n";

	s += "set output '1R_fi.png'\n";
	s += "plot '" + way + "L0.txt' with lines lc 'red' lw 2 title 'Lfi', ";
	s += "'" + way + "NL0.txt' with lines lc 'blue' lw 2 title 'NLfi'\n";

	s += "set output '1R_x.png'\n";
	s += "plot '" + way + "L1.txt' with lines lc 'red' lw 2 title 'Lx', ";
	s += "'" + way + "NL1.txt' with lines lc 'blue' lw 2 title 'NLx'\n";

	s += "set output '1R_dfi.png'\n";
	s += "plot '" + way + "L2.txt' with lines lc 'red' lw 2 title 'Ldfi', ";
	s += "'" + way + "NL2.txt' with lines lc 'blue' lw 2 title 'NLdfi'\n";

	s += "set output '1R_dx.png'\n";
	s += "plot '" + way + "L3.txt' with lines lc 'red' lw 2 title 'Ldx', ";
	s += "'" + way + "NL3.txt' with lines lc 'blue' lw 2 title 'NLdx'\n";

	s += "set output '1R_u.png'\n";
	s += "plot '" + way + "Lcontrol.txt' with lines lc 'red' lw 2 title 'Lu', ";
	s += "'" + way + "NLcontrol.txt' with lines lc 'blue' lw 2 title 'NLu'\n";
	PLOT(s);

}

void GetImageLinObserver(ContinuosSystemParametrs& SP, Vector<double>& x_start)
{
	int n = 8;
	SystemState state;
	double eps = pow(10, -10), h = pow(10, -3);
	int MaxIt = 20000;

	state.t = 0, state.x = x_start;
	ODE4 solve(state.t, x_start, h, eps, 3);
	solve.SetStructParam(SP);
	solve.WriteToFileWithControl(P.tmax, "LOBS", MaxIt);

	string s;
	s = setting;
	s += "set grid\n";

	s += "set output '1L_OBS_fi.png'\n";
	s += "plot '" + way + "LOBS4.txt' with lines lc 'red' lw 2 title 'Lfi', ";
	s += "'" + way + "LOBS0.txt' with lines lc 'blue' lw 2 title 'LObsfi'\n";

	s += "set output '1L_OBS_x.png'\n";
	s += "plot '" + way + "LOBS5.txt' with lines lc 'red' lw 2 title 'Lx', ";
	s += "'" + way + "LOBS1.txt' with lines lc 'blue' lw 2 title 'LObsx'\n";

	s += "set output '1L_OBS_dfi.png'\n";
	s += "plot '" + way + "LOBS6.txt' with lines lc 'red' lw 2 title 'Ldfi', ";
	s += "'" + way + "LOBS2.txt' with lines lc 'blue' lw 2 title 'LObsdfi'\n";

	s += "set output '1L_OBS_dx.png'\n";
	s += "plot '" + way + "LOBS7.txt' with lines lc 'red' lw 2 title 'Ldx', ";
	s += "'" + way + "LOBS3.txt' with lines lc 'blue' lw 2 title 'LObsdx'\n";

	s += "set output '1L_OBS_u.png'\n";
	s += "plot '" + way + "LOBScontrol.txt' with lines lc 'green' lw 2 title 'Lu'\n";
	PLOT(s);

}

void GetImageNonLinObserver(ContinuosSystemParametrs& SP, Vector<double>& pNon, Vector<double>& x_start)
{
	int n = 8;
	SystemState state;
	double eps = pow(10, -10), h = pow(10, -3), uVal;
	int MaxIt = 20000;

	state.t = 0, state.x = x_start;
	ODE4 solve(state.t, x_start, h, eps, 4);
	solve.SetStructParam(SP);
	solve.SetVecParam(pNon);
	solve.WriteToFileWithControl(P.tmax, "NLOBS", MaxIt);
	
	string s;
	s = setting;
	s += "set grid\n";

	s += "set output '1NL_OBS_fi.png'\n";
	s += "plot '" + way + "NLOBS4.txt' with lines lc 'red' lw 2 title 'NLfi', ";
	s += "'" + way + "NLOBS0.txt' with lines lc 'blue' lw 2 title 'NLObsfi'\n";

	s += "set output '1NL_OBS_x.png'\n";
	s += "plot '" + way + "NLOBS5.txt' with lines lc 'red' lw 2 title 'NLx', ";
	s += "'" + way + "NLOBS1.txt' with lines lc 'blue' lw 2 title 'NLObsx'\n";

	s += "set output '1NL_OBS_dfi.png'\n";
	s += "plot '" + way + "NLOBS6.txt' with lines lc 'red' lw 2 title 'NLdfi', ";
	s += "'" + way + "NLOBS2.txt' with lines lc 'blue' lw 2 title 'NLObsdfi'\n";

	s += "set output '1NL_OBS_dx.png'\n";
	s += "plot '" + way + "NLOBS7.txt' with lines lc 'red' lw 2 title 'NLdx', ";
	s += "'" + way + "NLOBS3.txt' with lines lc 'blue' lw 2 title 'NLObsdx'\n";

	s += "set output '1NL_OBS_u.png'\n";
	s += "plot '" + way + "NLOBScontrol.txt' with lines lc 'green' lw 2 title 'NLu'\n";
	PLOT(s);

}
