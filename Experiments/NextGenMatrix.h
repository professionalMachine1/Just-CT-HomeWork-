#pragma once

#ifndef _NEXTGENMATRIX_H
#define  _NEXTGENMATRIX_H

#include <initializer_list>
#include <iostream>
#include <math.h>

template <class ValType>
class Vector
{
protected:
	ValType* pVector;
	int size_;
public:
	//������������
	Vector(int _size = 1);
	Vector(const Vector& Vec);
	Vector(const ValType& Val, int _size);
	Vector(std::initializer_list<ValType> in_list);

	//�������� ������������
	Vector& operator= (const Vector& Vec);
	Vector& operator= (std::initializer_list<ValType> in_list);

	//�������� � ���������
	Vector  operator+(const Vector& Vec);
	Vector  operator-(const Vector& Vec);
	ValType  operator*(const Vector& Vec);
	Vector operator+=(const Vector& Vec);
	Vector operator-=(const Vector& Vec);

	//�������� � �������
	Vector  operator*(const ValType& Val);
	Vector  operator/(const ValType& Val);

	//������ �������
	int size() const { return size_; }

	//������ �������� �� �������
	//return (ind < 0) || (ind >= size) ? ValType() : pVector[i];
	ValType& operator[](const int& ind) const { return pVector[ind]; }

	//���������� ��������
	bool operator==(const Vector& Vec) const;
	bool operator!=(const Vector& Vec) const;
	bool operator>(const Vector& Vec) const;
	bool operator<(const Vector& Vec) const;
	bool operator>=(const Vector& Vec) const;
	bool operator<=(const Vector& Vec) const;

	//��������� � ���������
	Vector& operator++();
	Vector operator++(int);
	Vector& operator--();
	Vector operator--(int);

	//������� �����-������
	friend std::istream& operator >>(std::istream& is, Vector& Vec)
	{
		for (int i = 0; i < Vec.size_; i++)
			is >> Vec.pVector[i];
		return is;
	}
	friend std::ostream& operator<< (std::ostream& os, const Vector& Vec)
	{
		for (int i = 0; i < Vec.size_; i++)
			os << Vec.pVector[i] << "   ";
		os << "\n";
		return os;
	}

	~Vector();

	//����� ��� ������ ���������� ��������� �������� �������
	//Vector<ValType> operator*(const ValType& Val, const Vector<ValType>& Vec)
	//int MaxValue(const Vector<ValType>& Vec)
	//int MinValue(const Vector<ValType>& Vec)
	//Vector<ValType> abs(const Vector<ValType>& Vec)
	//ValType operator*(const ValType& Val1, const ValType& Val2)

};

/*

	������������

*/

template <class ValType>
Vector<ValType>::Vector(int _size)
{
	size_ = _size;
	pVector = new ValType[size_];
	for (int i = 0; i < size_; i++)
		pVector[i] = ValType();
}

template <class ValType>
Vector<ValType>::Vector(const ValType& Val, int _size)
{
	size_ = _size;
	pVector = new ValType[size_];
	for (int i = 0; i < size_; i++)
		pVector[i] = Val;
}

template <class ValType>
Vector<ValType>::Vector(const Vector& Vec)
{
	size_ = Vec.size_;
	pVector = new ValType[size_];
	for (int i = 0; i < size_; i++)
		pVector[i] = Vec.pVector[i];
}

template <typename ValType>
Vector<ValType>::Vector(std::initializer_list<ValType> in_list)
{
	size_ = in_list.size();
	pVector = new ValType[size_];

	auto it = in_list.begin();
	for (int i = 0; i < size_; i++)
		pVector[i] = *it, it++;
}

/*

	�������� ������������

*/

template <class ValType>
Vector<ValType>& Vector<ValType>::operator= (const Vector<ValType>& Vec)
{
	if (this != &Vec)
	{
		if (size_ != Vec.size_)
		{
			delete[]pVector;
			size_ = Vec.size_;
			pVector = new ValType[size_];
		}
		for (int i = 0; i < size_; i++)
			pVector[i] = Vec.pVector[i];
	}
	return *this;
}

template <class ValType>
Vector<ValType>& Vector<ValType>::operator= (std::initializer_list<ValType> in_list)
{
	if (size_ != in_list.size())
	{
		delete[]pVector;
		size_ = in_list.size();
		pVector = new ValType[size_];
	}

	auto it = in_list.begin();
	for (int i = 0; i < size_; i++)
		pVector[i] = *it, it++;
	return *this;
}

/*

	�������� � ���������

*/

template <class ValType>
Vector<ValType>  Vector<ValType>::operator+(const Vector<ValType>& Vec)
{
	//int n = size < Vec.size ? size : Vec.size;
	Vector<ValType> res(size_);
	for (int i = 0; i < size_; i++)
		res.pVector[i] = pVector[i] + Vec.pVector[i];
	return res;
}

template <class ValType>
Vector<ValType>  Vector<ValType>::operator-(const Vector<ValType>& Vec)
{
	Vector<ValType> res(size_);
	for (int i = 0; i < size_; i++)
		res.pVector[i] = pVector[i] - Vec.pVector[i];
	return res;
}

template <class ValType>
ValType  Vector<ValType>::operator*(const Vector<ValType>& Vec)
{
	ValType res = ValType();
	for (int i = 0; i < size_; i++)
		res += pVector[i] * Vec.pVector[i];
	return res;
}

template <class ValType>
Vector<ValType> Vector<ValType>::operator+=(const Vector<ValType>& Vec)
{
	Vector<ValType> res = *this;
	for (int i = 0; i < size_; i++)
		pVector[i] += Vec.pVector[i];
	return res;
}

template <class ValType>
Vector<ValType> Vector<ValType>::operator-=(const Vector<ValType>& Vec)
{
	Vector<ValType> res = *this;
	for (int i = 0; i < size_; i++)
		pVector[i] -= Vec.pVector[i];
	return res;
}

/*

	�������� � �������

*/

template <class ValType>
Vector<ValType>  Vector<ValType>::operator*(const ValType& Val)
{
	Vector<ValType> res(size_);
	for (int i = 0; i < size_; i++)
		res.pVector[i] = pVector[i] * Val;
	return res;
}

template <class ValType>
Vector<ValType>  Vector<ValType>::operator/(const ValType& Val)
{
	Vector<ValType> res(size_);
	for (int i = 0; i < size_; i++)
		res.pVector[i] = pVector[i] / Val;
	return res;
}

template <class ValType>
ValType operator*(const ValType& Val1, const ValType& Val2)
{
	return Val1 * Val2;
}

template <class ValType>
Vector<ValType> operator*(const ValType& Val, const Vector<ValType>& Vec)
{
	Vector<ValType> res(Vec.size());
	for (int i = 0; i < Vec.size(); i++)
		res[i] = Vec[i] * Val;
	return res;
}

/*

	���������� ��������

*/

template <class ValType>
bool Vector<ValType>::operator==(const Vector<ValType>& Vec) const
{
	if (size_ != Vec.size_)
		return false;
	for (int i = 0; i < size_; i++)
		if (pVector[i] != Vec.pVector[i])
			return false;
	return true;
}

template <class ValType>
bool Vector<ValType>::operator!=(const Vector<ValType>& Vec) const
{
	return !(*this == Vec);
}

template <class ValType>
bool Vector<ValType>::operator>(const Vector<ValType>& Vec) const
{
	if (size_ != Vec.size_)
		return false;
	for (int i = 0; i < size_; i++)
		if (pVector[i] <= Vec.pVector[i])
			return false;
	return true;
}

template <class ValType>
bool Vector<ValType>::operator<(const Vector<ValType>& Vec) const
{
	if (size_ != Vec.size_)
		return false;
	for (int i = 0; i < size_; i++)
		if (pVector[i] >= Vec.pVector[i])
			return false;
	return true;
}

template <class ValType>
bool Vector<ValType>::operator>=(const Vector<ValType>& Vec) const
{
	if (size_ != Vec.size_)
		return false;
	for (int i = 0; i < size_; i++)
		if (pVector[i] < Vec.pVector[i])
			return false;
	return true;
}

template <class ValType>
bool Vector<ValType>::operator<=(const Vector<ValType>& Vec) const
{
	if (size_ != Vec.size_)
		return false;
	for (int i = 0; i < size_; i++)
		if (pVector[i] > Vec.pVector[i])
			return false;
	return true;
}

/*

	��������� � ���������

*/

template <class ValType>
Vector<ValType>& Vector<ValType>::operator++()
{
	for (int i = 0; i < size_; i++)
		++pVector[i];
	return *this;
}

template <class ValType>
Vector<ValType> Vector<ValType>::operator++(int)
{
	Vector<ValType> res = *this;
	++* this;
	return res;
}

template <class ValType>
Vector<ValType>& Vector<ValType>::operator--()
{
	for (int i = 0; i < size_; i++)
		--pVector[i];
	return *this;
}

template <class ValType>
Vector<ValType> Vector<ValType>::operator--(int)
{
	Vector<ValType> res = *this;
	--* this;
	return res;
}

/*

	�������� ��������

*/

template <class ValType>
int MaxValue(const Vector<ValType>& Vec)
{
	ValType Max = Vec[0];
	int IndMax = 0;
	for (int i = 1; i < Vec.size(); i++)
		if (Vec[i] > Max)
		{
			Max = Vec[i];
			IndMax = i;
		}
	return IndMax;
}

template <class ValType>
int MinValue(const Vector<ValType>& Vec)
{
	ValType Min = Vec[0];
	int IndMin = 0;
	for (int i = 1; i < Vec.size(); i++)
		if (Vec[i] < Min)
		{
			Min = Vec[i];
			IndMin = i;
		}
	return IndMin;
}

template <class ValType>
Vector<ValType> abs(const Vector<ValType>& Vec)
{
	Vector<ValType> res(Vec.size());
	for (int i = 0; i < Vec.size(); i++)
		res[i] = fabs(Vec[i]);
	return res;
}

/*

	����������

*/

template <class ValType>
Vector<ValType>::~Vector()
{
	delete[]pVector;
}

/*



	����� ������



*/

template <class ValType>
class Matrix : public Vector<Vector<ValType> >
{
public:
	//������������ 
	Matrix(int _size = 1);
	Matrix(const Matrix& mt);
	Matrix(int rows, int cells);
	Matrix(const Vector<Vector<ValType> >& mt);
	Matrix(const ValType& Val, int rows, int cells);
	Matrix(std::initializer_list<std::initializer_list<ValType>> in_list);

	//�������� ������������
	Matrix& operator= (const Matrix& mt);
	Matrix& operator= (std::initializer_list<std::initializer_list<ValType>> in_list);

	//�������� � ���������
	Matrix operator*(const Matrix& mt);
	Matrix operator+(const Matrix& mt);
	Matrix operator-(const Matrix& mt);
	Matrix operator+=(const Matrix& mt);
	Matrix operator-=(const Matrix& mt);
	Matrix operator*=(const Matrix& mt);

	//������ �� ��������� �� �������
	Vector<ValType> GetCell(const int& ind);
	void PutCell(const int& ind, const Vector<ValType>& Vec);

	//�������� � �������
	Matrix  operator*(const ValType& Val);
	Matrix  operator/(const ValType& Val);

	//���������� ��������
	bool operator==(const Matrix& mt) const;
	bool operator!=(const Matrix& mt) const;
	bool operator>(const Matrix& mt) const;
	bool operator<(const Matrix& mt) const;
	bool operator>=(const Matrix& mt) const;
	bool operator<=(const Matrix& mt) const;

	//��������� � ���������
	Matrix& operator++();
	Matrix operator++(int);
	Matrix& operator--();
	Matrix operator--(int);

	friend std::istream& operator >> (std::istream& is, Matrix& mt)
	{
		for (int i = 0; i < mt.size_; i++)
			is >> mt.pVector[i];
		return is;
	}
	friend std::ostream& operator << (std::ostream& os, const Matrix& mt)
	{
		os << std::endl;
		for (int i = 0; i < mt.size_; i++)
			os << mt.pVector[i];
		return os;
	}

	//��������� �������� ������� ��� ������
	//Vector<ValType> operator*(const Vector<ValType>& Vec, const Matrix<ValType>& mt)
	//Vector<ValType> operator*(const Matrix<ValType>& mt, const Vector<ValType>& Vec)
	//Matrix<ValType> Single(int size_)
	//Matrix<ValType> pow(const Matrix<ValType>& mt, const int& degree)
	//Matrix<ValType> InverseMatrix(const Matrix<ValType>& mt)
	//Matrix<ValType> TransposedMatrix(const Matrix<ValType>& mt)
	//ValType Determinant(const Matrix<ValType>& mt)
	//Matrix<ValType> Minor(const Matrix<ValType>& mt, int ntrow, int ntcell)
	//Matrix<ValType> CellRowMultiply(const Vector<ValType>& cell, const Vector<ValType>& row)
	//Matrix<ValType> operator*(const ValType& Val, const Matrix<ValType>& mt)

};

/*

	������������

*/

template <class ValType>
Matrix<ValType>::Matrix(int _size) : Vector<Vector<ValType>>(_size)
{
	for (int i = 0; i < Matrix::size_; i++)
		Matrix::pVector[i] = Vector<ValType>(Matrix::size_);
}

template <class ValType>
Matrix<ValType>::Matrix(const Matrix<ValType>& mt) : Vector<Vector<ValType>>(mt) { }

template <class ValType>
Matrix<ValType>::Matrix(int rows, int cells) : Vector<Vector<ValType>>(rows)
{
	for (int i = 0; i < Matrix::size_; i++)
		Matrix::pVector[i] = Vector<ValType>(cells);
}

template <class ValType>
Matrix<ValType>::Matrix(const Vector<Vector<ValType> >& mt) : Vector<Vector<ValType>>(mt) { }

template <class ValType>
Matrix<ValType>::Matrix(const ValType& Val, int rows, int cells) : Vector<Vector<ValType>>(rows)
{
	for (int i = 0; i < Matrix::size_; i++)
		Matrix::pVector[i] = Vector<ValType>(Val, cells);
}

template <class ValType>
Matrix<ValType>::Matrix(std::initializer_list<std::initializer_list<ValType>> in_list)
{
	Matrix::size_ = in_list.size();
	Matrix::pVector = new Vector<ValType>[Matrix::size_];

	auto it = in_list.begin();
	for (int i = 0; i < Matrix::size_; i++)
		Matrix::pVector[i] = *it, it++;
}

/*

	�������� ������������

*/

template <class ValType>
Matrix<ValType>& Matrix<ValType>::operator= (const Matrix& mt)
{
	if (this != &mt)
	{
		if (Matrix::size_ != mt.size_)
		{
			delete[]Matrix::pVector;
			Matrix::size_ = mt.size_;
			Matrix::pVector = new Vector<ValType>[Matrix::size_];
		}
		for (int i = 0; i < Matrix::size_; i++)
			Matrix::pVector[i] = mt.pVector[i];
	}
	return *this;
}

template <class ValType>
Matrix<ValType>& Matrix<ValType>::operator= (std::initializer_list<std::initializer_list<ValType>> in_list)
{
	if (Matrix::size_ != in_list.size())
	{
		delete[]Matrix::pVector;
		Matrix::size_ = in_list.size();
		Matrix::pVector = new Vector<ValType>[Matrix::size_];
	}

	auto it = in_list.begin();
	for (int i = 0; i < Matrix::size_; i++)
		Matrix::pVector[i] = *it, it++;

	return *this;
}

/*

	�������� � ���������

*/

template <class ValType>
Matrix<ValType>  Matrix<ValType>::operator*(const Matrix<ValType>& mt)
{
	int row = Matrix::size_, cell = mt.pVector[0].size(), n = mt.size_;
	Matrix<ValType> res(row, cell);
	for (int i = 0; i < row; i++)
		for (int j = 0; j < cell; j++)
			for (int k = 0; k < n; k++)
				res.pVector[i][j] += Matrix::pVector[i][k] * mt.pVector[k][j];
	return res;
}

template <class ValType>
Matrix<ValType>  Matrix<ValType>::operator+(const Matrix<ValType>& mt)
{
	return Vector<Vector<ValType>>::operator+(mt);
}

template <class ValType>
Matrix<ValType>  Matrix<ValType>::operator+=(const Matrix<ValType>& mt)
{
	return Vector<Vector<ValType>>::operator+=(mt);
}

template <class ValType>
Matrix<ValType>  Matrix<ValType>::operator-=(const Matrix<ValType>& mt)
{
	return Vector<Vector<ValType>>::operator-=(mt);
}

template <class ValType>
Matrix<ValType>  Matrix<ValType>::operator-(const Matrix<ValType>& mt)
{
	return Vector<Vector<ValType>>::operator-(mt);
}

template <class ValType>
Matrix<ValType>  Matrix<ValType>::operator*=(const Matrix<ValType>& mt)
{
	Matrix<ValType> res = *this;
	*this = *this * mt;
	return res;
}

/*

	������ �� ��������� �� �������

*/

template <class ValType>
Vector<ValType> Matrix<ValType>::GetCell(const int& ind)
{
	Vector<ValType>	res(Matrix::size_);
	for (int i = 0; i < Matrix::size_; i++)
		res[i] = Matrix::pVector[i][ind];
	return res;
}

template <class ValType>
void Matrix<ValType>::PutCell(const int& ind, const Vector<ValType>& Vec)
{
	for (int i = 0; i < Matrix::size_; i++)
		Matrix::pVector[i][ind] = Vec[i];
}

/*

	�������� � �������

*/

template <class ValType>
Matrix<ValType>  Matrix<ValType>::operator*(const ValType& Val)
{
	Matrix<ValType> res = *this;
	for (int i = 0; i < res.size_; i++)
		res[i] = Matrix::pVector[i] * Val;
	return res;
}

template <class ValType>
Matrix<ValType>  Matrix<ValType>::operator/(const ValType& Val)
{
	Matrix<ValType> res = *this;
	for (int i = 0; i < res.size_; i++)
		res[i] = Matrix::pVector[i] / Val;
	return res;
}

/*

	���������� ��������

*/

template <class ValType>
bool Matrix<ValType>::operator==(const Matrix<ValType>& mt) const
{
	return Vector<Vector<ValType>>::operator==(mt);
}

template <class ValType>
bool Matrix<ValType>::operator!=(const Matrix<ValType>& mt) const
{
	return Vector<Vector<ValType>>::operator!=(mt);
}

template <class ValType>
bool Matrix<ValType>::operator>(const Matrix<ValType>& mt) const
{
	return Vector<Vector<ValType>>::operator>(mt);
}

template <class ValType>
bool Matrix<ValType>::operator<(const Matrix<ValType>& mt) const
{
	return Vector<Vector<ValType>>::operator<(mt);
}

template <class ValType>
bool Matrix<ValType>::operator>=(const Matrix<ValType>& mt) const
{
	return Vector<Vector<ValType>>::operator>=(mt);
}

template <class ValType>
bool Matrix<ValType>::operator<=(const Matrix<ValType>& mt) const
{
	return Vector<Vector<ValType>>::operator<=(mt);
}

/*

	��������� � ���������

*/

template <class ValType>
Matrix<ValType>& Matrix<ValType>::operator++()
{
	*this = Vector<Vector<ValType>>::operator++();
	return *this;
}

template <class ValType>
Matrix<ValType> Matrix<ValType>::operator++(int)
{
	return Vector<Vector<ValType>>::operator++(0);
}

template <class ValType>
Matrix<ValType>& Matrix<ValType>::operator--()
{
	*this = Vector<Vector<ValType>>::operator--();
	return *this;
}

template <class ValType>
Matrix<ValType> Matrix<ValType>::operator--(int)
{
	return Vector<Vector<ValType>>::operator--(0);
}

/*

	�������

*/

template <class ValType>
Vector<ValType> operator*(const Vector<ValType>& Vec, const Matrix<ValType>& mt)
{
	int n = Vec.size(), cell = mt[0].size();
	Vector<ValType> res(cell);
	for (int i = 0; i < cell; i++)
		for (int j = 0; j < n; j++)
			res[i] += Vec[j] * mt[j][i];
	return res;
}

template <class ValType>
Vector<ValType> operator*(const Matrix<ValType>& mt, const Vector<ValType>& Vec)
{
	Vector<ValType> res(mt.size());
	for (int i = 0; i < mt.size(); i++)
		res[i] = mt[i] * Vec;
	return res;
}

template <class ValType>
Matrix<ValType> Single(int size_)
{
	Matrix<ValType> res(size_);
	for (int i = 0; i < res.size(); i++)
		res[i][i] = 1;
	return res;
}

template <class ValType>
Matrix<ValType> pow(const Matrix<ValType>& mt, const int& degree)
{
	Matrix<ValType> res = Single<ValType>(mt.size());
	for (int i = 1; i <= degree; i++)
		res *= mt;
	return res;
}

template <class ValType>
Matrix<ValType> InverseMatrix(const Matrix<ValType>& mt)
{
	Matrix<ValType> res = Single<ValType>(mt.size()), matr = mt;
	ValType sup;
	for (int i = 1; i < matr.size(); i++)
	{
		if (fabs(matr[i - 1][i - 1]) <= pow(10, -14))
			for (int j = i; j < matr.size(); j++)
			{
				matr[i - 1] = matr[i - 1] + matr[j];
				res[i - 1] = res[i - 1] + res[j];
			}
		for (int j = i; j < matr.size(); j++)
		{
			sup = matr[j][i - 1] / matr[i - 1][i - 1];
			matr[j] = matr[j] - matr[i - 1] * sup;
			res[j] = res[j] - res[i - 1] * sup;
		}
	}
	for (int i = matr.size() - 2; i >= 0; i--)
	{
		for (int j = i; j >= 0; j--)
		{
			sup = matr[j][i + 1] / matr[i + 1][i + 1];
			matr[j] = matr[j] - matr[i + 1] * sup;
			res[j] = res[j] - res[i + 1] * sup;
		}
	}
	for (int i = 0; i < matr.size(); i++)
		res[i] = res[i] / matr[i][i];
	return res;
}

template <class ValType>
Matrix<ValType> TransposedMatrix(const Matrix<ValType>& mt)
{
	int row = mt.size(), cell = mt[0].size();
	Matrix<ValType> res(cell, row);
	for (int i = 0; i < cell; i++)
		for (int j = 0; j < row; j++)
			res[i][j] = mt[j][i];
	return res;
}

template <class ValType>
ValType Determinant(const Matrix<ValType>& mt)
{
	ValType res = ValType();

	if (mt.size() == 1)
		return mt[0][0];
	if (mt.size() == 2)
		return mt[0][0] * mt[1][1] - mt[0][1] * mt[1][0];

	Matrix<ValType> SupMT(mt.size() - 1);
	for (int i = 0; i < mt.size(); i++)
	{
		SupMT = Minor(mt, 0, i);
		res += pow(-1, i) * Determinant(SupMT) * mt[0][i];
	}

	return res;
}

template <class ValType>
Matrix<ValType> Minor(const Matrix<ValType>& mt, int ntrow, int ntcell)
{
	int row = mt.size(), cell = mt[0].size();
	Matrix<ValType> res(row - 1, cell - 1);
	for (int i = 0, ires = 0; i < row; i++)
	{
		if (i != ntrow)
		{
			for (int j = 0, jres = 0; j < cell; j++)
				if (j != ntcell)
					res[ires][jres] = mt[i][j], jres++;
			ires++;
		}
	}
	return res;
}

template <class ValType>
Matrix<ValType> CellRowMultiply(const Vector<ValType>& cell, const Vector<ValType>& row)
{
	Matrix<ValType> res(cell.size(), row.size());
	for (int i = 0; i < cell.size(); i++)
		for (int j = 0; j < row.size(); j++)
			res[i][j] = cell[i] * row[j];
	return res;
}

template <class ValType>
Matrix<ValType> operator*(const ValType& Val, const Matrix<ValType>& mt)
{
	Matrix<ValType> res(mt.size());
	for (int i = 0; i < mt.size(); i++)
		res[i] = Val * mt[i];
	return res;
}

#endif