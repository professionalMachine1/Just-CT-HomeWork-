#pragma once

#ifndef _TMATRIX_H_
#define _TMATRIX_H_

#include <iostream>
#include <iomanip>

using namespace std;

template <class ValType>
class TVector
{
protected:
	ValType* pVector;
	int size;
public:
	TVector(int s = 10);
	TVector(const TVector& v);
	~TVector();
	int Size() { return size; }
	ValType& operator[](int ind);
	bool operator==(const TVector& v) const;
	bool operator!=(const TVector& v) const;
	TVector& operator=(const TVector& v);
	TVector& operator=(const ValType& num);


	TVector  operator+(const ValType& num);
	TVector  operator-(const ValType& num);
	TVector  operator*(const ValType& num);
	TVector  operator/(const ValType& num);


	TVector  operator+(const TVector& v);
	TVector  operator-(const TVector& v);
	ValType  operator*(const TVector& v);

	int Max();
	int Min();

	TVector Absolute(const TVector& v);

	friend istream& operator >> (istream& is, TVector& v)
	{
		for (int i = 0; i < v.size; i++)
			is >> v.pVector[i];
		return is;
	}
	friend ostream& operator<<(ostream& os, const TVector& v)
	{
		for (int i = 0; i < v.size; i++)
			os << v.pVector[i] << "   ";
		os << "\n";
		return os;
	}
};

template <class ValType>
TVector<ValType>::TVector(int s)
{
	size = s;
	pVector = new ValType[size];
	for (int i = 0; i < size; i++)
		pVector[i] = 0;
}

template <class ValType>
TVector<ValType>::TVector(const TVector<ValType>& v)
{
	size = v.size;
	pVector = new ValType[size];
	for (int i = 0; i < size; i++)
		pVector[i] = v.pVector[i];
}

template <class ValType>
TVector<ValType>::~TVector()
{
	delete[]pVector;
}

template <class ValType>
ValType& TVector<ValType>::operator[](int ind)
{
	return pVector[ind];
}

template <class ValType>
bool TVector<ValType>::operator==(const TVector<ValType>& v) const
{
	if (size != v.size)
		return false;
	for (int i = 0; i < size; i++)
		if (pVector[i] != v.pVector[i])
			return false;
	return true;
}

template <class ValType>
bool TVector<ValType>::operator!=(const TVector<ValType>& v) const
{
	if (size != v.size)
		return true;
	for (int i = 0; i < size; i++)
		if (pVector[i] != v.pVector[i])
			return true;
	return false;
}

template <class ValType>
TVector<ValType>& TVector<ValType>::operator=(const TVector<ValType>& v)
{
	if (this != &v)
	{
		if (size != v.size)
		{
			delete[]pVector;
			pVector = new ValType[v.size];
		}
		size = v.size;
		for (int i = 0; i < size; i++)
		{
			pVector[i] = v.pVector[i];
		}
	}
	return *this;
}

template <class ValType>
TVector<ValType>& TVector<ValType>::operator=(const ValType& num)
{
	for (int i = 0; i < size; i++)
		pVector[i] = num;
	return *this;
}

template <class ValType>
TVector<ValType>  TVector<ValType>::operator+(const ValType& num)
{
	TVector res(size);
	for (int i = 0; i < size; i++)
	{
		res.pVector[i] = pVector[i] + num;
	}
	return res;

}

template <class ValType>
TVector<ValType>  TVector<ValType>::operator-(const ValType& num)
{
	TVector res(size);
	for (int i = 0; i < size; i++)
	{
		res.pVector[i] = pVector[i] - num;
	}
	return res;
}

template <class ValType>
TVector<ValType>  TVector<ValType>::operator*(const ValType& num)
{
	TVector res(size);
	for (int i = 0; i < size; i++)
	{
		res.pVector[i] = pVector[i] * num;
	}
	return res;
}

template <class ValType>
TVector<ValType>  TVector<ValType>::operator/(const ValType& num)
{
	TVector res(size);
	if (num != 0)
	{
		for (int i = 0; i < size; i++)
		{
			res.pVector[i] = pVector[i] / num;
		}
	}
	else if (num == 0)
	{
		for (int i = 0; i < size; i++)
		{
			res.pVector[i] = 0;
		}
	}
	return res;
}

template <class ValType>
int TVector<ValType>::Max()
{
	int ind = -1;
	ValType max = -1 * (pow(2, sizeof(ValType) * 8.0 - 1) - 1);
	for (int i = 0; i < size; i++)
		if (pVector[i] > max)
		{
			max = pVector[i];
			ind = i;
		}
	return ind;
}

template <class ValType>
int TVector<ValType>::Min()
{
	int ind = -1;
	ValType min = (pow(2, sizeof(ValType) * 8.0 - 1) - 1);
	for (int i = 0; i < size; i++)
		if (pVector[i] < min)
		{
			min = pVector[i];
			ind = i;
		}
	return ind;
}

template <class ValType>
TVector<ValType> TVector<ValType>::Absolute(const TVector& v)
{
	TVector res(v.size);
	for (int i = 0; i < v.size; i++)
	{
		res[i] = abs(v[i]);
	}
	return res;
}

template <class ValType>
TVector<ValType>  TVector<ValType>::operator+(const TVector<ValType>& v)
{
	TVector res(size);
	for (int i = 0; i < size; i++)
	{
		res.pVector[i] = pVector[i] + v.pVector[i];
	}
	return res;
}

template <class ValType>
TVector<ValType>  TVector<ValType>::operator-(const TVector<ValType>& v)
{
	TVector res(size);
	for (int i = 0; i < size; i++)
	{
		res.pVector[i] = pVector[i] - v.pVector[i];
	}
	return res;
}

template <class ValType>
ValType  TVector<ValType>::operator*(const TVector<ValType>& v)
{
	ValType res = 0;
	for (int i = 0; i < size; i++)
	{
		res += pVector[i] * v.pVector[i];
	}
	return res;
}

//----------------------------------------------//

template <class ValType>
class TMatrix : public TVector< TVector<ValType> >
{
public:
	TMatrix(int s = 10);
	TMatrix(const TMatrix& mt);
	TMatrix(const TVector<TVector<ValType> >& mt);
	bool operator==(const TMatrix& mt) const;
	bool operator!=(const TMatrix& mt) const;
	TMatrix& operator= (const TMatrix& mt);
	TMatrix  operator+ (const TMatrix& mt);
	TMatrix  operator- (const TMatrix& mt);
	TMatrix  operator* (const TMatrix& mt);
	TMatrix operator* (const ValType& num);
	TMatrix operator+ (const ValType& num);
	TMatrix operator- (const ValType& num);
	TMatrix operator/ (const ValType& num);
	TVector<ValType>  operator* (const TVector<ValType>& v);

	TMatrix POW(const int m);

	void SingleMatrix();
	TMatrix ReverseMatrix();
	void Transponir();

	friend istream& operator >> (istream& is, TMatrix& mt)
	{
		for (int i = 0; i < mt.size; i++)
			is >> mt.pVector[i];
		return is;
	}
	friend ostream& operator << (ostream& os, const TMatrix& mt)
	{
		for (int i = 0; i < mt.size; i++)
			os << mt.pVector[i];
		return os;
	}
};

template <class ValType>
TMatrix<ValType>::TMatrix(int s) : TVector< TVector<ValType> >(s)
{
	for (int i = 0; i < s; i++)
	{
		TMatrix::pVector[i] = TVector<ValType>(s);
	}
}

template <class ValType>
TMatrix<ValType>::TMatrix(const TMatrix<ValType>& mt) :TVector<TVector<ValType> >(mt)
{

}

template <class ValType>
TMatrix<ValType>::TMatrix(const TVector<TVector<ValType> >& mt) :TVector<TVector<ValType> >(mt)
{

}

template <class ValType>
bool TMatrix<ValType>::operator==(const TMatrix<ValType>& mt) const
{
	return TVector<TVector<ValType> >::operator==(mt);
}

template <class ValType>
bool TMatrix<ValType>::operator!=(const TMatrix<ValType>& mt) const
{
	return TVector<TVector<ValType> >::operator!=(mt);
}

//----------------------//

template <class ValType>
TMatrix<ValType> TMatrix<ValType>::operator* (const ValType& num)
{
	TMatrix res(TMatrix::size);
	for (int i = 0; i < TMatrix::size; i++)
		res.pVector[i] = TMatrix::pVector[i] * num;
	return res;
}

template <class ValType>
TMatrix<ValType> TMatrix<ValType>::operator+ (const ValType& num)
{
	return TVector<TVector<ValType> >::operator+(num);
}

template <class ValType>
TMatrix<ValType> TMatrix<ValType>::operator- (const ValType& num)
{
	return TVector<TVector<ValType> >::operator-(num);
}

template <class ValType>
TMatrix<ValType> TMatrix<ValType>::operator/ (const ValType& num)
{
	TMatrix res(TMatrix::size);
	for (int i = 0; i < TMatrix::size; i++)
		res.pVector[i] = TMatrix::pVector[i] / num;
	return res;
}

//-----------------//

template <class ValType>
TMatrix<ValType>& TMatrix<ValType>::operator= (const TMatrix<ValType>& mt)
{
	if (this != &mt)
	{
		if (TMatrix::size != mt.size)
		{
			delete[]TMatrix::pVector;
			TMatrix::pVector = new TVector<ValType>[mt.size];
		}
		TMatrix::size = mt.size;
		for (int i = 0; i < TMatrix::size; i++)
		{
			TMatrix::pVector[i] = mt.pVector[i];
		}
	}
	return *this;
}

template <class ValType>
TMatrix<ValType>  TMatrix<ValType>::operator+ (const TMatrix<ValType>& mt)
{
	return TVector<TVector<ValType> >::operator+(mt);
}

template <class ValType>
TMatrix<ValType>  TMatrix<ValType>::operator- (const TMatrix<ValType>& mt)
{
	return TVector<TVector<ValType> >::operator-(mt);
}

template <class ValType>
TMatrix<ValType>  TMatrix<ValType>::operator* (const TMatrix<ValType>& mt)
{
	int count = 0;
	TMatrix res(TMatrix::size);
	for (int i = 0; i < TMatrix::size; i++)
	{
		for (int j = 0; j < TMatrix::size; j++)
		{
			for (int t = 0; t < TMatrix::size; t++)
			{
				res[i][j] += TMatrix::pVector[i][t] * mt.pVector[t][j];
			}
		}
	}
	return res;
}

//------------------//

template <class ValType>
TVector<ValType>  TMatrix<ValType>::operator* (const TVector<ValType>& v)
{
	TVector<ValType> res(TMatrix::size);
	for (int i = 0; i < TMatrix::size; i++)
	{
		res[i] = TMatrix::pVector[i] * v;
	}
	return res;
}

template <class ValType>
TMatrix<ValType> TMatrix<ValType>::POW(const int m)
{
	TMatrix<ValType> res(TMatrix::size);
	res.SingleMatrix();

	if (m >= 0)
		for (int i = 0; i < m; i++)
			res = (*this) * res;
	return res;

}

template <class ValType>
void TMatrix<ValType>::SingleMatrix()
{
	for (int i = 0; i < TMatrix::size; i++)
		for (int j = 0; j < TMatrix::size; j++)
			TMatrix::pVector[i][j] = 0;
	for (int i = 0; i < TMatrix::size; i++)
	{
		TMatrix::pVector[i][i] = 1;
	}
}

template <class ValType>
void TMatrix<ValType>::Transponir()
{
	TMatrix<ValType> sup = *this;
	
	for (int i = 0; i < sup.size; i++)
		for (int j = 0; j < sup.size; j++)
			TMatrix::pVector[i][j] = sup.pVector[j][i];

}

template <class ValType>
TMatrix<ValType> TMatrix<ValType>::ReverseMatrix()
{
	TMatrix res(TMatrix::size), matr = *this;
	res.SingleMatrix();
	ValType sup;
	for (int i = 1; i < matr.size; i++)
	{
		if (fabs(matr.pVector[i - 1][i - 1]) <= pow(10, -15))
			for (int j = i; j < matr.size; j++)
			{
				matr[i - 1] = matr[i - 1] + matr[j];
				res[i - 1] = res[i - 1] + res[j];
			}
		for (int j = i; j < matr.size; j++)
		{
			sup = matr.pVector[j][i - 1] / matr.pVector[i - 1][i - 1];
			matr[j] = matr[j] - matr[i - 1] * sup;
			res[j] = res[j] - res[i - 1] * sup;
		}
	}
	for (int i = matr.size - 2; i >= 0; i--)
	{
		for (int j = i; j >= 0; j--)
		{
			sup = matr.pVector[j][i + 1] / matr.pVector[i + 1][i + 1];
			matr[j] = matr[j] - matr[i + 1] * sup;
			res[j] = res[j] - res[i + 1] * sup;
		}
	}
	for (int i = 0; i < matr.size; i++)
	{
		res[i] = res[i] / matr[i][i];
		matr[i] = matr[i] / matr[i][i];
	}
	return res;
}

#endif
