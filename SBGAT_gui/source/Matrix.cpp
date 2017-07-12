#include <iostream>
#include <cmath>
#include <cassert>
#include "Vect.hpp"
#include "Matrix.hpp"

//Overwritten copy constructor
Matrix::Matrix(const Matrix& otherMatrix)
{

  mNoRows = otherMatrix.mNoRows;
  mNoCols = otherMatrix.mNoCols;
  mData = new double* [mNoRows];
  for(int i = 0; i < mNoRows; i++)
    {
      mData[i] = new double [mNoCols];
    }
  for(int i = 0; i < mNoRows; i++)
    {
      for(int j = 0; j < mNoCols; j++)
	{
	  mData[i][j] = otherMatrix.mData[i][j];
	}
    }
}


//Constructor for vectors of given length
//Allocates memory and initializes entries to zero
Matrix::Matrix(int NoRows, int NoCols)
{
  assert(NoRows>0);
  assert(NoCols>0);
  mNoRows = NoRows;
  mNoCols = NoCols;
  mData = new double* [mNoRows];
  for(int i = 0; i < NoRows; i++)
    {
      mData[i] = new double [mNoCols];
    }
  for(int i = 0; i < mNoRows; i++)
    {
      for(int j = 0; j < mNoCols; j++)
      {
    	  mData[i][j] = 0.0;
      }
    }
}


//Overwrite destructor to free memory
Matrix::~Matrix()
{
  for(int i = 0; i < mNoRows; i++)
    {
      delete[] mData[i];
    }
  delete[] mData;
}



//Get the No of Rows
int Matrix::GetNoRows() const
{
  return mNoRows;
}


//Get the No of Cols
int Matrix::GetNoCols() const
{
  return mNoCols;
}


//Calculate the Determinant
double Matrix::CalcDet() const
{
  assert(mNoRows == mNoCols);
  double det = 0.0;

  if(mNoRows == 1)
    det = mData[0][0];
  else
    {
      //For more than one entry
      for(int i_out = 0; i_out < mNoRows; i_out++)
	{
	  Matrix sub_matrix(mNoRows-1,mNoCols-1);
	  for(int i = 0; i < mNoRows-1; i++)
	    {
	      for(int j = 0; j < i_out; j++)
		{
		  sub_matrix(i,j) = mData[i+1][j];
		}
	      for(int j = i_out; j < mNoRows-1; j++)
		{
		  sub_matrix(i,j) = mData[i+1][j+1];
		}
	    }
	  double sub_matrix_det = sub_matrix.CalcDet();

	  det += pow(-1.0, i_out)*mData[0][i_out]*sub_matrix_det;
	}
    }
  return det;
}



//Overload Operators
//Overloading the square brackets
double& Matrix::operator()(int i, int j)
{
  assert(i > -1);
  assert(i < mNoRows);
  assert(j > -1);
  assert(j < mNoCols);
  return mData[i][j];
}


//Overloading the assignment operator
Matrix& Matrix::operator = (const Matrix& otherMatrix)
{
  assert(mNoRows == otherMatrix.mNoRows);
  assert(mNoCols == otherMatrix.mNoCols);

  for(int i = 0; i < mNoRows; i++)
    {
      for(int j = 0; j < mNoCols; j++)
	{
	  mData[i][j] = otherMatrix.mData[i][j];
	}
    }
  return *this;
}


//Overloading the unary + operator
Matrix Matrix::operator+() const
{
  Matrix mat(mNoRows,mNoCols);
  for(int i = 0; i < mNoRows; i++)
    {
      for(int j = 0; j < mNoCols; j++)
	{
	  mat(i,j) = mData[i][j];
	}
    }
  return mat;
}


//Overloading the unary - operator
Matrix Matrix::operator-() const
{
  Matrix mat(mNoRows,mNoCols);
  for(int i = 0; i < mNoRows; i++)
    {
      for(int j = 0; j < mNoCols; j++)
	{
	  mat(i,j) = -mData[i][j];
	}
    }
  return mat;
}



//Overloading the binary + operator
Matrix Matrix::operator+(const Matrix& m1) const
{
  assert(mNoRows == m1.mNoRows);
  assert(mNoCols == m1.mNoCols);
  Matrix mat(mNoRows, mNoCols);

  for(int i = 0; i < mNoRows; i++)
    {
      for(int j = 0; j < mNoCols; j++)
	{
	  mat(i,j) = mData[i][j] + m1.mData[i][j];
	}
    }
  return mat;
}



//Overloading the binary - operator
Matrix Matrix::operator-(const Matrix& m1) const
{
  assert(mNoRows == m1.mNoRows);
  assert(mNoCols == m1.mNoCols);
  Matrix mat(mNoRows, mNoCols);

  for(int i = 0; i < mNoRows; i++)
    {
      for(int j = 0; j < mNoCols; j++)
	{
	  mat(i,j) = mData[i][j] - m1.mData[i][j];
	}
    }
  return mat;
}



//Overloading scalar multiplication
Matrix Matrix::operator*(double a) const
{

  Matrix mat(mNoRows, mNoCols);
  for(int i = 0; i < mNoRows; i++)
    {
      for(int j = 0; j < mNoCols; j++)
	{
	  mat(i,j) = mData[i][j]*a;
	}
    }
  return mat;
}



//Overloading  division per scalar
Matrix Matrix::operator/(double a) const
{
  assert(a!=0);
  Matrix mat(mNoRows, mNoCols);
  for(int i = 0; i < mNoRows; i++)
    {
      for(int j = 0; j < mNoCols; j++)
	{
	  mat(i,j) = mData[i][j]/a;
	}
    }
  return mat;
}


//Overloading Matrix multiplied by a vector
Vect operator*(const Matrix& m, const Vect& v)
{
  int original_vect_size = v.GetSize();
  assert(m.GetNoCols() == original_vect_size);
  int new_vect_length = m.GetNoRows();
  Vect new_vect(new_vect_length);

  for(int i = 0; i < new_vect_length; i++)
    {
      for(int j = 0; j < original_vect_size; j++)
	{
	  new_vect[i] += m.mData[i][j]*v.Read(j);
	}
    }
  return new_vect;
}


//Overload Vect multiplied by matrix
Vect operator*(const Vect& v, const Matrix& m)
{
  int original_vect_size = v.GetSize();
  assert(m.GetNoRows()==original_vect_size);
  int new_vect_length = m.GetNoCols();
  Vect new_vect(new_vect_length);

  for(int i = 0; i < new_vect_length; i++)
    {
      for(int j = 0; j < original_vect_size; j++)
	{
	  new_vect[i] += v.Read(j)*m.mData[j][i];
	}
    }
  return new_vect;
}


// Scalar x Matrix
Matrix operator*(double a, const Matrix& m)
{
  Matrix mat(m.mNoRows, m.mNoCols);
  for(int i = 0; i < m.mNoRows; i++)
    {
      for(int j = 0; j < m.mNoCols; j++)
	{
	  mat(i,j) = a*m.mData[i][j];
	}
    }
  return mat;
}




//MATLAB style det function
double det(const Matrix& m)
{
  return m.CalcDet();
}


//MATLAB style disp function
void disp(const Matrix& m)
{
  for(int i = 0; i < m.mNoRows; i++)
    {
      std::cout << "[ ";

      for(int j = 0; j < m.mNoCols; j++)
	{
	  std::cout << m.mData[i][j] << " ";
	}
      std::cout << "]" << std::endl;
    }
}



//MATLAB style inv function
Matrix inv(const Matrix& m)
{
  assert(m.mNoRows == m.mNoCols);
  assert(m.mNoRows < 4);

  Matrix m_inv(m.mNoRows,m.mNoCols);

  switch(m.mNoRows)
    {
    case 2:
      m_inv(0,0) = m.mData[1][1];
      m_inv(0,1) = -m.mData[0][1];
      m_inv(1,0) = -m.mData[1][0];
      m_inv(1,1) = m.mData[0][0];

      m_inv = m_inv/det(m);
      break;

    case 3:

      m_inv(0,0) = m.mData[1][1]*m.mData[2][2] - m.mData[2][1]*m.mData[1][2];
      m_inv(0,1) = m.mData[0][2]*m.mData[2][1] - m.mData[2][2]*m.mData[0][1];
      m_inv(0,2) = m.mData[0][1]*m.mData[1][2] - m.mData[1][1]*m.mData[0][2];

      m_inv(1,0) = m.mData[1][2]*m.mData[2][0] - m.mData[2][2]*m.mData[1][0];
      m_inv(1,1) = m.mData[0][0]*m.mData[2][2] - m.mData[2][0]*m.mData[0][2];
      m_inv(1,2) = m.mData[0][2]*m.mData[1][0] - m.mData[1][2]*m.mData[0][0];

      m_inv(2,0) = m.mData[1][0]*m.mData[2][1] - m.mData[2][0]*m.mData[1][1];
      m_inv(2,1) = m.mData[0][1]*m.mData[2][0] - m.mData[2][1]*m.mData[0][0];
      m_inv(2,2) = m.mData[0][0]*m.mData[1][1] - m.mData[1][0]*m.mData[0][1];

      m_inv = m_inv/det(m);
      break;
    }

  return m_inv;

}


//MATLAB style reshape function
Matrix reshape(const Matrix& m, int NoRows, int NoCols)
{
  assert(m.mNoRows*m.mNoCols == NoRows*NoCols);
  int N = m.mNoRows*m.mNoCols;
  Vect rshp(N);

  int k = 0;
  for(int i = 0; i < m.mNoRows; i++)
    {
      for(int j = 0; j < m.mNoCols; j++)
	{
	  rshp[k] = m.mData[i][j];
	  k++;
	}
    }

  Matrix m1(NoRows,NoCols);
  k = 0;
  for(int i = 0; i < NoRows; i++)
    {
      for(int j = 0; j < NoCols; j++)
	{
	  m1(i,j) = rshp.Read(k);
	  k++;
	}
    }

  return m1;
}




Matrix reshape(const Vect& m, int NoRows, int NoCols)
{
  assert(length(m) == NoRows*NoCols);

  Matrix m1(NoRows,NoCols);
  int k = 0;
  for(int i = 0; i < NoRows; i++)
    {
      for(int j = 0; j < NoCols; j++)
	{
	  m1(i,j) = m.Read(k);
	  k++;
	}
    }

  return m1;
}




// MATLAB style transpose function
Matrix transpose(const Matrix&m)
{
  int NoRows, NoCols;
  NoCols = m.GetNoRows();
  NoRows = m.GetNoCols();
  assert(NoRows > 0);
  assert(NoCols > 0);

  Matrix mat(NoCols, NoRows);

  for(int i = 0; i < NoRows; i++)
    {
      for(int j = 0; j < NoCols; j++)
	{
	  mat(i,j) = m.mData[j][i];
	}
    }
  return mat;
}
