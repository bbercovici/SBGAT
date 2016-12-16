#include <iostream>
#include <cmath>
#include <cassert>
#include "Vect.hpp"
using namespace std;


// Overridden copy constructor
Vect::Vect(const Vect& otherVect)
{

  mSize = otherVect.mSize;
  mData = new double [mSize];
  for (int i = 0; i < mSize; i++)
    {
      mData[i] = otherVect.mData[i];
    }

}


// constructor for vector of a given size
// allocate memory and initializes entries to zero
Vect::Vect(int size)
{
  assert(size > 0);
  mSize = size;
  mData = new double [mSize];
  for(int i  = 0; i < mSize; i++)
    {
      mData[i] = 0.0;
    }
}


// destructor, freeing memory
Vect::~Vect()
{
  delete[] mData;
}


// get the size of the vector
int Vect::GetSize() const
{
  return mSize;
}



// Definition of function disp
void Vect::Disp() const
{

  cout << "["; 
  for (int i = 0; i < mSize; i++){
      cout  << " " << mData[i];
  }
  cout << "]" << endl;

}


// Definition of function norm
double Vect::GetNorm() const
{

  double mod = 0.0;

  for(int i = 0; i < mSize; i++){
    mod += mData[i]*mData[i];
  }

  return sqrt(mod);

}



double Vect::Read(int i) const
{
  assert(i > -1);
  assert(i < mSize);
  return mData[i];
}



// Overloading square brackets
double& Vect::operator[] (int i)
{
  assert( i > -1);
  assert(i < mSize);
  return mData[i];
}


//Overloading the assignmenet operator
Vect& Vect::operator=(const Vect& otherVect)
{
  assert(mSize == otherVect.mSize);
  for(int i = 0; i < mSize; i++)
    {
      mData[i] = otherVect.mData[i];
    }
  return *this;
}


//Overloading the unary + operator
Vect Vect::operator+() const
{
  Vect v(mSize);
  for(int i = 0; i < mSize; i++)
    {
      v[i] = mData[i];
    }
  return v;
}


//Overloading the unary - operator
Vect Vect::operator-() const
{
  Vect v(mSize);
  for(int i = 0; i < mSize; i++)
    {
      if(mData[i] == 0.0)
	v[i] = 0;
      else
	v[i] = -mData[i];
    }
  return v;
}


//Overloading the binary + operator
Vect Vect::operator+(const Vect& v1) const
{
  assert(mSize == v1.mSize);
  Vect v(mSize);
  for(int i = 0; i < mSize; i++)
    {
      v[i] = mData[i] + v1.mData[i];
    }
  return v;
}


//Overloading the binary - operator
Vect Vect::operator-(const Vect& v1) const
{
  assert(mSize == v1.mSize);
  Vect v(mSize);
  for(int i = 0; i < mSize; i++)
    {
      v[i] = mData[i] - v1.mData[i];
    }
  return v;
}


//Overloading the scalar product
Vect Vect::operator*(double a) const
{
  Vect v(mSize);
  for(int i = 0; i < mSize; i++)
    {
      v[i] = a*mData[i];
    }
  return v;
}



//Overloading the division for a scalar
Vect Vect::operator/(double a) const
{
  Vect v(mSize);
  for(int i = 0; i < mSize; i++)
    {
      v[i] = mData[i]/a;
    }
  return v;
}



//Overloading the dot product
Vect operator*(double a, const Vect& v1)
{

  Vect v(v1.mSize);
  for(int i = 0; i < v1.mSize; i++)
    {
      v[i] = a*v1.mData[i];
    }
  return v;
}



//MATLAB style friend for cross product
Vect cross(const Vect& v, const Vect& v2)
{
  assert(v.mSize == 3);
  assert(v2.mSize == 3);
  
  Vect cross(3);

  cross[0] = v.mData[1]*v2.mData[2] - v.mData[2]*v2.mData[1];
  cross[1] = v.mData[2]*v2.mData[0] - v.mData[0]*v2.mData[2];
  cross[2] = v.mData[0]*v2.mData[1] - v.mData[1]*v2.mData[0];

  return cross;
}


// MATLAB style friend to get the size of a vector
void disp(const Vect& v)
{
  return v.Disp();
}


//MATLAB style friend for dot product
double dot(const Vect& v, const Vect& v2)
{
  assert(v.mSize == v2.mSize);
  double prod = 0.0;

  for(int i = 0; i < v.mSize; i++)
    {
      prod += v.mData[i]*v2.mData[i];
    }

  return prod;
}



// MATLAB style friend to get the size of a vector
int length(const Vect& v)
{
  return v.mSize;
}

// MATLAB style friend to get the norm of a vector
double norm(const Vect& v)
{
  return v.GetNorm();
}

// MATLAB style friend to get the size of a vector
int size(const Vect& v)
{
  return v.mSize;
}
