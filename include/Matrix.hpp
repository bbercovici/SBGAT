#ifndef MATRIXHEADERDEF
#define MATRIXHEADERDEF
#include "Vect.hpp"

class Matrix
{
private:
  double** mData; //matrix entries
  int mNoRows, mNoCols; //dimensions

public:
  //Constructors & Destructor
  Matrix(const Matrix& otherMatrix);
  Matrix(int NoRows, int NoCols);
  ~Matrix();

  //Functions operating on private members
  int GetNoRows() const;
  int GetNoCols() const;
  double CalcDet() const;
 
  //Overload operators
  double& operator()(int i, int j); //operator 0-based indexing
  Matrix& operator=(const Matrix& otherMatrix);
  Matrix operator+() const; // unary +
  Matrix operator-() const; // unary -
  Matrix operator+(const Matrix& m1) const; // binary +
  Matrix operator-(const Matrix& m1) const; // binary -
  Matrix operator*(double a) const; //product w scalar
  Matrix operator/(double a) const; //division per scalar


  //Friend functions
  friend Vect operator*(const Matrix& m, const Vect& v);
  friend Vect operator*(const Vect&v, const Matrix& m);
  friend Matrix operator*(double a, const Matrix& m);
  friend double det(const Matrix& m);
  friend void disp(const Matrix& m);
  friend Matrix inv(const Matrix& m);
  friend Matrix reshape(const Matrix& m, int NoRows, int NoCols);
  friend Matrix reshape(const Vect& m, int NoRows, int NoCols);
  friend Matrix transpose(const Matrix& m);
};

//prototype signatures for friend functions
Vect operator*(const Matrix& m, const Vect& v);
Vect operator*(const Vect& v, const Matrix& m);
Matrix operator*(double a, const Matrix& m);
double det(const Matrix& m);
void disp(const Matrix& m);
Matrix inv(const Matrix& m);
Matrix reshape(const Matrix& m, int NoRows, int NoCols);
Matrix reshape(const Vect& m, int NoRows, int NoCols);
Matrix transpose(const Matrix& m);

#endif
