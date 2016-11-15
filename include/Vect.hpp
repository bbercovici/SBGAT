#ifndef VECTHEADERDEF
#define VECTHEADERDEF

class Vect
{
private:
  double* mData; // data stored in vector
  int mSize;   // size of vector

public:
  Vect(const Vect& otherVect);
  Vect(int size);
  ~Vect();

  int GetSize() const;
  void Disp() const;
  double GetNorm() const;
  double Read(int i) const;

  double& operator[](int i); // zero-based indexing
  
  Vect& operator=(const Vect& otherVect);
  Vect operator +() const; // unary +
  Vect operator -() const; // unary -
  Vect operator +(const Vect& v1) const; // binary +
  Vect operator -(const Vect& v1) const; // binary -
  Vect operator *(double a) const; // scalar multiplication
  Vect operator /(double a) const; // division for a scalar

  //friend functions
  friend Vect operator*(double a, const Vect& v1); // scalar multiplication
  friend Vect cross(const Vect& v, const Vect& v2);
  friend void disp(const Vect& v);
  friend double dot(const Vect& v, const Vect& v2);
  friend int length(const Vect& v);
  friend double norm(const Vect& v);
  friend int size(const Vect& v);
};

// Prototype signature of friend functions
Vect operator*(double a, const Vect& v1); // scalar multiplication
Vect cross(const Vect& v, const Vect& v2);
void disp(const Vect& v);
double dot(const Vect& v, const Vect& v2);
int length(const Vect& v);
double norm(const Vect& v);
int size(const Vect& v);

#endif
