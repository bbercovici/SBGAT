#ifndef ASTEROIDHEADERDEF
#define ASTEROIDHEADERDEF

#include <fstream>
#include "Matrix.hpp"
#include "Vect.hpp"

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

class Asteroid
{
private:

  double mGs;   // G*dens product

  int mNOV; //No. of vertices
  int mNOF; //No. of facets
  int mNOE; //No. of edges

  // Vertex coordinates
  double* mX;
  double* mY;
  double* mZ;


  // Triangle Vertex List
  int** mListTri;

  // Face Normal Unit Vectors
  double** mListN;

  // Face Dyads
  double** mF;

  // Edge Vertex List
  int** mListE;

  // Edge Dyads
  double** mE;


public:

  Asteroid(std::ifstream& GravityFile);
  Asteroid(vtkSmartPointer<vtkPolyData> input_polydata, double Gs);

  ~Asteroid();

  double GetGs() const;

  int GetNOV() const;
  int GetNOF() const;
  int GetNOE() const;

  Vect GetX() const;
  Vect GetY() const;
  Vect GetZ() const;

  Vect GetListTri() const;
  Vect GetListN() const;
  Vect GetF() const;

  Vect GetListE() const;
  Vect GetE() const;

  friend Vect PolyGrav(Vect& Xsc, Asteroid& Body);

};


#endif
