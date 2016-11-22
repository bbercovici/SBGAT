
// Adapted from the original work
// of Nicola Baresi, CSML, University of Colorado Boulder (2014)


#ifndef ASTEROIDHEADERDEF
#define ASTEROIDHEADERDEF

#include <fstream>
#include <iostream>

#include "Matrix.hpp"
#include "Vect.hpp"

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <cassert>
#include <cmath>
#include <set>
#include <map>
#include <fstream>
#include <armadillo>

class Asteroid {
private:

  double mGs;   // G*dens product

  unsigned int mNOV; //No. of vertices
  unsigned int mNOF; //No. of facets
  unsigned int mNOE; //No. of edges

  // Vertex coordinates
  double* mX;
  double* mY;
  double* mZ;


  // Triangle Vertex List
  int** mListTri;

  // Face Normal Unit Vectors
  double** mListN;

  // Surface Gravity Acceleration Vectors
  double** surface_grav;

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

  double * get_X()  ;
  double * get_Y()  ;
  double * get_Z()  ;

  Vect GetListTri() const;
  Vect GetListN() const;
  Vect GetF() const;

  double ** get_ListN() ;
  int ** get_ListTri() ;

  Vect GetListE() const;
  Vect GetE() const;

  double ** get_surface_grav();

  void setmGs(double mGs);

  friend Vect PolyGrav(Vect& Xsc, Asteroid& Body, bool zero_indexed);
  Vect PolyGrav(Vect& Xsc, bool zero_indexed);
  void compute_global_pgm();
  /**
  Saves asteroid connectivity (vertices + facets) to an obj file
  */
  void write_to_obj(std::string filename);

};


#endif
