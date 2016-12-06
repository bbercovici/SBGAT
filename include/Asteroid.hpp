
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

  double get_spin_rate() const;
  arma::vec get_spin_axis() const;

  void set_spin_rate(double spin_rate);
  void set_spin_axis(arma::vec & spin_axis);
  void set_density(double density);


  void setmGs(double mGs);

  friend Vect PolyGrav(Vect& Xsc, Asteroid& Body, bool zero_indexed);
  Vect PolyGrav(Vect& Xsc, bool zero_indexed);
  void compute_global_pgm();

  /**
  Saves asteroid connectivity (vertices + facets) to an obj file
  @param filename Path to file
  */
  void write_to_obj(std::string filename);

  /**
  Saves the current surface acceleration (from the polyhedron gravity model)
  to file
  @param filename Path to file
  */
  void write_surface_acceleration(std::string filename);

  /**
  Loads a previously computed surface acceleration (from the polyhedron gravity model)
  @param filename Path to file
  @return - 1 if PGM gravity file was consistent with asteroid shape model
          - 0 otherwise (no surface acceleration effectively loaded)
  */
  int load_surface_acceleration(std::string filename);


protected:
  arma::vec spin_axis;
  double spin_rate;



};


#endif
