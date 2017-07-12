
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
#include <vtkCellData.h>
#include <vtkTriangle.h>

#include <cassert>
#include <cmath>
#include <set>
#include <map>
#include <fstream>
#include <armadillo>
#include <array>

#include <ctime>
#include <cstdio>
#include <chrono>
#include <sys/time.h>

/**
Declaration of the Asteroid Class. Stores
  - the shape model of the asteroid in the form of a vtkSmartPointer<vtkPolydata>
  - the PGM (raw accelerations evaluated at the center of each facet)
  - the principal axes/inertias (TBD)
  - the other physical properties of the asteroid
    - density
    - spin rate
    - spin axis direction
*/

class Asteroid {
private:
// G*dens product
  double mGs;

  // Surface Gravity Acceleration Vectors
  double** surface_grav;





public:

  Asteroid(vtkSmartPointer<vtkPolyData> polydata, double Gs);

  ~Asteroid();

  /**
  Returns the G * density product
  @return G * density product
  */
  double GetGs() const;

  double ** get_surface_grav();

  /**
  Returns the current spin rate of the asteroid
  @return Spin rate (rad/s)
  */
  double get_spin_rate() const;

  /**
  Returns the coordinates of the spin axis expressed in the
  asteroid's bodyframe.
  @return Spin axis unit direction
  */
  arma::vec get_spin_axis() const;

  /**
  Sets the asteroid's spin rate
  @param spin_rate Asteroid's spin rate (rad/s)
  */
  void set_spin_rate(double spin_rate);

  /**
  Sets the asteroid's spin axis coordinates. The provided
  vector is automatically normalized
  @param spin_axis Spin axis direction. Must have a non-zero norm
  */
  void set_spin_axis(arma::vec & spin_axis);

  /**
  Sets the asteroid's constant density
  @param density Density of asteroid (kg/m^3)
  */
  void set_density(double density);

  /**
  Returns the asteroid's density
  @return Asteroid's density (kg/m^3)
  */
  double get_density();

  /**
  Legacy method. Sets the G * density method by effectively changing the
  density of the asteroid
  @param mGs G * density product of the asteroid
  */
  void setmGs(double mGs);


  /**
  Computes the acceleration of gravity created by the asteroid
  at the provided location
  @param Xsc Coordinates at which the pgm acceleration must be evaluated (m)
  */
  Vect PolyGrav(Vect& Xsc);

  /**
  Evaluates the PGM gravity acceleration at the center of each facet on the
  asteroid
  */
  void compute_global_pgm();

  /**
  Returns a pointer to the polydata
  @return Pointer to polydata
  */
  vtkSmartPointer<vtkPolyData> get_polydata();

  /**
  Computes the acceleration of gravity created by the asteroid
  at the provided location
  @param Xsc Barycentric Coordinates at which the pgm acceleration must be evaluated (m)
  @return Acceleration expressed in barycentric coordinates (m/s^2)
  */
  arma::vec polygrav_vtk(arma::vec & Xsc);



  /**
  Saves the current surface acceleration (from the polyhedron gravity model)
  to file
  @param filename Path to file
  */
  void write_surface_acceleration(std::string filename);

  /**
  Loads a previously computed surface acceleration (from the polyhedron gravity model)
  @param filename Path to file
  @return 1 if PGM gravity file was consistent with asteroid shape model, 0 otherwise (no surface acceleration will be loaded)
  */
  int load_surface_acceleration(std::string filename);

  std::vector<std::array<unsigned int, 3> > facet_vertices_ids;

protected:
  arma::vec spin_axis;
  double spin_rate;
  vtkSmartPointer<vtkPolyData> polydata;
  std::vector < std::set< std::set <unsigned int> > > facet_edge_point_ids;
  std::map < std::set <unsigned int> , std::set< unsigned int> > edge_point_facet_ids;
  std::vector<arma::mat> F_vector;
  std::map<std::set<unsigned int>, arma::mat> edge_dyads;
};


#endif
