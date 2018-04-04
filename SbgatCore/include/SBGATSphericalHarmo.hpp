/*=========================================================================

  Program:   Small Body Geophysical Analysis
  Module:    SBGATSphericalHarmo.hpp

  Derived class from VTK's vtkPolyDataAlgorithm by Benjamin Bercovici  

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class  SBGATSphericalHarmo
 * @author Benjamin Bercovici
 * @brief  Computes/evaluates the outer spherical harmonics expansion of the exterior gravity
 field around a constant density polyhedron
 *
 * @details  Computes/evaluates the outer spherical harmonics expansion of the exterior gravity
 * field around a constant density polyhedron. Normalized or non-normalized coefficients can be computed
 * Adapted from the works of Yu Takahashi and Siamak Hesar by Benjamin Bercovici, University of Colorado Boulder
 * for more details, see 
 * Werner, R. a. (1997). 
 * Spherical harmonic coefficients for the potential of a constant-density polyhedron. 
 * Computers & Geosciences, 
 * 23(10), 
 * 1071–1077. 
 * https://doi.org/10.1016/S0098-3004(97)00110-6
*/

#ifndef SBGATSphericalHarmo_h
#define SBGATSphericalHarmo_h

#include <vtkFiltersCoreModule.h> // For export macro
#include <vtkPolyDataAlgorithm.h>
#include <armadillo>

class VTKFILTERSCORE_EXPORT SBGATSphericalHarmo : public vtkPolyDataAlgorithm{
public:
  /**
   * Constructs with initial values of zero.
   */
  static SBGATSphericalHarmo *New();

  vtkTypeMacro(SBGATSphericalHarmo,vtkPolyDataAlgorithm);
  void PrintSelf(std::ostream& os, vtkIndent indent) override;
  void PrintHeader(std::ostream& os, vtkIndent indent) override;
  void PrintTrailer(std::ostream& os, vtkIndent indent) override;

  /**
  Sets degree of spherical harmonics expansion
  @param deg degree of spherical harmonics expansion
  */
  void SetDegree(const unsigned int deg){
    this -> degree = deg;
    this -> degreeSet = true;
  }

  /**
  Sets reference radius in spherical harmonics expansion
  @param ref_radius reference radius in spherical harmonics expansion
  */
  void SetReferenceRadius(const double ref_radius){
    this -> referenceRadius = ref_radius;
    this -> referenceRadiusSet = true;
  }

  /**
  Sets polyhedron density
  @param density bulk density of polyhedron
  */
  void SetDensity(const double density){
    this -> density = density;
    this -> densitySet = true;
  }


  /*
  Will return normalized coefficients next (default is normalized)
  */
  void IsNormalized(){
    this -> normalized = true;
  }

  /*
  Will return non-normalized coefficients next (default is normalized)
  */
  void IsNonNormalized(){
    this -> normalized = false;
  }

  /*
  Returns Cnm arrays of coefficients
  @return Cnm arrays of coefficients
  */
  arma::mat GetCnm() {this -> Update(); return this -> Cnm;}


  /*
  Returns Cnm arrays of coefficients
  @param Cnm arrays of coefficients
  */
  void GetCnm(arma::mat & C_nm) {this -> Update(); C_nm = this -> Cnm;}


  /*
  Returns Snm arrays of coefficients
  @return Snm arrays of coefficients
  */
  arma::mat GetSnm() {this -> Update(); return this -> Snm;}

  /*
  Returns Snm arrays of coefficients
  @param Snm arrays of coefficients
  */
  void GetSnm(arma::mat & S_nm) {this -> Update(); S_nm = this -> Snm;}

  /**
  Returns the acceleration due to gravity at the specified point
  @param pos position at which the acceleration must be evaluated, expressed in the same frame as 
  the one used to build the spherical harmonics expansion
  */
  arma::vec GetAcceleration(const arma::vec & pos);

  /**
  Sets the scale factor to 1, indicative that the polydata has its coordinates expressed in meters
  */
  void SetScaleMeters() { this -> scaleFactor = 1; this -> scaleFactorSet = true;}

  /**
  Sets the scale factor to 1000, indicative that the polydata has its coordinates expressed in kilometers
  */
  void SetScaleKiloMeters() { this -> scaleFactor = 1000; this -> scaleFactorSet = true;}


protected:
  SBGATSphericalHarmo();
  ~SBGATSphericalHarmo() override;

  int RequestData(vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) override;

  arma::mat Cnm;
  arma::mat Snm;

  double referenceRadius;
  double density;
  double totalMass;
  double scaleFactor;

  bool normalized;
  unsigned int degree;

  bool degreeSet;
  bool densitySet;
  bool referenceRadiusSet;
  bool scaleFactorSet;




private:
  SBGATSphericalHarmo(const SBGATSphericalHarmo&) = delete;
  void operator=(const SBGATSphericalHarmo&) = delete;
};

#endif


