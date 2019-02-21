
/*=========================================================================

  Program:   Small Body Geophysical Analysis
  Module:    SBGATSphericalHarmo.hpp

  Class derived from VTK's vtkPolyDataAlgorithm by Benjamin Bercovici  

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 @class  SBGATSphericalHarmo
 @author Benjamin Bercovici
 @author Jay McMahon
@date October 2018
 @brief  Computes/evaluates the outer spherical harmonics expansion of the exterior gravity
 field around a constant density polyhedron

 @details  Computes/evaluates the outer spherical harmonics expansion of the exterior gravity
field around a constant density polyhedron. Normalized or non-normalized coefficients can be computed.
The computed coefficients are completely independent of the mass and density of the considered object
as they are only a geometric construct, thanks to the constant-density assumption. This class will always use results expressed in `meters` as their distance unit (e.g accelerations in m/s^2, potentials in m^2/s^2,...) . Unit consistency is enforced through the use of the SetScaleMeters()
and SetScaleKiloMeters() method. 

Adapted from the works of Yu Takahashi and Siamak Hesar by Benjamin Bercovici, University of Colorado Boulder
for more details, see 
Werner, R. a. (1997). 
Spherical harmonic coefficients for the potential of a constant-density polyhedron. 
Computers & Geosciences, 
23(10), 
1071–1077. 
https://doi.org/10.1016/S0098-3004(97)00110-6
@copyright MIT License, Benjamin Bercovici and Jay McMahon
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
  Sets reference radius in spherical harmonics expansion. Must be consistent 
  with the units in which the shape coordinates are expressed
  @param ref_radius reference radius in spherical harmonics expansion
  */
  void SetReferenceRadius(const double ref_radius){
    this -> referenceRadius = ref_radius;
    this -> referenceRadiusSet = true;
  }

  /**
  Sets polyhedron density 
  @param density bulk density of polyhedron (kg/m^3)
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
  Returns Cnm array of coefficients ordered like in the exemple below where degree = 5

    1     0     0     0    0      0
  C_10  C_11  C_12  C_13  C_14  C_15
  C_20  C_21  C_22  C_23  C_24  C_25 
  C_30  C_31  C_32  C_33  C_34  C_35 
  C_40  C_41  C_42  C_43  C_44  C_45
  C_50  C_51  C_52  C_53  C_54  C_55


  @return Cnm array of coefficients
  */
  arma::mat GetCnm() {this -> Update(); return this -> Cnm;}


  /*
  Returns Cnm array of coefficients ordered like in the exemple below where degree = 5

    1     0     0     0    0      0
  C_10  C_11  C_12  C_13  C_14  C_15
  C_20  C_21  C_22  C_23  C_24  C_25 
  C_30  C_31  C_32  C_33  C_34  C_35 
  C_40  C_41  C_42  C_43  C_44  C_45
  C_50  C_51  C_52  C_53  C_54  C_55


  @param[out] Cnm array of coefficients
  */
  void GetCnm(arma::mat & C_nm) {this -> Update(); C_nm = this -> Cnm;}


  /*
  
  Returns Snm array of coefficients ordered like in the exemple below where degree = 5

    0     0     0     0    0      0
    0  S_11  S_12  S_13  S_14  S_15
    0  S_21  S_22  S_23  S_24  S_25 
    0  S_31  S_32  S_33  S_34  S_35 
    0  S_41  S_42  S_43  S_44  S_45
    0  S_51  S_52  S_53  S_54  S_55


  @return Snm array of coefficients
  */
  arma::mat GetSnm() {this -> Update(); return this -> Snm;}

  /*
  Returns Snm array of coefficients ordered like in the exemple below where degree = 5

    0     0     0     0    0      0
    0  S_11  S_12  S_13  S_14  S_15
    0  S_21  S_22  S_23  S_24  S_25 
    0  S_31  S_32  S_33  S_34  S_35 
    0  S_41  S_42  S_43  S_44  S_45
    0  S_51  S_52  S_53  S_54  S_55

  @param[out] Snm array of coefficients
  */
  void GetSnm(arma::mat & S_nm) {this -> Update(); S_nm = this -> Snm;}

  /**
  Returns the acceleration due to gravity at the specified point
  @param pos position at which the acceleration must be evaluated (meters)
  @return acceleration (m / s ^ 2)
  */
  arma::vec::fixed<3> GetAcceleration(const arma::vec::fixed<3> & pos);


  /** 
  Evaluates the gravity gradient matrix (the partial derivative of the spherical 
  harmonics acceleration with respect to the position vector) at the prescribed
  location. Note that this gravity gradient matrix is expressed in the body-fixed frame of the considered object.
  @param[in] pos position at which the gravity gradient matrix must be evaluated, expressed in the same frame/same unit L as 
  the shape used to build the spherical harmonics expansion.
  @param[out] dAccdPos container holding the gravity gradient matrix (1 / s ^ 2)
  */
  void GetGravityGradientMatrix(const arma::vec::fixed<3> & pos,
    arma::mat::fixed<3,3> & dAccdPos);


  /** 
  Evaluates the partial derivative of the spherical harmonics acceleration
  with respect to the spherical harmonics coefficients. The coefficients and the partials of the acceleration w/r to the coefficients
  are ordered like so:
  
  Cnm :
(C_00 == 1)  0     0    0     0     0
  C_10     C_11    0    0     0     0 
  C_20     C_21  C_22   0     0     0  
  C_30     C_31  C_32  C_33   0     0  
  C_40     C_41  C_42  C_43  C_44   0
  C_50     C_51  C_52  C_53  C_54  C_55
  
  Snm :
  0    0     0     0    0      0
  0  S_11    0     0    0      0
  0  S_21  S_22    0    0      0 
  0  S_31  S_32  S_33   0      0 
  0  S_41  S_42  S_43  S_44    0
  0  S_51  S_52  S_53  S_54  S_55

  so 

  partial_C = [ dA/C_00,dA/C_10,dA/C_11,dA/C_20,dA/C_21,dA/C_22,...]
  partial_S = [ dA/S_11,dA/C_21,dA/S_22,dA/S_31,dA/S_32,dA/S_33,...]

  @param[in] pos position at which the partial derivatives must be evaluated, expressed in the same frame/same unit as 
  the shape used to build the spherical harmonics expansion.
  @param[out] partial_C container holding the partial derivative of the acceleration 
  with respect to the Cnm spherical harmonic coefficients. If the degree/order is n, then there are (n+1) * (n + 2)/2 non-zero Cnm coefficients
  
  @param[out] partial_S container holding the partial derivative of the acceleration 
  with respect to the Snm spherical harmonic coefficients. If the degree/order is n, then there are (n+1) * (n + 2)/2 - n non-zero Snm coefficients

  */
  void GetPartialHarmonics(const arma::vec::fixed<3> & pos,
    arma::mat & partial_C, 
    arma::mat & partial_S);
  /**
  Sets the scale factor to 1, indicative that the polydata has its coordinates expressed in meters
  */
  void SetScaleMeters() { this -> scaleFactor = 1; this -> scaleFactorSet = true;}

  /**
  Sets the scale factor to 1000, indicative that the polydata has its coordinates expressed in kilometers
  */
  void SetScaleKiloMeters() { this -> scaleFactor = 1000; this -> scaleFactorSet = true;}


  /**
  Exports the computed spherical harmonics expansion to 
  a JSON file. The saved fields are:
  - facets == number of facets
  - vertices == number of vertices
  - totalMass : {value, unit}
  - density : {value, unit}
  - reference_radius : {value, unit}
  - normalized == true if the coefficients are normalized
  - degree == degree of the spherical expansion
  - Cnm_coefs - vector of coefficients triplets {n,m,Cnm}
  - Snm_coefs - vector of coefficients triplets {n,m,Snm}
  @param path JSON file where the spherical harmonics model will be saved

  */
  void SaveToJson(std::string path) const;

  /**
  Loads a previously computed spherical harmonics expansion
  from a JSON file. Will set the appropriate fields in the SBGATSphericalHarmo object to
  allow calls to other methods. 
  The loadable fields are:
  - facets == number of facets (not needed for evaluation)
  - vertices == number of vertices (not needed for evaluation)
  - totalMass : {value, unit}
  - density : {value, unit}
  - reference_radius : {value, unit}
  - normalized == true if the coefficients are normalized
  - degree == degree of the spherical expansion
  - Cnm_coefs - vector of coefficients triplets {n,m,Cnm}
  - Snm_coefs - vector of coefficients triplets {n,m,Snm}
  @param path JSON file storing the spherical harmonics model
  */
  void LoadFromJson(std::string path);

  /**
  Sets the Cnm coefficients. There is normally no need to use this method outside of Sbgat's tests
  @param[in] Cnm coefficients
  */
  void SetCnm(arma::mat Cnm){this ->Cnm =Cnm;}

  /**
  Sets the Snm coefficients. There is normally no need to use this method outside of Sbgat's tests
  @param[in] Snm coefficients
  */
  void SetSnm(arma::mat Snm){this ->Snm =Snm;}



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
  double scaleFactor = 1;

  bool normalized;
  unsigned int degree;

  int n_facets;
  int n_vertices;

  bool degreeSet;
  bool densitySet;
  bool referenceRadiusSet;
  bool scaleFactorSet;
  bool setFromJSON;

private:
  SBGATSphericalHarmo(const SBGATSphericalHarmo&) = delete;
  void operator=(const SBGATSphericalHarmo&) = delete;
};

#endif


