/*=========================================================================

  Program:   Small Body Geophysical Analysis
  Module:    SBGATPolyhedronGravityModel.hpp

  Class derived from VTK's vtkPolyDataAlgorithm by Benjamin Bercovici  

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/**
@file SBGATPolyhedronGravityModel.hpp
@class  SBGATPolyhedronGravityModel
@author Benjamin Bercovici
@date October 2018

@brief  Evaluation of potential, acceleration caused by a constant-density polyhedron
 @details Computes the potential, acceleration caused by a polyhedron
 of constant density by evaluating the so called Polyhedron Gravity Model as derived by Werner and Scheeres.
The input must be a topologically-closed polyhedron. This class will always use results expressed in `meters` as their distance unit (e.g accelerations in m/s^2, potentials in m^2/s^2,...) . Unit consistency is enforced through the use of the SetScaleMeters()
and SetScaleKiloMeters() method. 
See Werner, R. A., & Scheeres, D. J. (1997). Exterior gravitation of a polyhedron derived and compared with harmonic and mascon gravitation representations of asteroid 4769 Castalia. Celestial Mechanics and Dynamical Astronomy, 65(3), 313–344. https://doi.org/10.1007/BF00053511
for further details. Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
@copyright MIT License, Benjamin Bercovici and Jay McMahon
*/

#ifndef SBGATPolyhedronGravityModel_hpp
#define SBGATPolyhedronGravityModel_hpp

#include <vtkFiltersCoreModule.h> // For export macro
#include <vtkPolyDataAlgorithm.h>
#include <armadillo>
#include <SBGATMassProperties.hpp>


class VTKFILTERSCORE_EXPORT SBGATPolyhedronGravityModel : public SBGATMassProperties{
public:

  /**
   * Constructs with initial values of zero.
   */

  static SBGATPolyhedronGravityModel *New();

  vtkTypeMacro(SBGATPolyhedronGravityModel,vtkPolyDataAlgorithm);
  void PrintSelf(std::ostream& os, vtkIndent indent) override;
  void PrintHeader(std::ostream& os, vtkIndent indent) override;
  void PrintTrailer(std::ostream& os, vtkIndent indent) override;


  /**
  Evaluates the Polyhedron Gravity Model potential at the specified point assuming 
  a constant density
  @param point pointer to coordinates of queried point, expressed in the same frame as
  the polydata
  @return PGM potential evaluated at the queried point (m ^ 2/ s ^2)
  */
  double GetPotential(double const * point) const;

  /**
  Evaluates the Polyhedron Gravity Model potential at the specified point assuming 
  a constant density
  @param point coordinates of queried point, expressed in the same frame as
  the polydata
  @return PGM potential evaluated at the queried point (m ^ 2 / s ^2)
  */
  double GetPotential(const arma::vec::fixed<3> & point) const;

  /**
  Evaluates the Polyhedron Gravity Model potential and acceleration at the specified point assuming 
  a constant density
  @param point coordinates of queried point, expressed in the same frame as
  the polydata used to construct the PGM
  @param[out] potential PGM potential evaluated at the queried point (m ^ 2 / s ^2)
  @param[out] acc PGM acceleration evaluated at the queried point (m / s ^2)
  */
  void GetPotentialAcceleration(double const * point,double & potential, 
    arma::vec::fixed<3> & acc) const;


  /**
  Evaluates the Polyhedron Gravity Model potential and acceleration at the specified point assuming 
  a constant density
  @param point coordinates of queried point, expressed in the same frame as
  the polydata used to construct the PGM
  @param[out] potential PGM potential evaluated at the queried point (m ^ 2 / s ^2)
  @param[out] acc PGM acceleration evaluated at the queried point (m / s ^2)
  */
  void GetPotentialAcceleration(const arma::vec::fixed<3> & point,double & potential, 
    arma::vec::fixed<3> & acc) const;


/**
Evaluates the Polyhedron Gravity Model potential, acceleration and gravity gradient matrix at the specified point assuming 
a constant density
@param point coordinates of queried point, expressed in the same frame as
the polydata used to construct the PGM
@param[out] potential PGM potential evaluated at the queried point (m ^ 2 / s ^2)
@param[out] acc PGM acceleration evaluated at the queried point (m / s ^2)
@param[out] gravity_gradient_mat PGM gravity gradient matrix evaluated at the queried point (1 / s ^2)
*/
  void GetPotentialAccelerationGravityGradient(double const  * point,double & potential, 
    arma::vec::fixed<3> & acc,arma::mat::fixed<3,3> & gravity_gradient_mat) const ;


/**
Get the body-fixed acceleration at the center of the specified facet
@param f facet index 
@return body-fixed acceleration (m/s^2)
*/
  arma::vec::fixed<3> GetBodyFixedAccelerationf(const int & f) const;



/**
Evaluates the Polyhedron Gravity Model potential, acceleration and gravity gradient matrix at the specified point assuming 
a constant density
@param point coordinates of queried point, expressed in the same frame as
the polydata used to construct the PGM
@param[out] potential PGM potential evaluated at the queried point (m ^ 2 / s ^2)
@param[out] acc PGM acceleration evaluated at the queried point (m / s ^2)
@param[out] gravity_gradient_mat PGM gravity gradient matrix evaluated at the queried point (1 / s ^2)
*/
  void GetPotentialAccelerationGravityGradient(const arma::vec::fixed<3> & point,double & potential, 
    arma::vec::fixed<3> & acc,arma::mat::fixed<3,3> & gravity_gradient_mat) const ;



/**
Evaluates the Polyhedron Gravity Model gravity gradient matrix at the specified point assuming 
a constant density
@param point coordinates of queried point, expressed in the same frame as
the polydata used to construct the PGM
@param[out] gravity_gradient_mat PGM gravity gradient matrix evaluated at the queried point (1 / s ^2)
*/
  arma::mat::fixed<3,3> GetGravityGradient(const arma::vec::fixed<3> & point) const ;




  /**
  Evaluates the Polyhedron Gravity Model acceleration at the specified point assuming 
  a constant density
  @param point coordinates of queried point, expressed in the same frame as
  the polydata used to construct the PGM
  @return PGM acceleration evaluated at the queried point (m / s ^2)
  */
  arma::vec::fixed<3> GetAcceleration(const arma::vec::fixed<3> & point) const;

  /**
  Evaluates the Polyhedron Gravity Model acceleration at the specified point assuming 
  a constant density
  @param point pointer to coordinates of queried point, expressed in the same frame as
  the input polydata
  @return PGM acceleration evaluated at the queried point (m / s ^2)
  */
  arma::vec::fixed<3> GetAcceleration(double const * point) const;

  /** 
  Determines whether the provided point lies inside or outside the shape
  @param point coordinates of queried point, expressed in the same frame as
  the polydata
  @param tolerance
  @return true if the polydata contains the point, false otherwise
  */
  bool Contains(double const * point, double tol = 1e-8) const;

  /** 
  Determines whether the provided point lies inside or outside the shape
  @param point coordinates of queried point, expressed in the same frame as
  the polydata
  @param tolerance
  @return true if the polydata contains the point, false otherwise
  */
  bool Contains(const arma::vec::fixed<3> & point,double tol = 1e-8) const;

  /**
  Computes the slope as the center of the designated facet
  @param f facet index
  @return slope (rad)
  */
  double GetSlope(const int & f ) const;


  /**
  Get polyhedron mass. Must have set the density and call
  @return mass (kg)
  */
  double GetMass() const{
    return this -> density * this -> GetVolume() * std::pow(this -> scaleFactor,3);
  }

  /**
  Return the performance factor of the f-th facet at the specified position
  @param pos position of field point
  @param f facet index
  @return omega_f
  */
  double GetPerformanceFactor(const arma::vec::fixed<3> & pos, const int & f) const;

  /**
  Return the performance factor of the f-th facet at the specified position
  @param pos position of field point
  @param f facet index
  @return omega_f
  */
  double GetPerformanceFactor( const double * pos, const int & f) const;

  /**
  Return the wire potential of the e-th edge at the specified position
  @param pos position of field point
  @param e edge index
  @return L_e
  */
  double GetLe(const arma::vec::fixed<3> & pos, const int & e) const;

  /**
  Return the wire potential of the e-th edge at the specified position
  @param pos position of field point
  @param e edge index
  @return L_e
  */
  double GetLe( const double * pos, const int & e) const;



  /**
  Evaluates the Polyhedron Gravity Model at the surface of the specified surface elements in the provided shape
  @param[in] selected_shape shape for which the surface polyhedron gravity model must be computed
  @param[in] queried_elements vector of elements indices where the polyhedron gravity model should be evaluated
  @param[in] is_in_meters true if the shape coordinates were expressed in meters, false if they were expressed in kilometers
  @param[in] density shape bulk density in kg/m^3
  @param[in] omega fixed angular velocity of shape expressed in rad/s
  @param[out] slopes vector storing the gravitational slopes (degrees) evaluated at the center of each queried element 
  @param[out] inertial_potentials vector storing the inertial gravitational potentials (m^2/s^2) evaluated at the center of each queried element 
  @param[out] body_fixed_potentials vector storing the inertial gravitational potentials (m^2/s^2) evaluated at the center of each queried element 
  @param[out] inertial_acc_magnitudes vector storing the inertial gravitational acceleration magnitudes  (m/s^2) evaluated at the center of each queried element 
  @param[out] body_fixed_acc_magnitudes vector storing the body-fixed gravitational acceleration magnitudes (m/s^2) evaluated at the center of each queried element 
  */
  static void ComputeSurfacePGM(
    vtkSmartPointer<vtkPolyData> selected_shape,
    const std::vector<unsigned int> & queried_elements,
    bool is_in_meters,
    double density,
    const arma::vec::fixed<3> & omega,
    std::vector<double> & slopes,
    std::vector<double> & inertial_potentials,
    std::vector<double> & body_fixed_potentials,
    std::vector<double> & inertial_acc_magnitudes,
    std::vector<double> & body_fixed_acc_magnitudes);




  /**
  Saves the provided Polyhedron Gravity Model to a file
  @param[in] selected_shape shape for which the surface polyhedron gravity model must be computed
  @param[in] queried_elements shape indices of elements where the polyhedron gravity model should be evaluated
  @param[in] is_in_meters true if the shape coordinates were expressed in meters, false if they were expressed in kilometers
  @param[in] mass mass of shape model (kg)
  @param[in] omega fixed angular velocity of shape (rad/s)
  @param[in] slopes vector storing the gravitational slopes (degrees) evaluated at the center of each queried element 
  @param[in] inertial_potentials vector storing the inertial gravitational potentials (m^2/s^2) evaluated at the center of each queried element 
  @param[in] body_fixed_potentials vector storing the inertial gravitational potentials (m^2/s^2) evaluated at the center of each queried element 
  @param[in] inertial_acc_magnitudes vector storing the inertial gravitational acceleration magnitudes  (m/s^2) evaluated at the center of each queried element 
  @param[in] body_fixed_acc_magnitudes vector storing the body-fixed gravitational acceleration magnitudes (m/s^2) evaluated at the center of each queried element 
  @param[in] path save path (ex: "pgm_surface.json")
  */
  static void SaveSurfacePGM(vtkSmartPointer<vtkPolyData> selected_shape,
    const std::vector<unsigned int> & queried_elements,
    bool is_in_meters,
    const double & mass,
    const arma::vec::fixed<3> & omega,
    const std::vector<double> & slopes,
    const std::vector<double> & inertial_potentials,
    const std::vector<double> & body_fixed_potentials,
    const std::vector<double> & inertial_acc_magnitudes,
    const std::vector<double> & body_fixed_acc_magnitudes,
    std::string path);

  /**
  Loads a previously computed surface Polyhedron Gravity Model from a file
  @param[out] mass mass of shape model (kg)
  @param[out] omega fixed angular velocity of shape (rad/s)
  @param[out] slopes vector storing the gravitational slopes (degrees) evaluated at the center of each queried element 
  @param[out] inertial_potentials vector storing the inertial gravitational potentials (m^2/s^2) evaluated at the center of each queried element 
  @param[out] body_fixed_potentials vector storing the body-fixed gravitational potentials (m^2/s^2) evaluated at the center of each queried element 
  @param[out] inertial_acc_magnitudes vector storing the inertial gravitational acceleration magnitudes  (m/s^2) evaluated at the center of each queried element 
  @param[out] body_fixed_acc_magnitudes vector storing the body-fixed gravitational acceleration magnitudes (m/s^2) evaluated at the center of each queried element 
  @param[in] path load path (ex: "pgm_surface.json")
  */
  static void LoadSurfacePGM(double & mass,
    arma::vec::fixed<3> & omega,
    std::vector<double> & slopes,
    std::vector<double> & inertial_potentials,
    std::vector<double> & body_fixed_potentials,
    std::vector<double> & inertial_acc_magnitudes,
    std::vector<double> & body_fixed_acc_magnitudes,
    std::string path);


  /**
  Return the Xe^E vector holding the e-th edge dyadic factors (Le,r_ie_0^T,Ee^T)^T
  @param pos field point
  @param e edge index
  @return Xe^E vector holding the e-th edge dyadic factors
  */
  arma::vec::fixed<10> GetXe(const arma::vec::fixed<3> & pos,const int & e) const;


  /**
  Return the Xf^F vector holding the f-th facet dyadic factors (omega_f,r_if_0^T,Ff^T)^T
  @param pos field point
  @param f facet index
  @return Xf^F vector holding the f-th facet dyadic factors
  */
  arma::vec::fixed<10> GetXf(const arma::vec::fixed<3> & pos,const int & f) const;

  /**
  Return the parametrization of the designated edge dyad Ee. This dyad 
  is stored in a flattened double container holding nine values and ordered like so

  E = {
    {E[0], E[1],  E[2]},
    {E[3], E[4],  E[5]},
    {E[6], E[7],  E[8]}
  };
  Due to the symmetric nature of E, its parametrization is 
  {E[0],E[4],E[8],E[1],E[2],E[5]};

  @param e edge index 
  @return Ee dyad parametrization
  */  
  arma::vec::fixed<6> GetEeParam(const int & e) const;


  /**
  Return the parametrization of the designated facet dyad Ff. This dyad 
  is stored in a flattened double container holding nine values and ordered like so

  F = {
    {F[0], F[1],  F[2]},
    {F[3], F[4],  F[5]},
    {F[6], F[7],  F[8]}
  };
  Due to the symmetric nature of E, its parametrization is 
  {F[0],F[4],F[8],F[1],F[2],F[5]};

  @param f facet index 
  @return Ff dyad parametrization
  */  
  arma::vec::fixed<6> GetFfParam(const int & f) const;


  /**
  Return the contribution of a specific edge to the potential at a specified field point
  @param Xe vector of dyadic coefficients for edge e at the prescribed fieldpoint
  @return contribution to the potential of this specific edge at the prescribed fieldpoint
  */
  static double GetUe(const arma::vec::fixed<10> & Xe);


  /**
  Return the contribution of a specific facet to the potential at a specified field point
  @param Xf vector of dyadic coefficients for facet f at the prescribed fieldpoint
  @return contribution to the potential of this specific facet at the prescribed fieldpoint
  */
  static double GetUf(const arma::vec::fixed<10> & Xf);


  /**
  Return the contribution of a specific edge to the acceleration at a specified field point
  @param Xe vector of dyadic coefficients for edge e at the prescribed fieldpoint
  @return contribution to the potential of this specific edge at the prescribed fieldpoint
  */
  static arma::vec::fixed<3> GetAe(const arma::vec::fixed<10> & Xe);


  /**
  Return the contribution of a specific facet to the acceleration at a specified field point
  @param Xf vector of dyadic coefficients for facet f at the prescribed fieldpoint
  @return contribution to the potential of this specific facet at the prescribed fieldpoint
  */
  static arma::vec::fixed<3> GetAf(const arma::vec::fixed<10> & Xf);

  /**
  Set the angular velocity vector
  @param Omega angular velocity expressed in body-frame (rad/s)
  */
  void SetOmega(arma::vec::fixed<3> Omega){this -> Omega = Omega;}

  /**
  Return the angular velocity vector
  @return Omega angular velocity expressed in body-frame (rad/s)
  */
  arma::vec::fixed<3> GetOmega() const {return this -> Omega;}

protected:
  SBGATPolyhedronGravityModel();
  ~SBGATPolyhedronGravityModel() override;

  void Clear();

  int RequestData(vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) override;

  double ** facet_dyads;
  double ** edge_dyads;

  arma::vec::fixed<3> Omega;


private:
  SBGATPolyhedronGravityModel(const SBGATPolyhedronGravityModel&) = delete;
  void operator=(const SBGATPolyhedronGravityModel&) = delete;
};

#endif


