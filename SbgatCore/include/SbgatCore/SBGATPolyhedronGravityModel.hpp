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
#include "SBGATMassProperties.hpp"


class VTKFILTERSCORE_EXPORT SBGATPolyhedronGravityModel : public vtkPolyDataAlgorithm{
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
  Returns coordinates of point at center of facet
  @param f facet index
  @return coordinates of f-th facet center
  */
  arma::vec::fixed<3> GetFacetCenter(const int & f) const;


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
@param Omega angular velocity vector expressed in the body-frame
@return body-fixed acceleration (m/s^2)
*/
  arma::vec::fixed<3> GetBodyFixedAccelerationf(const int & f,const arma::vec::fixed<3> & Omega) const;



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
  Sets the scale factor to 1, indicative that the polydata has its coordinates expressed in meters
  */
  void SetScaleMeters() { this -> scaleFactor = 1; this -> scaleFactorSet = true;this -> is_in_meters = true;}

  /**
  Sets the scale factor to 1000, indicative that the polydata has its coordinates expressed in kilometers
  */
  void SetScaleKiloMeters() { this -> scaleFactor = 1000; this -> scaleFactorSet = true;this -> is_in_meters = false;}


  /**
  Computes the slope as the center of the designated facet
  @param f facet index
  @param Omega angular velocity vector (rad/s)
  @return slope (rad)
  */
  double GetSlope(const int & f ,const arma::vec::fixed<3> & Omega) const;


  /**
  Returns the shape's scale factor
  @return scale factor
  */
  double GetScaleFactor() const {return this -> scaleFactor;}

  /**
  Sets polyhedron density
  @param density bulk density of polyhedron (kg/m^3)
  */
  void SetDensity(const double density){
    this -> density = density;
    this -> densitySet = true;
  }

  /**
  Returns length of considered edge in meters
  @param e edge index
  @return edge length(m )
*/  
  double GetEdgeLength(const int & e) const;

  /**
  Get polyhedron density
  @return density bulk density of polyhedron (kg/m^3)
  */
  double GetDensity() const{
    return this -> density ;
  }

  /**
  Get polyhedron mass. Must have set the density and call
  @return mass (kg)
  */
  double GetMass() const{
    return this -> density * this -> mass_properties -> GetVolume() * std::pow(this -> scaleFactor,3);
  }


  /**
  Gets pointer to associated mass properties
  @return mass properties
  */
  vtkSmartPointer<SBGATMassProperties> GetMassProperties() const{return this -> mass_properties;}

  /**
  Returns the performance factor of the f-th facet at the specified position
  @param pos position of field point
  @param f facet index
  @return omega_f
  */
  double GetOmegaf(const arma::vec::fixed<3> & pos, const int & f) const;

  /**
  Returns the performance factor of the f-th facet at the specified position
  @param pos position of field point
  @param f facet index
  @return omega_f
  */
  double GetOmegaf( const double * pos, const int & f) const;

  /**
  Returns the wire potential of the e-th edge at the specified position
  @param pos position of field point
  @param e edge index
  @return L_e
  */
  double GetLe(const arma::vec::fixed<3> & pos, const int & e) const;

  /**
  Returns the wire potential of the e-th edge at the specified position
  @param pos position of field point
  @param e edge index
  @return L_e
  */
  double GetLe( const double * pos, const int & e) const;


  /**
  Returns coordinates of the three vertices forming a facet
  @param[in] f facet index
  @param[out] r0 first vertex coordinates
  @param[out] r1 first vertex coordinates
  @param[out] r2 first vertex coordinates
  */
  void GetVerticesInFacet(const int & f,double * r0,double * r1, double * r2) const;



  /**
  Returns coordinates of the two vertices forming an edge
  @param[in] e edge index
  @param[out] r0 first vertex coordinates
  @param[out] r1 first vertex coordinates
  */
  void GetVerticesOnEdge(const int & e,double * r0,double * r1) const;


  /**
  Returns normal of facet f
  @param[in] f facet index
  @param[out] n normal at facet
  */
  void GetFacetNormal(const int & f, double * n) const{n[0] = this -> facet_normals[f][0];n[1] = this -> facet_normals[f][1];n[2] = this -> facet_normals[f][2];}


  /**
  Returns the indices of the two facets adjacent to the specified edge
  @param[in] e edge index
  @param[out] f0 index of first facet
  @param[out] f1 index of second facet
  */
  void GetIndicesOfAdjacentFacets(const int & e,int & f0, int & f1) const;


  /**
  Returns the vector from the field point to the first vertex on the designated edge
  @param pos field point
  @param e edge index
  @return the vector from the field point to the first vertex on the designated edge
  */
  arma::vec::fixed<3> GetRe(const arma::vec::fixed<3> & pos,const int & e) const;


  /**
  Returns the vector from the field point to the first vertex on the designated edge
  @param pos field point
  @param e edge index
  @return the vector from the field point to the first vertex on the designated edge
  */
  arma::vec::fixed<3> GetRe(const double * pos,const int & e) const;



  /**
  Returns the vector from the field point to the first vertex on the designated facet
  @param pos field point
  @param f facet index
  @return the vector from the field point to the first vertex on the designated facet
  */
  arma::vec::fixed<3> GetRf(const arma::vec::fixed<3> & pos,const int & f) const;


  /**
  Returns the vector from the field point to the first vertex on the designated facet
  @param pos field point
  @param f facet index
  @return the vector from the field point to the first vertex on the designated facet
  */
  arma::vec::fixed<3> GetRf(const double * pos,const int & f) const;


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
  Returns the Xe^E vector holding the e-th edge dyadic factors (Le,r_ie_0^T,Ee^T)^T
  @param pos field point
  @param e edge index
  @return Xe^E vector holding the e-th edge dyadic factors
  */
  arma::vec::fixed<10> GetXe(const arma::vec::fixed<3> & pos,const int & e) const;


  /**
  Returns the Xf^F vector holding the f-th facet dyadic factors (omega_f,r_if_0^T,Ff^T)^T
  @param pos field point
  @param f facet index
  @return Xf^F vector holding the f-th facet dyadic factors
  */
  arma::vec::fixed<10> GetXf(const arma::vec::fixed<3> & pos,const int & f) const;

  /**
  Returns the parametrization of the designated edge dyad Ee. This dyad 
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
  Returns the parametrization of the designated facet dyad Ff. This dyad 
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
  Returns the contribution of a specific edge to the potential at a specified field point
  @param Xe vector of dyadic coefficients for edge e at the prescribed fieldpoint
  @return contribution to the potential of this specific edge at the prescribed fieldpoint
  */
  static double GetUe(const arma::vec::fixed<10> & Xe);


  /**
  Returns the contribution of a specific facet to the potential at a specified field point
  @param Xf vector of dyadic coefficients for facet f at the prescribed fieldpoint
  @return contribution to the potential of this specific facet at the prescribed fieldpoint
  */
  static double GetUf(const arma::vec::fixed<10> & Xf);


  /**
  Returns the contribution of a specific edge to the acceleration at a specified field point
  @param Xe vector of dyadic coefficients for edge e at the prescribed fieldpoint
  @return contribution to the potential of this specific edge at the prescribed fieldpoint
  */
  static arma::vec::fixed<3> GetAe(const arma::vec::fixed<10> & Xe);


  /**
  Returns the contribution of a specific facet to the acceleration at a specified field point
  @param Xf vector of dyadic coefficients for facet f at the prescribed fieldpoint
  @return contribution to the potential of this specific facet at the prescribed fieldpoint
  */
  static arma::vec::fixed<3> GetAf(const arma::vec::fixed<10> & Xf);





  /**
  Returns the indices of the vertices forming the prescribed edge
  @param[in] e edge index
  @param[out] v0 first vertex index
  @param[out] v1 second vertex index
  */
  void GetIndicesVerticesOnEdge(const int & e, int & v0,int & v1){v0 = this -> edges[e][0];v1 = this -> edges[e][1];}

  /**
  Returns the indices of the vertices forming the prescribed facet
  @param[in] f facet index
  @param[out] v0 first vertex index
  @param[out] v1 second vertex index
  @param[out] v2 third vertex index
  */
  void GetIndicesVerticesInFacet(const int & f, int & v0,int & v1,int & v2){v0 = this -> facets[f][0];v1 = this -> facets[f][1];v2 = this -> facets[f][2];}


  /**
  Returns the non-normalized facet normal at this facet
  @param f facet index
  @return non-normalized facet normal
  */
  arma::vec::fixed<3> GetNonNormalizedFacetNormal(const int & f) const;


protected:
  SBGATPolyhedronGravityModel();
  ~SBGATPolyhedronGravityModel() override;


  void Clear();


  int RequestData(vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) override;

  double ** facet_dyads;
  double ** edge_dyads;
  double ** facet_normals;
  double ** vertices;

  double scaleFactor = 1;
  double density;

  int ** edges;
  int ** facets;
  int ** edge_facets_ids;

  bool scaleFactorSet;
  bool densitySet;
  bool is_in_meters = true;

  int N_facets;
  int N_edges;


  vtkSmartPointer<SBGATMassProperties> mass_properties;


private:
  SBGATPolyhedronGravityModel(const SBGATPolyhedronGravityModel&) = delete;
  void operator=(const SBGATPolyhedronGravityModel&) = delete;
};

#endif


