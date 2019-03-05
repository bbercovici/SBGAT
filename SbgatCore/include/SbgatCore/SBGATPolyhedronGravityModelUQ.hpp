/**
@file SBGATPolyhedronGravityModelUQ.hpp
@class  SBGATPolyhedronGravityModelUQ
@author Benjamin Bercovici
@date January 2019

@brief  Evaluation of the formal uncertainty in the potential (variance), acceleration (covariance) caused by a constant-density polyhedron
 @details Computes the potential variance, acceleration covariance associated to the gravity deriving from the polyhedron
 of constant density assuming that the underlying shape vertices are outcomes of a Gaussian distribution 
 of known mean and covariance
The input must be a topologically-closed polyhedron.
See Werner, R. A., & Scheeres, D. J. (1997). Exterior gravitation of a polyhedron derived and compared with harmonic and mascon gravitation representations of asteroid 4769 Castalia. Celestial Mechanics and Dynamical Astronomy, 65(3), 313–344. https://doi.org/10.1007/BF00053511
for further details. Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
@copyright MIT License, Benjamin Bercovici and Jay McMahon
*/

#ifndef SBGATPolyhedronGravityModelUQ_hpp
#define SBGATPolyhedronGravityModelUQ_hpp

#include <armadillo>
#include "SBGATMassProperties.hpp"
#include "SBGATPolyhedronGravityModel.hpp"
#include "SBGATMassPropertiesUQ.hpp"

class SBGATPolyhedronGravityModelUQ : public SBGATMassPropertiesUQ {
public:


  /**
  Evaluates the Polyhedron Gravity Model potential variance at the specified point assuming 
  a constant density
  @param point pointer to coordinates of queried point, expressed in the same frame as
  the polydata
  @return PGM potential variance evaluated at the queried point (m ^ 4/ s ^4)
  */
  double GetVariancePotential(double const * point) const;

  /**
  Evaluates the Polyhedron Gravity Model potential variance at the specified point assuming 
  a constant density
  @param point coordinates of queried point, expressed in the same frame as
  the polydata
  @return PGM potential variance evaluated at the queried point (m ^ 4 / s ^4)
  */
  double GetVariancePotential(const arma::vec::fixed<3> & point) const;

  /**
  Evaluates the Polyhedron Gravity Model potential variance and acceleration covariance at the specified point assuming 
  a constant density
  @param point coordinates of queried point, expressed in the same frame as
  the polydata used to construct the PGM
  @param[out] potential_var PGM potential variance evaluated at the queried point (m ^ 4 / s ^4)
  @param[out] acc_cov PGM acceleration covariance evaluated at the queried point (m^2 / s ^4)
  */
  void GetVariancePotentialAccelerationCovariance(double const * point,double & potential_var, 
    arma::mat::fixed<3,3> & acc_cov) const;



  /**
  Return the variance of the slope evaluated at the center of the designated facet. This method is NOT thread safe
  @param f facet index
  @return variance in slope (rad^2)
  */
  double GetVarianceSlope(const int & f );



  /**
  Return the variance of the slope evaluated at the center of the designated facets. This method is NOT thread safe
  @param[out] slope_variances
  @param[in] facets indices of facets where to evaluate the slope variance
  */
  void GetVarianceSlopes(std::vector<double> & slope_variances,const std::vector<int> & facets);


  /**
  Evaluates the Polyhedron Gravity Model potential variance and acceleration covariance at the specified point assuming 
  a constant density
  @param point coordinates of queried point, expressed in the same frame as
  the polydata used to construct the PGM
  @param[out] potential_var PGM potential variance evaluated at the queried point (m ^ 4 / s ^4)
  @param[out] acc_cov PGM acceleration covariance evaluated at the queried point (m^2 / s ^4)
  */
  void GetVariancePotentialAccelerationCovariance(const arma::vec::fixed<3> & point,double & potential_var, 
    arma::mat::fixed<3,3> & acc_cov) const;


  /**
  Runs a finite-differencing based test of the implemented PGM partials
  @param input path to obj file used to test the partials
  @param tol relative tolerance
  */
  static void TestPartials(std::string input , double tol, bool shape_in_meters);

  /**
  Obtain the partial derivative of the potential at the prescribed location
  due to a infinitesimal variation in the shape's control points
  @param pos position where to evaluate the partial derivative
  @return partial derivative of the potential with respect to the variation in the shape's control points
  */

  arma::rowvec GetPartialUPartialC(const arma::vec::fixed<3> & pos) const;

  /**
  Obtain the partial derivative of the acceleration at the prescribed location
  due to a infinitesimal variation in the shape's control points
  @param pos position where to evaluate the partial derivative
  @return partial derivative of the acceleration with respect to the variation in the shape's control points
  */
  arma::mat GetPartialAPartialC(const arma::vec::fixed<3> & pos) const;



  /**
  Get covariance in acceleration arising from the uncertain shape
  @param point coordinates where to evaluate the covariance
  @return covariance of acceleration
  */
  arma::mat::fixed<3,3> GetCovarianceAcceleration(double const * point) const;

   /**
  Get covariance in acceleration arising from the uncertain shape
  @param point coordinates where to evaluate the covariance
  @return covariance of acceleration
  */
  arma::mat::fixed<3,3> GetCovarianceAcceleration(const arma::vec::fixed<3> & point) const;


  /**
  Runs a Monte Carlo on the shape and samples accelerations & potentials at the provided positions
  @param[in] path_to_shape path to reference shape
  @param[in] density small body density in kg/m^3
  @param[in] shape_in_meters true if reference shape has its units expressed in meters, false otherwise
  @param[in] P_CC square root of the covariance of the shape vertices coordinates. Must be of dimensions (3N_C * 3N_C) where N_C is the
  number of vertices in the reference shape
  @param[in] N_samples number of shape outcomes to draw
  @param[in] all_positions vector storing all the position where acceleration & potential must be sampled
  @param[in] output_dir path ending in "/" where to save shape-related monte-carlo data. Only used
  if last argument is larger than 0
  @param[in] N_saved_shapes number of shape outcomes to save. must be lesser or equal than N_samples
  @param[out] deviations holds N_samples 3*N_C column vectors storing the deviation applied to the coordinates of the 
  reference shape at every sample
  @param[out] all_accelerations holds N_samples vectors, each storing the acceleration evaluated at the specified points
  @param[out] all_potentials holds N_samples vectors, each storing the potential evaluated at the specified points
  */

  static void RunMCUQPotentialAccelerationInertial(std::string path_to_shape,
    const double & density,
    const bool & shape_in_meters,
    const arma::mat & C_CC,
    const unsigned int & N_samples,
    const std::vector<arma::vec::fixed<3> > & all_positions,
    std::string output_dir,
    int N_saved_shapes,
    std::vector<arma::vec> & deviations,
    std::vector<std::vector<arma::vec::fixed<3> >> & all_accelerations,
    std::vector < std::vector<double> > & all_potentials );



/**
  Runs a Monte Carlo on the shape and samples accelerations at the provided positions
  @param[in] path_to_shape path to reference shape
  @param[in] density small body density in kg/m^3
  @param[in] shape_in_meters true if reference shape has its units expressed in meters, false otherwise
  @param[in] P_CC square root of the covariance of the shape vertices coordinates. Must be of dimensions (3N_C * 3N_C) where N_C is the
  number of vertices in the reference shape
  @param[in] N_samples number of shape outcomes to draw
  @param[in] all_positions vector storing all the position where acceleration & potential must be sampled
  @param[in] output_dir path ending in "/" where to save shape-related monte-carlo data. Only used
  if last argument is larger than 0
  @param[in] N_saved_shapes number of shape outcomes to save. must be lesser or equal than N_samples
  @param[out] deviations holds N_samples 3*N_C column vectors storing the deviation applied to the coordinates of the 
  reference shape at every sample
  @param[out] all_accelerations holds N_samples vectors, each storing the acceleration evaluated at the specified points
  */

static void RunMCUQAccelerationInertial(std::string path_to_shape,
  const double & density,
  const bool & shape_in_meters,
  const arma::mat & C_CC,
  const unsigned int & N_samples,
  const std::vector<arma::vec::fixed<3> > & all_positions,
  std::string output_dir,
  int N_saved_shapes,
  std::vector<arma::vec> & deviations,
  std::vector<std::vector<arma::vec::fixed<3> >> & all_accelerations);



  /**
  Runs a Monte Carlo on the shape and samples the slopes at the provided facets
  @param[in] path_to_shape path to reference shape
  @param[in] density small body density in kg/m^3
  @param[in] Omega angular velocity of small body in kg/m^3, expressed in the small body frame
  @param[in] shape_in_meters true if reference shape has its units expressed in meters, false otherwise
  @param[in] C_CC square root of the covariance of the shape vertices coordinates. Must be of dimensions (3N_C * 3N_C) where N_C is the
  number of vertices in the reference shape
  @param[in] period_standard_deviation standard deviation of the rotation period in seconds
  @param[in] N_samples number of shape outcomes to draw
  @param[in] all_facets vector storing all the facet indices where the gravitational slopes must be sampled
  @param[in] output_dir path ending in "/" where to save shape-related monte-carlo data. Only used
  if last argument is larger than 0
  @param[in] N_saved_shapes number of shape outcomes to save. must be lesser or equal than N_samples
  @param[out] deviations holds N_samples 3*N_C column vectors storing the deviation applied to the coordinates of the 
  reference shape at every sample
  @param[out] period_errors holds N_samples of the error on the rotation period

  @param[out] all_slopes holds N_samples vectors, each storing the slopes evaluated at the specified facets
  */

  static void RunMCUQSlopes(std::string path_to_shape,
    const double & density,
    const arma::vec::fixed<3> & Omega,
    const bool & shape_in_meters,
    const arma::mat & C_CC,
    const double & period_standard_deviation,
    const unsigned int & N_samples,
    const std::vector<int > & all_facets,
    std::string output_dir,
    int N_saved_shapes,
    std::vector<arma::vec> & deviations,
    std::vector<double> & period_errors,
    std::vector < std::vector<double> > & all_slopes );




  /**
  Sets the standard deviation of the rotation period
  @param standard deviation of the rotation period (s)
  */
  void SetPeriodErrorStandardDeviation(double rotation_period_sd){
    this -> period_standard_deviation = rotation_period_sd;
  }


  /**
  Return the partial derivative of the slope at the center of facet f relative to 
  the angular velocity magnitude and shape vertices coordinates
  @param f facet index
  @return partials
  */
  arma::rowvec GetPartialSlopePartialwPartialC(const int & f) const;

  
  /**
  Applies prescribed deviation to all the N_vertices control points and updates model
  @param delta_C deviation (3 * N_vertices x 1)
  */  
  virtual void ApplyDeviation(const arma::vec & delta_C);

  /**
  Return the partial derivative of the gravitation slope at the center of face tf relative to 
  the shape vertices coordinates
  @param f facet index
  @return partial derivative of the slope at the center of facet f relative to the shape vertices coordinates
  */

  arma::rowvec GetPartialSlopePartialC(const int & f) const;

  /**
  Runs a Monte Carlo on the shape and samples inertial accelerations & potentials at the provided position
  @param[in] path_to_shape path to reference shape
  @param[in] density small body density in kg/m^3
  @param[in] shape_in_meters true if reference shape has its units expressed in meters, false otherwise
  @param[in] C_CC square root of the covariance of the shape vertices coordinates. Must be of dimensions (3N_C * 3N_C) where N_C is the
  number of vertices in the reference shape
  @param[in] N_samples number of shape outcomes to draw
  @param[in] position the position where acceleration & potential must be sampled
  @param[in] output_dir path ending in "/" where to save shape-related monte-carlo data. Only used
  if last argument is larger than 0
  @param[in] N_saved_shapes number of shape outcomes to save. must be lesser or equal than N_samples
  @param[out] deviations holds N_samples 3*N_C column vectors storing the deviation applied to the coordinates of the 
  reference shape at every sample
  @param[out] accelerations holds N_samples accelerations evaluated at the specified point
  @param[out] potentials holds N_samples potential evaluated at the specified point
  */
  static void RunMCUQPotentialAccelerationInertial(std::string path_to_shape,
    const double & density,
    const bool & shape_in_meters,
    const arma::mat & C_CC,
    const unsigned int & N_samples,
    const arma::vec::fixed<3> & position,
    std::string output_dir,
    int N_saved_shapes,
    std::vector<arma::vec> & deviations,
    std::vector<arma::vec::fixed<3> > & accelerations,
    std::vector<double> & potentials);


/**
  Runs a Monte Carlo on the shape and samples inertial accelerations at the provided position
  @param[in] path_to_shape path to reference shape
  @param[in] density small body density in kg/m^3
  @param[in] shape_in_meters true if reference shape has its units expressed in meters, false otherwise
  @param[in] C_CC square root of the covariance of the shape vertices coordinates. Must be of dimensions (3N_C * 3N_C) where N_C is the
  number of vertices in the reference shape
  @param[in] N_samples number of shape outcomes to draw
  @param[in] position the position where acceleration & potential must be sampled
  @param[in] output_dir path ending in "/" where to save shape-related monte-carlo data. Only used
  if last argument is larger than 0
  @param[in] N_saved_shapes number of shape outcomes to save. must be lesser or equal than N_samples
  @param[out] deviations holds N_samples 3*N_C column vectors storing the deviation applied to the coordinates of the 
  reference shape at every sample
  @param[out] accelerations holds N_samples accelerations evaluated at the specified point
  */
static void RunMCUQAccelerationInertial(std::string path_to_shape,
  const double & density,
  const bool & shape_in_meters,
  const arma::mat & C_CC,
  const unsigned int & N_samples,
  const arma::vec::fixed<3> & position,
  std::string output_dir,
  int N_saved_shapes,
  std::vector<arma::vec> & deviations,
  std::vector<arma::vec::fixed<3> > & accelerations);




  /**
  Runs a Monte Carlo on the shape and samples the gravitational slopes at the provided facets
  @param[in] path_to_shape path to reference shape
  @param[in] density small body density in kg/m^3
  @param[in] Omega angular velocity of small body in kg/m^3, expressed in the small body frame
  @param[in] shape_in_meters true if reference shape has its units expressed in meters, false otherwise
  @param[in] C_CC square root of the covariance of the shape vertices coordinates. Must be of dimensions (3N_C * 3N_C) where N_C is the
  number of vertices in the reference shape
  @param[in] period_standard_deviation standard deviation of the rotation period in seconds
  @param[in] N_samples number of shape outcomes to draw
  @param[in] facet index of the facet where the surface pgm must be sampled
  @param[in] output_dir path ending in "/" where to save shape-related monte-carlo data. Only used
  if last argument is larger than 0
  @param[in] N_saved_shapes number of shape outcomes to save. must be lesser or equal than N_samples
  @param[out] deviations holds N_samples 3*N_C column vectors storing the deviation applied to the coordinates of the 
  reference shape at every sample
  @param[out] period_errors holds N_samples of the error on the rotation period
  @param[out] slopes holds N_samples slopes evaluated at the specified facet
  */
  static void RunMCUQSlopes(std::string path_to_shape,
    const double & density,
    const arma::vec::fixed<3> & Omega,
    const bool & shape_in_meters,
    const arma::mat & C_CC,
    const double & period_standard_deviation,
    const unsigned int & N_samples,
    const int & facet,
    std::string output_dir,
    int N_saved_shapes,
    std::vector<arma::vec> & deviations,
    std::vector<double> & period_errors,
    std::vector<double> & slopes);





protected:

  arma::vec GetBe() const;



  /**
  Get partial derivative of the angular velocity vector relative to 1) the angular velocity magnitude 2) the shape vertices
  coordinates
  @return partial derivative of Omega relative to its magnitude shape vertices coordinates
  */
  arma::mat PartialOmegaPartialwC() const;


  /**
  Return the partial derivative of the body-fixed angular velocity and the vertices coordinates relative
  to the angular velocity magnitude and shape vertices coordinates
  @return partial
  */
  arma::sp_mat PartialOmegaCPartialwC() const;


/**
Return the partial derivative of the body-fixed acceleration at the center of facet f relative
to the shape coordinates
@param f facet index
@return partial derivative
*/
  arma::mat PartialBodyFixedAccelerationfPartialC(const int & f) const;



  /**
  Return the partial derivative of the body-fixed acceleration at the center of facet f with respect to 
  the angular velocity and vertices coordinates
  @param f facet index
  @return partial derivative
  */
  arma::mat PartialBodyFixedAccelerationfPartialOmegaC(const int & f) const;




  arma::mat::fixed<3,3> PartialBodyFixedAccelerationfPartialOmega(const int & f) const;




  /**
  Adds to the properly initialized vector the partial derivative of the sum of all Ue
  @param[in] pos position where to evaluate the partials
  @param[out] partial partial derivative being evaluated
  */
  void AddPartialSumUePartialC(const arma::vec::fixed<3> & pos,arma::rowvec & partial) const;

   /**
  Add to the properly initialized vector the partial derivative of the sum of all Uf
  @param[in] pos position where to evaluate the partials
  @param[out] partial partial derivative being evaluated
  */
  void AddPartialSumUfPartialC(const arma::vec::fixed<3> & pos,arma::rowvec & partial) const;


  /**
  Adds to the properly initialized vector the partial derivative of the sum of all Acce
  @param[in] pos position where to evaluate the partials
  @param[out] partial partial derivative being evaluated
  */
  void AddPartialSumAccePartialC(const arma::vec::fixed<3> & pos,arma::mat & partial) const;

   /**
  Add to the properly initialized vector the partial derivative of the sum of all Accf
  @param[in] pos position where to evaluate the partials
  @param[out] partial partial derivative being evaluated
  */
  void AddPartialSumAccfPartialC(const arma::vec::fixed<3> & pos,arma::mat & partial) const;


  /**
  Applies deviation to the coordinates of the vertices on the prescribed edge
  and updates the pgm
  @param delta_Ae deviation
  @param e edge index
  */
  void ApplyAeDeviation(arma::vec::fixed<6> delta_Ae,const int & e);


  /**
  Applies deviation to the coordinates of the vertices in the prescribed facet
  and updates the pgm
  @param delta_Tf deviation
  @param f facet index
  */
  void ApplyTfDeviation(arma::vec::fixed<9> delta_Tf,const int & f);


  /**
  Return the partial derivative of an individual edge contribution to the potential (Ue) 
  with respect to the Xe^E vector holding the e-th edge dyadic factors
  @param pos position where to evaluate the partial
  @param e edge index
  @return PartialUePartialXe (1x10)
  */
  arma::rowvec::fixed<10> PartialUePartialXe(const arma::vec::fixed<3> & pos,const int & e) const;


  /**
  Return the partial derivative of an individual facet contribution to the potential (Uf) 
  with respect to the Xf^F vector holding the f-th facet dyadic factors
  @param pos position where to evaluate the partial
  @param f facet index
  @return PartialUfPartialXf (1x10)

  */
  arma::rowvec::fixed<10> PartialUfPartialXf(const arma::vec::fixed<3> & pos,
    const int & f) const;




  /**
  Return the partial derivative of an individual edge contribution to the acceleration (Acce) 
  with respect to the Xe^E vector holding the e-th edge dyadic factors
  @param pos position where to evaluate the partial
  @param e edge index
  @return PartialAccePartialXe (3x10)
  */
  arma::mat::fixed<3,10> PartialAccePartialXe(const arma::vec::fixed<3> & pos,const int & e) const;


  /**
  Return the partial derivative of an individual facet contribution to the acceleration (Accf) 
  with respect to the Xf^F vector holding the f-th facet dyadic factors
  @param pos position where to evaluate the partial
  @param f facet index
  @return PartialAccfPartialXf (3x10)
  */
  arma::mat::fixed<3,10> PartialAccfPartialXf(const arma::vec::fixed<3> & pos,const int & f) const;



  /**
  Return the partial derivative of Xf^F, the vector holding the f-th facet dyadic factors, 
  with respect to the vertices coordiantes constitutive of the f-th triangle (Tf) 
  @param pos position where to evaluate the partial
  @param f facet index
  @return PartialXfPartialTf (10x9)
  */
  arma::mat::fixed<10,9> PartialXfPartialTf(const arma::vec::fixed<3> & pos, const int & f) const;


  /**
  Return the partial derivative of the performance factor omega_f
  with respect to the vertices coordiantes constitutive of the f-th triangle (Tf) 
  @param pos position where to evaluate the partial
  @param f facet index
  @return PartialOmegafPartialTf (1x9)
  */
  arma::rowvec::fixed<9> PartialOmegafPartialTf(const arma::vec::fixed<3> & pos,const int & f) const;



  /**
  Return the partial derivative of Z_f = (alpha_f,gamma_f)^T (as in wf = 2 * arctan2(Z_f) )
  with respect to the unit vectors from the field point to the facet vertices
  @param UnitRf 3 unit vectors stacked up
  @return PartialZfPartialUnitRf (2x9)
  */
  static arma::mat::fixed<2,9> PartialZfPartialUnitRf(const arma::vec::fixed<9> & UnitRf);


  /**
  Return the partial derivative of arctan2(Z_f) w/r to Z_f 
  with respect to the unit vectors from the field point to the facet vertices
  @param Zf 
  @return PartialAtan2PartialZf (1x2)
  */
  static arma::rowvec::fixed<2> PartialAtan2PartialZf(const arma::vec::fixed<2> & Zf);


  /**
  Return the partial derivative of arctan(y/x)
  @param xy input
  @return partial derivative
  */
  static arma::rowvec::fixed<2> PartialOmegafPartialXY(const arma::vec::fixed<2> & xy);


  /**
  Return the partial derivative of the facet dyad parametrization (Ff)
  with respect to the vertices coordinates constitutive of the f-th triangle (Tf) 
  @param f facet index
  @return PartialFfPartialTf (6x9)
  */
  arma::mat::fixed<6,9> PartialFfPartialTf(const int & f) const;



  /**
  Return the partial derivative of a normalized vector n relative to the non-normalized
  vector N such that n = N / || N ||
  @param non_normalized_V non-normalized vector used to produce the normalized vector
  @return PartialNormalizedVPartialNonNormalizedV (3x3)
  */
  static arma::mat::fixed<3,3> PartialNormalizedVPartialNonNormalizedV(const arma::vec::fixed<3> & non_normalized_V);



  /**
  Return the partial derivative of the f-th facet dyad parametrization with respect to the 
  normalized normal of the f-th facet
  @param nf facet normal
  @return PartialFfPartialnf (6x3)

  */
  static arma::mat::fixed<6,3> PartialFfPartialnf(const arma::vec::fixed<3> & nf);



  /**
  Return the partial derivative of the wire potential Le 
  with respect to the coordinates of the two vertices forming the edge (stacked in Ae)
  @param pos position where to evaluate the partial
  @param e edge index
  @return PartialLePartialAe (1x6)
  */
  arma::rowvec::fixed<6> PartialLePartialAe(const arma::vec::fixed<3> & pos,const int & e) const;


  /**
  Return the partial derivative of field-point to edge-point vector
  with respect to the coordinates of the two vertices forming the edge (stacked in Ae)
  @return PartialRadiusEePartialAe (3x6)
  */
  arma::mat::fixed<3,6> PartialRadiusEePartialAe() const;


  /**
  Return the partial derivative of field-point to facet-point vector
  with respect to the coordinates of the three vertices forming the facet (stacked in Tf)
  @return PartialRadiusFfPartialTf (3x9)
  */
  arma::mat::fixed<3,9> PartialRadiusFfPartialTf() const;


  /**
  Return the partial derivative of the parametrization of the Xe dyadic vector
  with respect to the coordinates of the edges points and adjacent facets points
  @param pos position where to evaluate the partial
  @param e edge index
  @return PartialXePartialBe (10x24)
  */
  arma::mat::fixed<10,24> PartialXePartialBe(const arma::vec::fixed<3> & pos,const int & e) const;

  /**
  Return the partial derivative of the edge length le
  with respect to the coordinates of the edges points
  @param e edge index
  @return PartialEdgeLengthPartialAe (10x24)
  */
  arma::rowvec::fixed<6> PartialEdgeLengthPartialAe(const int & e) const;

  /**
  Return the partial derivative of the (q,r) component of the Ee dyad with respect to the 
  with respect to the coordinates of the edges points and adjacent facets points
  @param e edge index
  @param q row index
  @param r col index
  @return PartialEqrPartialBe (1x24)
  */
  arma::rowvec::fixed<24> PartialEqrPartialBe(const int & e,const int & q,const int & r) const;


  /**
  Return the partial derivative of the f-th facet slope argument (u as in slope = arcos(-u))
  with respect to the angular velocity and the shape coordinates
  @param f facet index
  @param body_fixed_acc body-fixed acceleration at the center of facet f
  @return partial derivative
  */
  arma::rowvec PartialSlopeArgumentPartialOmegaC(const int & f,
    const arma::vec::fixed<3> & body_fixed_acc) const;


  /**
  Return the partial derivative of the Ee dyad parametrization with respect 
  to the coordinates of the edges points and adjacent facets points
  @param e edge index
  @return PartialEPartialBe (6x24)
  */
  arma::mat::fixed<6,24> PartialEPartialBe(const int & e) const;


  /**
  Return the connectivity table associated with vector Be
  @param e edge index
  @return connectivity table
  */
  arma::sp_mat  PartialBePartialC(const int & e) const;



  /**
  Return the partial derivative of the slope s == arcos(-u) relative to the slope argument u
  @param u input parameter
  @param partial derivative of slope with respect to u
  */
  static double PartialSlopePartialSlopeArgument(const double & u);


  /**
  Given a prescribed global deviation of all of the shape's N control points,
  applies it and returns the deviation in each of the Be's vector (one per edge in the shape)
  @param delta deviation in all of the shape's N control points (3 x N_vertices)
  @return deviation in all of the shape's Be vectors (24 x N_edges)
  */
  arma::vec ApplyAndGetBeDeviation(const arma::vec & delta);


  static void TestPartialUePartialXe(std::string input , double tol, bool shape_in_meters) ;
  static void TestPartialUfPartialXf(std::string input , double tol, bool shape_in_meters) ;
  static void TestPartialXfPartialTf(std::string input , double tol, bool shape_in_meters) ;
  static void TestPartialOmegafPartialTf(std::string input , double tol, bool shape_in_meters) ;
  static void TestPartialZfPartialUnitRf(std::string input , double tol, bool shape_in_meters);
  static void TestPartialFfPartialTf(std::string input , double tol, bool shape_in_meters) ;
  static void TestPartialNormalizedVPartialNonNormalizedV(std::string input , double tol, bool shape_in_meters) ;
  static void TestPartialAtan2PartialZf(std::string input , double tol, bool shape_in_meters);
  static void TestPartialNfPartialTf(std::string input , double tol, bool shape_in_meters) ;
  static void TestPartialFfPartialnf(std::string input , double tol, bool shape_in_meters) ;
  static void TestPartialFfPartialNonNormalizedNf(std::string input , double tol, bool shape_in_meters);
  static void TestPartialLePartialAe(std::string input , double tol, bool shape_in_meters) ;
  static void TestPartialEePartialAe(std::string input , double tol, bool shape_in_meters) ;
  static void TestPartialEePartialTf(std::string input , double tol, bool shape_in_meters) ;
  static void TestPartialXePartialBe(std::string input , double tol, bool shape_in_meters) ;
  static void TestPartialEdgeLengthPartialAe(std::string input , double tol, bool shape_in_meters) ;
  static void TestPartialEPartialBe(std::string input , double tol, bool shape_in_meters) ;
  static void TestPartialUfPartialTf(std::string input , double tol, bool shape_in_meters);
  static void TestPartialUePartialBe(std::string input , double tol, bool shape_in_meters);
  static void TestPartialUPartialC(std::string input , double tol, bool shape_in_meters);
  static void TestGetPartialAPartialC(std::string input , double tol, bool shape_in_meters);
  static void TestPartialUfPartialC(std::string input , double tol, bool shape_in_meters);
  static void TestPartialUePartialC(std::string input , double tol, bool shape_in_meters);
  static void TestAddPartialSumUePartialC(std::string input , double tol, bool shape_in_meters);
  static void TestAddPartialSumUfPartialC(std::string input , double tol, bool shape_in_meters);
  static void TestAddPartialSumAccfPartialC(std::string input , double tol, bool shape_in_meters);
  static void TestAddPartialSumAccePartialC(std::string input , double tol, bool shape_in_meters);
  static void TestPartialBePartialC(std::string input , double tol, bool shape_in_meters);
  static void TestGetPartialSlopePartialwPartialC(std::string input , double tol, bool shape_in_meters);
  static void TestPartialBodyFixedAccelerationfPartialC(std::string input , double tol, bool shape_in_meters);
  static void TestPartialBodyFixedAccelerationfPartialwC(std::string input , double tol, bool shape_in_meters);
  static void TestPartialOmegaPartialwC(std::string input , double tol, bool shape_in_meters);
  static void TestPartialSlopeArgumentPartialOmegaC(std::string input , double tol, bool shape_in_meters);

  double period_standard_deviation;

  arma::mat precomputed_partialGpartialC;
  arma::mat precomputed_partialSigmapartialC;




};

#endif


