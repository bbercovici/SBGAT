/**
@file SBGATMassPropertiesUQUQ.hpp
@class  SBGATMassPropertiesUQUQ
@author Benjamin Bercovici
@date January 2019

@brief  Evaluation of the formal uncertainty in the volume, center of mass, inertia tensor parametrization
from a topologically-closed, constant density polyhedron.

@copyright MIT License, Benjamin Bercovici and Jay McMahon
*/

#ifndef SBGATMassPropertiesUQUQ_hpp
#define SBGATMassPropertiesUQUQ_hpp

#include <armadillo>
#include "SBGATMassProperties.hpp"


class SBGATMassPropertiesUQ {
public:


  /**
  Sets the mass model associated to this uncertainty quantification run
  and instantiates some containers whose size is depending upon
  the number of vertices in the input shape
  @param[in] pgm pointer to valid polyhedron gravity model
  */
  void SetMassProperties(vtkSmartPointer<SBGATMassProperties> mass_prop);

  /**
  Get the mass properties model associated to this uncertainty quantification run
  @return pgm pointer to valid mass properties model
  */
  vtkSmartPointer<SBGATMassProperties> GetMassProperties(){return this -> mass_prop ;}

  

  /**
  Returns a square root of the covariance matrix.
  The covariance square root is expressed in the original shape's unit squared (that is, 
  meters or kilometers). This method will attempt extracting the square root through 
  a cholesky decomposition of the covariance(P_CC = R * R.T) . Users can also ask for a square root derived from a spectral
  decomposition of the covariance (slower but more stable), where R = U * sqrt(D') * U.T with U == orthogonal matrix of eigenvectors and
  D' == diagonal matrix of clamped eigenvalues (there may be a few negative eigenvalues, these are clamped to 0). In case either of the methods fail, the 
  method return the identity matrix
  @param use_cholesky if true will make the method extract the square root of the covariance matrix through a lower Cholesky 
  decomposition of the covariance. If false, will use a clamped spectral decomposition instead
  @return covariance square root. Its unit are the same as those of the input shape (hence, m^2 or km^2)
  */
  arma::mat GetCovarianceSquareRoot(bool use_cholesky = true) const;

  /**
  Runs a finite-differencing based test of the implemented PGM partials
  @param input path to obj file used to test the partials
  @param tol relative tolerance
  */
  static void TestPartials(std::string input , double tol = 1e-2);

  /**
  Populates the covariance components over the entire shape by looping over each vertex Ci and:
  - Setting each vertex covariance to P_CiCi = standard_dev ^2 * n_i * n_i^T (n_i == surface normal)
  - Setting each correlaton matrix to P_CiCj = standard_dev ^2 * n_i * n_j^T * exp(- || Ci - Cj || ^2 / correl_distance^2) (n_i == surface normal)
  @param standard_dev standard deviation in normal component (m)
  @param correl_distance one-sigma correlation distance between points (m)
  */
  void ComputeVerticesCovarianceGlobal(const double & standard_dev,const double & correl_distance);


  /**
  Adds to the shape vertices covariance the contribution from an uncertainty region centered at vertex $regioN_verticesenter_index.
  The only 3x3 block-component of the global covariance affected by this method is the one
  starting at (3 * regioN_verticesenter_index,3 * regioN_verticesenter_index). The other components are off-diagonal
  @param standard_dev standard deviation in normal component (m)
  @param correl_distance one-sigma correlation distance between points (m)
  */
  void AddUncertaintyRegionToCovariance(int regioN_verticesenter_index,const double & standard_dev,const double & correl_distance);


  /**
  Sets the block P_Cv0_Cv1 in the total shape covariance to the prescribed value P. 
  When v0 != v1, this function must be called twice to set the two symmetric blocks
  on both sides of the diagonal. The covariance is expressed in the original shape's unit squared (that is, 
  meters squared or kilometers squared)
  @param P covariance/correlation of Cv0 and Cv1
  @param v0 index of first vertex
  @param v1 index of second vertex
  */
  void SetCovarianceComponent(const arma::mat::fixed<3,3> & P,const int & v0, const int & v1);

  /**
  Saves the vertices covariance in ascii format to the prescribed path
  @param path path where to save the covariance 
  */
  void SaveVerticesCovariance(std::string path) const;


  /**
  Returns a copy of the vertices covariance. Its unit are the same as those of the input
  shape (hence, m^2 or km^2)
  @return P_CC
  */
  arma::mat GetVerticesCovariance() const;


  /**
  Saves the non-zeros partitions of the full shape covariance to 
  a json file. Only saves the symmetric part of the covariance (i.e the upper diagonal terms)
  @param path path where to save the non-zero covariance partitions
  */
  void SaveNonZeroVerticesCovariance(std::string path) const;

  /**
  Populate the vertices covariance P_CC with the partitions saved in the provided json file
  @param path path where to save the non-zero covariance partitions
  @return 1 if loading sucessful, 0 otherwise (leaves P_CC untouched)
  */
  int LoadVerticesCovarianceFromJson(std::string path);


  /**
  Applies prescribed deviation to all the N_vertices control points and updates model
  @param delta_C deviation (3 * N_vertices x 1)
  */  
  void ApplyDeviation(const arma::vec & delta_C);


/**
Returns the partial derivative of the shape's center of mass with respect to the shape's vertices coordinates
@return partial derivative of the shape's center of mass with respect to the shape's vertices
*/
  arma::mat GetPartialCOMPartialC() const;

/**
Get partial derivatives of volume, center of mass and inertia parametrization relative to 
the shape's vertices coordinates
@param[out] dVdC partial derivative of volume relative to vertices coordinates
@param[out] dCOMdC partial derivative of center-of-mass relative to vertices coordinates
@param[out] dIdC partial derivative of inertia tensor parametrization relative to vertices coordinates
*/
  void GetPartialAllInertiaPartialC(arma::rowvec & dVdC,arma::mat & dCOMdC, arma::mat & dIdC) const;


protected:

  /**
  Returns the partial derivative of the volume of the tetrahedron subtended by facet f
  with respect to the facet coordinates
  @param f facet index
  @return partial derivative of tetrahedron volume with respect to facet coordinates
  */
  arma::rowvec::fixed<9> PartialDeltaVfPartialTf(const int & f) const;

  /**
  Returns the partial derivative of the center of mass of the considered tetrahedron
  with respect to the facet coordinates
  @param f facet index
  @return partial derivative of center of mass with respect to facet coordinates
  */
  arma::mat::fixed<3,9> PartialDeltaCMfPartialTf(const int & f) const;

  /**
  Returns the partial derivative of the tetrahedron's inertia tensor parametrization
  relative to the facet coordinates
  @param f facet index
  @return partial derivative of the tetrahedron's inertia tensor parametrization
  relative to the facet coordinates 
  */
  arma::mat::fixed<6,9> PartialDeltaIfPartialTf(const int & f) const;


  /**
  Returns the partial derivative of a tetrahedron's  inertia-times-volume tensor parametrization
  with respect to the subtending facet's vertices coordinates
  @param f facet index
  @return partial derivative of the tetrahedron's inertia-times-volume tensor parametrization
  with respect to the subtending facet's vertices coordinates
  */
  arma::mat::fixed<6,9> PartialDeltaIOverDeltaVPartialTf(const int & f) const;

  /**
  Applies deviation to the coordinates of the vertices in the prescribed facet
  and updates the pgm
  @param delta_Tf deviation
  @param f facet index
  */
  void ApplyTfDeviation(arma::vec::fixed<9> delta_Tf,const int & f);

  
  static void TestPartialDeltaVfPartialTf(std::string input,double tol);
  static void TestPartialCOMPartialC(std::string input,double tol);
  static void TestPartialDeltaIOverDeltaVPartialTf(std::string input,double tol);
  static void TestPartialDeltaIfPartialTf(std::string input,double tol);


  vtkSmartPointer<SBGATMassProperties> mass_prop;

  arma::mat P_CC;
  arma::sp_mat P_CC_sparse;

};

#endif


