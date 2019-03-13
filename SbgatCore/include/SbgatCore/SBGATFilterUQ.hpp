/**
@file SBGATFilterUQUQ.hpp
@class  SBGATFilterUQUQ
@author Benjamin Bercovici
@date January 2019

@brief  virtual class serving as the base of dedicated uncertainty quantification classes
in SBGAT

@copyright MIT License, Benjamin Bercovici and Jay McMahon
*/

#ifndef SBGATFilterUQ_hpp
#define SBGATFilterUQ_hpp

#include <armadillo>
#include <SBGATFilter.hpp>

class SBGATFilterUQ {
public:

  /**
  Sets the model associated to this uncertainty quantification container
  @param[in] pgm pointer to valid SBGATFilter
  */
  virtual void SetModel(vtkSmartPointer<SBGATFilter> model) = 0;

  /**
  Get the associated to this uncertainty quantification container
  @returnpointer to valid SBGATFilter model
  */
  vtkSmartPointer<SBGATFilter> GetModel(){return this -> model ;}

  /**
  Return a square root of the covariance matrix.
  The covariance square root is expressed in the original shape's unit squared (that is, 
  meters or kilometers). This method will attempt extracting the square root through 
  a cholesky decomposition of the covariance(P_CC = R * R.T) . Users can also ask for a square root derived from a spectral
  decomposition of the covariance (slower but more stable), where R = U * sqrt(D') * U.T with U == orthogonal matrix of eigenvectors and
  D' == diagonal matrix of clamped eigenvalues (there may be a few negative eigenvalues, these are clamped to 0). In case either of the methods fail, the 
  method return the identity matrix
  @param method picked from {"chol","eigen"}. "chol" will attempt extracting a lower cholesky decomposition
  of the covariance matrix. "eigen" looks for a spectral decomposition of P = U * D * U^T and returns U * sqrt(D) * U^T 
  where any negative eigenvalue in D has been clamped to zero. 
  @return covariance square root. Its unit are the same as those of the input shape (hence, m^2 or km^2)
  */
  arma::mat GetCovarianceSquareRoot(std::string method = "chol") const;

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
  @param standard_dev standard deviation in radial component (m)
  @param correl_distance one-sigma correlation distance between points (m)
  */
  void AddRadialUncertaintyRegionToCovariance(int regioN_verticesenter_index,const double & standard_dev,const double & correl_distance);

  /**
  Will 'regularize' the shape covariance by computing a spectral 
  decomposition P_CC = H D H^T and clamp any negative eigenvalue in D to zero. 
  After calling this method, the covariance thus becomes
  P_CC = H clamped(D) H^T
  @param return 0 if no eigenvalue was clamped, 1 otherwise
  */
  int RegularizeCovariance();


  /**
  Adds to the shape vertices covariance the contribution from an uncertainty region centered at vertex $regioN_verticesenter_index.
  The only 3x3 block-component of the global covariance affected by this method is the one
  starting at (3 * regioN_verticesenter_index,3 * regioN_verticesenter_index). The other components are off-diagonal
  @param standard_dev standard deviation in normal component (m)
  @param correl_distance one-sigma correlation distance between points (m)
  */
  void AddNormalUncertaintyRegionToCovariance(int regioN_verticesenter_index,const double & standard_dev,const double & correl_distance);


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
Sets the vertices covariance to the provided one
@parma P_CC (3 N_vertices * 3 N_vertices) covariance matrix (m^2 or km^2)
*/
  void SetCovariance(const arma::mat & P_CC);


  /**
  Saves the vertices covariance in ascii format to the prescribed path
  @param path path where to save the covariance 
  */
  void SaveVerticesCovariance(std::string path) const;


  /**
  Return a copy of the vertices covariance. Its unit are the same as those of the input
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
  virtual void ApplyDeviation(const arma::vec & delta_C) = 0;

  /**
  Return the connectivity table associated with vector Tf, the parametrization 
  of the vertices coordinates forming the f-th facet
  @param f facet index
  @return connectivity table
  */
  arma::sp_mat  PartialTfPartialC(const int & f) const;


  /**
  Return the partial derivative of the non-normalized normal at the f-th facet 
  with respect to the vertices coordinates constitutive of the f-th triangle (Tf) 
  @param f facet index
  @return PartialNfPartialTf (3x9)
  */
  arma::mat::fixed<3,9> PartialNfPartialTf(const int & f) const;



  /**
  Takes a slice of the shape along a plane of normal e_axis (axis in {0,1,2})
  and saves the corresponding interesections (line segments) to a text file
  @param axis direction of plane normal (axis in {0,1,2})
  @param path savepath
  @param c distance from the plane to (0,0,0) measured along the plane normal vector
  */
  void TakeAndSaveSlice(int axis,std::string path, const double & c) const ;


  /**
  Return a square root of the provided covariance matrix.
  This method will attempt extracting the square root through 
  a cholesky decomposition of the covariance(P_CC = R * R.T) . Users can also ask for a square root derived from a spectral
  decomposition of the covariance (slower but more stable), where R = U * sqrt(D') * U.T with U == orthogonal matrix of eigenvectors and
  D' == diagonal matrix of clamped eigenvalues (there may be a few negative eigenvalues, these are clamped to 0). In case either of the methods fail, the 
  method return the identity matrix
  @param method picked from {"chol","eigen"}. "chol" will attempt extracting a lower cholesky decomposition
  of the covariance matrix. "eigen" looks for a spectral decomposition of P = U * D * U^T and returns U * sqrt(D) * U^T 
  where any negative eigenvalue in D has been clamped to zero. 
  @param P_CC symmetric positive semi-definite matrix
  @return covariance square root. 
  */
  static arma::mat GetCovarianceSquareRoot(arma::mat P_CC,std::string method);


  /**

  Returns the Kullback-Leibler divergence of the two gaussians N0 and N1, D_KL(N0 || N1)
  @param m0 mean of first gaussian
  @param m1 mean of second gaussian
  @param P0 covariance of second gaussian
  @param P1 covariance of second gaussian
  @return KL divergence 
  */
  static  double KLDivergence(const arma::vec & m0,const arma::mat & m1,const arma::mat & P0,const arma::mat & P1);

protected:


  /**
  Applies deviation to the coordinates of the vertices in the prescribed facet
  and updates the pgm
  @param delta_Tf deviation
  @param f facet index
  */
  virtual void ApplyTfDeviation(arma::vec::fixed<9> delta_Tf,const int & f) = 0;



  void SaveSlice(int axis, std::string path, const std::vector<std::vector<arma::vec> > & lines) const;
  void TakeSlice(int axis,std::vector<std::vector<arma::vec> > & lines,const double & c) const;

  vtkSmartPointer<SBGATFilter> model;

  arma::mat P_CC;
  arma::sp_mat P_CC_sparse;

};

#endif


