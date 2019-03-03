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
#include <SBGATFilterUQ.hpp>



class SBGATMassPropertiesUQ : public SBGATFilterUQ{
public:


  /**
  Runs a finite-differencing based test of the implemented PGM partials
  @param input path to obj file used to test the partials
  @param tol relative tolerance
  @param shape_in_meters true if tested shape has its coordinates expressed in meters
  */
  static void TestPartials(std::string input , double tol ,bool shape_in_meters);


/**
Returns the partial derivative of the shape's center of mass with respect to the shape's vertices coordinates
@return partial derivative of the shape's center of mass with respect to the shape's vertices
*/
  arma::mat GetPartialComPartialC() const;


  /**
  Returns the partial derivative of the 6 unique components of the inertia tensor {I(0,0),I(1,1),I(2,2),I(0,1),I(0,2),I(1,2)}
  with respect to the shape coordinates
  @return partial derivative
  */
  arma::mat GetPartialIPartialC() const;


  /**
  Returns the partial derivative of the MRP orienting the body-frame (B) to principal-frame (P) dcm (PB)
  with respect to the shape vertices coordinates
  @return partial derivative
  */
  arma::mat GetPartialSigmaPartialC() const;


  /**
  Returns the partial derivative of the volume
  with respect to the shape coordinates
  @return partial derivative
  */
  arma::rowvec GetPartialVolumePartialC() const;



/**
Get partial derivatives of volume, center of mass and inertia parametrization relative to 
the shape's vertices coordinates
@param[out] dVdC partial derivative of volume relative to vertices coordinates
@param[out] dCOMdC partial derivative of center-of-mass relative to vertices coordinates
@param[out] dIdC partial derivative of inertia tensor parametrization relative to vertices coordinates
*/
  void GetPartialAllInertiaPartialC(arma::rowvec & dVdC,arma::mat & dCOMdC, arma::mat & dIdC) const;


  /**
  Applies prescribed deviation to all the N_vertices control points and updates model
  @param delta_C deviation (3 * N_vertices x 1)
  */  
  virtual void ApplyDeviation(const arma::vec & delta_C);

protected:

/**
  Returns the partial derivative of the unit density moments relative to a change in the inertia tensor parametrization
  @return  partial derivative of the unit density moments relative to a change in the inertia tensor parametrization
  */
  arma::mat GetPartialUnitDensityMomentsPartialI() const;



 


  /**
  Returns the partial derivative of the MRP orienting the body-frame (B) to principal-frame (P) dcm (PB)
  with respect to the inertia tensor parametrization
  @return partial derivative
  */
  arma::mat::fixed<3,6> GetPartialSigmaPartialI() const;


/**
Returns the partial derivative of the (q,r) component of the contribution of the f-facet 
to the shape's inertia tensor, relative to the f-facet 
vertices coordinates
@param f facet index
@param q first index
@param r second index
@param Tf 9x1 vector holding coordinates of vertices in facet EXPRESSED IN METERS
@return partial derivative of the (q,r) component of the contribution of the f-facet 
to the shape's inertia tensor
*/
  arma::rowvec::fixed<9> PartialEqDeltaIfErPartialTf(const int & f, const int & q, const int & r,const arma::vec::fixed<9> & Tf) const;

  /**
  Returns the partial derivative of the shape's center-of-mass
  with respect to the f-th facet coordinates
  @param f facet index
  @return partial derivative of shape's center-of-mass with respect to facet coordinates
  */
  arma::mat::fixed<3,9> PartialDeltaComPartialTf(const int & f) const;

  /**
  Returns the partial derivative of e_q.T * DeltaIOverDeltaVfEr * e_r with respect to the f-facet 
  vertices coordinates
  @param e_q first 3x1 vector canonical unit vector
  @param e_r first 3x1 vector canonical unit vector
  @param Tf 9x1 vector holding coordinates of vertices in facet EXPRESSED IN METERS
  @return the partial derivative of e_q.T * DeltaIOverDeltaVfEr * e_r with respect to the f-facet 
  vertices coordinates
  */
  arma::rowvec::fixed<9> PartialEqDeltaIOverDeltaVfErPartialTf(const arma::vec::fixed<3> & e_q,const arma::vec::fixed<3> & e_r,
    const arma::vec::fixed<9> & Tf) const;





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
  @return partial derivative of center of mass with respect to facet coordinates
  */
  static arma::mat::fixed<3,9> PartialDeltaCMfPartialTf();

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
  virtual void ApplyTfDeviation(arma::vec::fixed<9> delta_Tf,const int & f);

  
  static void TestPartialDeltaVfPartialTf(std::string input,double tol,bool shape_in_meters);
  static void TestPartialDeltaIOverDeltaVPartialTf(std::string input,double tol,bool shape_in_meters);
  static void TestPartialDeltaIfPartialTf(std::string input,double tol,bool shape_in_meters);
  static void TestGetPartialAllInertiaPartialCVSStandalone(std::string input,double tol,bool shape_in_meters);
  static void TestPartialDeltaVPartialC(std::string input,double tol,bool shape_in_meters);
  static void TestGetPartialVolumePartialC(std::string input,double tol,bool shape_in_meters);
  static void TestGetPartialComPartialC(std::string input,double tol,bool shape_in_meters);
  static void TestGetPartialIPartialC(std::string input,double tol,bool shape_in_meters);
  static void TestGetPartialAllInertiaPartialC(std::string input,double tol,bool shape_in_meters) ;
  static void TestPartialEqDeltaIfErPartialTf(std::string input,double tol,bool shape_in_meters);
  static void TestGetPartialSigmaPartialC(std::string input,double tol,bool shape_in_meters);

  

};

#endif


