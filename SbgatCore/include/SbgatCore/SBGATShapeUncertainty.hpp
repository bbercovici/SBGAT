/*=========================================================================

  Program:   Small Body Geophysical Analysis
  Module:    SBGATShapeUncertainty.hpp

  Class derived from VTK's vtkPolyDataAlgorithm by Benjamin Bercovici  

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
\file SBGATShapeUncertainty.hpp
\class  SBGATShapeUncertainty
\author Benjamin Bercovici 
\author Jay McMahon
\brief  Computes statistics in inertia properties of considered shape
\details Computes the standard deviation and covariances in volume, center of mass, 
inertia tensor, inertia moments, principal dimensions and principal axes. Interfaces with ShapeUQLig.
Note that this class will always shift the origin of the provided shape's coordinates to the shape's barycenter
\copyright MIT License, Benjamin Bercovici and Jay McMahon
*/

#ifndef SBGATShapeUncertainty_h
#define SBGATShapeUncertainty_h

#include <vtkFiltersCoreModule.h> // For export macro
#include <vtkPolyDataAlgorithm.h>
#include <armadillo>


template <class PointType> class ShapeModelBezier;
class ControlPoint;

class VTKFILTERSCORE_EXPORT SBGATShapeUncertainty : public vtkPolyDataAlgorithm{
public:
  /**
   * Constructs with initial values of zero.
   */
  static SBGATShapeUncertainty *New();

  vtkTypeMacro(SBGATShapeUncertainty,vtkPolyDataAlgorithm);
  void PrintSelf(std::ostream& os, vtkIndent indent) override;
  void PrintHeader(std::ostream& os, vtkIndent indent) override;
  void PrintTrailer(std::ostream& os, vtkIndent indent) override;

  /**
  Computes the inertia statistics of the provided shape using a linearized approach
  Refer to "Inertia of An Uncertain Small Body Shape" by Benjamin Bercovici and Jay McMahon, Icarus, 2019
  for further details.
  @param l correlation distance (L). Governs correlation decay
  @param sigma standard deviation of normal uncertainty (L)
  */
  void ComputeInertiaStatistics(double l,double sigma);



  /**
  Computes the inertia statistics of the provided shape using a Monte-Carlo sampling of the control point deviations
  Refer to "Inertia of An Uncertain Small Body Shape" by Benjamin Bercovici and Jay McMahon, Icarus, 2019
  for further details
  @param N_samples number of samples in the monte carlo
  @param l correlation distance (L). Governs correlation decay
  @param sigma standard deviation of normal uncertainty (L)
  @param output_dir directory where to save "slices" taken from the different shapes generated in the monte carlo. Only 20 slices are generated.
  If argument is left to default "" string, then nothing is saved
  */
  void ComputeInertiaStatisticsMC(int N_samples,double l,double sigma,std::string output_dir = "");


  /**
  Return variance in volume
  @return variance in volume (L^6)
  */
  double GetVolumeVariance() {return this -> volume_variance;}


  /**
  Return covariance in center-of-mass
  @return covariance in center-of-mass (L^2)
  */
  const arma::mat::fixed<3,3> & GetCOMCovariance()  {return this -> center_of_mass_covariance;}


  /**
  Return covariance in parametrization of inertia parametrization covariance
  @return covariance in inertia parametrization covariance(-)
  */
  const arma::mat::fixed<6,6> & GetInertiaParametrizationCovariance()  {return this -> inertia_tensor_parametrization_covariance;}


  /**
  Return covariance in principal dimensions 
  @return covariance in principal dimensions (L^2)
  */
  const arma::mat::fixed<3,3> & GetPrincipalDimensionsCovariance() {return this -> principal_dimensions_covariance;}


  /** 
  Return covariance in principal moments (A,B,C,volume)
  @return covariance in principal moments (-)
  */
  const arma::mat::fixed<4,4> & GetPrincipalMomentsCovariance() {return this -> principal_moments_covariance;}


  /**
  Return covariance in MRP parametrization of principal axes orientation
  @return covariance in MRP parametrization
  */
  const arma::mat::fixed<3,3> & GetPrincipalAxesMRPCovariance() {return this -> mrp_covariance; }




   /**
  Return covariance in non-normalized inertia principal directions
  @return covariance in non-normalized inertia principal directions
  */
  const arma::mat::fixed<9,9> & GetEvectorsPrincipalAxesCovariance() {return this -> Evectors_covariance; }

  
   /**
  Return covariance in unit-norm inertia principal directions
  @return covariance in inertia principal directions
  */
  const arma::mat::fixed<9,9> & GetEigenvectorPrincipalAxesCovariance() {return this -> eigenvectors_covariance; }





protected:
  SBGATShapeUncertainty();
  ~SBGATShapeUncertainty() override;

  int RequestData(vtkInformation* request,vtkInformationVector** inputVector, vtkInformationVector* outputVector) override;

  double volume_variance;
  arma::mat::fixed<3,3> center_of_mass_covariance;
  arma::mat::fixed<6,6> inertia_tensor_parametrization_covariance;
  arma::mat::fixed<3,3> principal_dimensions_covariance;
  arma::mat::fixed<4,4> principal_moments_covariance;
  arma::mat::fixed<9,9> Evectors_covariance;
  arma::mat::fixed<9,9> eigenvectors_covariance;
  arma::mat::fixed<3,3> mrp_covariance;


  std::shared_ptr<ShapeModelBezier<ControlPoint>> shape_model_bezier;


private:
  SBGATShapeUncertainty(const SBGATShapeUncertainty&) = delete;
  void operator=(const SBGATShapeUncertainty&) = delete;
};

#endif


