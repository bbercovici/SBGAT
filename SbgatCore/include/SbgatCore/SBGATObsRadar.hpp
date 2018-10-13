
/*=========================================================================

  Program:   Small Body Geophysical Analysis
  Module:    SBGATObsRadar.hpp

  Derived class from VTK's vtkPolyDataAlgorithm by Benjamin Bercovici  

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
@file SBGATObsRadar.hpp
@class  SBGATObsRadar
@author Benjamin Bercovici
@author Jay McMahon
@date October 2018
@brief  Computes range/range-rate Doppler images over the surface of
provided small body
@details Computes range/range-rate Doppler images over the surface of
provided small body
@copyright MIT License, Benjamin Bercovici and Jay McMahon

  
*/

#ifndef SBGATObsRadar_H
#define SBGATObsRadar_H

#include <vtkFiltersCoreModule.h> // For export macro
#include <vtkPolyDataAlgorithm.h>

#include <vtkModifiedBSPTree.h>
#include <vtkImageData.h>

#include <armadillo>
#include <array>

#include <SBGATObs.hpp>


/**
Vector of vector of radar observations. Each observation is comprised of (range,range-rate,incidence)
This incidence is not strictly speaking a measured quantity but is used to penalize the binned 
observations, blurring out returns collected at a high incidence
*/
typedef typename std::vector<std::vector<std::array<double, 3> > > SBGATRadarObsSequence;


class VTKFILTERSCORE_EXPORT SBGATObsRadar : public SBGATObs{
public:
  /**
   * Constructs with initial values of zero.
   */
  static SBGATObsRadar *New();

  vtkTypeMacro(SBGATObsRadar,vtkPolyDataAlgorithm);
  void PrintSelf(std::ostream& os, vtkIndent indent) override;
  void PrintHeader(std::ostream& os, vtkIndent indent) override;
  void PrintTrailer(std::ostream& os, vtkIndent indent) override;

  /**
  Collects range/range-rate samples over the surface of the specified system of small bodies at a given time after the epoch. The radar source is positionned at 1e6 * l meters from the target's center of mass, where 
  l is a measure of the first object's diagonal. 
  @param measurements_sequence reference to MeasurementsSequence, holding collected range/range-rate measurements at each observation time
  @param time observation timestamp
  @param N minimum number of measurements to produce over the smallest facet in the shape. The number of samples for any other facet will be equal to N * facet_surface_area / smallest_facet_surface_area
  @param radar_dir unit direction towards radar from target in inertial frame
  @param positions_vec vector of positions of each target's center-of-mass expressed in the primary's body frame
  @param velocities_vec vector of inertial velocities of each target's center-of-mass expressed in the primary's body frame
  @param mrps_vec vector of MRPs defining the inertial-to-body DCM [BN]
  @param omegas_vec vector of angular velocities of each body expressed in the inertial frame
  @param penalize_incidence if true, each measurement will be weighed by the cos(incidence) angle between
  the sampled point and the radar squared. If false, all measurements (in view of the radar and not blocked) are weighed equally
  have the same weight
  */
  void CollectMeasurements(SBGATRadarObsSequence & measurements_sequence,  
  const double & time,
  const int & N,
  const arma::vec & radar_dir,
  const std::vector<arma::vec> & positions_vec,
  const std::vector<arma::vec> & velocities_vec,
  const std::vector<arma::vec> & mrps_vec,
  const std::vector<arma::vec> & omegas_vec,
  const bool & penalize_indicence);




  /**
  Bins the provided measurements sequence into a series of 2d-histogram
  Will throw an std::runtime_error exception if 
  - either of the provided bin sizes are invalid (i.e <= 0)
  - the provided bin sizes yield empty histogram dimensions
  @param measurements_sequence reference to MeasurementsSequence, holding collected range/range-rate measurements at each observation time
  @param r_bin range bin size (m)
  @param rr_bin range-rate bin size (m/s)
  */
  void BinObservations(
    const SBGATRadarObsSequence & measurements_sequence,
    const double & r_bin,
    const double & rr_bin);

  /**
  Save the binned radar images to PNGs in the prescribed folder.
  The images will be normalized by the largest value in the observation sequence
  @param savepath path to folder where images will be saved (ex: "output/")
  */
  void SaveImages(std::string savepath);

  
  /**
  Return vector holding the computed images
  @return vector of radar images
  */
  std::vector<vtkSmartPointer<vtkImageData>> GetImages() const;

  /**
  Clears images, if any
  */

  void ClearImages() { this -> images.clear();}


protected:
  SBGATObsRadar();
  ~SBGATObsRadar() override;



  /**
  Ray traces the facets in view to the radar and measurements sequence with range/range-rate data
  if sampled point was in view of the radar
  @param measurements_sequence Sequence of measurements to add new image data to
  @param facets_in_view reference to a vector of vector holding indices of (maybe) illuminated facets for all considered bodies
  @param radar_dir radar direction expressed in inertial frame
  @param BN_dcms_vec vector holding the DCMs orienting the body frame of each body w/r to inertial
  @param positions_vec vector holding the position vector of the CM of each body w/r to the primary
  @param velocities_vec vector holding the velocities vector of the CM of each body
  @param omega_vec vector holding the angular velocities vector of each body
  */
  void reverse_ray_trace(SBGATRadarObsSequence & measurements_sequence,
    const std::vector<std::vector<int> > & facets_in_view,
    const arma::vec & radar_dir,
    const int N,
    const bool penalize_indicence,
    const std::vector<arma::mat> & BN_dcms_vec,
    const std::vector<arma::vec> & positions_vec,
    const std::vector<arma::vec> & velocities_vec,
    const std::vector<arma::vec> & omega_vec);


  std::vector<vtkSmartPointer<vtkImageData>> images;
  
  double max_value;

  

private:
  SBGATObsRadar(const SBGATObsRadar&) = delete;
  void operator=(const SBGATObsRadar&) = delete;
};

#endif


