/** MIT License

Copyright (c) 2018 Benjamin Bercovici and Jay McMahon

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/*=========================================================================

  Program:   Small Body Geophysical Analysis
  Module:    SBGATObsLightcurve.hpp

  Derived class from VTK's vtkPolyDataAlgorithm by Benjamin Bercovici  

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class  SBGATObsRadar
 * @author Benjamin Bercovici
 * @brief  Computes range/range-rate Doppler images over the surface of
 provided small body
 *
 * @details Computes range/range-rate Doppler images over the surface of
 provided small body
 *
*/

#ifndef SBGATObsLightcurve_H
#define SBGATObsLightcurve_H

#include <vtkFiltersCoreModule.h> // For export macro

#include <vtkModifiedBSPTree.h>
#include <vtkImageData.h>

#include <armadillo>
#include <array>

#include <SBGATObs.hpp>


class VTKFILTERSCORE_EXPORT SBGATObsLightcurve : public SBGATObs{
public:


  static SBGATObsLightcurve *New();

  vtkTypeMacro(SBGATObsLightcurve,SBGATObs);
  void PrintSelf(std::ostream& os, vtkIndent indent) override;
  void PrintHeader(std::ostream& os, vtkIndent indent) override;
  void PrintTrailer(std::ostream& os, vtkIndent indent) override;

  /**
  Computes collected luminosity over the surface of the considered small bodies at the 
  specified time after epoch. Exposure is instantaneous. The luminosity is computed as the 
  number of sample points in view of both the sun and the observer at each time (the "hit count"), normalized by the largest hit count in the observation sequence.
  The positions of each body is specified. The attitude of each body is supposed to derive from a constant-spin rotation such that
  each of the dcms BN is equal to identity when dt == 0.
  @param measurements reference to a vector of std::arrays holding (times,luminosity)
  @param time measurement time
  @param N minimum number of measurements to produce over the smallest facet in the shape. The number of samples for any other facet will be equal to N * facet_surface_area / smallest_facet_surface_area
  @param sun_dir unit direction towards sun from target in inertial frame
  @param observer_dir unit direction towards observer from target in inertial frame
  @param positions_vec vector of positions of each target's center-of-mass expressed in the primary's body frame
  @param velocities_vec vector of inertial velocities of each target's center-of-mass expressed in the primary's body frame
  @param mrps_vec vector of MRPs defining the inertial-to-body DCM [BN]
  @param omegas_vec vector of angular velocities of each body expressed in the inertial frame
  @param penalize_incidence if true, each measurement will be weighed by the cos(incidence) angle between
  the sampled point and the radar squared. If false, all measurements (in view of the radar and not blocked) are weighed equally
  have the same weight

  */
  void CollectMeasurements(
    std::vector<std::array<double, 2> > & measurements,
    const double & time,
    const int & N,
    const arma::vec & sun_dir,
    const arma::vec & observer_dir,
    const std::vector<arma::vec> & positions_vec,
    const std::vector<arma::vec> & velocities_vec,
    const std::vector<arma::vec> & mrps_vec,
    const std::vector<arma::vec> & omegas_vec,
    const bool & penalize_indicence);


  /**
  Save the raw data points from the unnormalized (times,luminosity) time series to file
  @param measurements reference to a vector of std::arrays holding (times,luminosity)
  @param savepath path to file where raw lightcurve data will be saved (ex: "output/lightcurve.txt")
  */
  void SaveLightCurveData(const std::vector<std::array<double, 2> > & measurements, std::string savepath);


protected:
  SBGATObsLightcurve();
  ~SBGATObsLightcurve() override;


  /**
  Ray traces all of the facets in view to the sun/observer and increment measurement
  counter if in view
  @param measurements_temp reference to an std::array holding (times,luminosity)
  @param facets_in_view reference to a vector of vector holding indices of (maybe) illuminated facets for all considered bodies
  @param sun_dir sun direction expressed in inertial frame
  @param observer_dir observer direction expressed in inertial frame
  @param BN_dcms_vec vector holding the DCMs orienting the body frame of each body w/r to inertial
  @param positions_vec vector holding the position vector of the CM of each body w/r to the primary
  */
  void reverse_ray_trace(std::array<double, 2>  & measurements_temp,
    const std::vector<std::vector<int> > & facets_in_view,
    const arma::vec & sun_dir,
    const arma::vec & observer_dir,
    const int N,
    const bool penalize_indicence,
    const std::vector<arma::mat> & BN_dcms_vec,
    const std::vector<arma::vec> & positions_vec);

  

private:
  SBGATObsLightcurve(const SBGATObsLightcurve&) = delete;
  void operator=(const SBGATObsLightcurve&) = delete;
};

#endif


