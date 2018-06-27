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
#include <vtkPolyDataAlgorithm.h>

#include <vtkModifiedBSPTree.h>
#include <vtkImageData.h>

#include <armadillo>
#include <array>




class VTKFILTERSCORE_EXPORT SBGATObsLightcurve : public vtkPolyDataAlgorithm{
public:

  /**
   * Constructs with initial values of zero.
   */
  static SBGATObsLightcurve *New();

  vtkTypeMacro(SBGATObsLightcurve,vtkPolyDataAlgorithm);
  void PrintSelf(std::ostream& os, vtkIndent indent) override;
  void PrintHeader(std::ostream& os, vtkIndent indent) override;
  void PrintTrailer(std::ostream& os, vtkIndent indent) override;

  /**
  Computes collected luminosity over the surface of the small body at the 
  specified time after epoch. Exposure is instantaneous. The luminosity is computed as the 
  number of sample points in view of both the sun and the observer at each time (the "hit count"), normalized by the largest hit count in the observation sequence.
  @param measurements reference to a vector of std::arrays holding (times,luminosity)
  @param N maximum number of measurements to produce over the largest facet in the shape. The number of samples for any other facet will be equal to N * facet_surface_area / larget_facet_surface_area
  @param dt time since epoch (s)
  @param period rotation period of target (s)
  @param sun_pos unit direction of sun with respect to target in inertial frame
  @param observer_pos unit direction of observer with respect to target in inertial frame
  @param spin (unit vector) direction of target's spin vector expressed in the target's body frame
  @param penalize_incidence if true, each measurement will be weighed by the cos(incidence) angle between
  the sampled point and the observer TIMES the cos(incidence) the sampled point and the sun. If false, all accepted measurements (in view of the observer and not blocked) 
  */
  void CollectMeasurementsSimpleSpin(
    std::vector<std::array<double, 2> > & measurements,
    const int & N,
    const double & dt,
    const double & period,
    const arma::vec & sun_dir,
    const arma::vec & observer_dir,
    const arma::vec & spin,
    const bool & penalize_indicence);


  /**
  Sets the scale factor to 1, indicative that the polydata has its coordinates expressed in meters (default)
  */
  void SetScaleMeters() { this -> scaleFactor = 1; }

  /**
  Sets the scale factor to 1000, indicative that the polydata has its coordinates expressed in kilometers
  */
  void SetScaleKiloMeters() { this -> scaleFactor = 1000; }


  /**
  Save the raw data points from the unnormalized (times,luminosity) time series to file
  @param measurements reference to a vector of std::arrays holding (times,luminosity)
  @param savepath path to file where raw lightcurve data will be saved (ex: "output/lightcurve.txt")
  */
  void SaveLightCurveData(const std::vector<std::array<double, 2> > & measurements, std::string savepath);


protected:
  SBGATObsLightcurve();
  ~SBGATObsLightcurve() override;

  int RequestData(vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) override;

  /**
  Determines what facets may be illuminated by the sun. Does not check for interests, 
  only compares the outbound normal orientation to the sun/facet direction
  @param body_index index of considered body
  @param facets_in_view reference to vector holding indices of (maybe) illuminated facets
  */
  void check_facet_illumination(const unsigned int & body_index,
    std::vector<int> & facets_in_view,
    const arma::vec & target_to_sun_dir_body_frame,
    const arma::vec & target_to_observer_dir_body_frame,
    double & max_area);



  /**
  Ray traces the facets in view to the sun/observer and increment measurement
  counter if in view
  @param measurements_temp reference to an std::array holding (times,luminosity)
  @param body_index index of considered body
  @param facets_in_view reference to vector holding indices of (maybe) illuminated facets
  */
  void reverse_ray_trace(std::array<double, 2> & measurements,
    const unsigned int & body_index,
    const std::vector<int> & facets_in_view,
    const arma::vec & target_to_sun_dir_body_frame,
    const arma::vec & target_to_observer_dir_body_frame,
    const int N,
    const double max_area,
    const bool penalize_indicence);


  std::vector<vtkSmartPointer<vtkModifiedBSPTree>> bspTree_vec;
  
  std::vector<arma::vec> cm_vec;

  double scaleFactor = 1;
  double max_value;

private:
  SBGATObsLightcurve(const SBGATObsLightcurve&) = delete;
  void operator=(const SBGATObsLightcurve&) = delete;
};

#endif


