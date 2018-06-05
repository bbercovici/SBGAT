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
  specified time after epoch.
  @param measurements reference to a vector of std::arrays holding (times,luminosity)
  @param N maximum number of measurements to produce per illuminated facet
  @param dt time since epoch (s)
  @param period rotation period of target (s)
  @param target_pos position of target with respect to sun in inertial frame
  @param observer_pos position of observer with respect to sun in inertial frame
  @param spin (unit vector) direction of target's spin vector expressed in the target's body frame
  */
  void CollectMeasurementsSimpleSpin(
    std::vector<std::array<double, 2> > & measurements,
    const int & N,
    const double & dt,
    const double & period,
    const arma::vec & target_pos,
    const arma::vec & observer_pos,
    const arma::vec & spin);


  /**
  Sets the scale factor to 1, indicative that the polydata has its coordinates expressed in meters (default)
  */
  void SetScaleMeters() { this -> scaleFactor = 1; }

  /**
  Sets the scale factor to 1000, indicative that the polydata has its coordinates expressed in kilometers
  */
  void SetScaleKiloMeters() { this -> scaleFactor = 1000; }



protected:
  SBGATObsLightcurve();
  ~SBGATObsLightcurve() override;

  int RequestData(vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) override;


  vtkSmartPointer<vtkModifiedBSPTree> bspTree;
  
  arma::vec center_of_mass;

  double scaleFactor = 1;
  double max_value;

private:
  SBGATObsLightcurve(const SBGATObsLightcurve&) = delete;
  void operator=(const SBGATObsLightcurve&) = delete;
};

#endif


