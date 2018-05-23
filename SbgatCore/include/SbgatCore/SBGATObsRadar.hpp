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
 * @class  SBGATObsRadar
 * @author Benjamin Bercovici
 * @brief  Computes range/range-rate Doppler images over the surface of
 provided small body
 *
 * @details Computes range/range-rate Doppler images over the surface of
 provided small body
 *
*/

#ifndef SBGATObsRadar_H
#define SBGATObsRadar_H

#include <vtkFiltersCoreModule.h> // For export macro
#include <vtkPolyDataAlgorithm.h>
#include <armadillo>
#include <array>

#include <vtkModifiedBSPTree.h>

class VTKFILTERSCORE_EXPORT SBGATObsRadar : public vtkPolyDataAlgorithm{
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
  Collects range/range-rate samples over the surface of the small body at the 
  specified time after epoch.
  The radar source is positionned at 1e6 * l meters from the target's center of mass, where 
  l is a measure of the object's diagonal
  @param measurements reference to std::vector holding collected range/range-rate measurements
  @param N maximum number of measurements to produce per illuminated facet
  @param dt time since epoch (hours)
  @param period rotation period of target (hours)
  @param dir (unit vector) direction of radar-to-target vector expressed in the target's body frame at epoch
  @param spin (unit vector) direction of target's spin vector expressed in the target's body frame
  */
  void CollectMeasurementsSimpleSpin(
    std::vector<std::array<double, 2> > & measurements,
    const int & N,
    const double & dt,
    const double & period,
    const arma::vec & dir,
    const arma::vec & spin);

  /**
  Collects range/range-rate samples over the surface of the small body at the 
  specified time after epoch
  @param measurements reference to std::vector holding collected range/range-rate measurements
  @param r_bin range bin size (m)
  @param rr_bin range-rate bin size (m/s)
  @param savepath path where resulting image is to be saved
  */
  void SaveImage(
    std::vector<std::array<double, 2> > & measurements,
    const double & r_bin,
    const double & rr_bin,
    std::string savepath);

  /**
  Sets the scale factor to 1, indicative that the polydata has its coordinates expressed in meters (default)
  */
  void SetScaleMeters() { this -> scaleFactor = 1; }

  /**
  Sets the scale factor to 1000, indicative that the polydata has its coordinates expressed in kilometers
  */
  void SetScaleKiloMeters() { this -> scaleFactor = 1000; }
 

protected:
  SBGATObsRadar();
  ~SBGATObsRadar() override;

  int RequestData(vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) override;


  vtkSmartPointer<vtkModifiedBSPTree> bspTree;
  arma::vec center_of_mass;

  double scaleFactor = 1;

private:
  SBGATObsRadar(const SBGATObsRadar&) = delete;
  void operator=(const SBGATObsRadar&) = delete;
};

#endif


