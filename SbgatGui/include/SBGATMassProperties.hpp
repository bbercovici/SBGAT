/*=========================================================================

  Program:   Small Body Geophysical Analysis
  Module:    SBGATMassProperties.h

  Derived class from VTK's vtkPolyDataAlgorithm by Benjamin Bercovici  

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class  SBGATMassProperties
 * @author Benjamin Bercovici
 * @brief  compute volume, area, shape index, center of mass,
 * inertia tensor and principal axes of a polyhedral mesh of constant density
 *
 * SBGATMassProperties estimates the volume, the surface area, and the
 * normalized shape index, center of mass (assuming constant density) and inertia tensor
 * of a triangle mesh. See "Inertia of Any Polyhedron" by Anthony R. Dobrovolskis, Icarus 124, 698–704 (1996) Article No. 0243
 * for further details
 *
*/

#ifndef SBGATMassProperties_h
#define SBGATMassProperties_h

#include <vtkFiltersCoreModule.h> // For export macro
#include <vtkPolyDataAlgorithm.h>
#include <armadillo>

class VTKFILTERSCORE_EXPORT SBGATMassProperties : public vtkPolyDataAlgorithm{
public:
  /**
   * Constructs with initial values of zero.
   */
  static SBGATMassProperties *New();

  vtkTypeMacro(SBGATMassProperties,vtkPolyDataAlgorithm);
  void PrintSelf(std::ostream& os, vtkIndent indent) override;
  void PrintHeader(std::ostream& os, vtkIndent indent) override;
  void PrintTrailer(std::ostream& os, vtkIndent indent) override;

  /**
   * Compute and return the volume.
   */
  double GetVolume() {this->Update(); return this->Volume;}

  /**
   * Compute and return the projected volume.
   * Typically you should compare this volume to the value returned by GetVolume
   * if you get an error (GetVolume()-GetVolumeProjected())*10000 that is greater
   * than GetVolume() this should identify a problem:
   * * Either the polydata is not closed
   * * Or the polydata contains triangle that are flipped
   */
  double GetVolumeProjected() {this->Update(); return this->VolumeProjected;}

  /**
   * Compute and return the volume projected on to each axis aligned plane.
   */
  double GetVolumeX() {this->Update(); return this->VolumeX;}
  double GetVolumeY() {this->Update(); return this->VolumeY;}
  double GetVolumeZ() {this->Update(); return this->VolumeZ;}

  /**
   * Compute and return the weighting factors for the maximum unit
   * normal component (MUNC).
   */
  double GetKx() {this->Update(); return this->Kx;}
  double GetKy() {this->Update(); return this->Ky;}
  double GetKz() {this->Update(); return this->Kz;}

  /**
   * Compute and return the area.
   */
  double GetSurfaceArea() {this->Update(); return this->SurfaceArea;}

  /**
   * Compute and return the min cell area.
   */
  double GetMinCellArea() {this->Update(); return this->MinCellArea;}

  /**
   * Compute and return the max cell area.
   */
  double GetMaxCellArea() {this->Update(); return this->MaxCellArea;}

  /**
   * Compute and return the normalized shape index. This characterizes the
   * deviation of the shape of an object from a sphere. A sphere's NSI
   * is one. This number is always >= 1.0.
   */
  double GetNormalizedShapeIndex()
  {this->Update(); return this->NormalizedShapeIndex;}


  /**
  * Compute and return the coordinates of the center of mass
  * evaluated in the frame of origin assuming a constant density distribution
  * across the shape
  */
  arma::vec::fixed<3> GetCenterOfMass(){
    this -> Update(); return this -> center_of_mass;}

    /**
  * Compute and return the coordinates of the center of mass
  * evaluated in the frame of origin assuming a constant density distribution
  * across the shape
  */
    void GetCenterOfMass(double * com){
      this -> Update(); 
      com[0] = this -> center_of_mass(0);
      com[1] = this -> center_of_mass(1);
      com[2] = this -> center_of_mass(2);
    }


  /**
  * Compute and return the inertia tensor
  * at the barycenter, evaluated in the frame of origin assuming a constant density distribution
  * across the shape
  */
    arma::mat::fixed<3,3> GetInertiaTensor(){
      this -> Update(); return this -> inertia_tensor;}

  /**
  * Compute and return the principal axes of the inertia tensor
  */
      arma::mat::fixed<3,3> GetPrincipalAxes(){
        this -> Update(); return this -> principal_axes;}


  /**
  * Compute and return the inertia moments assuming a constant density distribution
  * across the shape
  */
        arma::vec::fixed<3> GetInertiaMoments(){
          this -> Update(); return arma::eig_sym(this -> inertia_tensor);}


  /**
  * Compute and return the bounding box (xmin,xmax,ymin,ymax,zmin,zmax)
  */
          double * GetBoundingBox(){
            this -> Update(); return this -> bounds;
          }

        protected:
          SBGATMassProperties();
          ~SBGATMassProperties() override;

          int RequestData(vtkInformation* request,
            vtkInformationVector** inputVector,
            vtkInformationVector* outputVector) override;


          arma::vec::fixed<3> center_of_mass;
          arma::mat::fixed<3,3> inertia_tensor;
          arma::mat::fixed<3,3> principal_axes;

          double  SurfaceArea;
          double  MinCellArea;
          double  MaxCellArea;
          double  Volume;
  double  VolumeProjected; // == Projected area of triangles * average z values
  double  VolumeX;
  double  VolumeY;
  double  VolumeZ;
  double  Kx;
  double  Ky;
  double  Kz;
  double  NormalizedShapeIndex;
  double bounds[6];
private:
  SBGATMassProperties(const SBGATMassProperties&) = delete;
  void operator=(const SBGATMassProperties&) = delete;
};

#endif


