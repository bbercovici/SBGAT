/*=========================================================================

  Program:   Small Body Geophysical Analysis
  Module:    SBGATMassProperties.hpp

  Class derived from VTK's vtkPolyDataAlgorithm by Benjamin Bercovici  

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
\file SBGATMassProperties.hpp
\class  SBGATMassProperties
\author Benjamin Bercovici 
\author Jay McMahon
\brief  Computes volume, area, shape index, center of mass,
inertia tensor and principal axes of a polyhedral mesh of constant density
\details Computes the volume, the surface area, and the
normalized shape index, center of mass and inertia tensor of a topologically-closed, constant-density polyhedron.
This class will always use results expressed in `meters` as their distance unit (e.g center-of-mass coordinates in meters, volume in m^3,...) . Unit consistency is enforced through the use of the SetScaleMeters()
and SetScaleKiloMeters() method. 

See "Inertia of Any Polyhedron" by Anthony R. Dobrovolskis, Icarus 124, 698–704 (1996) Article No. 0243
for further details.  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
\copyright MIT License, Benjamin Bercovici and Jay McMahon
*/

#ifndef SBGATMassProperties_h
#define SBGATMassProperties_h

#include <vtkFiltersCoreModule.h> // For export macro
#include <vtkPolyDataAlgorithm.h>
#include <armadillo>
#include <SBGATFilter.hpp>

class VTKFILTERSCORE_EXPORT SBGATMassProperties : public SBGATFilter{
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
   * Compute and return the volume (m^3)
   */
  double GetVolume() { return this->Volume;}

  /**
   * Compute and return the projected volume.
   * Typically you should compare this volume to the value returned by GetVolume
   * if you get an error (GetVolume()-GetVolumeProjected())*10000 that is greater
   * than GetVolume() this should identify a problem:
   * * Either the polydata is not closed
   * * Or the polydata contains triangle that are flipped
   */
  double GetVolumeProjected() { return this->VolumeProjected;}

  /**
   * Compute and return the volume projected on to each axis aligned plane.
   */
  double GetVolumeX() { return this->VolumeX;}
  double GetVolumeY() { return this->VolumeY;}
  double GetVolumeZ() { return this->VolumeZ;}

  /**
   * Compute and return the weighting factors for the maximum unit
   * normal component (MUNC).
   */
  double GetKx() { return this->Kx;}
  double GetKy() { return this->Ky;}
  double GetKz() { return this->Kz;}

  /**
   * Compute and return the area in m^2
   */
  double GetSurfaceArea() { return this->SurfaceArea;
  }

  /**
   * Compute and return the min cell area in m^2
   */
  double GetMinCellArea() { return this->MinCellArea;
  }

  /**
   * Compute and return the max cell area in m^2
   */
  double GetMaxCellArea() { return this->MaxCellArea;
  }


  /**
  Checks whether the polydata is topologically closed or open
  If closed, the sum of the oriented surface area should be equal to zero
  */
  bool CheckClosed(){ return this -> IsClosed;}

  /**
   * Compute and return the normalized shape index. This characterizes the
   * deviation of the shape of an object from a sphere. A sphere's NSI
   * is one. This number is always >= 1.0.
   */
  double GetNormalizedShapeIndex()
  { return this->NormalizedShapeIndex;
  }


  /**
  * Compute and return the coordinates of the center of mass (m)
  * evaluated in the frame of origin assuming a constant density distribution
  * across the shape
  */
  const arma::vec::fixed<3> & GetCenterOfMass(){
   return this -> center_of_mass;
 }

  /**
  * Compute and return the coordinates of the center of mass (m)
  * evaluated in the frame of origin assuming a constant density distribution
  * across the shape
  */
 void GetCenterOfMass(double * com){

  com[0] = this -> center_of_mass(0);
  com[1] = this -> center_of_mass(1);
  com[2] = this -> center_of_mass(2);
}

  /**
  * Compute and return the dimensionless inertia tensor
  * evaluated in the frame of origin assuming a constant density distribution
  * across the shape. The normalization applied to the inertia tensor is I_norm = I / (mass * r_avg ^ 2) where r_avg = cbrt(3/4*Volume/pi)
  */
arma::mat::fixed<3,3> GetNormalizedInertiaTensor(){
 return this -> inertia_tensor;
}

  /**
  * Compute and return the dimensionless inertia tensor
  * evaluated in the frame of origin assuming a constant density distribution
  * across the shape. The normalization applied to the inertia tensor is I_norm = I / (rho) where rho is the density
  */
arma::mat::fixed<3,3> GetUnitDensityInertiaTensor(){
 return unit_density_inertia_tensor;
}


  /**
  * Compute and return the dcm orienting the principal axes of the small body relative to 
  the body coordinates frame. That is, denoting P the principal frame and B the frame in which the
  coordinates of the body are currently expressed, this method returns [PB]
  @return [PB] direction cosine matrix
  */
arma::mat::fixed<3,3> GetPrincipalAxes(){
 return this -> principal_axes;
}

  /**
  Computes and returns the principal dimensions (m) of the ellipsoid associated with the inertia tensor 
  tensor, sorted from the longest (smallest inertia) to shortest (largest inertia)
  @return principal dimensions associated with inertia tensor (m)
  */
arma::vec::fixed<3> GetPrincipalDimensions(){
 return this -> principal_dimensions;
}


  /**
  * Compute and return the normalized inertia moments assuming uniform density distribution
  * across the shape, sorted from the smallest inertia to the largest.
  * The normalization applied to the inertia tensor is I_norm = I / (mass * r_avg ^ 2) where r_avg = cbrt(3/4*Volume/pi)
  */
arma::vec::fixed<3> GetNormalizedInertiaMoments(){
 return normalized_principal_moments;
}

  /**
  * Compute and return the inertia moments assuming uniform unit density distribution
  * across the shape, sorted from the smallest inertia to the largest.
  */
arma::vec::fixed<3> GetUnitDensityInertiaMoments(){
 return unit_density_principal_moments;
}

  /**
  Returns the average radius of the shape (that is, the radius of a sphere occupying the same volume) (m)
  */
double GetAverageRadius(){
 return this -> r_avg;
}




    /**
    Computes the mass properties of the provided shape and saves the results to a JSON file
    @param shape pointer to considered shape
    @param path savepath (ex: "mass_properties.json")
    */
static void ComputeAndSaveMassProperties(vtkSmartPointer<vtkPolyData> shape,std::string path);


  /**
  Save the computed mass properties to a JSON file
  @param path savepath (ex: "mass_properties.json")

  */
void SaveMassProperties(std::string path) const ;

/**
  Returns signed contribution to total volume of tetrahedron subtended by facet
  f (m^3)
  @param f facet index
  @return signed volume of tetrahedron subtended by facet
  */
double GetDeltaV(const int & f) const;

  /**
  Returns coordinates of the tetrahedron's center-of-mass (m)
  @param f facet index
  @return coordinates of tetrahedron
  */
arma::vec::fixed<3> GetDeltaCM(const int & f) const ;

  /**
  Returns the unit-density tetrahedron's inertia tensor divided by tetrahedron's signed volume (m^2)
  @param f facet index
  @return tetrahedron's inertia tensor divided by tetrahedron's signed volume
  */
arma::mat::fixed<3,3> GetDeltaIOverDeltaV(const int & f) const ;

/**
Returns the parametrization of the the unit-density tetrahedron's inertia tensor
@param f facet index
@return parametrization of the tetrahedron's inertia tensor
*/
arma::vec::fixed<6> GetDeltaIf(const int & f) const;


protected:
  SBGATMassProperties();
  ~SBGATMassProperties() override;

  int RequestData(vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) override;

  

  arma::vec::fixed<3> center_of_mass;
  arma::mat::fixed<3,3> inertia_tensor;
  arma::mat::fixed<3,3> principal_axes;

  arma::vec::fixed<3> normalized_principal_moments;
  arma::vec::fixed<3> unit_density_principal_moments;
  arma::mat::fixed<3,3> unit_density_inertia_tensor;
  arma::vec::fixed<3> principal_dimensions;


  double  SurfaceArea;
  double  MinCellArea;
  double  MaxCellArea;
  double  Volume;
  double  VolumeProjected; 
  double  VolumeX;
  double  VolumeY;
  double  VolumeZ;
  double  Kx;
  double  Ky;
  double  Kz;
  double  NormalizedShapeIndex;
  double r_avg;
  bool IsClosed;



private:
  SBGATMassProperties(const SBGATMassProperties&) = delete;
  void operator=(const SBGATMassProperties&) = delete;
};

#endif


