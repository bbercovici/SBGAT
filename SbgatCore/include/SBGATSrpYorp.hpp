/*=========================================================================

  Program:   Small Body Geophysical Analysis
  Module:    SBGATSrpYorp.hpp

  Derived class from VTK's vtkPolyDataAlgorithm by Benjamin Bercovici  

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class  SBGATSrpYorp
 * @author Benjamin Bercovici
 * @brief  Computation of Fourier decomposition of force/torques caused by SRP over a shape model
 *
 * @details Wrapper around YORPLib to allow VTK's Polydata to be utilized as shape containers.
 see Scheeres, D. J. (2007). `The dynamical evolution of uniformly rotating asteroids subject to YORP`. Icarus, 188(2), 430–450. https://doi.org/10.1016/j.icarus.2006.12.015
 for more details
*/

#ifndef SBGATSrpYorp_HEADER
#define SBGATSrpYorp_HEADER

#include <vtkFiltersCoreModule.h> // For export macro
#include <vtkPolyDataAlgorithm.h>
#include <armadillo>

class VTKFILTERSCORE_EXPORT SBGATSrpYorp : public vtkPolyDataAlgorithm{
public:

  /**
   * Constructs with initial values of zero.
   */

  static SBGATSrpYorp *New();

  vtkTypeMacro(SBGATSrpYorp,vtkPolyDataAlgorithm);
  void PrintSelf(std::ostream& os, vtkIndent indent) override;
  void PrintHeader(std::ostream& os, vtkIndent indent) override;
  void PrintTrailer(std::ostream& os, vtkIndent indent) override;


  /**
  Setter to the longitude step. default is 1 deg
  @param lambdaDel longitude step (degrees)
  */
  void set_lambdaDel(double lambdaDel);

  /**
  Setter to the declination step. default is 1 deg
  @param deltaDel declination step (degrees)
  */
  void set_deltaDel(double deltaDel);

  /**
  Setter to the maximum order of the Fourier decomposition. default is 2
  @param maxFourier maximum decomposition order
  */  
  void set_maxFourier(int maxFourier);


  /**
  Setter to numrefine. default is 5
  @param numrefine
  */  
  void set_numrefine(int numrefine);

  

  /**
  Setter to number of voxel per axis. default is 40
  @param numVox voxel pex axis
  */  
  void set_numVox(int numVox);

  /**
  Setter to reflectivity
  @param rho
  */
  void set_rho(double rho);

   /**
  Setter to specular component
  @param spec
  */
  void set_spec(double spec);



  /**
  Setter to number of ray "bounces" to be accounted for. default is 3
  @param howManyBounces number of ray reflections to keep track of
  */  
  void set_howManyBounces(int howManyBounces);

  /**
  Setter to directory where results are stored
  @param outputFileBaseName save directory
  */  
  void set_outputFileBaseName(std::string outputFileBaseName);



protected:
  SBGATSrpYorp();
  ~SBGATSrpYorp() override;

  int RequestData(vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) override;

  std::string outputFileBaseName;
  int numVox;
  double rho, spec;
  double lambdaDel;
  double deltaDel;
  double maxFourier;
  int howManyBounces;
  int numrefine;

private:
  SBGATSrpYorp(const SBGATSrpYorp&) = delete;
  void operator=(const SBGATSrpYorp&) = delete;
};

#endif


