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
 * @details to be completed
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
  Setter to the longitude step
  @param lambdaDel longitude step (degrees)
  */
  void set_lambdaDel(double lambdaDel);

  /**
  Setter to the declination step
  @param deltaDel declination step (degrees)
  */
  void set_deltaDel(double deltaDel);

  /**
  Setter to the maximum order of the Fourier decomposition
  @param maxFourier maximum decomposition order
  */  
  void set_maxFourier(int maxFourier);


  /**
  Setter to numrefine
  @param numrefine
  */  
  void set_numrefine(int numrefine);

  

  /**
  Setter voxel grid size
  @param maxFourier maximum decomposition order
  */  
  void set_numVox(int numVox);



  /**
  Setter to number of ray "bounces" to be accounted for 
  @param howManyBounces number of ray reflections to keep track of
  */  
  void set_howManyBounces(int howManyBounces);

  /**
  Setter to path where results are stored
  @param outputFileBaseName save path
  */  
  void set_outputFileBaseName(std::string outputFileBaseName);



protected:
  SBGATSrpYorp();
  ~SBGATSrpYorp() override;

  int RequestData(vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) override;

  std::string objFileName, opticalFileName, outputFileBaseName, inFileName;
  int numVox, opticalFlag;
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


