/*=========================================================================

  Program:   Small Body Geophysical Analysis
  Module:    SBGATChildTest.hpp

  Class derived from VTK's vtkPolyDataAlgorithm by Benjamin Bercovici  

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
\file SBGATChildTest.hpp
\class  SBGATChildTest
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

#ifndef SBGATChildTest_h
#define SBGATChildTest_h


#include <vtkFiltersCoreModule.h> // For export macro
#include <vtkPolyDataAlgorithm.h>
#include "SBGATFilter.hpp"

#include <armadillo>

class VTKFILTERSCORE_EXPORT SBGATChildTest : public SBGATFilter{
public:
  /**
   * Constructs with initial values of zero.
   */
  static SBGATChildTest *New();

  vtkTypeMacro(SBGATChildTest,vtkPolyDataAlgorithm);
  void PrintSelf(std::ostream& os, vtkIndent indent) override;
  void PrintHeader(std::ostream& os, vtkIndent indent) override;
  void PrintTrailer(std::ostream& os, vtkIndent indent) override;

 
protected:
  SBGATChildTest();
  ~SBGATChildTest() override;

  int RequestData(vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector);


private:
  SBGATChildTest(const SBGATChildTest&) = delete;
  void operator=(const SBGATChildTest&) = delete;
};

#endif


