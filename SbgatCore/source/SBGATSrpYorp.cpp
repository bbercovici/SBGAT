/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SBGATSrpYorp.cxx

  Derived class from VTK's vtkPolyDataAlgorithm by Benjamin Bercovici  
  
  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "SBGATSrpYorp.hpp"

#include <vtkObjectFactory.h>
#include <vtkCell.h>
#include <vtkDataObject.h>
#include <vtkIdList.h>
#include <vtkMath.h>
#include <vtkSmartPointer.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <RigidBodyKinematics.hpp>


// Should prefix those or wrap them within a namespace
#include <Body.h>
#include <SRPModel.h>

vtkStandardNewMacro(SBGATSrpYorp);

//----------------------------------------------------------------------------
// Constructs with initial 0 values.
SBGATSrpYorp::SBGATSrpYorp(){

  this -> numVox = 40;
  this -> rho = 1;
  this -> spec = 0;
  this -> lambdaDel = 1;
  this -> deltaDel = 1;
  this -> maxFourier = 2;
  this -> howManyBounces = 3;
  this -> numrefine = 5;


  this->SetNumberOfOutputPorts(0);
}

//----------------------------------------------------------------------------
// Destroy any allocated memory.
SBGATSrpYorp::~SBGATSrpYorp(){
}




void SBGATSrpYorp::set_lambdaDel(double lambdaDel){
  this -> lambdaDel = lambdaDel;
}



void SBGATSrpYorp::set_deltaDel(double deltaDel){
  this -> deltaDel = deltaDel;
}


void SBGATSrpYorp::set_maxFourier(int maxFourier){
  this -> maxFourier = maxFourier;
}



void SBGATSrpYorp::set_howManyBounces(int howManyBounces){
  this -> howManyBounces = howManyBounces;
}


void SBGATSrpYorp::set_numrefine(int numrefine){
  this -> numrefine = numrefine;
}


void SBGATSrpYorp::set_outputFileBaseName(std::string outputFileBaseName){
  this -> outputFileBaseName = outputFileBaseName;
}


void SBGATSrpYorp::set_numVox(int numVox){
  this -> numVox = numVox;
}



void SBGATSrpYorp::set_rho(double rho){
  this -> rho = rho;
}

void SBGATSrpYorp::set_spec(double spec){
  this -> rho = spec;
}

//----------------------------------------------------------------------------
// Description:
// This method measures volume, surface area, normalized shape index, the surface area, the
//  normalized shape index, center of mass (assuming constant density) and inertia tensor
// (assuming constant density) of a triangle mesh.

int SBGATSrpYorp::RequestData(
  vtkInformation* vtkNotUsed( request ),
  vtkInformationVector** inputVector,
  vtkInformationVector* vtkNotUsed( outputVector )){
  vtkInformation *inInfo =
  inputVector[0]->GetInformationObject(0);


  if (outputFileBaseName.size() == 0){
    throw std::runtime_error("no output file was specified");
  }

  // call ExecuteData
  vtkPolyData * input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkIdType cellId, numCells, numPts, numIds;

  numCells = input->GetNumberOfCells();
  numPts = input->GetNumberOfPoints();
  if (numCells < 1 || numPts < 1){
    vtkErrorMacro( << "No data to measure...!");
    return 1;
  }

  vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
  ptIds->Allocate(VTK_CELL_SIZE);

  // The facets are browsed to created the facet dyads
  for (cellId=0; cellId < numCells; cellId++){
    if ( input->GetCellType(cellId) != VTK_TRIANGLE){
      vtkWarningMacro(<< "Input data type must be VTK_TRIANGLE not "
        << input->GetCellType(cellId));
      continue;
    }

    input->GetCellPoints(cellId,ptIds);
    numIds = ptIds->GetNumberOfIds();
    assert(numIds == 3);

  } 
  
  std::vector<std::vector<double> > vertices(numPts);
  std::vector<std::vector<int> > facets(numCells);


  // The vertex coordinates are extracted
  #pragma omp parallel for
  for(int i = 0; i < numPts; ++i) {

    double verts[3];
    input -> GetPoint(i,verts);
    std::vector<double> coords = {verts[0], verts[1], verts[2]};
    vertices[i] = coords;
  }

  #pragma omp parallel for
  for(int i = 0; i < numCells; ++i) {
    
    vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
    ptIds -> Allocate(VTK_CELL_SIZE);

    input -> GetCellPoints(i,ptIds);

    std::vector<int> facet = {ptIds -> GetId(0),ptIds -> GetId(1),ptIds -> GetId(2)};

    facets[i] = facet;

  }


  Body targetObj(vertices,facets,this -> rho,this -> spec);

  // Determine limits for voxel grid
  double xmax = 1.01 * targetObj.getMaxDim(1);
  double ymax = 1.01 * targetObj.getMaxDim(2);
  double zmax = 1.01 * targetObj.getMaxDim(3);

  // Fill out Voxel Grid
  targetObj.setVoxelGrid(xmax, ymax, zmax, this -> numVox);

  // Compute and save SRP Fourier coefficients
  SRPModel targetSRP(this -> lambdaDel, 
    this -> deltaDel, 
    this -> maxFourier, 
    &targetObj, 
    this -> howManyBounces, 
    this -> numrefine);

  targetSRP.writeSRPCoeffsFile(this -> outputFileBaseName, 2*(90.0/this -> deltaDel) + 1);


  return 1;
}

void SBGATSrpYorp::PrintHeader(ostream& os, vtkIndent indent) {

}
void SBGATSrpYorp::PrintTrailer(ostream& os, vtkIndent indent) {

}


//----------------------------------------------------------------------------
void SBGATSrpYorp::PrintSelf(std::ostream& os, vtkIndent indent){


}
