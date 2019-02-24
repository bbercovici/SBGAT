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

  Class derived from VTK's vtkPolyDataAlgorithm by Benjamin Bercovici  

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "SBGATFilter.hpp"

#include <vtkObjectFactory.h>
#include <vtkCell.h>
#include <vtkDataObject.h>
#include <vtkIdList.h>
#include <vtkMath.h>
#include <vtkSmartPointer.h>
#include <vtkCleanPolyData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPolyDataNormals.h>
#include <vtkFloatArray.h>
#include <vtkExtractEdges.h>
#include <vtkLine.h>
#include <vtkCellData.h>


#include <json.hpp>
vtkStandardNewMacro(SBGATFilter);

//----------------------------------------------------------------------------
// Constructs with initial 0 values.
SBGATFilter::SBGATFilter(){

  this -> N_facets = 0;
  this -> N_edges = 0;
  this -> N_vertices = 0;
  

  this->SetNumberOfOutputPorts(0);
}

//----------------------------------------------------------------------------
// Destroy any allocated memory.
SBGATFilter::~SBGATFilter(){
  this -> Clear();

}

//----------------------------------------------------------------------------
// Description:
// This method measures volume, surface area, normalized shape index, the surface area, the
//  normalized shape index, center of mass (assuming constant density) and inertia tensor
// (assuming constant density) of a triangle mesh.

int SBGATFilter::RequestData(
  vtkInformation* vtkNotUsed( request ),
  vtkInformationVector** inputVector,
  vtkInformationVector* vtkNotUsed( outputVector )){
  vtkInformation *inInfo =
  inputVector[0]->GetInformationObject(0);

  // Any data previously owned is erased
  this -> Clear();

 // call ExecuteData
  vtkPolyData * input_unclean = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkSmartPointer<vtkCleanPolyData> cleaner =
  vtkSmartPointer<vtkCleanPolyData>::New();
  cleaner -> SetInputData (input_unclean);
  cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
  cleaner -> Update();

  vtkPolyData * input = cleaner -> GetOutput();

  input -> GetBounds(this -> bounds);

  arma::vec::fixed<3> bbox_min = {this -> scaleFactor * this -> bounds[0],this -> scaleFactor * this -> bounds[2],this -> scaleFactor * this -> bounds[4]};
  arma::vec::fixed<3> bbox_max = {this -> scaleFactor * this -> bounds[1],this -> scaleFactor * this -> bounds[3],this -> scaleFactor * this -> bounds[5]};
 
  vtkIdType cellId, numCells, numPts, numIds;

  numCells = input -> GetNumberOfCells();
  numPts = input -> GetNumberOfPoints();
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
    // Generate normals
  vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();

  normalGenerator->SetInputData(input);
  normalGenerator->ComputePointNormalsOff();
  normalGenerator->ComputeCellNormalsOn();
  normalGenerator->Update();

  vtkPolyData * input_with_normals = normalGenerator -> GetOutput();
  vtkFloatArray * normals =  vtkFloatArray::SafeDownCast(input_with_normals->GetCellData()->GetArray("Normals"));
    
  // Required by vtkPolyData::GetPointCells 
  input -> BuildLinks();

  
  // The vertex coordinates are extracted
  this -> vertices = new double * [input -> GetNumberOfPoints()];

  #pragma omp parallel for
  for(int i = 0; i < input -> GetNumberOfPoints(); ++i) {
    this -> vertices[i] = new double[3];
    input -> GetPoint(i,this -> vertices[i]);
  }

  // The facet are created
  this -> facet_normals = new double * [numCells];
  this -> facets = new int * [numCells];


  #pragma omp parallel for
  for(int i = 0; i < numCells; ++i) {
    this -> facet_normals[i] = new double[3];
    this -> facets[i] = new int[3];

    vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
    ptIds -> Allocate(VTK_CELL_SIZE);

    input -> GetCellPoints(i,ptIds);

    double normal[3];
    normals -> GetTuple(i,normal);
    
    this -> facet_normals[i][0] = normal[0];
    this -> facet_normals[i][1] = normal[1];
    this -> facet_normals[i][2] = normal[2];
    this -> facets[i][0] = ptIds -> GetId(0);
    this -> facets[i][1] = ptIds -> GetId(1);
    this -> facets[i][2] = ptIds -> GetId(2);

  }

  // The edges are extracted
  vtkSmartPointer<vtkExtractEdges> extractEdges = vtkSmartPointer<vtkExtractEdges>::New();
  extractEdges -> SetInputData(input);
  extractEdges -> Update();

  unsigned int edge_count = extractEdges -> GetOutput() -> GetNumberOfCells();
  std::vector<std::array<vtkIdType,4>> edge_points_ids_facet_ids(edge_count);

  // Should get rid of this map and use a vector instead
  // This loop cannot be parallelized since GetPointCells is 
  // not thread-safe
  for(vtkIdType i = 0; i < edge_count; i++){

    vtkSmartPointer<vtkLine> line = vtkLine::SafeDownCast(extractEdges->GetOutput()->GetCell(i));
    edge_points_ids_facet_ids[i][0] = line->GetPointIds()->GetId(0);
    edge_points_ids_facet_ids[i][1] = line->GetPointIds()->GetId(1);

    // We need to find the facets forming this edge
    vtkSmartPointer<vtkIdList> facet_ids_point_1 = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> facet_ids_point_2 = vtkSmartPointer<vtkIdList>::New();

    input -> GetPointCells  ( edge_points_ids_facet_ids[i][0],facet_ids_point_1 );
    input -> GetPointCells  ( edge_points_ids_facet_ids[i][1],facet_ids_point_2 );



    // Now, we find the two facet indices showing up in both facet_ids_point_1 and facet_ids_point_2
    facet_ids_point_1 -> IntersectWith  ( facet_ids_point_2 ) ;

    if (facet_ids_point_1 -> GetNumberOfIds() != 2){
      throw(std::runtime_error("In SBGATPolyhedronGravityModel.cpp: the intersection of the facet id lists should have exactly 2 items, not " + std::to_string(facet_ids_point_1 -> GetNumberOfIds())));
    }

    edge_points_ids_facet_ids[i][2] = facet_ids_point_1->GetId(0);
    edge_points_ids_facet_ids[i][3] = facet_ids_point_1->GetId(1);
    
  }


  // The edges dyads are created
  this -> edges = new int * [edge_points_ids_facet_ids.size()];
  this -> edge_facets_ids = new int * [edge_points_ids_facet_ids.size()];


  #pragma omp parallel for
  for(unsigned int i = 0; i < edge_points_ids_facet_ids.size(); ++i) {
    this -> edges[i] = new int[2];
    this -> edge_facets_ids[i] = new int[2];

    unsigned int p0_index = edge_points_ids_facet_ids[i][0];
    unsigned int p1_index = edge_points_ids_facet_ids[i][1];
    unsigned int fA_index = edge_points_ids_facet_ids[i][2];
    unsigned int fB_index = edge_points_ids_facet_ids[i][3];

    this -> edges[i][0] = p0_index;
    this -> edges[i][1] = p1_index;

    this -> edge_facets_ids[i][0] = fA_index;
    this -> edge_facets_ids[i][1] = fB_index;


  }

  this -> N_edges = edge_count;
  this -> N_facets = numCells;
  this -> N_vertices = input -> GetNumberOfPoints();

  // Check that the Euler characteristic == 2
  assert (input -> GetNumberOfPoints() - edge_count + numCells == 2);

  return 1;
}



void SBGATFilter::PrintHeader(ostream& os, vtkIndent indent) {

}
void SBGATFilter::PrintTrailer(ostream& os, vtkIndent indent) {

}


//----------------------------------------------------------------------------
void SBGATFilter::PrintSelf(std::ostream& os, vtkIndent indent){

  vtkPolyData *input = vtkPolyData::SafeDownCast(this->GetInput(0));
  if (!input){
    return;
  }
  
}


void SBGATFilter::Clear(){



  if (this -> N_facets > 0){

  //Vertices
    for(int i = 0; i < this -> N_vertices; ++i) {
      delete[] this -> vertices[i];   
    }
    delete[] this -> vertices;

 
  //Facets
    for(int i = 0; i < this -> N_facets; ++i) {
      delete[] this -> facets[i];   
    }
    delete[] this -> facets;

  //Facet normals
    for(int i = 0; i < this -> N_facets; ++i) {
      delete[] this -> facet_normals[i];   
    }
    delete[] this -> facet_normals;


  //Edges 
    for(int i = 0; i < this -> N_edges; ++i) {
      delete[] this -> edges[i];   
    }
    delete[] this -> edges;

  // Edge facets ids

    for (int i = 0; i < this -> N_edges; ++i){
      delete[] this -> edge_facets_ids[i];
    }
    delete[] this -> edge_facets_ids;
  }

  this -> N_edges = 0;
  this -> N_facets = 0;
  this -> N_vertices = 0;


}


arma::vec::fixed<3> SBGATFilter::GetNonNormalizedFacetNormal(const int & f) const{

  double r0[3], r1[3], r2[3];

  this -> GetVerticesInFacet(f,r0,r1,r2);

  arma::vec::fixed<3> r0_arma = {r0[0],r0[1],r0[2]};
  arma::vec::fixed<3> r1_arma = {r1[0],r1[1],r1[2]};
  arma::vec::fixed<3> r2_arma = {r2[0],r2[1],r2[2]};
  return arma::cross(r1_arma - r0_arma,r2_arma - r1_arma);

}


void SBGATFilter::GetVerticesInFacet(const int & f,double * r0,double * r1, double * r2) const{

  r0[0] = this -> vertices[this -> facets[f][0]][0];
  r0[1] = this -> vertices[this -> facets[f][0]][1];
  r0[2] = this -> vertices[this -> facets[f][0]][2];

  r1[0] = this -> vertices[this -> facets[f][1]][0];
  r1[1] = this -> vertices[this -> facets[f][1]][1];
  r1[2] = this -> vertices[this -> facets[f][1]][2];

  r2[0] = this -> vertices[this -> facets[f][2]][0];
  r2[1] = this -> vertices[this -> facets[f][2]][1];
  r2[2] = this -> vertices[this -> facets[f][2]][2];


}


void SBGATFilter::GetVerticesOnEdge(const int & e,double * r0,double * r1) const{

  r0[0] = this -> vertices[this -> edges[e][0]][0];
  r0[1] = this -> vertices[this -> edges[e][0]][1];
  r0[2] = this -> vertices[this -> edges[e][0]][2];

  r1[0] = this -> vertices[this -> edges[e][1]][0];
  r1[1] = this -> vertices[this -> edges[e][1]][1];
  r1[2] = this -> vertices[this -> edges[e][1]][2];

}


void SBGATFilter::GetIndicesOfAdjacentFacets(const int & e,int & f0, int & f1) const{

  f0 = this -> edge_facets_ids[e][0];
  f1 = this -> edge_facets_ids[e][1];

}

arma::vec::fixed<3> SBGATFilter::GetFacetCenter(const int & f) const{
  
  double r0[3];
  double r1[3];
  double r2[3];
  this -> GetVerticesInFacet(f,r0,r1,r2);

  return 1./3 * arma::vec({r0[0] + r1[0] + r2[0], r0[1] + r1[1] + r2[1], r0[2] + r1[2] + r2[2] } );


}


