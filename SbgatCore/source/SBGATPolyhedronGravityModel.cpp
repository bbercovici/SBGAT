/** MIT License

Copyright (c) 2018 Benjamin Bercovici

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

  Program:   Visualization Toolkit
  Module:    SBGATMassProperties.cxx

  Derived class from VTK's vtkPolyDataAlgorithm by Benjamin Bercovici  

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "SBGATPolyhedronGravityModel.hpp"

#include <vtkObjectFactory.h>
#include <vtkCell.h>
#include <vtkDataObject.h>
#include <vtkIdList.h>
#include <vtkMath.h>
#include <vtkSmartPointer.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPolyDataNormals.h>
#include <RigidBodyKinematics.hpp>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkExtractEdges.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkLine.h>
#include <set>
#include <vtkMath.h>

vtkStandardNewMacro(SBGATPolyhedronGravityModel);

//----------------------------------------------------------------------------
// Constructs with initial 0 values.
SBGATPolyhedronGravityModel::SBGATPolyhedronGravityModel(){

	this -> N_facets = 0;
	this -> N_edges = 0;

	
	this->SetNumberOfOutputPorts(0);
}

//----------------------------------------------------------------------------
// Destroy any allocated memory.
SBGATPolyhedronGravityModel::~SBGATPolyhedronGravityModel(){

	this -> Clear();


}


//----------------------------------------------------------------------------
// Description:
// This method computes the internal variables required for the 
// polyhedron gravity model evaluation

int SBGATPolyhedronGravityModel::RequestData(
	vtkInformation* vtkNotUsed( request ),
	vtkInformationVector** inputVector,
	vtkInformationVector* vtkNotUsed( outputVector )){
	vtkInformation *inInfo =
	inputVector[0]->GetInformationObject(0);



	if (!(this -> densitySet && this -> scaleFactorSet)){
		throw(std::runtime_error("Trying to evaluate polyhedron gravity model although the density and scale factor may have not been properly set"));
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

    // Generate normals
	vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();

	normalGenerator->SetInputData(input);
	normalGenerator->ComputePointNormalsOff();
	normalGenerator->ComputeCellNormalsOn();
	normalGenerator->Update();

	vtkPolyData * input_with_normals = normalGenerator -> GetOutput();
	
	// Required by vtkPolyData::GetPointCells	
	input -> BuildLinks();

	vtkFloatArray * normals =  vtkFloatArray::SafeDownCast(input_with_normals->GetCellData()->GetArray("Normals"));
	
	// Any data previously owned is erased
	this -> Clear();

	// The vertex coordinates are extracted
	this -> vertices = new double * [input -> GetNumberOfPoints()];

	#pragma omp parallel for
	for(int i = 0; i < input -> GetNumberOfPoints(); ++i) {
		this -> vertices[i] = new double[3];
		input -> GetPoint(i,this -> vertices[i]);
	}

	// The facet dyads are created
	this -> facet_dyads = new double * [numCells];
	this -> facet_normals = new double * [numCells];
	this -> facets = new int * [numCells];


	#pragma omp parallel for
	for(int i = 0; i < numCells; ++i) {
		this -> facet_dyads[i] = new double[9];
		this -> facet_normals[i] = new double[3];
		this -> facets[i] = new int[3];

		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
		ptIds -> Allocate(VTK_CELL_SIZE);

		input -> GetCellPoints(i,ptIds);


		double normal[3];
		normals -> GetTuple(i,normal);
		this -> facet_dyads[i][0] = normal[0] * normal[0];
		this -> facet_dyads[i][1] = normal[0] * normal[1];
		this -> facet_dyads[i][2] = normal[0] * normal[2];
		this -> facet_dyads[i][3] = normal[1] * normal[0];
		this -> facet_dyads[i][4] = normal[1] * normal[1];
		this -> facet_dyads[i][5] = normal[1] * normal[2];
		this -> facet_dyads[i][6] = normal[2] * normal[0];
		this -> facet_dyads[i][7] = normal[2] * normal[1];
		this -> facet_dyads[i][8] = normal[2] * normal[2];
		this -> facet_normals[i][0] = normal[0];
		this -> facet_normals[i][1] = normal[1];
		this -> facet_normals[i][2] = normal[2];
		this -> facets[i][0] = ptIds -> GetId(0);
		this -> facets[i][1] = ptIds -> GetId(1);
		this -> facets[i][2] = ptIds -> GetId(2);

	}

	// The edges are extracted

	vtkSmartPointer<vtkExtractEdges> extractEdges = 
	vtkSmartPointer<vtkExtractEdges>::New();
	extractEdges->SetInputData(input);
	extractEdges->Update();


	unsigned int edge_count = extractEdges->GetOutput()->GetNumberOfCells();
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

		input -> GetPointCells	(	edge_points_ids_facet_ids[i][0],facet_ids_point_1 );
		input -> GetPointCells	(	edge_points_ids_facet_ids[i][1],facet_ids_point_2 );

		// Now, we find the two indices showing up in both facet_ids_point_1 and facet_ids_point_2
		facet_ids_point_1 -> IntersectWith	(	facet_ids_point_2	)	;

		if (facet_ids_point_1 -> GetNumberOfIds() != 2){
			throw(std::runtime_error("In SBGATPolyhedronGravityModel.cpp: the intersection of the facet id lists should have exactly 2 items, not " + std::to_string(facet_ids_point_2 -> GetNumberOfIds())));
		}
		
		edge_points_ids_facet_ids[i][2] = facet_ids_point_1->GetId(0);
		edge_points_ids_facet_ids[i][3] = facet_ids_point_1->GetId(1);
		
	}


	// The edges dyads are created
	this -> edge_dyads = new double * [edge_points_ids_facet_ids.size()];
	this -> edges = new int * [edge_points_ids_facet_ids.size()];

	#pragma omp parallel for
	for(unsigned int i = 0; i < edge_points_ids_facet_ids.size(); ++i) {
		this -> edge_dyads[i] = new double[9];
		this -> edges[i] = new int[2];

		unsigned int p0_index = edge_points_ids_facet_ids[i][0];
		unsigned int p1_index = edge_points_ids_facet_ids[i][1];
		unsigned int fA_index = edge_points_ids_facet_ids[i][2];
		unsigned int fB_index = edge_points_ids_facet_ids[i][3];

		double nA[3];
		double nB[3];

		double p0[3];
		double p1[3];

		input -> GetPoint(p0_index,p0);
		input -> GetPoint(p1_index,p1);

		nA[0] = this -> facet_normals[fA_index][0];
		nA[1] = this -> facet_normals[fA_index][1];
		nA[2] = this -> facet_normals[fA_index][2];

		nB[0] = this -> facet_normals[fB_index][0];
		nB[1] = this -> facet_normals[fB_index][1];
		nB[2] = this -> facet_normals[fB_index][2];

		double edge_dir[3];
		vtkMath::Cross(nA,nB,edge_dir);
		vtkMath::Normalize(edge_dir);
		double p1_m_p0[3];
		vtkMath::Subtract(p1,p0,p1_m_p0);
		
		if (vtkMath::Dot(p1_m_p0,edge_dir) < 0){
			vtkMath::MultiplyScalar(edge_dir,-1.);
		} 

		double edge_normal_A_to_B[3];
		double edge_normal_B_to_A[3];

		vtkMath::Cross(nA,edge_dir,edge_normal_A_to_B);
		vtkMath::Cross(nB,edge_dir,edge_normal_B_to_A);
		vtkMath::MultiplyScalar(edge_normal_A_to_B,-1.);

		double dyad_A[3][3];
		double dyad_B[3][3];


		vtkMath::Outer(nA,edge_normal_A_to_B,dyad_A);
		vtkMath::Outer(nB,edge_normal_B_to_A,dyad_B);

		this -> edge_dyads[i][0] = dyad_A[0][0] + dyad_B[0][0] ;
		this -> edge_dyads[i][1] = dyad_A[0][1] + dyad_B[0][1] ;
		this -> edge_dyads[i][2] = dyad_A[0][2] + dyad_B[0][2] ;
		this -> edge_dyads[i][3] = dyad_A[1][0] + dyad_B[1][0] ;
		this -> edge_dyads[i][4] = dyad_A[1][1] + dyad_B[1][1] ;
		this -> edge_dyads[i][5] = dyad_A[1][2] + dyad_B[1][2] ;
		this -> edge_dyads[i][6] = dyad_A[2][0] + dyad_B[2][0] ;
		this -> edge_dyads[i][7] = dyad_A[2][1] + dyad_B[2][1] ;
		this -> edge_dyads[i][8] = dyad_A[2][2] + dyad_B[2][2] ;


		this -> edges[i][0] = p0_index;
		this -> edges[i][1] = p1_index;


	}

	this -> N_edges = edge_count;
	this -> N_facets = numCells;
	this -> N_facets = numCells;


	this -> mass_properties = vtkSmartPointer<SBGATMassProperties>::New();
	this -> mass_properties -> SetInputData(input);
	this -> mass_properties -> Update();

	// Check that the Euler characteristic == 2
	assert (input -> GetNumberOfPoints() - edge_count + numCells == 2);

	// Check that the shape is topologically closed
	assert(this -> mass_properties -> CheckClosed());

	return 1;
}




double SBGATPolyhedronGravityModel::GetPotential(arma::vec point){
	return this -> GetPotential(point.colptr(0));
}


double SBGATPolyhedronGravityModel::GetPotential(double * point) {

	double potential = 0;

	// Facet loop
	#pragma omp parallel for reduction(+:potential)
	for (vtkIdType facet_index = 0; facet_index < this -> N_facets; ++ facet_index) {


		double * r0 = this -> vertices[this -> facets[facet_index][0]];
		double * r1 = this -> vertices[this -> facets[facet_index][1]];
		double * r2 = this -> vertices[this -> facets[facet_index][2]];


		double r0m[3];
		double r1m[3];
		double r2m[3];

		vtkMath::Subtract(r0,point,r0m);
		vtkMath::Subtract(r1,point,r1m);
		vtkMath::Subtract(r2,point,r2m);

		double R0 = vtkMath::Norm(r0m);
		double R1 = vtkMath::Norm(r1m);
		double R2 = vtkMath::Norm(r2m);

		double r1m_cross_r2m[3];

		vtkMath::Cross(r1m,r2m,r1m_cross_r2m);

		double wf = 2 * std::atan2(vtkMath::Dot(r0m,r1m_cross_r2m),R0 * R1 * R2 + 
			R0 * vtkMath::Dot(r1m,r2m) + R1 * vtkMath::Dot(r0m,r2m) + R2 * vtkMath::Dot(r0m,r1m));

		double * F = this -> facet_dyads[facet_index];

		double a[3] = {
			F[0] * r0m[0] + F[1] * r0m[1] +  F[2] * r0m[2],
			F[3] * r0m[0] + F[4] * r0m[1] +  F[5] * r0m[2],
			F[6] * r0m[0] + F[7] * r0m[1] +  F[8] * r0m[2]
		};
		

		potential += - wf * vtkMath::Dot(r0m,a);

	}

	// Edge loop
	#pragma omp parallel for reduction(+:potential)
	for (unsigned int edge_index = 0; edge_index < this -> N_edges; ++ edge_index) {

		double * r0 = this -> vertices[this -> edges[edge_index][0]];
		double * r1 = this -> vertices[this -> edges[edge_index][1]];

		
		double r0m[3];
		double r1m[3];
		double rem[3];

		vtkMath::Subtract(r0,point,r0m);
		vtkMath::Subtract(r1,point,r1m);
		vtkMath::Subtract(r1m,r0m,rem);

		double R0 = vtkMath::Norm(r0m);
		double R1 = vtkMath::Norm(r1m);
		double Re = vtkMath::Norm(rem);

		double Le = std::log((R0 + R1 + Re) / (R0 + R1 - Re));

		double * E = this -> edge_dyads[edge_index];

		double a[3] = {
			E[0] * r0m[0] + E[1] * r0m[1] +  E[2] * r0m[2],
			E[3] * r0m[0] + E[4] * r0m[1] +  E[5] * r0m[2],
			E[6] * r0m[0] + E[7] * r0m[1] +  E[8] * r0m[2]
		};


		potential += Le * vtkMath::Dot(r0m,a);


	}

	potential *= 0.5 * arma::datum::G / std::pow(this -> scaleFactor,3) * this ->density;

	return potential;

}


bool SBGATPolyhedronGravityModel::Contains(double * point, double tol ) {

	double laplacian = 0;

	// Facet loop
	#pragma omp parallel for reduction(+:laplacian)
	for (vtkIdType facet_index = 0; facet_index < this -> N_facets; ++ facet_index) {

		double * r0 = this -> vertices[this -> facets[facet_index][0]];
		double * r1 = this -> vertices[this -> facets[facet_index][1]];
		double * r2 = this -> vertices[this -> facets[facet_index][2]];

		double r0m[3];
		double r1m[3];
		double r2m[3];

		vtkMath::Subtract(r0,point,r0m);
		vtkMath::Subtract(r1,point,r1m);
		vtkMath::Subtract(r2,point,r2m);

		double R0 = vtkMath::Norm(r0m);
		double R1 = vtkMath::Norm(r1m);
		double R2 = vtkMath::Norm(r2m);

		double r1m_cross_r2m[3];

		vtkMath::Cross(r1m,r2m,r1m_cross_r2m);

		double wf = 2 * std::atan2(vtkMath::Dot(r0m,r1m_cross_r2m),R0 * R1 * R2 + 
			R0 * vtkMath::Dot(r1m,r2m) + R1 * vtkMath::Dot(r0m,r2m) + R2 * vtkMath::Dot(r0m,r1m));

		laplacian += wf;

	}

	if (std::abs(laplacian) / (4 * arma::datum::pi) < tol) {
		return false;
	}
	else {
		return true;
	}

}



arma::vec SBGATPolyhedronGravityModel::GetAcceleration(arma::vec point){
	return this-> GetAcceleration(  point.colptr(0));
}

arma::vec SBGATPolyhedronGravityModel::GetAcceleration(double * point) {

	double acc_x = 0;
	double acc_y = 0;
	double acc_z = 0;


	// Facet loop
	#pragma omp parallel for reduction(+:acc_x,acc_y,acc_z)
	for (vtkIdType facet_index = 0; facet_index < this -> N_facets; ++ facet_index) {

		double * r0 = this -> vertices[this -> facets[facet_index][0]];
		double * r1 = this -> vertices[this -> facets[facet_index][1]];
		double * r2 = this -> vertices[this -> facets[facet_index][2]];

		double r0m[3];
		double r1m[3];
		double r2m[3];

		vtkMath::Subtract(r0,point,r0m);
		vtkMath::Subtract(r1,point,r1m);
		vtkMath::Subtract(r2,point,r2m);

		double R0 = vtkMath::Norm(r0m);
		double R1 = vtkMath::Norm(r1m);
		double R2 = vtkMath::Norm(r2m);

		double r1m_cross_r2m[3];

		vtkMath::Cross(r1m,r2m,r1m_cross_r2m);

		double wf = 2 * std::atan2(vtkMath::Dot(r0m,r1m_cross_r2m),R0 * R1 * R2 + 
			R0 * vtkMath::Dot(r1m,r2m) + R1 * vtkMath::Dot(r0m,r2m) + R2 * vtkMath::Dot(r0m,r1m));

		double * F = this -> facet_dyads[facet_index];

		acc_x += wf *( F[0] * r0m[0] + F[1] * r0m[1] +  F[2] * r0m[2]);
		acc_y += wf *( F[3] * r0m[0] + F[4] * r0m[1] +  F[5] * r0m[2]);
		acc_z += wf *( F[6] * r0m[0] + F[7] * r0m[1] +  F[8] * r0m[2]);


	}

	// Edge loop
	#pragma omp parallel for reduction(-:acc_x,acc_y,acc_z)
	for (unsigned int edge_index = 0; edge_index < this -> N_edges; ++ edge_index) {

		double * r0 = this -> vertices[this -> edges[edge_index][0]];
		double * r1 = this -> vertices[this -> edges[edge_index][1]];

		
		double r0m[3];
		double r1m[3];
		double rem[3];

		vtkMath::Subtract(r0,point,r0m);
		vtkMath::Subtract(r1,point,r1m);
		vtkMath::Subtract(r1m,r0m,rem);

		double R0 = vtkMath::Norm(r0m);
		double R1 = vtkMath::Norm(r1m);
		double Re = vtkMath::Norm(rem);

		double Le = std::log((R0 + R1 + Re) / (R0 + R1 - Re));

		double * E = this -> edge_dyads[edge_index];

		acc_x -= Le *( E[0] * r0m[0] + E[1] * r0m[1] +  E[2] * r0m[2]);
		acc_y -= Le *( E[3] * r0m[0] + E[4] * r0m[1] +  E[5] * r0m[2]);
		acc_z -= Le *( E[6] * r0m[0] + E[7] * r0m[1] +  E[8] * r0m[2]);

	}

	arma::vec acc = {acc_x,acc_y,acc_z};
	acc *= arma::datum::G  / std::pow(this -> scaleFactor,3)* this -> density;

	return acc;

}


void SBGATPolyhedronGravityModel::PrintHeader(ostream& os, vtkIndent indent) {

}
void SBGATPolyhedronGravityModel::PrintTrailer(ostream& os, vtkIndent indent) {

}


void SBGATPolyhedronGravityModel::Clear(){
	if (this -> N_facets > 0){
		int N_vertices = this -> N_edges - this -> N_facets + 2;

	//Vertices
		for(unsigned int i = 0; i < N_vertices; ++i) {
			delete[] this -> vertices[i];   
		}
		delete[] this -> vertices;

	//Facet dyads
		for(unsigned int i = 0; i < this -> N_facets; ++i) {
			delete[] this -> facet_dyads[i];   
		}
		delete[] this -> facet_dyads;

	//Facets
		for(unsigned int i = 0; i < this -> N_facets; ++i) {
			delete[] this -> facets[i];   
		}
		delete[] this -> facets;

	//Facet normals
		for(unsigned int i = 0; i < this -> N_facets; ++i) {
			delete[] this -> facet_normals[i];   
		}
		delete[] this -> facet_normals;


	//Edge dyads
		for(unsigned int i = 0; i < this -> N_edges; ++i) {
			delete[] this -> edge_dyads[i];   
		}
		delete[] this -> edge_dyads;

	//Edges 
		for(unsigned int i = 0; i < this -> N_edges; ++i) {
			delete[] this -> edges[i];   
		}
		delete[] this -> edges;
	}
}



//----------------------------------------------------------------------------
void SBGATPolyhedronGravityModel::PrintSelf(std::ostream& os, vtkIndent indent){

	vtkPolyData *input = vtkPolyData::SafeDownCast(this->GetInput(0));
	if (!input)
	{
		return;
	}
}










