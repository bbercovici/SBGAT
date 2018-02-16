/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SBGATMassProperties.cxx

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
	
	this->SetNumberOfOutputPorts(0);
}

//----------------------------------------------------------------------------
// Destroy any allocated memory.
SBGATPolyhedronGravityModel::~SBGATPolyhedronGravityModel(){

	//Facet dyads
	for(unsigned int i = 0; i < this -> N_facets; ++i) {
		delete[] this -> facet_dyads[i];   
	}
	delete[] this -> facet_dyads;

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
	

	// The facet dyads are created
	this -> facet_dyads = new double * [numCells];
	this -> facet_normals = new double * [numCells];


	#pragma omp parallel for
	for(int i = 0; i < numCells; ++i) {
		this -> facet_dyads[i] = new double[9];
		this -> facet_normals[i] = new double[3];

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

	// #pragma omp parallel for
	for(int i = 0; i < edge_points_ids_facet_ids.size(); ++i) {
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

	this -> mass_properties = vtkSmartPointer<SBGATMassProperties>::New();
	this -> mass_properties -> SetInputData(input);
	this -> mass_properties -> Update();

	// Check that the Euler characteristic == 2
	assert (input -> GetNumberOfPoints() - edge_count + numCells == 2);

}


double SBGATPolyhedronGravityModel::ComputePgmPotential(double * point ,const double density,const double G) {

	double potential = 0;

	vtkPolyData * input = this -> GetPolyDataInput(0);

	// Facet loop
	#pragma omp parallel for reduction(+:potential)
	for (vtkIdType facet_index = 0; facet_index < this -> N_facets; ++ facet_index) {

		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
		ptIds->Allocate(VTK_CELL_SIZE);

		input->GetCellPoints(facet_index,ptIds);

		double r0[3];
		double r1[3];
		double r2[3];

		input -> GetPoint(ptIds -> GetId(0),r0);
		input -> GetPoint(ptIds -> GetId(1),r1);
		input -> GetPoint(ptIds -> GetId(2),r2);

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

		double r0[3];
		double r1[3];

		input -> GetPoint(this -> edges[edge_index][0],r0);
		input -> GetPoint(this -> edges[edge_index][1],r1);
		
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

	potential *= 0.5 * G * density;

	return potential;

}


bool SBGATPolyhedronGravityModel::Contains(double * point, double tol ) {

	double laplacian = 0;
	vtkPolyData * input = this -> GetPolyDataInput(0);

	// Facet loop
	// #pragma omp parallel for reduction(+:laplacian)
	for (vtkIdType facet_index = 0; facet_index < input -> GetNumberOfCells(); ++ facet_index) {

		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
		ptIds->Allocate(VTK_CELL_SIZE);

		input->GetCellPoints(facet_index,ptIds);

		double r0[3];
		double r1[3];
		double r2[3];

		input -> GetPoint(ptIds -> GetId(0),r0);
		input -> GetPoint(ptIds -> GetId(1),r1);
		input -> GetPoint(ptIds -> GetId(2),r2);

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

	if (std::abs(laplacian) < tol) {
		return false;
	}
	else {
		return true;
	}

}


arma::vec SBGATPolyhedronGravityModel::ComputePgmAcceleration(double * point ,const double density,const double G) {

	double acc_x = 0;
	double acc_y = 0;
	double acc_z = 0;


	vtkPolyData * input = this -> GetPolyDataInput(0);

	// Facet loop
	#pragma omp parallel for reduction(+:acc_x,acc_y,acc_z)
	for (vtkIdType facet_index = 0; facet_index < this -> N_facets; ++ facet_index) {

		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
		ptIds->Allocate(VTK_CELL_SIZE);

		input->GetCellPoints(facet_index,ptIds);

		double r0[3];
		double r1[3];
		double r2[3];

		input -> GetPoint(ptIds -> GetId(0),r0);
		input -> GetPoint(ptIds -> GetId(1),r1);
		input -> GetPoint(ptIds -> GetId(2),r2);

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

		double r0[3];
		double r1[3];

		input -> GetPoint(this -> edges[edge_index][0],r0);
		input -> GetPoint(this -> edges[edge_index][1],r1);
		
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
	acc *= G * density;

	return acc;

}


void SBGATPolyhedronGravityModel::PrintHeader(ostream& os, vtkIndent indent) {

}
void SBGATPolyhedronGravityModel::PrintTrailer(ostream& os, vtkIndent indent) {

}

//----------------------------------------------------------------------------
void SBGATPolyhedronGravityModel::PrintSelf(std::ostream& os, vtkIndent indent){

	vtkPolyData *input = vtkPolyData::SafeDownCast(this->GetInput(0));
	if (!input)
	{
		return;
	}
  //   os << "\tVolumeX: " << this->GetVolumeX () << "\n";
  //   os << "\tVolumeY: " << this->GetVolumeY () << "\n";
  //   os << "\tVolumeZ: " << this->GetVolumeZ () << "\n";
  //   os << "\tKx: " << this->GetKx () << "\n";
  //   os << "\tKy: " << this->GetKy () << "\n";
  //   os << "\tKz: " << this->GetKz () << "\n";
  //   os << "\tVolume:  " << this->GetVolume  () << "\n";
  // //os << indent << "Volume Projected:  " << this->GetVolumeProjected  () << "\n";
  // //os << indent << "Volume Error:  " <<
  // //  fabs(this->GetVolume() - this->GetVolumeProjected())   << "\n";
  //   os << "\tSurface Area: " << this->GetSurfaceArea () << "\n";
  //   os << "\tMin Cell Area: " << this->GetMinCellArea () << "\n";
  //   os << "\tMax Cell Area: " << this->GetMaxCellArea () << "\n";
  //   os << "\tNormalized Shape Index: "
  //   << this->GetNormalizedShapeIndex () << "\n";
}










