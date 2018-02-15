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
	for(unsigned int i = 0; i < sizeof(this -> facet_dyads) / sizeof(this -> facet_dyads[0]); ++i) {
		delete[] this -> facet_dyads[i];   
	}
	delete[] this -> facet_dyads;

	//Facet normals
	for(unsigned int i = 0; i < sizeof(this -> facet_normals) / sizeof(this -> facet_normals[0]); ++i) {
		delete[] this -> facet_normals[i];   
	}
	delete[] this -> facet_normals;


	//Edge dyads
	for(unsigned int i = 0; i < sizeof(this -> edge_dyads) / sizeof(this -> edge_dyads[0]); ++i) {
		delete[] this -> edge_dyads[i];   
	}
	delete[] this -> edge_dyads;

	//Edges 
	for(unsigned int i = 0; i < sizeof(this -> edges) / sizeof(this -> edges[0]); ++i) {
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
		
		edge_points_ids_facet_ids[i][2] = line->GetPointIds()->GetId(0);
		edge_points_ids_facet_ids[i][3] = line->GetPointIds()->GetId(1);
		
	}

	
	// The edges dyads are created
	this -> edge_dyads = new double * [edge_points_ids_facet_ids.size()];
	this -> edges = new int * [edge_points_ids_facet_ids.size()];

	#pragma omp parallel for
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


		SBGATPolyhedronGravityModel::ComputeEdgeDyad(this -> edge_dyads[i],nA,nB,p0,p1);
		this -> edges[i][0] = p0_index;
		this -> edges[i][1] = p1_index;


	}

	this -> mass_properties = vtkSmartPointer<SBGATMassProperties>::New();
	this -> mass_properties -> SetInputData(input);
	this -> mass_properties -> Update();

	// Check that the Euler characteristic == 2
	assert (input -> GetNumberOfPoints() - edge_count + numCells == 2);

}

void SBGATPolyhedronGravityModel::ComputeEdgeDyad(double * edge_dyad,
	double * nA, double * nB,double * p0, double * p1) {

	
	arma::vec nA_arma = {nA[0],nA[1],nA[2]};
	arma::vec nB_arma = {nB[0],nB[1],nB[2]};

	arma::vec p1_arma = {p0[0],p0[1],p0[2]};
	arma::vec p2_arma = {p1[0],p1[1],p1[2]};

	arma::vec edge_dir = arma::normalise(arma::cross(nA_arma,nB_arma));

	if (arma::dot(edge_dir,p2_arma - p1_arma) < 0){
		edge_dir = - edge_dir;
	}

	arma::vec edge_normal_A_to_B = - arma::normalise(arma::cross(nA_arma,edge_dir));
	arma::vec edge_normal_B_to_A = arma::normalise(arma::cross(nB_arma,edge_dir));

	arma::mat edge_dyad_arma = nA_arma * edge_normal_A_to_B.t() + nB_arma * edge_normal_B_to_A.t();

	edge_dyad[0] = edge_dyad_arma(0,0);
	edge_dyad[1] = edge_dyad_arma(0,1);
	edge_dyad[2] = edge_dyad_arma(0,2);
	edge_dyad[3] = edge_dyad_arma(1,0);
	edge_dyad[4] = edge_dyad_arma(1,1);
	edge_dyad[5] = edge_dyad_arma(1,2);
	edge_dyad[6] = edge_dyad_arma(2,0);
	edge_dyad[7] = edge_dyad_arma(2,1);
	edge_dyad[8] = edge_dyad_arma(2,2);

}


double SBGATPolyhedronGravityModel::ComputePgmPotential(double * point ,const double density) {

	double potential = 0;

	vtkPolyData * input = this -> GetPolyDataInput(0);

	// Facet loop
	#pragma omp parallel for reduction(+:potential)
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

		r0m[0] = r0[0] - point[0];
		r0m[1] = r0[1] - point[1];
		r0m[2] = r0[2] - point[2];

		r1m[0] = r1[0] - point[0];
		r1m[1] = r1[1] - point[1];
		r1m[2] = r1[2] - point[2];

		r2m[0] = r2[0] - point[0];
		r2m[1] = r2[1] - point[1];
		r2m[2] = r2[2] - point[2];


		double R1 = std::sqrt( r0m[0] * r0m[0]
			+ r0m[1] * r0m[1]
			+ r0m[2] * r0m[2]       );

		double R2 = std::sqrt( r1m[0] * r1m[0]
			+ r1m[1] * r1m[1]
			+ r1m[2] * r1m[2]      );


		double R3 = std::sqrt( r2m[0] * r2m[0]
			+ r2m[1] * r2m[1]
			+ r2m[2] * r2m[2]      );

		double r2_cross_r3_0 = r1m[1] * r2m[2] - r1m[2] * r2m[1];
		double r2_cross_r3_1 = r2m[0] * r1m[2] - r2m[2] * r1m[0];
		double r2_cross_r3_2 = r1m[0] * r2m[1] - r1m[1] * r2m[0];


		double wf = 2 * std::atan2(
			r0m[0] * r2_cross_r3_0 + r0m[1] * r2_cross_r3_1 + r0m[2] * r2_cross_r3_2,

			R1 * R2 * R3 + R1 * (r1m[0] * r2m[0] + r1m[1] * r2m[1]  + r1m[2] * r2m[2] )
			+ R2 * (r2m[0] * r0m[0] + r2m[1] * r0m[1] + r2m[2] * r0m[2])
			+ R3 * (r0m[0] * r1m[0] + r0m[1] * r1m[1] + r0m[2] * r1m[2]));


		double * F = this -> facet_dyads[facet_index];

		double ax = wf * (F[0] * r0m[0] + F[1] * r0m[1] +  F[2] * r0m[2]);
		double ay = wf * (F[3] * r0m[0] + F[4] * r0m[1] +  F[5] * r0m[2]);
		double az = wf * (F[6] * r0m[0] + F[7] * r0m[1] +  F[8] * r0m[2]);

		potential += - (ax * r0m[0] + ay * r0m[1] + az * r0m[2]);

	}

	// Edge loop
	#pragma omp parallel for reduction(+:potential)
	for (unsigned int edge_index = 0; edge_index < sizeof(this -> edge_dyads) / sizeof(this -> edge_dyads[0]); ++ edge_index) {


		double r0[3];
		double r1[3];

		input -> GetPoint(this -> edges[edge_index][0],r0);
		input -> GetPoint(this -> edges[edge_index][1],r1);

		double r0m[3];
		double r1m[3];

		r0m[0] = r0[0] - point[0];
		r0m[1] = r0[1] - point[1];
		r0m[2] = r0[2] - point[2];

		r1m[0] = r1[0] - point[0];
		r1m[1] = r1[1] - point[1];
		r1m[2] = r1[2] - point[2];

		double R1 = std::sqrt( r0m[0] * r0m[0]
			+ r0m[1] * r0m[1]
			+ r0m[2] * r0m[2]       );

		double R2 = std::sqrt( r1m[0] * r1m[0]
			+ r1m[1] * r1m[1]
			+ r1m[2] * r1m[2]      );

		double Re = std::sqrt( (r1m[0] - r0m[0]) * (r1m[0] - r0m[0])
			+ (r1m[1] - r0m[1]) * (r1m[1] - r0m[1])
			+ (r1m[2] - r0m[2]) * (r1m[2] - r0m[2])      );


		double Le = std::log((R1 + R2 + Re) / (R1 + R2 - Re));


		double * E = this -> edge_dyads[edge_index];

		double ax = Le * (E[0] * r0m[0] + E[1] * r0m[1] +  E[2] * r0m[2]);
		double ay = Le * (E[3] * r0m[0] + E[4] * r0m[1] +  E[5] * r0m[2]);
		double az = Le * (E[6] * r0m[0] + E[7] * r0m[1] +  E[8] * r0m[2]);

		potential += (ax * r0m[0] + ay * r0m[1] + az * r0m[2]);

	}

	potential *= 0.5 * arma::datum::G * density;

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

		r0m[0] = r0[0] - point[0];
		r0m[1] = r0[1] - point[1];
		r0m[2] = r0[2] - point[2];

		r1m[0] = r1[0] - point[0];
		r1m[1] = r1[1] - point[1];
		r1m[2] = r1[2] - point[2];

		r2m[0] = r2[0] - point[0];
		r2m[1] = r2[1] - point[1];
		r2m[2] = r2[2] - point[2];


		double R0 = std::sqrt( r0m[0] * r0m[0]
			+ r0m[1] * r0m[1]
			+ r0m[2] * r0m[2]       );

		double R1 = std::sqrt( r1m[0] * r1m[0]
			+ r1m[1] * r1m[1]
			+ r1m[2] * r1m[2]      );


		double R2 = std::sqrt( r2m[0] * r2m[0]
			+ r2m[1] * r2m[1]
			+ r2m[2] * r2m[2]      );

		double r1_cross_r2_0 = r1m[1] * r2m[2] - r1m[2] * r2m[1];
		double r1_cross_r2_1 = r2m[0] * r1m[2] - r2m[2] * r1m[0];
		double r1_cross_r2_2 = r1m[0] * r2m[1] - r1m[1] * r2m[0];


		double wf = 2 * std::atan2(
			r0m[0] * r1_cross_r2_0 + r0m[1] * r1_cross_r2_1 + r0m[2] * r1_cross_r2_2,

			R0 * R1 * R2 + R0 * (r1m[0] * r2m[0] + r1m[1] * r2m[1]  + r1m[2] * r2m[2] )
			+ R1 * (r2m[0] * r0m[0] + r2m[1] * r0m[1] + r2m[2] * r0m[2])
			+ R2 * (r0m[0] * r1m[0] + r0m[1] * r1m[1] + r0m[2] * r1m[2]));

		laplacian += wf;

	}

	if (std::abs(laplacian) < tol) {
		return false;
	}
	else {
		return true;
	}

}


	// arma::vec DynamicAnalyses::pgm_acceleration(double * point , const double mu) const {

	// 	double ax = 0;
	// 	double ay = 0;
	// 	double az = 0;

	// // Facet loop
	// #pragma omp parallel for reduction(+:ax,ay,az) if (USE_OMP_DYNAMIC_ANALYSIS)
	// 	for (unsigned int facet_index = 0; facet_index < this -> shape_model -> get_NFacets(); ++ facet_index) {

	// 		std::vector<std::shared_ptr<Vertex > > * vertices = this -> shape_model -> get_facets() -> at(facet_index) -> get_vertices();

	// 		const double * r0 =  vertices -> at(0) -> get_coordinates() -> colptr(0);
	// 		const double * r1 =  vertices -> at(1) -> get_coordinates() -> colptr(0);
	// 		const double * r2 =  vertices -> at(2) -> get_coordinates() -> colptr(0);

	// 		double r0m[3];
	// 		double r1m[3];
	// 		double r2m[3];

	// 		r0m[0] = r0[0] - point[0];
	// 		r0m[1] = r0[1] - point[1];
	// 		r0m[2] = r0[2] - point[2];

	// 		r1m[0] = r1[0] - point[0];
	// 		r1m[1] = r1[1] - point[1];
	// 		r1m[2] = r1[2] - point[2];

	// 		r2m[0] = r2[0] - point[0];
	// 		r2m[1] = r2[1] - point[1];
	// 		r2m[2] = r2[2] - point[2];


	// 		double R1 = std::sqrt( r0m[0] * r0m[0]
	// 			+ r0m[1] * r0m[1]
	// 			+ r0m[2] * r0m[2]       );

	// 		double R2 = std::sqrt( r1m[0] * r1m[0]
	// 			+ r1m[1] * r1m[1]
	// 			+ r1m[2] * r1m[2]      );


	// 		double R3 = std::sqrt( r2m[0] * r2m[0]
	// 			+ r2m[1] * r2m[1]
	// 			+ r2m[2] * r2m[2]      );

	// 		double r2_cross_r3_0 = r1m[1] * r2m[2] - r1m[2] * r2m[1];
	// 		double r2_cross_r3_1 = r2m[0] * r1m[2] - r2m[2] * r1m[0];
	// 		double r2_cross_r3_2 = r1m[0] * r2m[1] - r1m[1] * r2m[0];


	// 		double wf = 2 * std::atan2(
	// 			r0m[0] * r2_cross_r3_0 + r0m[1] * r2_cross_r3_1 + r0m[2] * r2_cross_r3_2,

	// 			R1 * R2 * R3 + R1 * (r1m[0] * r2m[0] + r1m[1] * r2m[1]  + r1m[2] * r2m[2] )
	// 			+ R2 * (r2m[0] * r0m[0] + r2m[1] * r0m[1] + r2m[2] * r0m[2])
	// 			+ R3 * (r0m[0] * r1m[0] + r0m[1] * r1m[1] + r0m[2] * r1m[2]));


	// 		arma::mat * Fdyad = this -> shape_model -> get_facets() -> at(facet_index) -> get_facet_dyad();

	// 		double * F_col_0 = Fdyad -> colptr(0);
	// 		double * F_col_1 = Fdyad -> colptr(1);
	// 		double * F_col_2 = Fdyad -> colptr(2);

	// 		ax += wf * (F_col_0[0] * r0m[0] + F_col_1[0] * r0m[1] +  F_col_2[0] * r0m[2]);
	// 		ay += wf * (F_col_0[1] * r0m[0] + F_col_1[1] * r0m[1] +  F_col_2[1] * r0m[2]);
	// 		az += wf * (F_col_0[2] * r0m[0] + F_col_1[2] * r0m[1] +  F_col_2[2] * r0m[2]);

	// 	}


	// // Edge loop
	// #pragma omp parallel for reduction(-:ax,ay,az) if (USE_OMP_DYNAMIC_ANALYSIS)
	// 	for (unsigned int edge_index = 0; edge_index < this -> shape_model -> get_NEdges(); ++ edge_index) {

	// 		const double * r0 =  this -> shape_model -> get_edges() -> at(edge_index) -> get_v0() -> get_coordinates() -> colptr(0);
	// 		const double * r1 =  this -> shape_model -> get_edges() -> at(edge_index) -> get_v1() -> get_coordinates() -> colptr(0);

	// 		double r0m[3];
	// 		double r1m[3];

	// 		r0m[0] = r0[0] - point[0];
	// 		r0m[1] = r0[1] - point[1];
	// 		r0m[2] = r0[2] - point[2];

	// 		r1m[0] = r1[0] - point[0];
	// 		r1m[1] = r1[1] - point[1];
	// 		r1m[2] = r1[2] - point[2];


	// 		double R1 = std::sqrt( r0m[0] * r0m[0]
	// 			+ r0m[1] * r0m[1]
	// 			+ r0m[2] * r0m[2]       );

	// 		double R2 = std::sqrt( r1m[0] * r1m[0]
	// 			+ r1m[1] * r1m[1]
	// 			+ r1m[2] * r1m[2]      );

	// 		double Re = std::sqrt( (r1m[0] - r0m[0]) * (r1m[0] - r0m[0])
	// 			+ (r1m[1] - r0m[1]) * (r1m[1] - r0m[1])
	// 			+ (r1m[2] - r0m[2]) * (r1m[2] - r0m[2])      );


	// 		double Le = std::log((R1 + R2 + Re) / (R1 + R2 - Re));


	// 		arma::mat * Edyad = this -> shape_model -> get_edges() -> at(edge_index) -> get_edge_dyad();

	// 		double * E_col_0 = Edyad -> colptr(0);
	// 		double * E_col_1 = Edyad -> colptr(1);
	// 		double * E_col_2 = Edyad -> colptr(2);

	// 		ax -= Le * (E_col_0[0] * r0m[0] + E_col_1[0] * r0m[1] +  E_col_2[0] * r0m[2]);
	// 		ay -= Le * (E_col_0[1] * r0m[0] + E_col_1[1] * r0m[1] +  E_col_2[1] * r0m[2]);
	// 		az -= Le * (E_col_0[2] * r0m[0] + E_col_1[2] * r0m[1] +  E_col_2[2] * r0m[2]);

	// 	}

	// 	arma::vec acceleration = {ax, ay, az};
	// 	acceleration = acceleration * mu / this -> shape_model -> get_volume();

	// 	return acceleration;

	// }



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










