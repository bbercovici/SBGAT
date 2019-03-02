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
#include <vtkCleanPolyData.h>
#include <json.hpp>

#include <vtkLine.h>
#include <set>
#include <vtkMath.h>
#include <array>

#pragma omp declare reduction( + : arma::vec::fixed<3> : omp_out += omp_in ) \
initializer( omp_priv = arma::zeros<arma::vec>(3))

#pragma omp declare reduction( + : arma::mat::fixed<3,3> : omp_out += omp_in ) \
initializer( omp_priv = arma::zeros<arma::mat>(3,3))

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
	inputVector[0] -> GetInformationObject(0);

	if (!(this -> densitySet && this -> scaleFactorSet)){
		throw(std::runtime_error("Trying to evaluate polyhedron gravity model although the density and scale factor have not been properly set"));
	}

  	// call ExecuteData
	vtkPolyData * input_unclean = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkSmartPointer<vtkCleanPolyData> cleaner =
	vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner -> SetInputData (input_unclean);
	cleaner -> SetOutputPointsPrecision	( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
	cleaner -> Update();

	vtkPolyData * input = cleaner -> GetOutput();

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

	normalGenerator -> SetInputData(input);
	normalGenerator -> ComputePointNormalsOff();
	normalGenerator -> ComputeCellNormalsOn();
	normalGenerator -> Update();

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

		input -> GetPointCells	(	edge_points_ids_facet_ids[i][0],facet_ids_point_1 );
		input -> GetPointCells	(	edge_points_ids_facet_ids[i][1],facet_ids_point_2 );



		// Now, we find the two facet indices showing up in both facet_ids_point_1 and facet_ids_point_2
		facet_ids_point_1 -> IntersectWith	(	facet_ids_point_2	)	;

		if (facet_ids_point_1 -> GetNumberOfIds() != 2){
			throw(std::runtime_error("In SBGATPolyhedronGravityModel.cpp: the intersection of the facet id lists should have exactly 2 items, not " + std::to_string(facet_ids_point_1 -> GetNumberOfIds())));
		}

		edge_points_ids_facet_ids[i][2] = facet_ids_point_1->GetId(0);
		edge_points_ids_facet_ids[i][3] = facet_ids_point_1->GetId(1);
		
	}


	// The edges dyads are created
	this -> edge_dyads = new double * [edge_points_ids_facet_ids.size()];
	this -> edges = new int * [edge_points_ids_facet_ids.size()];
	this -> edge_facets_ids = new int * [edge_points_ids_facet_ids.size()];


	#pragma omp parallel for
	for(unsigned int i = 0; i < edge_points_ids_facet_ids.size(); ++i) {
		this -> edge_dyads[i] = new double[9];
		this -> edges[i] = new int[2];
		this -> edge_facets_ids[i] = new int[2];

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

		this -> edge_facets_ids[i][0] = fA_index;
		this -> edge_facets_ids[i][1] = fB_index;


	}

	this -> N_edges = edge_count;
	this -> N_facets = numCells;

	this -> mass_properties = vtkSmartPointer<SBGATMassProperties>::New();
	this -> mass_properties -> SetInputData(input);
	if (this -> is_in_meters){
		this -> mass_properties -> SetScaleMeters();
	}
	else{
		this -> mass_properties -> SetScaleKiloMeters();
	}

	this -> mass_properties -> Update();

	// Check that the Euler characteristic == 2
	assert (input -> GetNumberOfPoints() - edge_count + numCells == 2);

	// Check that the shape is topologically closed
	assert(this -> mass_properties -> CheckClosed());

	return 1;
}


double SBGATPolyhedronGravityModel::GetPotential(const arma::vec::fixed<3> & point) const{
	return this -> GetPotential(point.colptr(0));
}


double SBGATPolyhedronGravityModel::GetPotential(double const * point) const{

	double potential = 0;

	// Facet loop
	#pragma omp parallel for reduction(+:potential)
	for (vtkIdType facet_index = 0; facet_index < this -> N_facets; ++ facet_index) {

		potential += this -> GetUf(this -> GetXf(point,facet_index));

	}

	// Edge loop
	#pragma omp parallel for reduction(+:potential)
	for (int edge_index = 0; edge_index < this -> N_edges; ++ edge_index) {

		potential += this -> GetUe(this -> GetXe(point,edge_index));

	}

	return potential * 0.5 * arma::datum::G * this -> density ;

}

bool SBGATPolyhedronGravityModel::Contains(const arma::vec::fixed<3> & point,double tol) const{
	return this -> Contains(point.colptr(0),tol);
}

bool SBGATPolyhedronGravityModel::Contains(double const * point, double tol ) const{

	double laplacian = 0;
	double point_scaled[3] = {point[0],point[1],point[2]};
	vtkMath::MultiplyScalar(point_scaled,1./this -> scaleFactor);



	// Facet loop
	#pragma omp parallel for reduction(+:laplacian)
	for (vtkIdType facet_index = 0; facet_index < this -> N_facets; ++ facet_index) {
		laplacian += this -> GetOmegaf(point_scaled, facet_index);
	}

	if (std::abs(laplacian) / (4 * arma::datum::pi) < tol) {
		return false;
	}
	else {
		return true;
	}

}


double SBGATPolyhedronGravityModel::GetEdgeLength(const int & e) const{

	return std::sqrt(vtkMath::Distance2BetweenPoints(this -> vertices[this -> edges[e][0]],this -> vertices[this -> edges[e][1]])) * this -> scaleFactor;
}


arma::vec::fixed<3> SBGATPolyhedronGravityModel::GetAcceleration(const arma::vec::fixed<3> & point) const{
	return this-> GetAcceleration(  point.colptr(0));
}

arma::vec::fixed<3> SBGATPolyhedronGravityModel::GetAcceleration(double const * point) const{

	double point_scaled[3] = {point[0],point[1],point[2]};
	vtkMath::MultiplyScalar(point_scaled,1./this -> scaleFactor);

	arma::vec::fixed<3> acc = arma::zeros<arma::vec>(3);

	// Facet loop
	#pragma omp parallel for reduction(+:acc)
	for (vtkIdType facet_index = 0; facet_index < this -> N_facets; ++ facet_index) {

		double * r0 = this -> vertices[this -> facets[facet_index][0]];

		double r0m[3];
		
		vtkMath::Subtract(r0,point_scaled,r0m);

		double wf = this -> GetOmegaf( point, facet_index);
		double * F = this -> facet_dyads[facet_index];

		acc(0) += wf *( F[0] * r0m[0] + F[1] * r0m[1] +  F[2] * r0m[2]);
		acc(1) += wf *( F[3] * r0m[0] + F[4] * r0m[1] +  F[5] * r0m[2]);
		acc(2) += wf *( F[6] * r0m[0] + F[7] * r0m[1] +  F[8] * r0m[2]);

	}

	// Edge loop
	#pragma omp parallel for reduction(+:acc)
	for (int edge_index = 0; edge_index < this -> N_edges; ++ edge_index) {

		double * r0 = this -> vertices[this -> edges[edge_index][0]];
		
		double r0m[3];
		
		vtkMath::Subtract(r0,point_scaled,r0m);
		
		double Le = this -> GetLe( point, edge_index);

		double * E = this -> edge_dyads[edge_index];

		acc(0) += - Le *( E[0] * r0m[0] + E[1] * r0m[1] +  E[2] * r0m[2]);
		acc(1) += - Le *( E[3] * r0m[0] + E[4] * r0m[1] +  E[5] * r0m[2]);
		acc(2) += - Le *( E[6] * r0m[0] + E[7] * r0m[1] +  E[8] * r0m[2]);

	}

	return acc * arma::datum::G * this -> density * this -> scaleFactor;

}


void SBGATPolyhedronGravityModel::GetPotentialAcceleration(const arma::vec::fixed<3> & point,double & potential, 
	arma::vec::fixed<3> & acc) const{

	this -> GetPotentialAcceleration(point.colptr(0),potential, acc);

}



void SBGATPolyhedronGravityModel::GetPotentialAcceleration(double const  * point,double & potential, 
	arma::vec::fixed<3> & acc) const {

	double pot = 0;
	double acc_x = 0;
	double acc_y = 0;
	double acc_z = 0;

	double point_scaled[3] = {point[0],point[1],point[2]};
	vtkMath::MultiplyScalar(point_scaled,1./this -> scaleFactor);

	// Facet loop
	#pragma omp parallel for reduction(+:acc_x,acc_y,acc_z,pot)
	for (vtkIdType facet_index = 0; facet_index < this -> N_facets; ++ facet_index) {

		double * r0 = this -> vertices[this -> facets[facet_index][0]];
		
		double r0m[3];

		vtkMath::Subtract(r0,point_scaled,r0m);
		
		double wf = this -> GetOmegaf( point, facet_index);

		double * F = this -> facet_dyads[facet_index];

		double a[3] = {
			F[0] * r0m[0] + F[1] * r0m[1] +  F[2] * r0m[2],
			F[3] * r0m[0] + F[4] * r0m[1] +  F[5] * r0m[2],
			F[6] * r0m[0] + F[7] * r0m[1] +  F[8] * r0m[2]
		};

		acc_x += wf * a[0];
		acc_y += wf * a[1];
		acc_z += wf * a[2];

		pot += - wf * vtkMath::Dot(r0m,a);


	}

	// Edge loop
	#pragma omp parallel for reduction(-:acc_x,acc_y,acc_z,pot)
	for (int edge_index = 0; edge_index < this -> N_edges; ++ edge_index) {

		double * r0 = this -> vertices[this -> edges[edge_index][0]];

		
		double r0m[3];
		
		vtkMath::Subtract(r0,point_scaled,r0m);
		
		
		double Le = this -> GetLe( point, edge_index);

		double * E = this -> edge_dyads[edge_index];

		double a[3] = {
			E[0] * r0m[0] + E[1] * r0m[1] +  E[2] * r0m[2],
			E[3] * r0m[0] + E[4] * r0m[1] +  E[5] * r0m[2],
			E[6] * r0m[0] + E[7] * r0m[1] +  E[8] * r0m[2]
		};


		pot += Le * vtkMath::Dot(r0m,a);

		acc_x -= Le * a[0];
		acc_y -= Le * a[1];
		acc_z -= Le * a[2];

	}

	acc(0) = acc_x;
	acc(1) = acc_y;
	acc(2) = acc_z;

	acc *= arma::datum::G  * this -> density * this -> scaleFactor ;
	pot *= 0.5 * arma::datum::G * this -> density * this -> scaleFactor* this -> scaleFactor;

	potential = pot;

}

void SBGATPolyhedronGravityModel::GetPotentialAccelerationGravityGradient(const arma::vec::fixed<3> & point,double & potential, 
	arma::vec::fixed<3> & acc,arma::mat::fixed<3,3> & gravity_gradient_mat) const{

	this -> GetPotentialAccelerationGravityGradient(point.colptr(0),potential, acc,gravity_gradient_mat);

}

void SBGATPolyhedronGravityModel::GetPotentialAccelerationGravityGradient(double const  * point,double & potential, 
	arma::vec::fixed<3> & acc,arma::mat::fixed<3,3> & gravity_gradient_mat) const{

	double point_scaled[3] = {point[0],point[1],point[2]};
	vtkMath::MultiplyScalar(point_scaled,1./this -> scaleFactor);




	double pot = 0;
	double acc_x = 0;
	double acc_y = 0;
	double acc_z = 0;
	double grav_mat_acc_xx,grav_mat_acc_xy,grav_mat_acc_xz;
	double grav_mat_acc_yx,grav_mat_acc_yy,grav_mat_acc_yz;
	double grav_mat_acc_zx,grav_mat_acc_zy,grav_mat_acc_zz;

	grav_mat_acc_xx = 0;
	grav_mat_acc_xy = 0;
	grav_mat_acc_xz = 0;
	grav_mat_acc_yx = 0;
	grav_mat_acc_yy = 0;
	grav_mat_acc_yz = 0;
	grav_mat_acc_zx = 0;
	grav_mat_acc_zy = 0;
	grav_mat_acc_zz = 0;

	gravity_gradient_mat = arma::zeros<arma::mat>(3,3);

	// Facet loop
	#pragma omp parallel for reduction(+:acc_x,acc_y,acc_z), reduction(-:pot,grav_mat_acc_xx,grav_mat_acc_xy,grav_mat_acc_xz,grav_mat_acc_yx,grav_mat_acc_yy,grav_mat_acc_yz,grav_mat_acc_zx,grav_mat_acc_zy,grav_mat_acc_zz)
	for (vtkIdType facet_index = 0; facet_index < this -> N_facets; ++ facet_index) {

		double * r0 = this -> vertices[this -> facets[facet_index][0]];
		
		double r0m[3];

		vtkMath::Subtract(r0,point_scaled,r0m);

		double wf = this -> GetOmegaf( point, facet_index);

		double * F = this -> facet_dyads[facet_index];

		double a[3] = {
			F[0] * r0m[0] + F[1] * r0m[1] +  F[2] * r0m[2],
			F[3] * r0m[0] + F[4] * r0m[1] +  F[5] * r0m[2],
			F[6] * r0m[0] + F[7] * r0m[1] +  F[8] * r0m[2]
		};

		arma::mat::fixed<3,3> F_arma = {
			{F[0],F[1],F[2]},
			{F[3],F[4],F[5]},
			{F[6],F[7],F[8]}
		};

		acc_x += wf * a[0];
		acc_y += wf * a[1];
		acc_z += wf * a[2];

		pot -= wf * vtkMath::Dot(r0m,a);


		grav_mat_acc_xx -= F_arma(0,0) * wf;
		grav_mat_acc_yx -= F_arma(1,0) * wf;
		grav_mat_acc_zx -= F_arma(2,0) * wf;

		grav_mat_acc_xy -= F_arma(0,1) * wf;
		grav_mat_acc_yy -= F_arma(1,1) * wf;
		grav_mat_acc_zy -= F_arma(2,1) * wf;

		grav_mat_acc_xz -= F_arma(0,2) * wf;
		grav_mat_acc_yz -= F_arma(1,2) * wf;
		grav_mat_acc_zz -= F_arma(2,2) * wf;


	}

	// Edge loop
	#pragma omp parallel for reduction(-:acc_x,acc_y,acc_z), reduction(+:pot,grav_mat_acc_xz,grav_mat_acc_xx,grav_mat_acc_xy,grav_mat_acc_yx,grav_mat_acc_yy,grav_mat_acc_yz,grav_mat_acc_zx,grav_mat_acc_zy,grav_mat_acc_zz)
	for (int edge_index = 0; edge_index < this -> N_edges; ++ edge_index) {

		double * r0 = this -> vertices[this -> edges[edge_index][0]];


		double r0m[3];

		vtkMath::Subtract(r0,point_scaled,r0m);

		double Le = this -> GetLe( point, edge_index);

		double * E = this -> edge_dyads[edge_index];

		double a[3] = {
			E[0] * r0m[0] + E[1] * r0m[1] +  E[2] * r0m[2],
			E[3] * r0m[0] + E[4] * r0m[1] +  E[5] * r0m[2],
			E[6] * r0m[0] + E[7] * r0m[1] +  E[8] * r0m[2]
		};

		arma::mat::fixed<3,3> E_arma = {
			{E[0],E[1],E[2]},
			{E[3],E[4],E[5]},
			{E[6],E[7],E[8]}
		};

		pot += Le * vtkMath::Dot(r0m,a);

		acc_x -= Le * a[0];
		acc_y -= Le * a[1];
		acc_z -= Le * a[2];

		grav_mat_acc_xx += E_arma(0,0) * Le;
		grav_mat_acc_xy += E_arma(0,1) * Le;
		grav_mat_acc_xz += E_arma(0,2) * Le;

		grav_mat_acc_yx += E_arma(1,0) * Le;
		grav_mat_acc_yy += E_arma(1,1) * Le;
		grav_mat_acc_yz += E_arma(1,2) * Le;

		grav_mat_acc_zx += E_arma(2,0) * Le;
		grav_mat_acc_zy += E_arma(2,1) * Le;
		grav_mat_acc_zz += E_arma(2,2) * Le;

		

	}

	acc(0) = acc_x;
	acc(1) = acc_y;
	acc(2) = acc_z;

	gravity_gradient_mat = {
		{grav_mat_acc_xx,grav_mat_acc_xy,grav_mat_acc_xz},
		{grav_mat_acc_yx,grav_mat_acc_yy,grav_mat_acc_yz},
		{grav_mat_acc_zx,grav_mat_acc_zy,grav_mat_acc_zz}
	};

	acc *= arma::datum::G  * this -> density * this -> scaleFactor ;
	pot *= 0.5 * arma::datum::G * this -> density *  this -> scaleFactor *  this -> scaleFactor;
	gravity_gradient_mat *= arma::datum::G  * this -> density;

	potential = pot;

}


void SBGATPolyhedronGravityModel::PrintHeader(ostream& os, vtkIndent indent) {

}
void SBGATPolyhedronGravityModel::PrintTrailer(ostream& os, vtkIndent indent) {

}


void SBGATPolyhedronGravityModel::Clear(){
	if (this -> N_facets > 0){
		int N_vertices = this -> N_edges - this -> N_facets + 2;

	//Vertices
		for(int i = 0; i < N_vertices; ++i) {
			delete[] this -> vertices[i];   
		}
		delete[] this -> vertices;

	//Facet dyads
		for(int i = 0; i < this -> N_facets; ++i) {
			delete[] this -> facet_dyads[i];   
		}
		delete[] this -> facet_dyads;

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


	//Edge dyads
		for(int i = 0; i < this -> N_edges; ++i) {
			delete[] this -> edge_dyads[i];   
		}
		delete[] this -> edge_dyads;

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
}



//----------------------------------------------------------------------------
void SBGATPolyhedronGravityModel::PrintSelf(std::ostream& os, vtkIndent indent){

	vtkPolyData *input = vtkPolyData::SafeDownCast(this->GetInput(0));
	if (!input)
	{
		return;
	}
}









void SBGATPolyhedronGravityModel::ComputeSurfacePGM(
	vtkSmartPointer<vtkPolyData> selected_shape,
	const std::vector<unsigned int> & queried_elements,
	bool is_in_meters,
	double density,
	const arma::vec::fixed<3> & omega,
	std::vector<double> & slopes,
	std::vector<double> & inertial_potentials,
	std::vector<double> & body_fixed_potentials,
	std::vector<double> & inertial_acc_magnitudes,
	std::vector<double> & body_fixed_acc_magnitudes){

	vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
	pgm_filter -> SetInputData(selected_shape);
	pgm_filter -> SetDensity(density);

	vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
	mass_prop -> SetInputData(selected_shape);

	if (is_in_meters){
		pgm_filter -> SetScaleMeters();
		mass_prop -> SetScaleMeters();
	}
	else{
		pgm_filter -> SetScaleKiloMeters();
		mass_prop -> SetScaleKiloMeters();
	}

	pgm_filter -> Update();
	mass_prop -> Update();

	slopes.clear();
	inertial_potentials.clear();
	body_fixed_potentials.clear();
	inertial_acc_magnitudes.clear();
	body_fixed_acc_magnitudes.clear();
	
	for (int i = 0; i < static_cast<int>(selected_shape -> GetNumberOfCells()); ++i){

		slopes.push_back(std::numeric_limits<double>::quiet_NaN());
		inertial_potentials.push_back(std::numeric_limits<double>::quiet_NaN());
		body_fixed_potentials.push_back(std::numeric_limits<double>::quiet_NaN());
		inertial_acc_magnitudes.push_back(std::numeric_limits<double>::quiet_NaN());
		body_fixed_acc_magnitudes.push_back(std::numeric_limits<double>::quiet_NaN());

	}


	arma::vec::fixed<3> com = mass_prop -> GetCenterOfMass();

    // The queried facets are browsed and their surface PGM evaluated
    #pragma omp parallel for
	for ( int e = 0; e <  queried_elements.size(); ++e){

		int el = queried_elements[e];

		unsigned int cellId = queried_elements.at(el);

		double p0[3];
		double p1[3];
		double p2[3];

		pgm_filter -> GetVerticesInFacet(el,p0,p1, p2);


		arma::vec::fixed<3> p0_arma = {p0[0],p0[1],p0[2]};
		arma::vec::fixed<3> p1_arma = {p1[0],p1[1],p1[2]};
		arma::vec::fixed<3> p2_arma = {p2[0],p2[1],p2[2]};
		arma::vec::fixed<3> facet_center = 1./3 * (p0_arma + p1_arma + p2_arma);
		arma::vec::fixed<3> normal = arma::normalise(arma::cross(p1_arma - p0_arma,p2_arma - p0_arma));

		double potential,slope;
		arma::vec::fixed<3> acc,acc_body_fixed;

		pgm_filter -> GetPotentialAcceleration(facet_center * pgm_filter -> GetScaleFactor(),potential,acc);
		acc_body_fixed = acc - arma::cross(omega,arma::cross(omega,facet_center * pgm_filter -> GetScaleFactor() - com));

		slope = std::acos(arma::dot(-arma::normalise(acc_body_fixed),normal)) * 180./arma::datum::pi;

		slopes[cellId] = slope;
		inertial_potentials[cellId] = potential;
		body_fixed_potentials[cellId] = potential + 0.5 * arma::dot(RBK::tilde(omega) * (facet_center * pgm_filter -> GetScaleFactor() - com),
			RBK::tilde(omega) * (facet_center * pgm_filter -> GetScaleFactor() - com));
		inertial_acc_magnitudes[cellId] = arma::norm(acc);
		body_fixed_acc_magnitudes[cellId] = arma::norm(acc_body_fixed);

	}	

}

void SBGATPolyhedronGravityModel::SaveSurfacePGM(vtkSmartPointer<vtkPolyData> selected_shape,
	const std::vector<unsigned int> & queried_elements,
	bool is_in_meters,
	const double & mass,
	const arma::vec::fixed<3> & omega,
	const std::vector<double> & slopes,
	const std::vector<double> & inertial_potentials,
	const std::vector<double> & body_fixed_potentials,
	const std::vector<double> & inertial_acc_magnitudes,
	const std::vector<double> & body_fixed_acc_magnitudes,
	std::string path){


  // The surface gravity model is saved to a JSON file
  // The JSON fieds are:
  // - facets == number of facets of the input shape
  // - vertices == number of vertices of the input shape
  // - omega : { value : {omega_x,omega_y,omega_z},unit}
  // - slopes - vector of slopes : {index, value, unit}
  // - inertial_potentials - vector of inertial potentials : {index, value, unit}
  // - body_fixed_potentials - vector of body-fixed potentials : {index, value, unit}
  // - inertial_acc_magnitudes - vector of inertial acc_magnitudes : {index, value, unit}
  // - body_fixed_acc_magnitudes - vector of body-fixed acc_magnitudes : {index, value, unit}


	nlohmann::json surface_pgm_json;
	nlohmann::json omega_json = {
		{"value",{omega(0),omega(1),omega(2)}},
		{"unit","rad/s"}
	};


	surface_pgm_json["facets"] = selected_shape -> GetNumberOfCells();
	surface_pgm_json["vertices"] = selected_shape -> GetNumberOfPoints();

	std::string distance_unit,potential_unit,acceleration_unit;
	distance_unit = "m";
	potential_unit = "m^2/s^2";
	acceleration_unit = "m/s^2";

	

	nlohmann::json slopes_json,
	inertial_potentials_json,
	body_fixed_potentials_json,
	inertial_acc_magnitudes_json,
	body_fixed_acc_magnitudes_json;

	for (int i = 0; i < queried_elements.size(); ++i){
		int index = queried_elements[i];
		
		nlohmann::json slope = { 
			{"index", index}, 
			{"value", slopes[i]},
			{"unit","deg"} 
		};

		nlohmann::json inertial_potential = { 
			{"index", index}, 
			{"value", inertial_potentials[i]},
			{"unit",potential_unit} 
		};

		nlohmann::json body_fixed_potential = { 
			{"index", index}, 
			{"value", body_fixed_potentials[i]},
			{"unit",potential_unit} 
		};


		nlohmann::json inertial_acc_magnitude = { 
			{"index", index}, 
			{"value", inertial_acc_magnitudes[i]},
			{"unit",acceleration_unit} 
		};
		nlohmann::json body_fixed_acc_magnitude = { 
			{"index", index}, 
			{"value", body_fixed_acc_magnitudes[i]},
			{"unit",acceleration_unit} 
		};


		slopes_json.push_back(slope);
		inertial_potentials_json.push_back(inertial_potential);
		body_fixed_potentials_json.push_back(body_fixed_potential);
		inertial_acc_magnitudes_json.push_back(inertial_acc_magnitude);
		body_fixed_acc_magnitudes_json.push_back(body_fixed_acc_magnitude);

	}

	surface_pgm_json["mass"] = mass;
	surface_pgm_json["omega"] = omega_json;
	surface_pgm_json["slopes"] = slopes_json;
	surface_pgm_json["inertial_acc_magnitudes"] = inertial_acc_magnitudes_json;
	surface_pgm_json["body_fixed_acc_magnitudes"] = body_fixed_acc_magnitudes_json;
	surface_pgm_json["inertial_potentials"] = inertial_potentials_json;
	surface_pgm_json["body_fixed_potentials"] = body_fixed_potentials_json;



	std::ofstream o(path);
	o << std::setw(4) << surface_pgm_json << std::endl;


}



void SBGATPolyhedronGravityModel::LoadSurfacePGM(double & mass,
	arma::vec::fixed<3> & omega,
	std::vector<double> & slopes,
	std::vector<double> & inertial_potentials,
	std::vector<double> & body_fixed_potentials,
	std::vector<double> & inertial_acc_magnitudes,
	std::vector<double> & body_fixed_acc_magnitudes,
	std::string path){

  // The surface gravity model is saved to a JSON file
  // The JSON fieds are:
  // - facets == number of facets of the input shape
  // - vertices == number of vertices of the input shape
  // - omega : { value : {omega_x,omega_y,omega_z},unit}
  // - slopes - vector of slopes : {index, value, unit}
  // - inertial_potentials - vector of inertial potentials : {index, value, unit}
  // - body_fixed_potentials - vector of body-fixed potentials : {index, value, unit}
  // - inertial_acc_magnitudes - vector of inertial acc_magnitudes : {index, value, unit}
  // - body_fixed_acc_magnitudes - vector of body-fixed acc_magnitudes : {index, value, unit}

	// The JSON container is created
	nlohmann::json surface_pgm_json;

  // The file is loaded into the container
	std::ifstream i(path);  
	i >> surface_pgm_json;


	slopes.resize(static_cast<int>(surface_pgm_json["facets"]));
	inertial_potentials.resize(static_cast<int>(surface_pgm_json["facets"]));
	body_fixed_potentials.resize(static_cast<int>(surface_pgm_json["facets"]));

	inertial_acc_magnitudes.resize(static_cast<int>(surface_pgm_json["facets"]));
	body_fixed_acc_magnitudes.resize(static_cast<int>(surface_pgm_json["facets"]));


	for (int k = 0; k < slopes.size(); ++k){
		slopes[k] = std::numeric_limits<double>::quiet_NaN();
		inertial_potentials[k] = std::numeric_limits<double>::quiet_NaN();
		body_fixed_potentials[k] = std::numeric_limits<double>::quiet_NaN();
		inertial_acc_magnitudes[k] = std::numeric_limits<double>::quiet_NaN();
		body_fixed_acc_magnitudes[k] = std::numeric_limits<double>::quiet_NaN();		
	}


	nlohmann::json omega_json;

	try{
		omega_json = surface_pgm_json.at("omega");
		omega(0) = omega_json["value"][0];
		omega(1) = omega_json["value"][1];
		omega(2) = omega_json["value"][2];

	}
	catch (nlohmann::detail::parse_error & e){
		throw(std::runtime_error("Error loading angular velocity in SBGATPolyhedronGravityModel::LoadSurfacePGM. Can't find field `omega`"));
	}

	try{
		mass = surface_pgm_json.at("mass");
	}
	catch (nlohmann::detail::parse_error & e){
		throw(std::runtime_error("Error loading mass in SBGATPolyhedronGravityModel::LoadSurfacePGM. Can't find field `mass`"));
	}


	nlohmann::json slopes_json;
	
	try{
		slopes_json = surface_pgm_json.at("slopes");
		for (auto slope : slopes_json){
			slopes[slope["index"]] = slope["value"];
		}
	}
	catch (nlohmann::detail::parse_error & e){
		throw(std::runtime_error("Error loading slopes in SBGATPolyhedronGravityModel::LoadSurfacePGM. Can't find field `slopes`"));
	}

	nlohmann::json inertial_potentials_json;
	try{
		inertial_potentials_json = surface_pgm_json.at("inertial_potentials");
		for (auto inertial_potential : inertial_potentials_json){
			inertial_potentials[inertial_potential["index"]] = inertial_potential["value"];
			
		}
	}
	catch (nlohmann::detail::parse_error & e){
		throw(std::runtime_error("Error loading inertial potentials in SBGATPolyhedronGravityModel::LoadSurfacePGM. Can't find field `inertial_potentials`"));
	}

	nlohmann::json body_fixed_potentials_json;
	try{
		body_fixed_potentials_json = surface_pgm_json.at("body_fixed_potentials");
		for (auto body_fixed_potential : body_fixed_potentials_json){
			body_fixed_potentials[body_fixed_potential["index"]] = body_fixed_potential["value"];
			
		}
	}
	catch (nlohmann::detail::parse_error & e){
		throw(std::runtime_error("Error loading inertial potentials in SBGATPolyhedronGravityModel::LoadSurfacePGM. Can't find field `body_fixed_potentials`"));
	}

	nlohmann::json inertial_acc_json;

	try{
		inertial_acc_json = surface_pgm_json.at("inertial_acc_magnitudes");
		for (auto acc : inertial_acc_json){
			inertial_acc_magnitudes[acc["index"]] = acc["value"];
			
		}

	}
	catch (nlohmann::detail::parse_error & e){
		throw(std::runtime_error("Error loading inertial accelerations in SBGATPolyhedronGravityModel::LoadSurfacePGM. Can't find field `inertial_acc_magnitudes`"));
	}


	nlohmann::json body_fixed_acc_magnitudes_json;
	try{
		body_fixed_acc_magnitudes_json = surface_pgm_json.at("body_fixed_acc_magnitudes");
		for (auto acc : body_fixed_acc_magnitudes_json){
			body_fixed_acc_magnitudes[acc["index"]] = acc["value"];
			
		}
	}	
	catch (nlohmann::detail::parse_error & e){
		throw(std::runtime_error("Error loading body-fixed accelerations in SBGATPolyhedronGravityModel::LoadSurfacePGM. Can't find field `body_fixed_acc_magnitudes`"));
	}


}

double SBGATPolyhedronGravityModel::GetOmegaf(const arma::vec::fixed<3> & pos, const int & f) const{
	return this -> GetOmegaf(pos.colptr(0),f);
}

double SBGATPolyhedronGravityModel::GetOmegaf( const double * pos, const int & f) const{

	const double * r0 = this -> vertices[this -> facets[f][0]];
	const double * r1 = this -> vertices[this -> facets[f][1]];
	const double * r2 = this -> vertices[this -> facets[f][2]];

	double r0m[3];
	double r1m[3];
	double r2m[3];

	double pos_scaled[3] = {pos[0],pos[1],pos[2]};
	vtkMath::MultiplyScalar(pos_scaled,1./this -> GetScaleFactor());

	vtkMath::Subtract(r0,pos_scaled,r0m);
	vtkMath::Subtract(r1,pos_scaled,r1m);
	vtkMath::Subtract(r2,pos_scaled,r2m);

	double R0 = vtkMath::Norm(r0m);
	double R1 = vtkMath::Norm(r1m);
	double R2 = vtkMath::Norm(r2m);

	double r1m_cross_r2m[3];

	vtkMath::Cross(r1m,r2m,r1m_cross_r2m);

	return 2 * std::atan2(vtkMath::Dot(r0m,r1m_cross_r2m),R0 * R1 * R2 + 
		R0 * vtkMath::Dot(r1m,r2m) + R1 * vtkMath::Dot(r0m,r2m) + R2 * vtkMath::Dot(r0m,r1m));


}

double SBGATPolyhedronGravityModel::GetLe(const arma::vec::fixed<3> & pos, const int & e) const{
	return this -> GetLe(pos.colptr(0),e);

}

double SBGATPolyhedronGravityModel::GetLe( const double * pos, const int & e) const{

	double * r0 = this -> vertices[this -> edges[e][0]];
	double * r1 = this -> vertices[this -> edges[e][1]];


	double r0m[3];
	double r1m[3];
	double rem[3];
	double pos_scaled[3] = {pos[0],pos[1],pos[2]};
	vtkMath::MultiplyScalar(pos_scaled,1./this -> GetScaleFactor());

	vtkMath::Subtract(r0,pos_scaled,r0m);
	vtkMath::Subtract(r1,pos_scaled,r1m);
	vtkMath::Subtract(r1m,r0m,rem);

	double R0 = vtkMath::Norm(r0m);
	double R1 = vtkMath::Norm(r1m);
	double Re = vtkMath::Norm(rem);

	return std::log((R0 + R1 + Re) / (R0 + R1 - Re));

}


arma::vec::fixed<3> SBGATPolyhedronGravityModel::GetRe(const arma::vec::fixed<3> & pos,const int & e) const{
	
	return this -> GetRe(pos.colptr(0),e);
}

arma::vec::fixed<3> SBGATPolyhedronGravityModel::GetRe(const double * pos,const int & e) const{

	
	double re[3];
	double pos_scaled[3] = {pos[0],pos[1],pos[2]};
	vtkMath::MultiplyScalar(pos_scaled,1./this -> GetScaleFactor());


	vtkMath::Subtract(this -> vertices[this -> edges[e][0]],pos_scaled,re);

	return this -> GetScaleFactor() * arma::vec({re[0],re[1],re[2]});
}

arma::vec::fixed<10> SBGATPolyhedronGravityModel::GetXe(const arma::vec::fixed<3> & pos,const int & e) const{

	arma::vec::fixed<10> Xe;

	Xe(0) = this -> GetLe(pos,e);
	Xe.subvec(1,3) = this -> GetRe(pos,e);
	Xe.subvec(4,9) = this -> GetEeParam(e);
	return Xe;

}





arma::vec::fixed<6> SBGATPolyhedronGravityModel::GetEeParam(const int & e) const{


	const double * E = this -> edge_dyads[e];

	

	return {E[0],E[4],E[8],E[1],E[2],E[5]};

}

double SBGATPolyhedronGravityModel::GetUe(const arma::vec::fixed<10> & Xe){

   // {E[0],E[4],E[8],E[1],E[2],E[5]};

	const arma::vec::fixed<6> & E_vec = Xe.subvec(4,9);
	const arma::vec::fixed<3> & r_ei_0= Xe.subvec(1,3);
	const double & Le = Xe(0);

	arma::mat::fixed<3,3> Ee = {
		{E_vec[0],E_vec[3],E_vec[4]},
		{E_vec[3],E_vec[1],E_vec[5]},
		{E_vec[4],E_vec[5],E_vec[2]}
	};

	return Le * arma::dot(r_ei_0,Ee * r_ei_0);

}

arma::vec::fixed<3> SBGATPolyhedronGravityModel::GetAe(const arma::vec::fixed<10> & Xe){


	const arma::vec::fixed<6> & E_vec = Xe.subvec(4,9);
	const arma::vec::fixed<3> & r_ei_0= Xe.subvec(1,3);
	const double & Le = Xe(0);

	arma::mat::fixed<3,3> Ee = {
		{E_vec[0],E_vec[3],E_vec[4]},
		{E_vec[3],E_vec[1],E_vec[5]},
		{E_vec[4],E_vec[5],E_vec[2]}
	};

	return - Le * Ee * r_ei_0;

}





arma::vec::fixed<3> SBGATPolyhedronGravityModel::GetRf(const arma::vec::fixed<3> & pos,const int & f) const{
	
	return this -> GetRf(pos.colptr(0),f);
}

arma::vec::fixed<3> SBGATPolyhedronGravityModel::GetRf(const double * pos,const int & f) const{

	
	double rf[3];
	double pos_scaled[3] = {pos[0],pos[1],pos[2]};
	vtkMath::MultiplyScalar(pos_scaled,1./this -> GetScaleFactor());

	vtkMath::Subtract(this -> vertices[this -> facets[f][0]],pos_scaled,rf);

	return arma::vec({rf[0],rf[1],rf[2]}) * this -> GetScaleFactor();
}

arma::vec::fixed<10> SBGATPolyhedronGravityModel::GetXf(const arma::vec::fixed<3> & pos,const int & f) const{

	arma::vec::fixed<10> Xf;

	Xf(0) = this -> GetOmegaf(pos,f);
	Xf.subvec(1,3) = this -> GetRf(pos,f);
	Xf.subvec(4,9) = this -> GetFfParam(f);
	return Xf;

}





arma::vec::fixed<6> SBGATPolyhedronGravityModel::GetFfParam(const int & f) const{


	//  this -> facet_dyads[i][0] = normal[0] * normal[0];
	// 	this -> facet_dyads[i][1] = normal[0] * normal[1];
	// 	this -> facet_dyads[i][2] = normal[0] * normal[2];
	// 	this -> facet_dyads[i][3] = normal[1] * normal[0];
	// 	this -> facet_dyads[i][4] = normal[1] * normal[1];
	// 	this -> facet_dyads[i][5] = normal[1] * normal[2];
	// 	this -> facet_dyads[i][6] = normal[2] * normal[0];
	// 	this -> facet_dyads[i][7] = normal[2] * normal[1];
	// 	this -> facet_dyads[i][8] = normal[2] * normal[2];



	const double * F = this -> facet_dyads[f];

	return {F[0],F[4],F[8],F[1],F[2],F[5]};

}

double SBGATPolyhedronGravityModel::GetUf(const arma::vec::fixed<10> & Xf){

   // {F[0],F[4],F[8],F[1],F[2],F[5]};

	const arma::vec::fixed<6> & F_vec = Xf.subvec(4,9);
	const arma::vec::fixed<3> & r_fi_0= Xf.subvec(1,3);
	const double & omega_f = Xf(0);

	arma::mat::fixed<3,3> F = {
		{F_vec[0],F_vec[3],F_vec[4]},
		{F_vec[3],F_vec[1],F_vec[5]},
		{F_vec[4],F_vec[5],F_vec[2]}
	};

	return - omega_f * arma::dot(r_fi_0,F * r_fi_0);


}

arma::vec::fixed<3> SBGATPolyhedronGravityModel::GetAf(const arma::vec::fixed<10> & Xf){

   // {F[0],F[4],F[8],F[1],F[2],F[5]};

	const arma::vec::fixed<6> & F_vec = Xf.subvec(4,9);
	const arma::vec::fixed<3> & r_fi_0= Xf.subvec(1,3);
	const double & omega_f = Xf(0);

	arma::mat::fixed<3,3> F = {
		{F_vec[0],F_vec[3],F_vec[4]},
		{F_vec[3],F_vec[1],F_vec[5]},
		{F_vec[4],F_vec[5],F_vec[2]}
	};

	return omega_f * F * r_fi_0;


}

arma::vec::fixed<3> SBGATPolyhedronGravityModel::GetNonNormalizedFacetNormal(const int & f) const{

	double r0[3], r1[3], r2[3];

	this -> GetVerticesInFacet(f,r0,r1,r2);

	arma::vec::fixed<3> r0_arma = {r0[0],r0[1],r0[2]};
	arma::vec::fixed<3> r1_arma = {r1[0],r1[1],r1[2]};
	arma::vec::fixed<3> r2_arma = {r2[0],r2[1],r2[2]};
	return arma::cross(r1_arma - r0_arma,r2_arma - r1_arma) * std::pow(this -> scaleFactor,2);

}


void SBGATPolyhedronGravityModel::GetVerticesInFacet(const int & f,double * r0,double * r1, double * r2) const{

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


void SBGATPolyhedronGravityModel::GetVerticesOnEdge(const int & e,double * r0,double * r1) const{

	r0[0] = this -> vertices[this -> edges[e][0]][0];
	r0[1] = this -> vertices[this -> edges[e][0]][1];
	r0[2] = this -> vertices[this -> edges[e][0]][2];

	r1[0] = this -> vertices[this -> edges[e][1]][0];
	r1[1] = this -> vertices[this -> edges[e][1]][1];
	r1[2] = this -> vertices[this -> edges[e][1]][2];

}


void SBGATPolyhedronGravityModel::GetIndicesOfAdjacentFacets(const int & e,int & f0, int & f1) const{

	f0 = this -> edge_facets_ids[e][0];
	f1 = this -> edge_facets_ids[e][1];

}

arma::vec::fixed<3> SBGATPolyhedronGravityModel::GetFacetCenter(const int & f) const{
	
	double r0[3];
	double r1[3];
	double r2[3];
	this -> GetVerticesInFacet(f,r0,r1,r2);

	return 1./3 * this -> scaleFactor * arma::vec({r0[0] + r1[0] + r2[0], r0[1] + r1[1] + r2[1], r0[2] + r1[2] + r2[2] } );


}


arma::vec::fixed<3> SBGATPolyhedronGravityModel::GetBodyFixedAccelerationf(const int & f,const arma::vec::fixed<3> & Omega) const{

	return (this -> GetAcceleration(this -> GetFacetCenter(f)) - RBK::tilde(Omega) * RBK::tilde(Omega) * (this -> GetFacetCenter(f) - this -> mass_properties -> GetCenterOfMass()));


}


double SBGATPolyhedronGravityModel::GetSlope(const int & f ,const arma::vec::fixed<3> & Omega) const{

	return std::acos(arma::dot(-arma::normalise(this -> GetBodyFixedAccelerationf(f,Omega)),arma::normalise(this -> GetNonNormalizedFacetNormal(f))));
}

arma::mat::fixed<3,3> SBGATPolyhedronGravityModel::GetGravityGradient(const arma::vec::fixed<3> & point) const {



	double point_scaled[3] = {point[0],point[1],point[2]};
	vtkMath::MultiplyScalar(point_scaled,1./this -> scaleFactor);
	arma::mat::fixed<3,3> gravity_gradient = arma::zeros<arma::mat>(3,3);

	// Facet loop
	#pragma omp parallel for reduction(+:gravity_gradient)
	for (vtkIdType facet_index = 0; facet_index < this -> N_facets; ++ facet_index) {

		double * r0 = this -> vertices[this -> facets[facet_index][0]];
		
		double r0m[3];

		vtkMath::Subtract(r0,point_scaled,r0m);

		double wf = this -> GetOmegaf( point, facet_index);

		double * F = this -> facet_dyads[facet_index];

		
		arma::mat::fixed<3,3> F_arma = {
			{F[0],F[1],F[2]},
			{F[3],F[4],F[5]},
			{F[6],F[7],F[8]}
		};

		gravity_gradient += (-w_f * F_arma);


	}

	// Edge loop
	#pragma omp parallel for reduction(+:gravity_gradient)
	for (int edge_index = 0; edge_index < this -> N_edges; ++ edge_index) {

		double * r0 = this -> vertices[this -> edges[edge_index][0]];

		double r0m[3];

		vtkMath::Subtract(r0,point_scaled,r0m);

		double Le = this -> GetLe( point, edge_index);

		double * E = this -> edge_dyads[edge_index];

		arma::mat::fixed<3,3> E_arma = {
			{E[0],E[1],E[2]},
			{E[3],E[4],E[5]},
			{E[6],E[7],E[8]}
		};


		gravity_gradient += Le * E_arma;

		

	}

	return gravity_gradient * arma::datum::G  * this -> density;

}



