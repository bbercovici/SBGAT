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

	//Free each sub-array
	for(int i = 0; i < sizeof(this -> facet_dyads) / sizeof(this -> facet_dyads[0]); ++i) {
		delete[] this -> facet_dyads[i];   
	}
    //Free the array of pointers
	delete[] this -> facet_dyads;
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

	input = normalGenerator -> GetOutput();
	// Required by vtkPolyData::GetPointCell	
	input -> BuildLinks();

	vtkFloatArray * normals =  vtkFloatArray::SafeDownCast(input->GetCellData()->GetArray("Normals"));
	
	auto start = std::chrono::system_clock::now();

	// The facet dyads are created
	this -> facet_dyads = new double * [numCells];
	#pragma omp parallel for
	for(int i = 0; i < numCells; ++i) {
		this -> facet_dyads[i] = new double[9];
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
	}


	// The edge dyads are created

	vtkSmartPointer<vtkExtractEdges> extractEdges = 
	vtkSmartPointer<vtkExtractEdges>::New();
	extractEdges->SetInputData(input);
	extractEdges->Update();

	// The output of extractEdges may contain duplicates. For this reason, 
	// it is looped through
	std::map<std::set<vtkIdType>,std::set<vtkIdType> > unique_edges;
	
	for(vtkIdType i = 0; i < extractEdges->GetOutput()->GetNumberOfCells(); i++){

		vtkSmartPointer<vtkLine> line = vtkLine::SafeDownCast(extractEdges->GetOutput()->GetCell(i));
		std::set<vtkIdType> edge_points_ids;
		edge_points_ids.insert(line->GetPointIds()->GetId(0));
		edge_points_ids.insert(line->GetPointIds()->GetId(1));

		if (unique_edges.find(edge_points_ids) == unique_edges.end()){
			
			// This is a new unique edge
			// We need to find the facets forming this edge

			vtkSmartPointer<vtkIdList> facet_ids_point_1 = vtkSmartPointer<vtkIdList>::New();
			vtkSmartPointer<vtkIdList> facet_ids_point_2 = vtkSmartPointer<vtkIdList>::New();

			vtkIdType p1 = *edge_points_ids.begin() ;
			vtkIdType p2 = *(++edge_points_ids.begin()) ;

			std::cout << p1 << " " << p2 << std::endl;

			input -> GetPointCells	(	p1,facet_ids_point_1 );
			input -> GetPointCells	(	p2,facet_ids_point_2 );

			std::cout << " list 1 has " << facet_ids_point_1 -> GetNumberOfIds() << " facet ids " << std::endl;
			for (unsigned int j = 0; j < facet_ids_point_1 -> GetNumberOfIds(); ++j){
				std::cout << facet_ids_point_1 -> GetId(j) << std::endl;
			}
			std::cout << " first 2 list has " << facet_ids_point_2 -> GetNumberOfIds() << " facet ids " << std::endl;
			for (unsigned int j = 0; j < facet_ids_point_2 -> GetNumberOfIds(); ++j){
				std::cout << facet_ids_point_2 -> GetId(j) << std::endl;
			}



			// Now, we find the two indices showing up in both facet_ids_point_1 and facet_ids_point_2
			facet_ids_point_1 -> IntersectWith	(	facet_ids_point_2	)	;

			if (facet_ids_point_2 -> GetNumberOfIds() != 2){
				throw(std::runtime_error("In SBGATPolyhedronGravityModel.cpp: The intersection of the facet id lists should have exactly 2 items, not " + std::to_string(facet_ids_point_2 -> GetNumberOfIds())));
			}

			std::set<vtkIdType> edge_facet_ids;
			edge_facet_ids.insert(facet_ids_point_1-> GetId(0));
			edge_facet_ids.insert(facet_ids_point_1-> GetId(1));

			unique_edges[edge_points_ids] = edge_facet_ids;

		}

	}

	std::cout << "There are " << unique_edges.size() << " unique edges\n";

	



	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end-start;

	std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";



  // this -> edge_dyads;
}


void SBGATPolyhedronGravityModel::get_facet_dyad(unsigned int facet_index,double * facet_dyad) const{
	for (unsigned int i = 0; i < 9; ++i){

		facet_dyad[i] = this  -> facet_dyads[facet_index][i];
	}

}



bool SBGATPolyhedronGravityModel::GetCellNormals(vtkPolyData* polydata){

	vtkFloatArray* normalDataFloat =
	vtkFloatArray::SafeDownCast(polydata->GetCellData()->GetArray("Normals"));
	if(normalDataFloat){
		std::cout << "normals were found\n";
	}
	return normalDataFloat;
	
}



	// double DynamicAnalyses::pgm_potential(double * point ,const double mu) const {

	// 	double potential = 0;

	// // Facet loop
	// #pragma omp parallel for reduction(+:potential) if (USE_OMP_DYNAMIC_ANALYSIS)
	// 	for (unsigned int facet_index = 0; facet_index < this -> shape_model -> get_NFacets(); ++ facet_index) {

	// 		std::vector<std::shared_ptr<Vertex > > * vertices = this -> shape_model -> get_facets() -> at(facet_index) -> get_vertices();

	// 		const double * r1 =  vertices -> at(0) -> get_coordinates() -> colptr(0);
	// 		const double * r2 =  vertices -> at(1) -> get_coordinates() -> colptr(0);
	// 		const double * r3 =  vertices -> at(2) -> get_coordinates() -> colptr(0);

	// 		double r1m[3];
	// 		double r2m[3];
	// 		double r3m[3];

	// 		r1m[0] = r1[0] - point[0];
	// 		r1m[1] = r1[1] - point[1];
	// 		r1m[2] = r1[2] - point[2];

	// 		r2m[0] = r2[0] - point[0];
	// 		r2m[1] = r2[1] - point[1];
	// 		r2m[2] = r2[2] - point[2];

	// 		r3m[0] = r3[0] - point[0];
	// 		r3m[1] = r3[1] - point[1];
	// 		r3m[2] = r3[2] - point[2];


	// 		double R1 = std::sqrt( r1m[0] * r1m[0]
	// 			+ r1m[1] * r1m[1]
	// 			+ r1m[2] * r1m[2]       );

	// 		double R2 = std::sqrt( r2m[0] * r2m[0]
	// 			+ r2m[1] * r2m[1]
	// 			+ r2m[2] * r2m[2]      );


	// 		double R3 = std::sqrt( r3m[0] * r3m[0]
	// 			+ r3m[1] * r3m[1]
	// 			+ r3m[2] * r3m[2]      );

	// 		double r2_cross_r3_0 = r2m[1] * r3m[2] - r2m[2] * r3m[1];
	// 		double r2_cross_r3_1 = r3m[0] * r2m[2] - r3m[2] * r2m[0];
	// 		double r2_cross_r3_2 = r2m[0] * r3m[1] - r2m[1] * r3m[0];


	// 		double wf = 2 * std::atan2(
	// 			r1m[0] * r2_cross_r3_0 + r1m[1] * r2_cross_r3_1 + r1m[2] * r2_cross_r3_2,

	// 			R1 * R2 * R3 + R1 * (r2m[0] * r3m[0] + r2m[1] * r3m[1]  + r2m[2] * r3m[2] )
	// 			+ R2 * (r3m[0] * r1m[0] + r3m[1] * r1m[1] + r3m[2] * r1m[2])
	// 			+ R3 * (r1m[0] * r2m[0] + r1m[1] * r2m[1] + r1m[2] * r2m[2]));


	// 		arma::mat * Fdyad = this -> shape_model -> get_facets() -> at(facet_index) -> get_facet_dyad();

	// 		double * F_col_0 = Fdyad -> colptr(0);
	// 		double * F_col_1 = Fdyad -> colptr(1);
	// 		double * F_col_2 = Fdyad -> colptr(2);

	// 		double ax = wf * (F_col_0[0] * r1m[0] + F_col_1[0] * r1m[1] +  F_col_2[0] * r1m[2]);
	// 		double ay = wf * (F_col_0[1] * r1m[0] + F_col_1[1] * r1m[1] +  F_col_2[1] * r1m[2]);
	// 		double az = wf * (F_col_0[2] * r1m[0] + F_col_1[2] * r1m[1] +  F_col_2[2] * r1m[2]);

	// 		potential += - (ax * r1m[0] + ay * r1m[1] + az * r1m[2]);

	// 	}

	// // Edge loop
	// #pragma omp parallel for reduction(+:potential) if (USE_OMP_DYNAMIC_ANALYSIS)
	// 	for (unsigned int edge_index = 0; edge_index < this -> shape_model -> get_NEdges(); ++ edge_index) {

	// 		const double * r1 =  this -> shape_model -> get_edges() -> at(edge_index) -> get_v0() -> get_coordinates() -> colptr(0);
	// 		const double * r2 =  this -> shape_model -> get_edges() -> at(edge_index) -> get_v1() -> get_coordinates() -> colptr(0);

	// 		double r1m[3];
	// 		double r2m[3];

	// 		r1m[0] = r1[0] - point[0];
	// 		r1m[1] = r1[1] - point[1];
	// 		r1m[2] = r1[2] - point[2];

	// 		r2m[0] = r2[0] - point[0];
	// 		r2m[1] = r2[1] - point[1];
	// 		r2m[2] = r2[2] - point[2];


	// 		double R1 = std::sqrt( r1m[0] * r1m[0]
	// 			+ r1m[1] * r1m[1]
	// 			+ r1m[2] * r1m[2]       );

	// 		double R2 = std::sqrt( r2m[0] * r2m[0]
	// 			+ r2m[1] * r2m[1]
	// 			+ r2m[2] * r2m[2]      );

	// 		double Re = std::sqrt( (r2m[0] - r1m[0]) * (r2m[0] - r1m[0])
	// 			+ (r2m[1] - r1m[1]) * (r2m[1] - r1m[1])
	// 			+ (r2m[2] - r1m[2]) * (r2m[2] - r1m[2])      );


	// 		double Le = std::log((R1 + R2 + Re) / (R1 + R2 - Re));


	// 		arma::mat * Edyad = this -> shape_model -> get_edges() -> at(edge_index) -> get_edge_dyad();

	// 		double * E_col_0 = Edyad -> colptr(0);
	// 		double * E_col_1 = Edyad -> colptr(1);
	// 		double * E_col_2 = Edyad -> colptr(2);

	// 		double ax = Le * (E_col_0[0] * r1m[0] + E_col_1[0] * r1m[1] +  E_col_2[0] * r1m[2]);
	// 		double ay = Le * (E_col_0[1] * r1m[0] + E_col_1[1] * r1m[1] +  E_col_2[1] * r1m[2]);
	// 		double az = Le * (E_col_0[2] * r1m[0] + E_col_1[2] * r1m[1] +  E_col_2[2] * r1m[2]);

	// 		potential += (ax * r1m[0] + ay * r1m[1] + az * r1m[2]);

	// 	}
	// 	potential *= 0.5 * mu / this -> shape_model -> get_volume();


	// 	return potential;

	// }

	// arma::vec DynamicAnalyses::pgm_acceleration(double * point , const double mu) const {

	// 	double ax = 0;
	// 	double ay = 0;
	// 	double az = 0;

	// // Facet loop
	// #pragma omp parallel for reduction(+:ax,ay,az) if (USE_OMP_DYNAMIC_ANALYSIS)
	// 	for (unsigned int facet_index = 0; facet_index < this -> shape_model -> get_NFacets(); ++ facet_index) {

	// 		std::vector<std::shared_ptr<Vertex > > * vertices = this -> shape_model -> get_facets() -> at(facet_index) -> get_vertices();

	// 		const double * r1 =  vertices -> at(0) -> get_coordinates() -> colptr(0);
	// 		const double * r2 =  vertices -> at(1) -> get_coordinates() -> colptr(0);
	// 		const double * r3 =  vertices -> at(2) -> get_coordinates() -> colptr(0);

	// 		double r1m[3];
	// 		double r2m[3];
	// 		double r3m[3];

	// 		r1m[0] = r1[0] - point[0];
	// 		r1m[1] = r1[1] - point[1];
	// 		r1m[2] = r1[2] - point[2];

	// 		r2m[0] = r2[0] - point[0];
	// 		r2m[1] = r2[1] - point[1];
	// 		r2m[2] = r2[2] - point[2];

	// 		r3m[0] = r3[0] - point[0];
	// 		r3m[1] = r3[1] - point[1];
	// 		r3m[2] = r3[2] - point[2];


	// 		double R1 = std::sqrt( r1m[0] * r1m[0]
	// 			+ r1m[1] * r1m[1]
	// 			+ r1m[2] * r1m[2]       );

	// 		double R2 = std::sqrt( r2m[0] * r2m[0]
	// 			+ r2m[1] * r2m[1]
	// 			+ r2m[2] * r2m[2]      );


	// 		double R3 = std::sqrt( r3m[0] * r3m[0]
	// 			+ r3m[1] * r3m[1]
	// 			+ r3m[2] * r3m[2]      );

	// 		double r2_cross_r3_0 = r2m[1] * r3m[2] - r2m[2] * r3m[1];
	// 		double r2_cross_r3_1 = r3m[0] * r2m[2] - r3m[2] * r2m[0];
	// 		double r2_cross_r3_2 = r2m[0] * r3m[1] - r2m[1] * r3m[0];


	// 		double wf = 2 * std::atan2(
	// 			r1m[0] * r2_cross_r3_0 + r1m[1] * r2_cross_r3_1 + r1m[2] * r2_cross_r3_2,

	// 			R1 * R2 * R3 + R1 * (r2m[0] * r3m[0] + r2m[1] * r3m[1]  + r2m[2] * r3m[2] )
	// 			+ R2 * (r3m[0] * r1m[0] + r3m[1] * r1m[1] + r3m[2] * r1m[2])
	// 			+ R3 * (r1m[0] * r2m[0] + r1m[1] * r2m[1] + r1m[2] * r2m[2]));


	// 		arma::mat * Fdyad = this -> shape_model -> get_facets() -> at(facet_index) -> get_facet_dyad();

	// 		double * F_col_0 = Fdyad -> colptr(0);
	// 		double * F_col_1 = Fdyad -> colptr(1);
	// 		double * F_col_2 = Fdyad -> colptr(2);

	// 		ax += wf * (F_col_0[0] * r1m[0] + F_col_1[0] * r1m[1] +  F_col_2[0] * r1m[2]);
	// 		ay += wf * (F_col_0[1] * r1m[0] + F_col_1[1] * r1m[1] +  F_col_2[1] * r1m[2]);
	// 		az += wf * (F_col_0[2] * r1m[0] + F_col_1[2] * r1m[1] +  F_col_2[2] * r1m[2]);

	// 	}


	// // Edge loop
	// #pragma omp parallel for reduction(-:ax,ay,az) if (USE_OMP_DYNAMIC_ANALYSIS)
	// 	for (unsigned int edge_index = 0; edge_index < this -> shape_model -> get_NEdges(); ++ edge_index) {

	// 		const double * r1 =  this -> shape_model -> get_edges() -> at(edge_index) -> get_v0() -> get_coordinates() -> colptr(0);
	// 		const double * r2 =  this -> shape_model -> get_edges() -> at(edge_index) -> get_v1() -> get_coordinates() -> colptr(0);

	// 		double r1m[3];
	// 		double r2m[3];

	// 		r1m[0] = r1[0] - point[0];
	// 		r1m[1] = r1[1] - point[1];
	// 		r1m[2] = r1[2] - point[2];

	// 		r2m[0] = r2[0] - point[0];
	// 		r2m[1] = r2[1] - point[1];
	// 		r2m[2] = r2[2] - point[2];


	// 		double R1 = std::sqrt( r1m[0] * r1m[0]
	// 			+ r1m[1] * r1m[1]
	// 			+ r1m[2] * r1m[2]       );

	// 		double R2 = std::sqrt( r2m[0] * r2m[0]
	// 			+ r2m[1] * r2m[1]
	// 			+ r2m[2] * r2m[2]      );

	// 		double Re = std::sqrt( (r2m[0] - r1m[0]) * (r2m[0] - r1m[0])
	// 			+ (r2m[1] - r1m[1]) * (r2m[1] - r1m[1])
	// 			+ (r2m[2] - r1m[2]) * (r2m[2] - r1m[2])      );


	// 		double Le = std::log((R1 + R2 + Re) / (R1 + R2 - Re));


	// 		arma::mat * Edyad = this -> shape_model -> get_edges() -> at(edge_index) -> get_edge_dyad();

	// 		double * E_col_0 = Edyad -> colptr(0);
	// 		double * E_col_1 = Edyad -> colptr(1);
	// 		double * E_col_2 = Edyad -> colptr(2);

	// 		ax -= Le * (E_col_0[0] * r1m[0] + E_col_1[0] * r1m[1] +  E_col_2[0] * r1m[2]);
	// 		ay -= Le * (E_col_0[1] * r1m[0] + E_col_1[1] * r1m[1] +  E_col_2[1] * r1m[2]);
	// 		az -= Le * (E_col_0[2] * r1m[0] + E_col_1[2] * r1m[1] +  E_col_2[2] * r1m[2]);

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










