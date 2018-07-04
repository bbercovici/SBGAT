
#include <SBGATObs.hpp>


void SBGATObs::find_max_facet_surface_area(){

	double max_area = -1;
	vtkIdType numCells, numIds;

	vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
	ptIds -> Allocate(VTK_CELL_SIZE);

	for (auto input : this -> polydata_vec){

		numCells = input -> GetNumberOfCells();


		for (vtkIdType cellId=0; cellId < numCells; cellId++){

			if ( input->GetCellType(cellId) != VTK_TRIANGLE){
				throw(std::runtime_error("Input data type must be VTK_TRIANGLE not " + std::to_string(input->GetCellType(cellId))));
			}

			input -> GetCellPoints(cellId,ptIds);
			numIds = ptIds -> GetNumberOfIds();
			assert(numIds == 3);

			double p0[3];
			double p1[3];
			double p2[3];

			input->GetPoint(ptIds->GetId(0), p0);
			input->GetPoint(ptIds->GetId(1), p1);
			input->GetPoint(ptIds->GetId(2), p2);

			double e0[3];
			double e1[3];

			vtkMath::Subtract(p1,p0,e0);
			vtkMath::Subtract(p2,p0,e1);

			double n[3];

			vtkMath::Cross(e0,e1,n);
			max_area = std::max(max_area,vtkMath::Norm(n)/2);

		}
	}

	this-> max_area = max_area;

}

void SBGATObs::prefind_facets_inview (
	
	std::vector<int> & facets_in_view,
	const unsigned int & body_index,
	const std::vector<arma::vec> & dir_to_check_vec,
	const std::vector<arma::mat> & BN_dcms_vec,
	const std::vector<arma::vec> & positions_vec){

	vtkIdType numCells, numIds;

	vtkPolyData * input = this -> polydata_vec[body_index];

	vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
	ptIds -> Allocate(VTK_CELL_SIZE);
	numCells = input -> GetNumberOfCells();

	for (vtkIdType cellId=0; cellId < numCells; cellId++){

		if ( input->GetCellType(cellId) != VTK_TRIANGLE){
			throw(std::runtime_error("Input data type must be VTK_TRIANGLE not " + std::to_string(input->GetCellType(cellId))));
		}

		input -> GetCellPoints(cellId,ptIds);
		numIds = ptIds -> GetNumberOfIds();
		assert(numIds == 3);

		double p0[3];
		double p1[3];
		double p2[3];

		input->GetPoint(ptIds->GetId(0), p0);
		input->GetPoint(ptIds->GetId(1), p1);
		input->GetPoint(ptIds->GetId(2), p2);

		double e0[3];
		double e1[3];

		vtkMath::Subtract(p1,p0,e0);
		vtkMath::Subtract(p2,p0,e1);

		double n[3];

		vtkMath::Cross(e0,e1,n);
		vtkMath::Normalize(n);

    // The normal vector is oriented using the body's dcm
		arma::vec n_body_frame = {n[0],n[1],n[2]};

    // Check if all the provided directions point above this facet's horizon
		bool in_view = true;

		for(auto dir : dir_to_check_vec){
			if (arma::dot(n_body_frame,BN_dcms_vec[body_index] * dir) < 0){
				in_view = false;
				break;
			}
		}

		if (in_view){
			facets_in_view.push_back(cellId);
		}

	}
}


bool SBGATObs::check_line_for_intersect(const int & origin_body_index,
	const arma::vec & start_point_origin_body,
	const arma::vec & end_point_inertial,
	const std::vector<arma::mat> & BN_dcms_vec,
	const std::vector<arma::vec> & positions_vec,
	const double & tol) const{

	for (int considered_body_index = 0; considered_body_index < this -> number_of_bodies; ++considered_body_index){

		// The KD tree of the considered body is aligned with the axes of the B frame 
		// the center of mass of the considered body is not necessarily at (0,0,0) so this offset must be accounted for

		arma::vec start_point_inertial = BN_dcms_vec[origin_body_index].t() * (start_point_origin_body - this -> center_of_mass_vec[origin_body_index]) + positions_vec[origin_body_index];
		arma::vec start_point_considered = BN_dcms_vec[considered_body_index] *  (start_point_inertial - positions_vec[considered_body_index]) + this -> center_of_mass_vec[considered_body_index];
		arma::vec end_point_considered = BN_dcms_vec[considered_body_index] *  (end_point_inertial - positions_vec[considered_body_index]) + this -> center_of_mass_vec[considered_body_index];

		vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
		vtkSmartPointer<vtkPoints> verts = vtkSmartPointer<vtkPoints>::New();

		 // Checking if the considered line is intercepted by any of the considered bodies
		this -> bspTree_vec[considered_body_index] -> IntersectWithLine(start_point_considered.colptr(0), 
			end_point_considered.colptr(0), tol, verts, cellIds);

		 if (verts -> GetNumberOfPoints() > 0){
		 	return true;
		 }

	}

	return false;

}









