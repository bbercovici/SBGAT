
#include <SBGATObs.hpp>

void SBGATObs::prefind_facets_inview(
	const unsigned int & body_index,
	std::vector<int> & facets_in_view,
	const std::vector<arma::vec> & dir_to_check_vec,
	double & max_area,
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
		max_area = std::max(max_area,vtkMath::Norm(n)/2);
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
