
#include <SBGATObs.hpp>
#include <SBGATMassProperties.hpp>


vtkStandardNewMacro(SBGATObs);


//----------------------------------------------------------------------------
// Constructs with initial 0 values.
SBGATObs::SBGATObs(){

	this -> SetNumberOfOutputPorts(0);
	this -> SetNumberOfInputPorts(1);

}


SBGATObs::~SBGATObs(){

}



int SBGATObs::RequestData(
	vtkInformation* vtkNotUsed( request ),
	vtkInformationVector** inputVector,
	vtkInformationVector* vtkNotUsed( outputVector )){

	this -> bspTree_vec.clear();
	this -> polydata_vec.clear();

	vtkSmartPointer<SBGATMassProperties> mass_filter = vtkSmartPointer<SBGATMassProperties>::New();

  // Processing the primary
	vtkInformation *inInfo0 = inputVector[0]->GetInformationObject(0);
	vtkPolyData * primary = vtkPolyData::SafeDownCast(inInfo0->Get(vtkDataObject::DATA_OBJECT()));


	vtkSmartPointer<vtkModifiedBSPTree> tree = vtkSmartPointer<vtkModifiedBSPTree>::New();
	tree -> SetDataSet(primary);
	tree -> BuildLocator();

	this -> bspTree_vec.push_back(tree);
	this -> polydata_vec.push_back(primary);

	mass_filter -> SetInputData(primary);
	mass_filter -> Update();
	this -> center_of_mass_vec.push_back(mass_filter -> GetCenterOfMass());
	
  	// Processing the secondary, if any
	if(inputVector[0] -> GetNumberOfInformationObjects() == 2){

		vtkInformation *inInfo1 = inputVector[0]->GetInformationObject(1);
		vtkPolyData * secondary =  vtkPolyData::SafeDownCast(inInfo1->Get(vtkDataObject::DATA_OBJECT()));
		vtkSmartPointer<vtkModifiedBSPTree> tree = vtkSmartPointer<vtkModifiedBSPTree>::New();
		tree -> SetDataSet(secondary);
		tree -> BuildLocator();
		this -> bspTree_vec.push_back(tree);
		this -> polydata_vec.push_back(secondary);

		mass_filter -> SetInputData(secondary);
		mass_filter -> Update();
		this -> center_of_mass_vec.push_back(mass_filter -> GetCenterOfMass());

	}


	this -> number_of_bodies = polydata_vec.size();

 // The surface area of the largest facet amongst all considered shapes is found
	this -> find_min_facet_surface_area();

	return 1;
}


void SBGATObs::find_min_facet_surface_area(){

	double min_area = std::numeric_limits<double>::infinity();
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
			min_area = std::min(min_area,vtkMath::Norm(n)/2);

		}
	}

	this-> min_area = min_area;

}



int SBGATObs::FillInputPortInformation( int port, vtkInformation* info ){
	if ( port == 0 ){
		info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
		info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), true);
		return 1;

	}
	
	return 0;
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



void SBGATObs::PrintHeader(ostream& os, vtkIndent indent) {

}
void SBGATObs::PrintTrailer(ostream& os, vtkIndent indent) {

}




//----------------------------------------------------------------------------
void SBGATObs::PrintSelf(std::ostream& os, vtkIndent indent){

	vtkPolyData *input = vtkPolyData::SafeDownCast(this->GetInput(0));
	if (!input)
	{
		return;
	}

}








