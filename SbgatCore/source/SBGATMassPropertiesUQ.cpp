#include <SBGATMassPropertiesUQ.hpp>
#include <RigidBodyKinematics.hpp>
#include <vtkCleanPolyData.h>
#include <vtkOBJReader.h>
#include <SBGATPolyhedronGravityModel.hpp>
#include <SBGATObjWriter.hpp>

#pragma omp declare reduction( + : arma::mat : omp_out += omp_in ) \
initializer( omp_priv = arma::zeros<arma::mat>(omp_orig.n_rows,omp_orig.n_cols))


#pragma omp declare reduction( + : arma::rowvec : omp_out += omp_in ) \
initializer( omp_priv = arma::zeros<arma::rowvec>(omp_orig.n_cols))


arma::rowvec::fixed<9> SBGATMassPropertiesUQ::PartialDeltaVfPartialTf(const int & f) const{

	double r0[3];
	double r1[3];
	double r2[3];

	arma::rowvec::fixed<9> partial;
	this -> model -> GetVerticesInFacet(f,r0,r1,r2);
	arma::vec::fixed<3> C0 = {r0[0],r0[1],r0[2]};
	arma::vec::fixed<3> C1 = {r1[0],r1[1],r1[2]};
	arma::vec::fixed<3> C2 = {r2[0],r2[1],r2[2]};

	partial.subvec(0,2) = arma::cross(C1,C2).t();
	partial.subvec(3,5) = - C0.t() * RBK::tilde(C2);
	partial.subvec(6,8) = C0.t() * RBK::tilde(C1);

	return 1./6 * partial;

	
}

arma::mat::fixed<3,9> SBGATMassPropertiesUQ::PartialDeltaCMfPartialTf(){

	arma::mat::fixed<3,9> triple_I;
	triple_I.cols(0,2) = arma::eye<arma::mat>(3,3);
	triple_I.cols(3,5) = arma::eye<arma::mat>(3,3);
	triple_I.cols(6,8) = arma::eye<arma::mat>(3,3);

	return  1./4 * triple_I ;

}

void SBGATMassPropertiesUQ::TestPartials(std::string input,double tol,bool shape_in_meters){

	std::cout << "\tRunning SBGATMassPropertiesUQ::TestPartials on " << input << std::endl;

	SBGATMassPropertiesUQ::TestPartialDeltaVfPartialTf(input,tol,shape_in_meters);
	SBGATMassPropertiesUQ::TestPartialDeltaIfPartialTf(input,tol,shape_in_meters);
	
	SBGATMassPropertiesUQ::TestGetPartialVolumePartialC(input,tol,shape_in_meters);
	SBGATMassPropertiesUQ::TestGetPartialComPartialC(input,tol,shape_in_meters);

	SBGATMassPropertiesUQ::TestPartialEqDeltaIfErPartialTf(input,tol,shape_in_meters);
	SBGATMassPropertiesUQ::TestGetPartialIPartialC(input,tol,shape_in_meters);

	SBGATMassPropertiesUQ::TestGetPartialAllInertiaPartialC(input,tol,shape_in_meters);

	SBGATMassPropertiesUQ::TestGetPartialSigmaPartialC(input,tol,shape_in_meters);

}



void SBGATMassPropertiesUQ::TestPartialDeltaVfPartialTf(std::string input,double tol,bool shape_in_meters){

	std::cout << "\t In TestPartialDeltaVfPartialTf ... ";
	int successes = 0;
	arma::arma_rng::set_seed(0);
	int N = 1000;


	// Reading
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader -> SetFileName(input.c_str());
	reader -> Update(); 

	// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cleaner =
	vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner -> SetInputConnection (reader -> GetOutputPort());
	cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
	cleaner -> Update();

	#pragma omp parallel for reduction(+:successes)
	
	for (int i = 0; i < N ; ++i){

		vtkSmartPointer<vtkPolyData> input_shape = vtkSmartPointer<vtkPolyData>::New();

		input_shape -> DeepCopy(cleaner -> GetOutput());

	// Creating the PGM dyads
		vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
		mass_prop -> SetInputData(input_shape);
		mass_prop -> SetDensity(1970); 
		if (shape_in_meters) 
			mass_prop -> SetScaleMeters();
		else
			mass_prop -> SetScaleKiloMeters();

		mass_prop -> Update();

		SBGATMassPropertiesUQ shape_uq;
		shape_uq.SetModel(mass_prop);
		shape_uq.PrecomputeMassPropertiesPartials();


		int N_facets = vtkPolyData::SafeDownCast(mass_prop -> GetInput()) -> GetNumberOfCells();
		
		arma::ivec f_vec = arma::randi<arma::ivec>(1,arma::distr_param(0,N_facets - 1));
		int f = f_vec(0);



	// Nominal Volume
		double volume = mass_prop -> GetDeltaV(f);

	// Deviation
		arma::vec::fixed<9> delta_Tf = 1e-3 * arma::randn<arma::vec>(9) / mass_prop -> GetScaleFactor();

	// Linear dUf
		double dV_lin = arma::dot(shape_uq.PartialDeltaVfPartialTf(f) , mass_prop -> GetScaleFactor()* delta_Tf);

	// Apply Tf deviation
		shape_uq. ApplyTfDeviation(delta_Tf,f);

	// Perturbed Uf
		double volume_p = mass_prop -> GetDeltaV(f);
		
	// Non-linear dUf
		double dV = volume_p - volume;


		if(std::abs(dV - dV_lin) / std::abs(dV_lin) < tol){
			++successes;
		}



	}

	std::cout << "\t Passed TestPartialDeltaVfPartialTf with " << double(successes)/N * 100 << " \% of successes. \n";


}


void SBGATMassPropertiesUQ::TestGetPartialVolumePartialC(std::string input,double tol,bool shape_in_meters){

	std::cout << "\t In TestGetPartialVolumePartialC ... ";
	int successes = 0;
	arma::arma_rng::set_seed(0);
	int N = 1000;
	#pragma omp parallel for reduction(+:successes)
	
	for (int i = 0; i < N ; ++i){

		// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(input.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner =
		vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();


	// Creating the PGM dyads
		vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
		mass_prop -> SetInputConnection(cleaner -> GetOutputPort());
		mass_prop -> SetDensity(1970); 
		if (shape_in_meters) 
			mass_prop -> SetScaleMeters();
		else
			mass_prop -> SetScaleKiloMeters();

		mass_prop -> Update();

		SBGATMassPropertiesUQ shape_uq;
		shape_uq.SetModel(mass_prop);
		shape_uq.PrecomputeMassPropertiesPartials();


	// Nominal Volume
		double volume = mass_prop -> GetVolume();

	// Deviation
		arma::vec delta_C = 1e-3 * arma::randn<arma::vec>(3 * mass_prop -> GetN_vertices()  ) / mass_prop -> GetScaleFactor();
		
	// Linear dUf
		double dV_lin = mass_prop -> GetScaleFactor() * arma::dot(shape_uq.GetPartialVolumePartialC() , delta_C);

	// Apply Tf deviation
		shape_uq. ApplyDeviation(delta_C);

	// Perturbed Uf
		double volume_p = mass_prop -> GetVolume();
		
	// Non-linear dUf
		double dV = volume_p - volume;

		if(std::abs(dV - dV_lin) / std::abs(dV_lin) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestGetPartialVolumePartialC with " << double(successes)/N * 100 << " \% of successes. \n";


}



void SBGATMassPropertiesUQ::TestGetPartialComPartialC(std::string input,double tol,bool shape_in_meters){


	std::cout << "\t In TestGetPartialComPartialC ... ";
	int successes = 0;
	arma::arma_rng::set_seed(0);
	int N = 1000;

	#pragma omp parallel for reduction(+:successes)
	
	for (int i = 0; i < N ; ++i){

		// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(input.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner =
		vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();


	// Creating the PGM dyads
		vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
		mass_prop -> SetInputConnection(cleaner -> GetOutputPort());
		mass_prop -> SetDensity(1970); 
		if (shape_in_meters) 
			mass_prop -> SetScaleMeters();
		else
			mass_prop -> SetScaleKiloMeters();

		mass_prop -> Update();

		SBGATMassPropertiesUQ shape_uq;
		shape_uq.SetModel(mass_prop);
		shape_uq.PrecomputeMassPropertiesPartials();

	// Nominal Volume
		arma::vec::fixed<3> com = mass_prop -> GetCenterOfMass();

	// Deviation
		arma::vec delta_C = 1e-3 * arma::randn<arma::vec>(3 * mass_prop -> GetN_vertices()  ) / mass_prop -> GetScaleFactor();

	// Linear dUf
		arma::vec::fixed<3> dcom_lin = shape_uq.GetPartialComPartialC() * delta_C * mass_prop -> GetScaleFactor() ;

	// Apply Tf deviation
		shape_uq. ApplyDeviation(delta_C);

	// Perturbed Uf
		arma::vec::fixed<3> com_p = mass_prop -> GetCenterOfMass();
		
	// Non-linear dUf
		arma::vec::fixed<3> dcom = com_p - com;

		if(arma::norm(dcom - dcom_lin) / arma::norm(dcom_lin) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestGetPartialComPartialC with " << double(successes)/N * 100 << " \% of successes. \n";

}



void SBGATMassPropertiesUQ::ApplyTfDeviation(arma::vec::fixed<9> delta_Tf,const int & f){
	
	int v0_index,v1_index,v2_index;
	this -> model -> GetIndicesVerticesInFacet(f,v0_index,v1_index,v2_index);

	double r0[3],r1[3],r2[3];
	vtkPolyData * polydata = vtkPolyData::SafeDownCast(this -> model -> GetInput());

	polydata -> GetPoint(v0_index,r0);
	polydata -> GetPoint(v1_index,r1);
	polydata -> GetPoint(v2_index,r2);


	r0[0] += delta_Tf(0);
	r0[1] += delta_Tf(1);
	r0[2] += delta_Tf(2);

	r1[0] += delta_Tf(3);
	r1[1] += delta_Tf(4);
	r1[2] += delta_Tf(5);

	r2[0] += delta_Tf(6);
	r2[1] += delta_Tf(7);
	r2[2] += delta_Tf(8);

	polydata -> GetPoints() -> SetPoint(v0_index,r0);
	polydata -> GetPoints() -> SetPoint(v1_index,r1);
	polydata -> GetPoints() -> SetPoint(v2_index,r2);

	polydata -> GetPoints() -> Modified();

	polydata -> Modified();
	this -> model -> Modified();

	this -> model -> Update();

}




void SBGATMassPropertiesUQ::TestGetPartialIPartialC(std::string input,double tol,bool shape_in_meters) {


	std::cout << "\t In TestGetPartialIPartialC ... ";
	int successes = 0;
	arma::arma_rng::set_seed(0);
	int N = 1000;
	#pragma omp parallel for reduction(+:successes)
	
	for (int i = 0; i < N ; ++i){

		// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(input.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner =
		vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();


	// Creating the PGM dyads
		vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
		mass_prop -> SetInputConnection(cleaner -> GetOutputPort());
		mass_prop -> SetDensity(1970);
		if (shape_in_meters) 
			mass_prop -> SetScaleMeters();
		else
			mass_prop -> SetScaleKiloMeters();

		mass_prop -> Update();

		SBGATMassPropertiesUQ shape_uq;
		shape_uq.SetModel(mass_prop);
		shape_uq.PrecomputeMassPropertiesPartials();

	// Nominal 
		arma::mat deltaI_mat = mass_prop -> GetUnitDensityInertiaTensor() ;
		arma::vec::fixed<6> deltaI = {deltaI_mat(0,0),deltaI_mat(1,1),deltaI_mat(2,2),deltaI_mat(0,1),deltaI_mat(0,2),deltaI_mat(1,2)};

	// Deviation
		arma::vec delta_C = 1e-3 * arma::randn<arma::vec>(3 * mass_prop -> GetN_vertices()  ) / mass_prop -> GetScaleFactor();

	// Linear 
		arma::vec::fixed<6> ddeltaI_lin = mass_prop -> GetScaleFactor() * shape_uq.GetPartialIPartialC() * delta_C;

	// Apply deviation
		shape_uq. ApplyDeviation(delta_C);

	// Nominal 
		arma::mat deltaI_mat_p = mass_prop -> GetUnitDensityInertiaTensor() ;
		arma::vec::fixed<6> deltaI_p = {deltaI_mat_p(0,0),deltaI_mat_p(1,1),deltaI_mat_p(2,2),deltaI_mat_p(0,1),deltaI_mat_p(0,2),deltaI_mat_p(1,2)};

	// Non-linear 
		arma::vec::fixed<6> ddeltaI = deltaI_p - deltaI;

		if(arma::norm(arma::vectorise(ddeltaI - ddeltaI_lin)) / arma::norm(arma::vectorise(ddeltaI_lin)) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestGetPartialIPartialC with " << double(successes)/N * 100 << " \% of successes. \n";


}


void SBGATMassPropertiesUQ::TestGetPartialSigmaPartialC(std::string input,double tol,bool shape_in_meters) {

	std::cout << "\t In TestGetPartialSigmaPartialC ... ";
	int successes = 0;
	arma::arma_rng::set_seed(0);
	int N = 1000;
	#pragma omp parallel for reduction(+:successes)
	
	for (int i = 0; i < N ; ++i){

		// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(input.c_str());
		reader -> Update(); 

		// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner =
		vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();


	// Creating the PGM dyads
		vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
		mass_prop -> SetInputConnection(cleaner -> GetOutputPort());
		mass_prop -> SetDensity(1970);
		if (shape_in_meters) 
			mass_prop -> SetScaleMeters();
		else
			mass_prop -> SetScaleKiloMeters();

		mass_prop -> Update();

		SBGATMassPropertiesUQ shape_uq;
		shape_uq.SetModel(mass_prop);
		shape_uq.PrecomputeMassPropertiesPartials();

	// Nominal 
		arma::vec::fixed<3> sigma = RBK::dcm_to_mrp(mass_prop -> GetPrincipalAxes()) ;

	// Deviation
		arma::vec delta_C = 1e-3 * arma::randn<arma::vec>(3 * mass_prop -> GetN_vertices()  ) / mass_prop -> GetScaleFactor();

	// Linear 
		arma::vec::fixed<3> dsigma_lin = mass_prop -> GetScaleFactor() * shape_uq.GetPartialSigmaPartialC() * delta_C;

	// Apply deviation
		shape_uq. ApplyDeviation(delta_C);

	// Perturbed 
		arma::vec::fixed<3> sigma_p = RBK::dcm_to_mrp(mass_prop -> GetPrincipalAxes()) ;
		
	// Non-linear 
		arma::vec::fixed<3> dsigma = RBK::dcm_to_mrp(RBK::mrp_to_dcm(sigma_p) * RBK::mrp_to_dcm(sigma).t());

		if(arma::norm(RBK::dcm_to_mrp(RBK::mrp_to_dcm(dsigma) * RBK::mrp_to_dcm(dsigma_lin).t())) / arma::norm(dsigma_lin) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestGetPartialSigmaPartialC with " << double(successes)/N * 100 << " \% of successes. \n";


}


void SBGATMassPropertiesUQ::TestGetPartialAllInertiaPartialC(std::string input,double tol,bool shape_in_meters) {


	std::cout << "\t In TestGetPartialAllInertiaPartialC ... ";
	int successes = 0;
	arma::arma_rng::set_seed(0);
	int N = 1000;

	#pragma omp parallel for reduction(+:successes) 
	for (int i = 0; i < N ; ++i){

		// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(input.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner =
		vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();


	// Creating the PGM dyads
		vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
		mass_prop -> SetInputConnection(cleaner -> GetOutputPort());
		mass_prop -> SetDensity(1970);
		if (shape_in_meters) 
			mass_prop -> SetScaleMeters();
		else
			mass_prop -> SetScaleKiloMeters();

		mass_prop -> Update();

		SBGATMassPropertiesUQ shape_uq;
		shape_uq.SetModel(mass_prop);
		shape_uq.PrecomputeMassPropertiesPartials();

	// Nominal 
		double V = mass_prop -> GetVolume();
		arma::vec::fixed<3> com = mass_prop -> GetCenterOfMass();

		arma::mat deltaI_mat = mass_prop -> GetUnitDensityInertiaTensor() ;
		arma::vec::fixed<6> deltaI = {deltaI_mat(0,0),deltaI_mat(1,1),deltaI_mat(2,2),deltaI_mat(0,1),deltaI_mat(0,2),deltaI_mat(1,2)};

	// Deviation
		arma::vec delta_C = 1e-3 * arma::randn<arma::vec>(3 * mass_prop -> GetN_vertices()  ) / mass_prop -> GetScaleFactor();

	// Linear 

		double dV_lin = arma::dot(shape_uq.GetPartialVolumePartialC(),mass_prop -> GetScaleFactor() * delta_C);
		arma::vec::fixed<3> dcom_lin = shape_uq.GetPartialComPartialC() * mass_prop -> GetScaleFactor() * delta_C;
		arma::vec::fixed<6> ddeltaI_lin = mass_prop -> GetScaleFactor() * shape_uq.GetPartialIPartialC() * delta_C;

	// Apply deviation
		shape_uq. ApplyDeviation(delta_C);

		double V_p = mass_prop -> GetVolume();
		arma::vec::fixed<3> com_p = mass_prop -> GetCenterOfMass();


		double dV =V_p - V;


		arma::vec::fixed<3> dcom = com_p - com;


		arma::mat deltaI_mat_p = mass_prop -> GetUnitDensityInertiaTensor() ;
		arma::vec::fixed<6> deltaI_p = {deltaI_mat_p(0,0),deltaI_mat_p(1,1),deltaI_mat_p(2,2),deltaI_mat_p(0,1),deltaI_mat_p(0,2),deltaI_mat_p(1,2)};



	// Non-linear 
		arma::vec::fixed<6> ddeltaI = deltaI_p - deltaI;

		if(arma::abs(arma::vectorise(ddeltaI - ddeltaI_lin)).max() / arma::norm(arma::vectorise(ddeltaI_lin)) < tol){
			


			if (std::abs(dV - dV_lin)/std::abs(dV_lin) < tol){

				if(arma::abs(dcom - dcom_lin).max() / arma::norm(dcom_lin) < tol){
					++successes;

				}
			}
		}

	}

	std::cout << "\t Passed TestGetPartialAllInertiaPartialC with " << double(successes)/N * 100 << " \% of successes \n";


}


void SBGATMassPropertiesUQ::TestPartialDeltaIfPartialTf(std::string input,double tol,bool shape_in_meters) {


	std::cout << "\t In TestPartialDeltaIfPartialTf ... ";
	int successes = 0;
	arma::arma_rng::set_seed(0);
	int N = 1000;
	#pragma omp parallel for reduction(+:successes)
	
	for (int i = 0; i < N ; ++i){

		// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(input.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner =
		vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();


	// Creating the PGM dyads
		vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
		mass_prop -> SetInputConnection(cleaner -> GetOutputPort());
		mass_prop -> SetDensity(1970); 
		if (shape_in_meters) 
			mass_prop -> SetScaleMeters();
		else
			mass_prop -> SetScaleKiloMeters();

		mass_prop -> Update();

		SBGATMassPropertiesUQ shape_uq;
		shape_uq.SetModel(mass_prop);
		shape_uq.PrecomputeMassPropertiesPartials();

		int N_facets = vtkPolyData::SafeDownCast(mass_prop -> GetInput()) -> GetNumberOfCells();
		arma::ivec f_vec = arma::randi<arma::ivec>(1,arma::distr_param(0,N_facets - 1));
		int f = f_vec(0);
	// Nominal 
		arma::mat deltaIf_mat = mass_prop -> GetDeltaIOverDeltaV(f) * mass_prop -> GetDeltaV(f) ;
		arma::vec::fixed<6> deltaIf = {deltaIf_mat(0,0),deltaIf_mat(1,1),deltaIf_mat(2,2),deltaIf_mat(0,1),deltaIf_mat(0,2),deltaIf_mat(1,2)};

	// Deviation
		arma::vec::fixed<9> delta_Tf = 1e-3 * arma::randn<arma::vec>(9) / mass_prop -> GetScaleFactor();

	// Linear 
		arma::vec::fixed<6> ddeltaIf_lin = mass_prop -> GetScaleFactor() * shape_uq.PartialDeltaIfPartialTf(f) * delta_Tf;

	// Apply Tf deviation
		shape_uq. ApplyTfDeviation(delta_Tf,f);

	// Nominal 
		arma::mat deltaIf_mat_p = mass_prop -> GetDeltaIOverDeltaV(f) * mass_prop -> GetDeltaV(f) ;
		arma::vec::fixed<6> deltaIf_p = {deltaIf_mat_p(0,0),deltaIf_mat_p(1,1),deltaIf_mat_p(2,2),deltaIf_mat_p(0,1),deltaIf_mat_p(0,2),deltaIf_mat_p(1,2)};

	// Non-linear 
		arma::vec::fixed<6> ddeltaIf = deltaIf_p - deltaIf;

		if(arma::norm(arma::vectorise(ddeltaIf - ddeltaIf_lin)) / arma::norm(arma::vectorise(ddeltaIf_lin)) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestPartialDeltaIfPartialTf with " << double(successes)/N * 100 << " \% of successes. \n";

}



arma::mat::fixed<3,9> SBGATMassPropertiesUQ::PartialDeltaComPartialTf(const int & f) const{

	SBGATMassProperties * mass_model = nullptr;

	mass_model = SBGATMassProperties::SafeDownCast(this -> model);
	if (!mass_model){
		mass_model = SBGATPolyhedronGravityModel::SafeDownCast(this -> model);
	}


	return (1./mass_model -> GetVolume() * ((mass_model -> GetDeltaCM(f) 
		- mass_model -> GetCenterOfMass()) * this -> PartialDeltaVfPartialTf(f)  
	+ mass_model -> GetDeltaV(f) * SBGATMassPropertiesUQ::PartialDeltaCMfPartialTf()));

}


arma::rowvec::fixed<9> SBGATMassPropertiesUQ::PartialEqDeltaIfErPartialTf(const int & f,const int & q, const int & r,const arma::vec::fixed<9> & Tf) const{

	SBGATMassProperties * mass_model = nullptr;

	mass_model = SBGATMassProperties::SafeDownCast(this -> model);
	if (!mass_model){
		mass_model = SBGATPolyhedronGravityModel::SafeDownCast(this -> model);
	}

	arma::vec::fixed<3> e_q = arma::zeros<arma::vec >(3);
	arma::vec::fixed<3> e_r = arma::zeros<arma::vec >(3);
	
	e_q(q) = 1;
	e_r(r) = 1;

	arma::rowvec::fixed<9> partial = (arma::dot(e_q,
		mass_model -> GetDeltaIOverDeltaV(f) * e_r) * this -> PartialDeltaVfPartialTf(f) 
	+ mass_model -> GetDeltaV(f) * this -> PartialEqDeltaIOverDeltaVfErPartialTf(e_q,e_r,Tf));

	return partial;



}



arma::rowvec::fixed<9> SBGATMassPropertiesUQ::PartialEqDeltaIOverDeltaVfErPartialTf(const arma::vec::fixed<3> & e_q,const arma::vec::fixed<3> & e_r,
	const arma::vec::fixed<9> & Tf) const {

	arma::rowvec::fixed<9> partial;

	arma::mat::fixed<3,9> A,A_0,A_1,A_2;

	A.cols(0,2) = arma::eye<arma::mat>(3,3);
	A.cols(3,5) = arma::eye<arma::mat>(3,3);
	A.cols(6,8) = arma::eye<arma::mat>(3,3);

	A_0.cols(0,2) = arma::eye<arma::mat>(3,3);
	A_0.cols(3,5) = arma::zeros<arma::mat>(3,3);
	A_0.cols(6,8) = arma::zeros<arma::mat>(3,3);

	A_1.cols(0,2) = arma::zeros<arma::mat>(3,3);
	A_1.cols(3,5) = arma::eye<arma::mat>(3,3);
	A_1.cols(6,8) = arma::zeros<arma::mat>(3,3);

	A_2.cols(0,2) = arma::zeros<arma::mat>(3,3);
	A_2.cols(3,5) = arma::zeros<arma::mat>(3,3);
	A_2.cols(6,8) = arma::eye<arma::mat>(3,3);


	partial = Tf.t() * (A.t() * (RBK::tilde(e_q) * RBK::tilde(e_r) + RBK::tilde(e_r) * RBK::tilde(e_q)) * A
		+ A_0.t() * (RBK::tilde(e_q) * RBK::tilde(e_r) + RBK::tilde(e_r) * RBK::tilde(e_q)) * A_0
		+ A_1.t() * (RBK::tilde(e_q) * RBK::tilde(e_r) + RBK::tilde(e_r) * RBK::tilde(e_q)) * A_1
		+ A_2.t() * (RBK::tilde(e_q) * RBK::tilde(e_r) + RBK::tilde(e_r) * RBK::tilde(e_q)) * A_2);



	return -1./20 * partial;


}


arma::mat::fixed<6,9> SBGATMassPropertiesUQ::PartialDeltaIfPartialTf(const int & f) const{

	arma::mat::fixed<6,9> partial;

	double r0[3];
	double r1[3];
	double r2[3];

	this -> model -> GetVerticesInFacet(f,r0,r1,r2);

	arma::vec::fixed<9> Tf = arma::vec({r0[0],r0[1],r0[2],r1[0],r1[1],r1[2],r2[0],r2[1],r2[2]});

	partial.row(0) = this -> PartialEqDeltaIfErPartialTf(f,0,0,Tf);
	partial.row(1) = this -> PartialEqDeltaIfErPartialTf(f,1,1,Tf);
	partial.row(2) = this -> PartialEqDeltaIfErPartialTf(f,2,2,Tf);
	partial.row(3) = this -> PartialEqDeltaIfErPartialTf(f,0,1,Tf);
	partial.row(4) = this -> PartialEqDeltaIfErPartialTf(f,0,2,Tf);
	partial.row(5) = this -> PartialEqDeltaIfErPartialTf(f,1,2,Tf);

	return partial ;
}


void SBGATMassPropertiesUQ::TestPartialEqDeltaIfErPartialTf(std::string input,double tol,bool shape_in_meters){

	std::cout << "\t In TestPartialEqDeltaIfErPartialTf ... ";
	int successes = 0;
	arma::arma_rng::set_seed(0);
	int N = 1000;


	// Reading
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader -> SetFileName(input.c_str());
	reader -> Update(); 

	// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cleaner =
	vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner -> SetInputConnection (reader -> GetOutputPort());
	cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
	cleaner -> Update();



	#pragma omp parallel for reduction(+:successes)
	
	for (int i = 0; i < N ; ++i){

		vtkSmartPointer<vtkPolyData> input_shape = vtkSmartPointer<vtkPolyData>::New();

		input_shape -> DeepCopy(cleaner -> GetOutput());

	// Creating the PGM dyads
		vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
		mass_prop -> SetInputData(input_shape);
		mass_prop -> SetDensity(1970); 
		if (shape_in_meters) 
			mass_prop -> SetScaleMeters();
		else
			mass_prop -> SetScaleKiloMeters();

		mass_prop -> Update();

		SBGATMassPropertiesUQ shape_uq;
		shape_uq.SetModel(mass_prop);
		shape_uq.PrecomputeMassPropertiesPartials();

		int N_facets = vtkPolyData::SafeDownCast(mass_prop -> GetInput()) -> GetNumberOfCells();
		
		arma::ivec f_vec = arma::randi<arma::ivec>(1,arma::distr_param(0,N_facets - 1));
		int f = f_vec(0);

		arma::ivec indices_vec = arma::randi<arma::ivec>(2,arma::distr_param(0,2));
		
		arma::vec::fixed<3> e_q = arma::zeros<arma::vec>(3);
		arma::vec::fixed<3> e_r = arma::zeros<arma::vec>(3);

		e_q(indices_vec(0)) = 1;
		e_r(indices_vec(1)) = 1;


	// Nominal 
		double Iqf = arma::dot(e_q,mass_prop -> GetDeltaIOverDeltaV(f) * mass_prop -> GetDeltaV(f) * e_r);

	// Deviation
		arma::vec::fixed<9> delta_Tf = 1e-3 * arma::randn<arma::vec>(9) / mass_prop -> GetScaleFactor();

	// Linear dUf

		double r0[3];
		double r1[3];
		double r2[3];

		mass_prop -> GetVerticesInFacet(f,r0,r1,r2);

		arma::vec::fixed<9> Tf = arma::vec({r0[0],r0[1],r0[2],r1[0],r1[1],r1[2],r2[0],r2[1],r2[2]});

		double dIqf_lin = arma::dot(shape_uq.PartialEqDeltaIfErPartialTf(f,indices_vec(0),indices_vec(1), Tf) , mass_prop -> GetScaleFactor()* delta_Tf);

	// Apply Tf deviation
		shape_uq. ApplyTfDeviation(delta_Tf,f);

	// Perturbed 
		double Iqf_p = arma::dot(e_q,mass_prop -> GetDeltaIOverDeltaV(f) * mass_prop -> GetDeltaV(f) * e_r);
		
		
	// Non-linear dUf
		double dIqf = Iqf_p - Iqf;


		if(std::abs(dIqf - dIqf_lin) / std::abs(dIqf_lin) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestPartialEqDeltaIfErPartialTf with " << double(successes)/N * 100 << " \% of successes. \n";


}






void SBGATMassPropertiesUQ::ApplyDeviation(const arma::vec & delta_C){

	vtkPolyData * polydata = vtkPolyData::SafeDownCast(this -> model -> GetInput());
	int N_C = polydata -> GetNumberOfPoints();
	assert(3 * N_C == delta_C.n_rows);

	double r[3];

	for (int i = 0; i < N_C; ++i){

		polydata -> GetPoint(i,r);

		r[0] += delta_C(3 * i);
		r[1] += delta_C(3 * i + 1);
		r[2] += delta_C(3 * i + 2);

		polydata -> GetPoints() -> SetPoint(i,r);
	}

	polydata -> GetPoints() -> Modified();
	polydata -> Modified();

	this -> model -> Modified();
	SBGATMassProperties::SafeDownCast(this -> model) -> Update();

}


arma::mat::fixed<3,6> SBGATMassPropertiesUQ::PartialSigmaPartialI() const{
	SBGATMassProperties * mass_model = nullptr;

	mass_model = SBGATMassProperties::SafeDownCast(this -> model);
	
	if (!mass_model){
		mass_model = SBGATPolyhedronGravityModel::SafeDownCast(this -> model);
	}

	
	arma::mat::fixed<3,3> BP = mass_model -> GetPrincipalAxes().t();

	arma::vec::fixed<3> moments = mass_model -> GetUnitDensityInertiaMoments();


	arma::mat::fixed<3,6> W1,W2,W3;
	W1 = W2 = W3 = arma::zeros<arma::mat>(3,6);

	W1(0,0) = 1;
	W1(1,3) = 1;
	W1(2,4) = 1;

	W2(0,3) = 1;
	W2(1,1) = 1;
	W2(2,5) = 1;

	W3(0,4) = 1;
	W3(1,5) = 1;
	W3(2,2) = 1;


	arma::mat::fixed<9,3> H = arma::zeros<arma::mat>(9,3);
	arma::mat::fixed<9,6> V = arma::zeros<arma::mat>(9,6);
	
	arma::mat::fixed<3,3> D = arma::diagmat(moments);

	int index = 0;
	for (int q = 0; q < 3; ++q){
		
		arma::vec::fixed<3> e_q = arma::zeros<arma::vec>(3);
		e_q(q) = 1;

		arma::vec::fixed<3> f_q = BP * e_q;
		
		for (int r = 0; r < 3; ++r){
			
			arma::vec::fixed<3> e_r= arma::zeros<arma::vec>(3);
			e_r(r) = 1;

			arma::vec::fixed<3> f_r = BP * e_r;

			double delta_rq;
			
			if (r==q){
				delta_rq = 1;
			}
			else{
				delta_rq = 0;
			}

			arma::mat::fixed<3,6> Fq = arma::zeros<arma::mat>(3,6);
			
			Fq.row(0) = f_q.t() * W1;
			Fq.row(1) = f_q.t() * W2;
			Fq.row(2) = f_q.t() * W3;



			arma::rowvec::fixed<6> J_rq = f_r.t() * Fq;

			H.row(index) = e_r.t() * ( D * RBK::tilde(e_q) - RBK::tilde(D * e_q));
			V.row(index) = (J_rq - delta_rq * e_r.t() * this -> precomputed_partialUnitDensityMomentsPartialI)/4;
			++index;

		}
	}

	return arma::inv(H.t() * H) * H.t() * V;


}




arma::mat::fixed<3,6>  SBGATMassPropertiesUQ::PartialUnitDensityMomentsPartialI() const{
	

	SBGATMassProperties * mass_model = nullptr;

	mass_model = SBGATMassProperties::SafeDownCast(this -> model);
	if (!mass_model){
		mass_model = SBGATPolyhedronGravityModel::SafeDownCast(this -> model);
	}

	arma::mat::fixed<3,3> eigvec = mass_model -> GetPrincipalAxes().t();


	arma::mat::fixed<3,6> W1,W2,W3;
	W1 = W2 = W3 = arma::zeros<arma::mat>(3,6);

	W1(0,0) = 1;
	W1(1,3) = 1;
	W1(2,4) = 1;

	W2(0,3) = 1;
	W2(1,1) = 1;
	W2(2,5) = 1;

	W3(0,4) = 1;
	W3(1,5) = 1;
	W3(2,2) = 1;

	arma::mat::fixed<3,6> U1,U2,U3;

	U1.row(0) = eigvec.col(0).t() * W1;
	U1.row(1) = eigvec.col(0).t() * W2;
	U1.row(2) = eigvec.col(0).t() * W3;

	U2.row(0) = eigvec.col(1).t() * W1;
	U2.row(1) = eigvec.col(1).t() * W2;
	U2.row(2) = eigvec.col(1).t() * W3;

	U3.row(0) = eigvec.col(2).t() * W1;
	U3.row(1) = eigvec.col(2).t() * W2;
	U3.row(2) = eigvec.col(2).t() * W3;

	arma::mat::fixed<3,6> dMdI = arma::zeros<arma::mat>(3,6);
	dMdI.row(0) = eigvec.col(0).t() * U1;
	dMdI.row(1) = eigvec.col(1).t() * U2;
	dMdI.row(2) = eigvec.col(2).t() * U3;
	return dMdI;



}


 /**
  Evaluates the partial of the volume, center of mass and mrp orienting the principal axes
  relative to the vertices coordinates and stores the computed partials in designated containers
  */
void SBGATMassPropertiesUQ::PrecomputeMassPropertiesPartials(){


	arma::rowvec dVdC = arma::zeros<arma::rowvec>(3 * this -> model -> GetN_vertices());
	arma::mat dCOMdC = arma::zeros<arma::mat>(3,3 * this -> model -> GetN_vertices());
	arma::mat dIdC = arma::zeros<arma::mat>(6,3 * this -> model -> GetN_vertices());



	#pragma omp parallel for reduction(+:dVdC) reduction(+:dCOMdC,dIdC)
	for (int f = 0; f < this -> model -> GetN_facets(); ++f){

		arma::sp_mat connect_table = this -> PartialTfPartialC(f);

		dVdC += this -> PartialDeltaVfPartialTf(f) * connect_table;
		dCOMdC += this -> PartialDeltaComPartialTf(f) * connect_table;
		dIdC += this -> PartialDeltaIfPartialTf(f) * connect_table;

	}

	this -> precomputed_partialVpartialC.clear();
	this -> precomputed_partialGpartialC.clear();
	this -> precomputed_partialSigmapartialC.clear();
	this -> precomputed_partialIpartialC.clear();

	this -> precomputed_partialVpartialC = dVdC;
	this -> precomputed_partialGpartialC = dCOMdC;
	this -> precomputed_partialIpartialC = dIdC;
	this -> precomputed_partialUnitDensityMomentsPartialI = this -> PartialUnitDensityMomentsPartialI();
	this -> precomputed_partialSigmapartialI = this -> PartialSigmaPartialI();
	this -> precomputed_partialSigmapartialC = this -> precomputed_partialSigmapartialI * precomputed_partialIpartialC;


}


void SBGATMassPropertiesUQ::RunMCUQVolumeCOMInertia(std::string path_to_shape,
	const double & density,
	const bool & shape_in_meters,
	const arma::mat & C_CC,
	const unsigned int & N_samples,
	std::string output_dir,
	int N_saved_shapes,
	arma::mat & deviations,
	arma::vec & all_volumes,
	arma::mat &  all_com,
	arma::mat & all_inertia){



	// Reading
	vtkSmartPointer<vtkOBJReader> reader_mc = vtkSmartPointer<vtkOBJReader>::New();
	reader_mc -> SetFileName(path_to_shape.c_str());
	reader_mc -> Update(); 

		// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cleaner_mc = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner_mc -> SetInputConnection (reader_mc -> GetOutputPort());
	cleaner_mc -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );

	cleaner_mc -> Update();




	if (all_volumes.size() != N_samples){
		all_volumes.clear();
		all_volumes = arma::vec (N_samples);
	}
	if (deviations.size() != N_samples){
		deviations.clear();
		deviations = arma::mat(3 * cleaner_mc -> GetOutput() -> GetNumberOfPoints(), N_samples);
	}
	if (all_com.size() != N_samples){
		all_com.clear();
		all_com = arma::mat(3,N_samples);
	}
	if (all_inertia.size() != N_samples){
		all_inertia.clear();
		all_inertia = arma::mat(6,N_samples);
	}
	


	for (unsigned int i = 0; i < N_samples ; ++i){
		deviations.col(i) = C_CC * arma::randn<arma::vec>(3 * cleaner_mc -> GetOutput() -> GetNumberOfPoints());
	}


	#pragma omp parallel for
	for (unsigned int i = 0; i < N_samples ; ++i){

		vtkSmartPointer<vtkPolyData> shape_copy = vtkSmartPointer<vtkPolyData>::New();
		shape_copy -> DeepCopy(cleaner_mc -> GetOutput());

		// Creating the PGM dyads
		vtkSmartPointer<SBGATMassProperties> mass_filter_mc = vtkSmartPointer<SBGATMassProperties>::New();
		mass_filter_mc -> SetInputData(shape_copy);
		mass_filter_mc -> SetDensity(density); 
		
		if (shape_in_meters){
			mass_filter_mc -> SetScaleMeters();
		}
		else{
			mass_filter_mc -> SetScaleKiloMeters();
		}

		SBGATMassPropertiesUQ shape_uq_mc;
		shape_uq_mc. SetModel(mass_filter_mc);
		shape_uq_mc. ApplyDeviation(deviations.col(i));


		all_volumes(i) = mass_filter_mc -> GetVolume();
		all_com.col(i) = mass_filter_mc -> GetCenterOfMass();

		arma::mat I = mass_filter_mc -> GetUnitDensityInertiaTensor() ;
		arma::vec::fixed<6> I_param = {I(0,0),I(1,1),I(2,2),I(0,1),I(0,2),I(1,2)};
		all_inertia.col(i) = I_param;


		if (i < N_saved_shapes){
			shape_uq_mc. TakeAndSaveSlice(0,output_dir + "slice_x_" + std::to_string(i) + ".txt",0);
			shape_uq_mc. TakeAndSaveSlice(1,output_dir + "slice_y_" + std::to_string(i) + ".txt",0);
			shape_uq_mc. TakeAndSaveSlice(2,output_dir + "slice_z_" + std::to_string(i) + ".txt",0);
			vtkSmartPointer<SBGATObjWriter> writer = vtkSmartPointer<SBGATObjWriter>::New();

			writer -> SetInputData(mass_filter_mc -> GetInput());

			writer -> SetFileName(std::string({output_dir + "mc_shape_" + std::to_string(i) + ".obj"}).c_str());
			writer -> Update();	


		}
	}


}


void SBGATMassPropertiesUQ::GetMassPropertiesUncertainties(double & volume_variance,
	arma::mat::fixed<3,3> & cov_com,
	arma::mat::fixed<6,6> & cov_I,
	arma::mat::fixed<3,3> & cov_sigma) const{



	const arma::rowvec & partialVolumePartialC = this -> GetPartialVolumePartialC();
	const arma::mat & partialIPartialC = this -> GetPartialIPartialC();
	const arma::mat & partialComPartialC = this -> GetPartialComPartialC();
	const arma::mat & partialSigmaPartialC = this -> GetPartialSigmaPartialC();

	
	volume_variance = arma::dot(partialVolumePartialC.t(),this -> P_CC * partialVolumePartialC.t());
	cov_com = partialComPartialC * this -> P_CC * partialComPartialC.t();
	cov_I = partialIPartialC * this -> P_CC * partialIPartialC.t();


}





