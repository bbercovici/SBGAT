#include <SBGATMassPropertiesUQ.hpp>
#include <RigidBodyKinematics.hpp>
#include <vtkCleanPolyData.h>
#include <vtkOBJReader.h>


#pragma omp declare reduction( + : arma::mat : omp_out += omp_in ) \
initializer( omp_priv = arma::zeros<arma::mat>(omp_orig.n_rows,omp_orig.n_cols))


#pragma omp declare reduction( + : arma::rowvec : omp_out += omp_in ) \
initializer( omp_priv = arma::zeros<arma::rowvec>(omp_orig.n_cols))


void SBGATMassPropertiesUQ::SetMassProperties(vtkSmartPointer<SBGATMassProperties> mass_prop){
	this -> mass_prop = mass_prop;
}

arma::rowvec::fixed<9> SBGATMassPropertiesUQ::PartialDeltaVfPartialTf(const int & f) const{

	double r0[3];
	double r1[3];
	double r2[3];

	arma::rowvec::fixed<9> partial;
	this -> mass_prop -> GetVerticesInFacet(f,r0,r1,r2);
	arma::vec::fixed<3> C0 = {r0[0],r0[1],r0[2]};
	arma::vec::fixed<3> C1 = {r1[0],r1[1],r1[2]};
	arma::vec::fixed<3> C2 = {r2[0],r2[1],r2[2]};

	partial.subvec(0,2) = arma::cross(C1,C2).t();
	partial.subvec(3,5) = - C0.t() * RBK::tilde(C2);
	partial.subvec(6,8) = C0.t() * RBK::tilde(C1);

	return 1./6 * partial * std::pow(this -> mass_prop -> GetScaleFactor(),2);

	
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
	SBGATMassPropertiesUQ::TestGetPartialIPartialC(input,tol,shape_in_meters);

	SBGATMassPropertiesUQ::TestGetPartialAllInertiaPartialC(input,tol,shape_in_meters);
	SBGATMassPropertiesUQ::TestGetPartialAllInertiaPartialCVSStandalone(input,tol,shape_in_meters);

}



void SBGATMassPropertiesUQ::TestGetPartialAllInertiaPartialCVSStandalone(std::string input,double tol,bool shape_in_meters){
	std::cout << "\t In TestGetPartialAllInertiaPartialCVSStandalone ... \n";


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
	mass_prop -> SetScaleMeters();
	mass_prop -> Update();

	SBGATMassPropertiesUQ shape_uq;
	shape_uq.SetMassProperties(mass_prop);

	
	arma::rowvec dVdC_standalone = shape_uq.GetPartialVolumePartialC();
	arma::mat dComdC_standalone = shape_uq.GetPartialComPartialC();
	arma::mat dIdC_standalone = shape_uq.GetPartialIPartialC();



	arma::rowvec dVdC ;
	arma::mat dComdC;
	arma::mat dIdC;

	shape_uq.GetPartialAllInertiaPartialC(dVdC,dComdC,dIdC);

	dIdC_standalone.cols(0,5).print("dIdC_standalone: \n");
	dIdC_standalone.cols(0,5).print("dIdC: \n");

	arma::mat error = dIdC_standalone - dIdC;

	arma::vec error_index = arma::index_max(error,1);

	std::cout << dIdC_standalone.col(error_index(0)) << " | " << dIdC.col(error_index(0)) << std::endl;
	std::cout << dIdC_standalone.col(error_index(1)) << " | " << dIdC.col(error_index(1)) << std::endl;
	std::cout << dIdC_standalone.col(error_index(2)) << " | " << dIdC.col(error_index(2)) << std::endl;
	std::cout << dIdC_standalone.col(error_index(3)) << " | " << dIdC.col(error_index(3)) << std::endl;
	std::cout << dIdC_standalone.col(error_index(4)) << " | " << dIdC.col(error_index(4)) << std::endl;

	std::cout << dIdC_standalone.col(error_index(5)) << " | " << dIdC.col(error_index(5)) << std::endl;

	
	
	std::cout << arma::abs(dIdC_standalone - dIdC).max() << std::endl;

	assert(arma::abs(dVdC_standalone - dVdC).max() < 1e-10);
	assert(arma::abs(dComdC_standalone - dComdC).max() < 1e-10);
	assert(arma::abs(dIdC_standalone - dIdC).max() < 1e-10);


	std::cout << "\t Done running TestGetPartialAllInertiaPartialCVSStandalone ... \n";

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
	shape_uq.SetMassProperties(mass_prop);



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
		shape_uq.SetMassProperties(mass_prop);

		int N_facets = vtkPolyData::SafeDownCast(mass_prop -> GetInput()) -> GetNumberOfCells();
		
		arma::ivec f_vec = arma::randi<arma::ivec>(1,arma::distr_param(0,N_facets - 1));
		int f = f_vec(0);



	// Nominal Volume
		double volume = shape_uq.GetMassProperties() -> GetDeltaV(f);

	// Deviation
		arma::vec::fixed<9> delta_Tf = 1e-3 * arma::randn<arma::vec>(9) / mass_prop -> GetScaleFactor();

	// Linear dUf
		double dV_lin = arma::dot(shape_uq.PartialDeltaVfPartialTf(f) , mass_prop -> GetScaleFactor()* delta_Tf);

	// Apply Tf deviation
		shape_uq. ApplyTfDeviation(delta_Tf,f);

	// Perturbed Uf
		double volume_p = shape_uq.GetMassProperties() -> GetDeltaV(f);
		
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
		shape_uq.SetMassProperties(mass_prop);

	// Nominal Volume
		double volume = shape_uq.GetMassProperties() -> GetVolume();

	// Deviation
		arma::vec delta_C = 1e-3 * arma::randn<arma::vec>(3 * mass_prop -> GetN_vertices()  ) / mass_prop -> GetScaleFactor();
		
	// Linear dUf
		double dV_lin = mass_prop -> GetScaleFactor() * arma::dot(shape_uq.GetPartialVolumePartialC() , delta_C);

	// Apply Tf deviation
		shape_uq. ApplyDeviation(delta_C);

	// Perturbed Uf
		double volume_p = shape_uq.GetMassProperties() -> GetVolume();
		
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
		shape_uq.SetMassProperties(mass_prop);

	// Nominal Volume
		arma::vec::fixed<3> com = shape_uq.GetMassProperties() -> GetCenterOfMass();

	// Deviation
		arma::vec delta_C = 1e-3 * arma::randn<arma::vec>(3 * mass_prop -> GetN_vertices()  ) / mass_prop -> GetScaleFactor();

	// Linear dUf
		arma::vec::fixed<3> dcom_lin = shape_uq.GetPartialComPartialC() * delta_C * mass_prop -> GetScaleFactor() ;

	// Apply Tf deviation
		shape_uq. ApplyDeviation(delta_C);

	// Perturbed Uf
		arma::vec::fixed<3> com_p = shape_uq.GetMassProperties() -> GetCenterOfMass();
		
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
	this -> mass_prop -> GetIndicesVerticesInFacet(f,v0_index,v1_index,v2_index);

	double r0[3],r1[3],r2[3];
	vtkPolyData * polydata = vtkPolyData::SafeDownCast(this -> mass_prop -> GetInput());

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
	this -> mass_prop -> Modified();

	this -> mass_prop -> Update();

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
		shape_uq.SetMassProperties(mass_prop);

	// Nominal 
		arma::mat deltaI_mat = shape_uq.GetMassProperties() -> GetUnitDensityInertiaTensor() ;
		arma::vec::fixed<6> deltaI = {deltaI_mat(0,0),deltaI_mat(1,1),deltaI_mat(2,2),deltaI_mat(0,1),deltaI_mat(0,2),deltaI_mat(1,2)};

	// Deviation
		arma::vec delta_C = 1e-3 * arma::randn<arma::vec>(3 * mass_prop -> GetN_vertices()  ) / mass_prop -> GetScaleFactor();

	// Linear 
		arma::vec::fixed<6> ddeltaI_lin = mass_prop -> GetScaleFactor() * shape_uq.GetPartialIPartialC() * delta_C;

	// Apply deviation
		shape_uq. ApplyDeviation(delta_C);

	// Nominal 
		arma::mat deltaI_mat_p = shape_uq.GetMassProperties() -> GetUnitDensityInertiaTensor() ;
		arma::vec::fixed<6> deltaI_p = {deltaI_mat_p(0,0),deltaI_mat_p(1,1),deltaI_mat_p(2,2),deltaI_mat_p(0,1),deltaI_mat_p(0,2),deltaI_mat_p(1,2)};

	// Non-linear 
		arma::vec::fixed<6> ddeltaI = deltaI_p - deltaI;

		if(arma::norm(arma::vectorise(ddeltaI - ddeltaI_lin)) / arma::norm(arma::vectorise(ddeltaI_lin)) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestGetPartialIPartialC with " << double(successes)/N * 100 << " \% of successes. \n";


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
		shape_uq.SetMassProperties(mass_prop);

	// Nominal 
		double V = shape_uq.GetMassProperties() -> GetVolume();
		arma::vec::fixed<3> com = shape_uq.GetMassProperties() -> GetCenterOfMass();

		arma::mat deltaI_mat = shape_uq.GetMassProperties() -> GetUnitDensityInertiaTensor() ;
		arma::vec::fixed<6> deltaI = {deltaI_mat(0,0),deltaI_mat(1,1),deltaI_mat(2,2),deltaI_mat(0,1),deltaI_mat(0,2),deltaI_mat(1,2)};

	// Deviation
		arma::vec delta_C = 1e-3 * arma::randn<arma::vec>(3 * mass_prop -> GetN_vertices()  ) / mass_prop -> GetScaleFactor();

	// Linear 

		arma::rowvec dVdC;

		arma::mat dIdC,dComdC;

		shape_uq.GetPartialAllInertiaPartialC(dVdC,dComdC,dIdC);


		double dV_lin = arma::dot(dVdC,mass_prop -> GetScaleFactor() * delta_C);
		arma::vec::fixed<3> dcom_lin = dComdC * mass_prop -> GetScaleFactor() * delta_C;
		arma::vec::fixed<6> ddeltaI_lin = mass_prop -> GetScaleFactor() * dIdC * delta_C;

	// Apply deviation
		shape_uq. ApplyDeviation(delta_C);

		double V_p = shape_uq.GetMassProperties() -> GetVolume();
		arma::vec::fixed<3> com_p = shape_uq.GetMassProperties() -> GetCenterOfMass();


		double dV =V_p - V;


		arma::vec::fixed<3> dcom = com_p - com;


		arma::mat deltaI_mat_p = shape_uq.GetMassProperties() -> GetUnitDensityInertiaTensor() ;
		arma::vec::fixed<6> deltaI_p = {deltaI_mat_p(0,0),deltaI_mat_p(1,1),deltaI_mat_p(2,2),deltaI_mat_p(0,1),deltaI_mat_p(0,2),deltaI_mat_p(1,2)};



	// Non-linear 
		arma::vec::fixed<6> ddeltaI = deltaI_p - deltaI;

		if(arma::norm(arma::vectorise(ddeltaI - ddeltaI_lin)) / arma::norm(arma::vectorise(ddeltaI_lin)) < tol){
			
			if (std::abs(dV - dV_lin)/std::abs(dV_lin) < tol){

				if(arma::norm(arma::vectorise(dcom - dcom_lin)) / arma::norm(dcom_lin) < tol){
					++successes;

				}
			}
		}

	}

	std::cout << "\t Passed TestGetPartialAllInertiaPartialC with " << double(successes)/N * 100 << " \% of successes. \n";


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
		shape_uq.SetMassProperties(mass_prop);

		int N_facets = vtkPolyData::SafeDownCast(mass_prop -> GetInput()) -> GetNumberOfCells();
		arma::ivec f_vec = arma::randi<arma::ivec>(1,arma::distr_param(0,N_facets - 1));
		int f = f_vec(0);
	// Nominal 
		arma::mat deltaIf_mat = shape_uq.GetMassProperties() -> GetDeltaIOverDeltaV(f) * shape_uq.GetMassProperties() -> GetDeltaV(f) ;
		arma::vec::fixed<6> deltaIf = {deltaIf_mat(0,0),deltaIf_mat(1,1),deltaIf_mat(2,2),deltaIf_mat(0,1),deltaIf_mat(0,2),deltaIf_mat(1,2)};

	// Deviation
		arma::vec::fixed<9> delta_Tf = 1e-3 * arma::randn<arma::vec>(9) / mass_prop -> GetScaleFactor();

	// Linear 
		arma::vec::fixed<6> ddeltaIf_lin = mass_prop -> GetScaleFactor() * shape_uq.PartialDeltaIfPartialTf(f) * delta_Tf;

	// Apply Tf deviation
		shape_uq. ApplyTfDeviation(delta_Tf,f);

	// Nominal 
		arma::mat deltaIf_mat_p = shape_uq.GetMassProperties() -> GetDeltaIOverDeltaV(f) * shape_uq.GetMassProperties() -> GetDeltaV(f) ;
		arma::vec::fixed<6> deltaIf_p = {deltaIf_mat_p(0,0),deltaIf_mat_p(1,1),deltaIf_mat_p(2,2),deltaIf_mat_p(0,1),deltaIf_mat_p(0,2),deltaIf_mat_p(1,2)};

	// Non-linear 
		arma::vec::fixed<6> ddeltaIf = deltaIf_p - deltaIf;

		if(arma::norm(arma::vectorise(ddeltaIf - ddeltaIf_lin)) / arma::norm(arma::vectorise(ddeltaIf_lin)) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestPartialDeltaIfPartialTf with " << double(successes)/N * 100 << " \% of successes. \n";


}


arma::mat SBGATMassPropertiesUQ::GetPartialComPartialC() const{

	arma::mat partial = arma::zeros<arma::mat>(3,3 * this -> mass_prop -> GetN_vertices());


	#pragma omp parallel for reduction(+:partial)
	for (int f = 0; f < this -> mass_prop -> GetN_facets(); ++f){

		partial += this -> PartialDeltaComPartialTf(f) * this -> mass_prop -> PartialTfPartialC(f);

	}

	return partial ;

}


void SBGATMassPropertiesUQ::GetPartialAllInertiaPartialC(arma::rowvec & dVdC,arma::mat & dCOMdC,
	arma::mat & dIdC) const{


	dVdC = arma::zeros<arma::rowvec>(3 * this -> mass_prop -> GetN_vertices());
	dCOMdC = arma::zeros<arma::mat>(3,3 * this -> mass_prop -> GetN_vertices());
	dIdC = arma::zeros<arma::mat>(6,3 * this -> mass_prop -> GetN_vertices());

	#pragma omp parallel for reduction(+:dVdC) reduction(+:dCOMdC,dIdC)
	for (int f = 0; f < this -> mass_prop -> GetN_facets(); ++f){

		arma::sp_mat connect_table = this -> mass_prop -> PartialTfPartialC(f);

		dVdC += this -> PartialDeltaVfPartialTf(f) * connect_table;
		dCOMdC += this -> PartialDeltaComPartialTf(f) * connect_table;
		dIdC += this -> PartialDeltaIfPartialTf(f) * connect_table;

	}

}


arma::rowvec SBGATMassPropertiesUQ::GetPartialVolumePartialC() const{
	
	arma::rowvec partial = arma::zeros<arma::rowvec>(3 * this -> mass_prop -> GetN_vertices());


	#pragma omp parallel for reduction(+:partial) 
	for (int f = 0; f < this -> mass_prop -> GetN_facets(); ++f){

		partial += this -> PartialDeltaVfPartialTf(f) * this -> mass_prop -> PartialTfPartialC(f);

	}

	return partial;

}


arma::mat::fixed<3,9> SBGATMassPropertiesUQ::PartialDeltaComPartialTf(const int & f) const{


	return (1./this -> mass_prop -> GetVolume() * ((this -> mass_prop -> GetDeltaCM(f) 
		- this -> mass_prop -> GetCenterOfMass()) * this -> PartialDeltaVfPartialTf(f)  
	+ this -> mass_prop -> GetDeltaV(f) * SBGATMassPropertiesUQ::PartialDeltaCMfPartialTf()));

}


arma::rowvec::fixed<9> SBGATMassPropertiesUQ::PartialEqDeltaIfErPartialTf(const int & f,const int & q, const int & r,const arma::vec::fixed<9> & Tf) const{


	arma::vec::fixed<3> e_q = arma::zeros<arma::vec >(3);
	arma::vec::fixed<3> e_r = arma::zeros<arma::vec >(3);
	
	e_q(q) = 1;
	e_r(r) = 1;

	arma::rowvec::fixed<9> partial = (arma::dot(e_q,this -> mass_prop -> GetDeltaIOverDeltaV(f) * e_r) * this -> PartialDeltaVfPartialTf(f) 
		+ this -> mass_prop -> GetDeltaV(f) * this -> PartialEqDeltaIOverDeltaVfErPartialTf(e_q,e_r,Tf));

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



	return -1./20 * partial * this -> mass_prop -> GetScaleFactor();


}


arma::mat::fixed<6,9> SBGATMassPropertiesUQ::PartialDeltaIfPartialTf(const int & f) const{

	arma::mat::fixed<6,9> partial;

	double r0[3];
	double r1[3];
	double r2[3];

	this -> mass_prop -> GetVerticesInFacet(f,r0,r1,r2);

	arma::vec::fixed<9> Tf = {r0[0],r0[1],r0[2],r1[0],r1[1],r1[2],r2[0],r2[1],r2[2]};

	partial.row(0) = this -> PartialEqDeltaIfErPartialTf(f,0,0,Tf);
	partial.row(1) = this -> PartialEqDeltaIfErPartialTf(f,1,1,Tf);
	partial.row(2) = this -> PartialEqDeltaIfErPartialTf(f,2,2,Tf);
	partial.row(3) = this -> PartialEqDeltaIfErPartialTf(f,0,1,Tf);
	partial.row(4) = this -> PartialEqDeltaIfErPartialTf(f,0,2,Tf);
	partial.row(5) = this -> PartialEqDeltaIfErPartialTf(f,1,2,Tf);

	return partial ;
}



arma::mat SBGATMassPropertiesUQ::GetPartialIPartialC() const{

	arma::mat partial = arma::zeros<arma::mat>(6,3 * this -> mass_prop -> GetN_vertices());

	#pragma omp parallel for
	for (int f = 0; f < this -> mass_prop -> GetN_facets(); ++f){

		partial += this -> PartialDeltaIfPartialTf(f) * this -> mass_prop ->  PartialTfPartialC(f);

	}

	return partial * std::pow(this -> mass_prop -> GetScaleFactor(),4);

}





void SBGATMassPropertiesUQ::ApplyDeviation(const arma::vec & delta_C){

	vtkPolyData * polydata = vtkPolyData::SafeDownCast(this -> mass_prop -> GetInput());
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

	this -> mass_prop -> Modified();
	this -> mass_prop -> Update();

}






