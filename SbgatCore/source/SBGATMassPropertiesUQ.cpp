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

	return 1./6 * partial;

	
}

arma::mat::fixed<3,9> SBGATMassPropertiesUQ::PartialDeltaCMfPartialTf(const int & f) const{

	arma::mat::fixed<3,9> triple_I;
	triple_I.cols(0,2) = arma::eye<arma::mat>(3,3);
	triple_I.cols(3,5) = arma::eye<arma::mat>(3,3);
	triple_I.cols(6,8) = arma::eye<arma::mat>(3,3);

	return  1./4 * triple_I ;

}

void SBGATMassPropertiesUQ::TestPartials(std::string input,double tol){

	SBGATMassPropertiesUQ::TestPartialDeltaVfPartialTf(input,tol);
	SBGATMassPropertiesUQ::TestPartialDeltaIfPartialTf(input,tol);

	SBGATMassPropertiesUQ::TestPartialCOMPartialC(input,tol);

}

void SBGATMassPropertiesUQ::TestPartialDeltaVfPartialTf(std::string input,double tol){

	std::cout << "\t In TestPartialDeltaVfPartialTf ... ";
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
		mass_prop -> SetScaleMeters();
		mass_prop -> Update();

		SBGATMassPropertiesUQ shape_uq;
		shape_uq.SetMassProperties(mass_prop);

		int N_facets = vtkPolyData::SafeDownCast(mass_prop -> GetInput()) -> GetNumberOfCells();
		int f = arma::randi(arma::distr_param(0,N_facets - 1));

	// Nominal Volume
		double volume = shape_uq.GetMassProperties() -> GetDeltaV(f);

	// Deviation
		arma::vec::fixed<9> delta_Tf = 1e-3 * arma::randn<arma::vec>(9);

	// Linear dUf
		double dV_lin = arma::dot(shape_uq.PartialDeltaVfPartialTf(f) , delta_Tf);

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

void SBGATMassPropertiesUQ::TestPartialDeltaIfPartialTf(std::string input,double tol){

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
		mass_prop -> SetScaleMeters();
		mass_prop -> Update();

		SBGATMassPropertiesUQ shape_uq;
		shape_uq.SetMassProperties(mass_prop);

		int f = arma::randi(arma::distr_param(0,mass_prop -> GetN_facets() - 1));

	// Nominal 
		arma::vec::fixed<6> DeltaI = shape_uq.GetMassProperties() -> GetDeltaIf(f);

	// Deviation
		arma::vec::fixed<9> delta_Tf = 1e-3 * arma::randn<arma::vec>(9);

	// Linear
		arma::vec::fixed<6> dI_lin = shape_uq.PartialDeltaIfPartialTf(f) * delta_Tf;

	// Apply Tf deviation
		shape_uq. ApplyTfDeviation(delta_Tf,f);

	// Perturbed
		arma::vec::fixed<6> DeltaI_p = shape_uq.GetMassProperties() -> GetDeltaIf(f);
		
	// Non-linear
		arma::vec::fixed<6> dI = DeltaI_p - DeltaI;

		if(arma::norm(dI - dI_lin) / arma::norm(dI_lin) < tol){
			++successes;
		}



	}

	std::cout << "\t Passed TestPartialDeltaIfPartialTf with " << double(successes)/N * 100 << " \% of successes. \n";






}


void SBGATMassPropertiesUQ::TestPartialCOMPartialC(std::string input,double tol){


	std::cout << "\t In TestPartialCOMPartialC ... ";
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
		mass_prop -> SetScaleMeters();
		mass_prop -> Update();

		SBGATMassPropertiesUQ shape_uq;
		shape_uq.SetMassProperties(mass_prop);

	// Nominal Volume
		arma::vec::fixed<3> com = shape_uq.GetMassProperties() -> GetCenterOfMass();

	// Deviation
		arma::vec delta_C = 1e-3 * arma::randn<arma::vec>(3 * mass_prop -> GetN_vertices()  );

	// Linear dUf
		arma::vec::fixed<3> dcom_lin = shape_uq.GetPartialCOMPartialC() * delta_C;

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

	std::cout << "\t Passed TestPartialCOMPartialC with " << double(successes)/N * 100 << " \% of successes. \n";

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


arma::mat::fixed<6,9> SBGATMassPropertiesUQ::PartialDeltaIOverDeltaVPartialTf(const int & f) const {



}

void SBGATMassPropertiesUQ::TestPartialDeltaIOverDeltaVPartialTf(std::string input,double tol) {


	std::cout << "\t In TestPartialDeltaIOverDeltaVPartialTf ... ";
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
		mass_prop -> SetScaleMeters();
		mass_prop -> Update();

		SBGATMassPropertiesUQ shape_uq;
		shape_uq.SetMassProperties(mass_prop);

		int N_facets = vtkPolyData::SafeDownCast(mass_prop -> GetInput()) -> GetNumberOfCells();
		int f = arma::randi(arma::distr_param(0,N_facets - 1));

	// // Nominal 
	// 	arma::mat::fixed<3,3> deltaITimesVol = shape_uq.GetMassProperties() -> GetDeltaIOverDeltaV(f);

	// // Deviation
	// 	arma::vec::fixed<9> delta_Tf = 1e-3 * arma::randn<arma::vec>(9);

	// // Linear 
	// 	arma::mat::fixed<3,3> ddeltaITimesVol_lin = shape_uq.PartialDeltaIOverDeltaVPartialTf(f) * delta_Tf;

	// // Apply Tf deviation
	// 	shape_uq. ApplyTfDeviation(delta_Tf,f);

	// // Perturbed 
	// 	arma::mat::fixed<3,3> deltaITimesVol_p = shape_uq.GetMassProperties() -> GetDeltaIOverDeltaV(f);
		
	// // Non-linear 
	// 	arma::mat::fixed<3,3> ddeltaITimesVol = deltaITimesVol_p - deltaITimesVol;

	// 	if(arma::norm(arma::vectorise(dcom - dcom_lin)) / arma::norm(arma::vectorise(dcom_lin)) < tol){
	// 		++successes;
	// 	}

	}

	std::cout << "\t Passed TestPartialDeltaIOverDeltaVPartialTf with " << double(successes)/N * 100 << " \% of successes. \n";


}


arma::mat SBGATMassPropertiesUQ::GetPartialCOMPartialC() const{

	arma::mat partial = arma::zeros<arma::mat>(3,3 * this -> mass_prop -> GetN_vertices());

	const arma::vec::fixed<3> & com = this -> mass_prop -> GetCenterOfMass();

	#pragma omp parallel for reduction(+:partial)
	for (int f = 0; f < this -> mass_prop -> GetN_facets(); ++f){

		double DV_f = this -> mass_prop -> GetDeltaV(f);
		arma::vec::fixed<3> DCom_f = this -> mass_prop -> GetDeltaCM(f);
		partial += ((DCom_f - com) * this -> PartialDeltaVfPartialTf(f)  + DV_f * this -> PartialDeltaCMfPartialTf(f)) * this -> mass_prop -> PartialTfPartialC(f);

	}

	return 1./this -> mass_prop -> GetVolume() * partial;

}


void SBGATMassPropertiesUQ::GetPartialAllInertiaPartialC(arma::rowvec & dVdC,arma::mat & dCOMdC,
	arma::mat & dIdC) const{
	const arma::vec::fixed<3> & com = this -> mass_prop -> GetCenterOfMass();
	double volume =  this -> mass_prop -> GetVolume();
	
	#pragma omp parallel for reduction(+:dVdC) reduction(+:dCOMdC) reduction(+:dIdC)
	for (int f = 0; f < this -> mass_prop -> GetN_facets(); ++f){

		arma::sp_mat connect_table = this -> mass_prop -> PartialTfPartialC(f);
		arma::rowvec dVdTf = this -> PartialDeltaVfPartialTf(f);
		double DV_f = this -> mass_prop -> GetDeltaV(f);
		arma::vec::fixed<3> DCom_f = this -> mass_prop -> GetDeltaCM(f);

		dVdC += dVdTf * connect_table;
		dCOMdC += 1./volume * ((DCom_f - com) * dVdTf  + DV_f * this -> PartialDeltaCMfPartialTf(f)) * connect_table;

		dIdC += this -> PartialDeltaIfPartialTf(f) * connect_table;

	}

}


arma::mat::fixed<6,9> SBGATMassPropertiesUQ::PartialDeltaIfPartialTf(const int & f) const{

	arma::mat::fixed<6,9> partial;
	throw(std::runtime_error("PartialDeltaIPartialTf is not yet implemented"));
	return partial;
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



